#include "deneb_timescheme_rosenbrock.h"

#include <iomanip>
#include <sstream>

#include "avocado.h"
#include "deneb_config_macro.h"
#include "deneb_contour.h"
#include "deneb_data.h"
#include "deneb_equation.h"
#include "deneb_limiter.h" 
#include "deneb_pressurefix.h" 
#include "deneb_artificial_viscosity.h" 
#include "deneb_saveload.h"

namespace deneb {
// -------------------- Rosenbrock R-K ------------------- //
TimeschemeRosenbrock::TimeschemeRosenbrock() : Timescheme(true){};
TimeschemeRosenbrock::~TimeschemeRosenbrock() {
  VecDestroy(&rhs_);
  VecDestroy(&delta_);
  for (auto vec : stage_solution_) VecDestroy(&vec);
  stage_solution_.clear();
};
void TimeschemeRosenbrock::BuildData(void) {
  const int& num_states = DENEB_EQUATION->GetNumStates();
  const int& num_bases = DENEB_DATA->GetNumBases();
  const int& num_cells = DENEB_DATA->GetNumCells();
  const int& num_global_cells = DENEB_DATA->GetNumGlobalCells();
  length_ = num_cells * num_states * num_bases;
  const int global_length = num_global_cells * num_states * num_bases;

  VecCreate(MPI_COMM_WORLD, &rhs_);
  VecSetSizes(rhs_, length_, global_length);
  VecSetFromOptions(rhs_);
  VecAssemblyBegin(rhs_);
  VecAssemblyEnd(rhs_);
  VecCreate(MPI_COMM_WORLD, &delta_);
  VecDuplicate(rhs_, &delta_);

  local_timestep_.resize(num_cells);
  computing_cost_ = 0.0;

  InitializeSystemMatrix(sysmat_);
  solver_.Initialize(sysmat_);

  stage_solution_.resize(num_stages_);
  for (int istage = 0; istage < num_stages_; istage++) {
    VecCreate(MPI_COMM_WORLD, &stage_solution_[istage]);
    VecDuplicate(delta_, &stage_solution_[istage]);
  }

  InitSolution();
}
void TimeschemeRosenbrock::Marching(void) {
  DENEB_LIMITER->Limiting(&solution_[0]); 
  DENEB_ARTIFICIAL_VISCOSITY->ComputeArtificialViscosity(
      &solution_[0]); 
  auto& config = AVOCADO_CONFIG;
  const std::string& dir = config->GetConfigValue(RETURN_DIR);
  DENEB_CONTOUR->FaceGrid(dir + RETURN_POST_DIR + "Face/Grid" +
                          std::to_string(iteration_) + ".plt");
  DENEB_CONTOUR->CellGrid(dir + RETURN_POST_DIR + "Grid" +
                          std::to_string(iteration_) + ".plt");
  DENEB_CONTOUR->FaceSolution(
      dir + RETURN_POST_DIR + "Face/Iter" + std::to_string(iteration_) + ".plt",
      &solution_[0], GetCurrentTime());
  DENEB_CONTOUR->CellSolution(
      dir + RETURN_POST_DIR + "Iter" + std::to_string(iteration_) + ".plt",
      &solution_[0], GetCurrentTime());

  bool is_stop, is_post, is_save;
  const double* solution = &solution_[0];
  double *rhs_ptr, *delta_ptr, *stage_solution_ptr;
  double error_norm;
  std::vector<double> rhs_dt(length_, 0.0);
  std::vector<double> pseudo_solution(length_, 0.0);
  std::vector<int> sub_iteration(num_stages_, 0);

  while (iteration_ < max_iteration_) {
    START_TIMER();
    // Computing time step
    unsigned __int64 time_step =
        ConvertTime(ComputeGlobalTimestep(solution, local_timestep_));
    is_stop = stop_.CheckTimeFinish(current_time_, time_step);
    is_post = post_.CheckTimeEvent(current_time_, time_step);
    is_save = save_.CheckTimeEvent(current_time_, time_step);

    const double t = GetCurrentTime();
    const double dt = ConvertTime(time_step);
    const double dt_inv = 1.0 / dt;

    // Updating solution
    DENEB_EQUATION->PreProcess(solution);

    DENEB_EQUATION->ComputeSystemMatrix(solution, sysmat_, t);
    MatShift(sysmat_, dt_inv / coeff_gamma_);

    const bool isRHSdt = DENEB_EQUATION->ComputeRHSdt(solution, &rhs_dt[0], t);

    for (int istage = 0; istage < num_stages_; istage++) {
      cblas_dcopy(length_, solution, 1, &pseudo_solution[0], 1);
      for (int jstage = 0; jstage < istage; jstage++) {
        VecGetArray(stage_solution_[jstage], &stage_solution_ptr);
        cblas_daxpy(length_, -coeff_a_[istage][jstage],
                    stage_solution_ptr, 1, &pseudo_solution[0], 1);
        VecRestoreArray(stage_solution_[jstage], &stage_solution_ptr);
      }

      DENEB_ARTIFICIAL_VISCOSITY->ComputeArtificialViscosity(
          &pseudo_solution[0]);                         
      DENEB_LIMITER->Limiting(&pseudo_solution[0]);     
      DENEB_PRESSUREFIX->Execute(&pseudo_solution[0]);  

      VecGetArray(rhs_, &rhs_ptr);
      DENEB_EQUATION->ComputeRHS(&pseudo_solution[0], rhs_ptr,
                                 t + coeff_alpha_[istage] * dt);
      for (int jstage = 0; jstage < istage; jstage++) {
        VecGetArray(stage_solution_[jstage], &stage_solution_ptr);
        cblas_daxpy(length_, coeff_c_[istage][jstage] * dt_inv,
                    stage_solution_ptr, 1, rhs_ptr, 1);
        VecRestoreArray(stage_solution_[jstage], &stage_solution_ptr);
      }

      if (isRHSdt)
        cblas_daxpy(length_, coeff_r_[istage] * dt, &rhs_dt[0], 1, rhs_ptr, 1);

      VecRestoreArray(rhs_, &rhs_ptr);
      sub_iteration[istage] =
          solver_.Solve(sysmat_, rhs_, stage_solution_[istage]);
    }
    VecCopy(stage_solution_[0], delta_);
    VecScale(delta_, coeff_m_[0]);
    for (int istage = 1; istage < num_stages_; istage++)
      VecAXPY(delta_, coeff_m_[istage], stage_solution_[istage]);

    // Check divergence
    VecNorm(delta_, NORM_2, &error_norm);
    double iteration_multiply = 1.0;
    for (auto& iter : sub_iteration) iteration_multiply *= iter;
    if ((!std::isnormal(error_norm)) || (std::abs(iteration_multiply)<1.0E-6)) {
      MASTER_MESSAGE("GMRES diverged!\n");
      is_stop = true;
    } else {
      VecGetArray(delta_, &delta_ptr);
      cblas_daxpy(length_, -1.0, delta_ptr, 1, &solution_[0], 1);
      VecRestoreArray(delta_, &delta_ptr);
      DENEB_LIMITER->Limiting(&solution_[0]);     
      DENEB_PRESSUREFIX->Execute(&solution_[0]);  

      // Updating time and iteration
      current_time_ += time_step;
      iteration_++;
    }

    // Measuring computing cost
    const double cost = STOP_TIMER();
    computing_cost_ += cost;

    // Printing message
    std::stringstream ss;
    ss << "Iter=" << std::setw(2) << iteration_;
    ss << " |  PhyT=" << std::scientific << std::setprecision(6)
       << GetCurrentTime();
    ss << " |  ComT=" << std::fixed << std::setprecision(2) << computing_cost_;
    ss << " |  ComT/iter=" << std::scientific << std::setprecision(4) << cost;
    ss << " |  GMRES subiter=(" << sub_iteration[0];
    for (int i = 1; i < num_stages_; i++) ss << ", " << sub_iteration[i];
    ss << ") |  delta norm=" << std::scientific << std::setprecision(3)
       << error_norm;
    MASTER_MESSAGE(ss.str() + "\n");

    // Checking interruption
    is_stop = is_stop || stop_.CheckIterationFinish(iteration_);
    is_post = is_post || post_.CheckIterationEvent(iteration_);
    is_save = is_save || save_.CheckIterationEvent(iteration_);
    
    if (is_post) {
      DENEB_CONTOUR->FaceSolution(dir + RETURN_POST_DIR + "Face/Iter" +
                                      std::to_string(iteration_) + ".plt",
                                  &solution_[0], GetCurrentTime());
      DENEB_CONTOUR->CellSolution(
          dir + RETURN_POST_DIR + "Iter" + std::to_string(iteration_) + ".plt",
          &solution_[0], GetCurrentTime());
    }
    if (is_save) {
      SaveLoad::SaveData data;
      data.iteration_ = iteration_;
      data.strandid_ = DENEB_CONTOUR->GetStrandID();
      data.time_ = current_time_;
      DENEB_SAVELOAD->Save(
          dir + RETURN_SAVE_DIR + "Iter" + std::to_string(iteration_) + ".SAVE",
          &solution_[0], data);
    }
    if (is_stop) break;
  }
  MASTER_MESSAGE("Computing cost: " + std::to_string(computing_cost_) + "\n");
  DENEB_CONTOUR->FaceSolution(
      dir + RETURN_POST_DIR + "Face/Iter" + std::to_string(iteration_) + ".plt",
      &solution_[0], GetCurrentTime());
  DENEB_CONTOUR->CellSolution(
      dir + RETURN_POST_DIR + "Iter" + std::to_string(iteration_) + ".plt",
      &solution_[0], GetCurrentTime());
}

// ------------------------ ROS3P ------------------------ //
// Lang-Verwer ROS3P (RO3-3)
TimeschemeROS3P::TimeschemeROS3P() : TimeschemeRosenbrock() {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeROS3P"));
  MASTER_MESSAGE("Implicit = " + std::string(is_implicit_ ? "true" : "false") +
                 "\n");

  // time order and the number of stages
  time_order_ = 3;
  num_stages_ = 3;

  coeff_a_.resize(num_stages_);
  coeff_c_.resize(num_stages_);

  // coefficients
  coeff_gamma_ = 0.5 + std::sqrt(3.0) / 6.0;
  coeff_a_[1] = {1.0 / coeff_gamma_};
  coeff_a_[2] = {1.0 / coeff_gamma_, 0.0};
  coeff_c_[1] = {-1.0 / pow(coeff_gamma_, 2.0)};
  coeff_c_[2] = {
      -(1.0 + (2.0 - 0.5 / coeff_gamma_) / coeff_gamma_) / coeff_gamma_,
      -(2.0 - 0.5 / coeff_gamma_) / coeff_gamma_};
  coeff_m_ = {(1.0 + (2.0 / 3.0 - 1.0 / (6.0 * coeff_gamma_)) / coeff_gamma_) /
                  coeff_gamma_,
              (2.0 / 3.0 - 1.0 / (6.0 * coeff_gamma_)) / coeff_gamma_,
              1.0 / (3.0 * coeff_gamma_)};

  coeff_r_ = {coeff_gamma_, coeff_gamma_ - 1.0, 0.5 - 2.0 * coeff_gamma_};
  coeff_alpha_ = {0.0, 1.0, 1.0};
}
TimeschemeROS3P::~TimeschemeROS3P(){};
void TimeschemeROS3P::BuildData(void) {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeROS3P::BuildData()"));
  TimeschemeRosenbrock::BuildData();
}
void TimeschemeROS3P::Marching(void) {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeROS3P::Marching()"));
  TimeschemeRosenbrock::Marching();
}

// ----------------------- RODAS3 ------------------------ //
// Hairer-Wanner RODAS3 (RO3-4)
TimeschemeRODAS3::TimeschemeRODAS3() : TimeschemeRosenbrock() {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeRODAS3"));
  MASTER_MESSAGE("Implicit = " + std::string(is_implicit_ ? "true" : "false") +
                 "\n");

  // time order and the number of stages
  time_order_ = 3;
  num_stages_ = 4;

  coeff_a_.resize(num_stages_);
  coeff_c_.resize(num_stages_);

  // coefficients
  coeff_gamma_ = 0.5;
  coeff_a_[1] = {0.0};
  coeff_a_[2] = {2.0, 0.0};
  coeff_a_[3] = {2.0, 0.0, 1.0};
  coeff_c_[1] = {4.0};
  coeff_c_[2] = {1.0, -1.0};
  coeff_c_[3] = {1.0, -1.0, -8.0 / 3.0};
  coeff_m_ = {2.0, 0.0, 1.0, 1.0};

  coeff_r_ = {0.5, 1.5, 0.0, 0.0};
  coeff_alpha_ = {0.0, 0.0, 1.0, 1.0};
}
TimeschemeRODAS3::~TimeschemeRODAS3(){};
void TimeschemeRODAS3::BuildData(void) {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeRODAS3::BuildData()"));
  TimeschemeRosenbrock::BuildData();
}
void TimeschemeRODAS3::Marching(void) {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeRODAS3::Marching()"));
  TimeschemeRosenbrock::Marching();
}

// ------------------------ ROS4 ------------------------- //
// Shampine ROS4 (RO4-4)
TimeschemeROS4::TimeschemeROS4() : TimeschemeRosenbrock() {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeROS4"));
  MASTER_MESSAGE("Implicit = " + std::string(is_implicit_ ? "true" : "false") +
                 "\n");

  // time order and the number of stages
  time_order_ = 4;
  num_stages_ = 4;

  coeff_a_.resize(num_stages_);
  coeff_c_.resize(num_stages_);

  // coefficients
  coeff_gamma_ = 0.5;
  coeff_a_[1] = {2.0};
  coeff_a_[2] = {48.0 / 25.0, 6.0 / 25.0};
  coeff_a_[3] = {48.0 / 25.0, 6.0 / 25.0, 0.0};
  coeff_c_[1] = {-8.0};
  coeff_c_[2] = {372.0 / 25.0, 12.0 / 5.0};
  coeff_c_[3] = {-112.0 / 125.0, -54.0 / 125.0, -2.0 / 5.0};
  coeff_m_ = {19.0 / 9.0, 0.5, 25.0 / 108.0, 125.0 / 108.0};

  coeff_r_ = {0.5, -1.5, 2.42, 0.116};
  coeff_alpha_ = {0.0, 1.0, 0.6, 0.6};
}
TimeschemeROS4::~TimeschemeROS4(){};
void TimeschemeROS4::BuildData(void) {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeROS4::BuildData()"));
  TimeschemeRosenbrock::BuildData();
}
void TimeschemeROS4::Marching(void) {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeROS4::Marching()"));
  TimeschemeRosenbrock::Marching();
}
// ----------------------- RODASP ------------------------ //
// Steinebach RODASP (RO4-6)
TimeschemeRODASP::TimeschemeRODASP() : TimeschemeRosenbrock() {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeRODASP"));
  MASTER_MESSAGE("Implicit = " + std::string(is_implicit_ ? "true" : "false") +
                 "\n");

  // time order and the number of stages
  time_order_ = 4;
  num_stages_ = 6;

  coeff_a_.resize(num_stages_);
  coeff_c_.resize(num_stages_);

  // coefficients
  coeff_gamma_ = 0.25;
  coeff_a_[1] = {3.0};
  coeff_a_[2] = {1.831036793486759, 0.4955183967433795};
  coeff_a_[3] = {2.304376582692669, -0.05249275245743001, -1.176798761832782};
  coeff_a_[4] = {-7.170454962423024, -4.741636671481785, -16.31002631330971,
                 -1.062004044111401};
  coeff_a_[5] = {-7.170454962423024, -4.741636671481785, -16.31002631330971,
                 -1.062004044111401, 1.0};
  coeff_c_[1] = {-12.0};
  coeff_c_[2] = {-8.791795173947035, -2.207865586973518};
  coeff_c_[3] = {10.81793056857153, 6.780270611428266, 19.53485944642410};
  coeff_c_[4] = {34.19095006749676, 15.49671153725963, 54.74760875964130,
                 14.16005392148534};
  coeff_c_[5] = {34.62605830930532, 15.30084976114473, 56.99955578662667,
                 18.40807009793095, -5.714285714285717};
  coeff_m_ = {-7.170454962423024,
              -4.741636671481785,
              -16.31002631330971,
              -1.062004044111401,
              1.0,
              1.0};

  coeff_r_ = {0.25, -0.5, -0.023504, -0.0362, 0.0, 0.0};
  coeff_alpha_ = {0.0, 0.75, 0.21, 0.63, 1.0, 1.0};
}
TimeschemeRODASP::~TimeschemeRODASP(){};
void TimeschemeRODASP::BuildData(void) {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeRODASP::BuildData()"));
  TimeschemeRosenbrock::BuildData();
}
void TimeschemeRODASP::Marching(void) {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeRODASP::Marching()"));
  TimeschemeRosenbrock::Marching();
}
// ----------------------- RODAS5 ------------------------ //
// Di Marzo RODAS5-Rod5_1 (RO5-8)
TimeschemeRODAS5::TimeschemeRODAS5() : TimeschemeRosenbrock() {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeRODAS5"));
  MASTER_MESSAGE("Implicit = " + std::string(is_implicit_ ? "true" : "false") +
                 "\n");

  // time order and the number of stages
  time_order_ = 5;
  num_stages_ = 8;

  coeff_a_.resize(num_stages_);
  coeff_c_.resize(num_stages_);

  // coefficients
  coeff_gamma_ = 0.19;
  coeff_a_[1] = {2.0};
  coeff_a_[2] = {3.040894194418781, 1.041747909077569};
  coeff_a_[3] = {2.576417536461461, 1.622083060776640, -0.9089668560264532};
  coeff_a_[4] = {2.760842080225597, 1.446624659844071, -0.3036980084553738,
                 0.2877498600325443};
  coeff_a_[5] = {-14.09640773051259, 6.925207756232704, -41.47510893210728,
                 2.343771018586405, 24.13215229196062};
  coeff_a_[6] = {-14.09640773051259, 6.925207756232704, -41.47510893210728,
                 2.343771018586405,  24.13215229196062, 1.0};
  coeff_a_[7] = {-14.09640773051259,
                 6.925207756232704,
                 -41.47510893210728,
                 2.343771018586405,
                 24.13215229196062,
                 1.0,
                 1.0};
  coeff_c_[1] = {-10.31323885133993};
  coeff_c_[2] = {-21.04823117650003, -7.234992135176716};
  coeff_c_[3] = {32.22751541853323, -4.943732386540191, 19.44922031041879};
  coeff_c_[4] = {-20.69865579590063, -8.816374604402768, 1.260436877740897,
                 -0.7495647613787146};
  coeff_c_[5] = {-46.22004352711257, -17.49534862857472, -289.6389582892057,
                 93.60855400400906, 318.3822534212147};
  coeff_c_[6] = {34.20013733472935, -14.15535402717690, 57.82335640988400,
                 25.83362985412365, 1.408950972071624,  -6.551835421242162};
  coeff_c_[7] = {42.57076742291101, -13.80770672017997, 93.98938432427124,
                 18.77919633714503, -31.58359187223370, -6.685968952921985,
                 -5.810979938412932};
  coeff_m_ = {-14.09640773051259,
              6.925207756232704,
              -41.47510893210728,
              2.343771018586405,
              24.13215229196062,
              1.0,
              1.0,
              1.0};

  coeff_r_ = {0.19,
              -0.182307922533373,
              -0.319231832186875,
              0.344982862472536,
              -0.377417564392090,
              0.0,
              0.0,
              0.0};
  coeff_alpha_ = {
      0.0, 0.38, 0.387850999832152, 0.483971893787382, 0.457047700881957, 1.0,
      1.0, 1.0};
}
TimeschemeRODAS5::~TimeschemeRODAS5(){};
void TimeschemeRODAS5::BuildData(void) {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeRODAS5::BuildData()"));
  TimeschemeRosenbrock::BuildData();
}
void TimeschemeRODAS5::Marching(void) {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeRODAS5::Marching()"));
  TimeschemeRosenbrock::Marching();
}
// ----------------------- ROW6A ------------------------- //
// Kaps and Wanner ROW6A (RO6-6)
TimeschemeROW6A::TimeschemeROW6A() : TimeschemeRosenbrock() {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeROW6A"));
  MASTER_MESSAGE("Implicit = " + std::string(is_implicit_ ? "true" : "false") +
                 "\n");

  // time order and the number of stages
  time_order_ = 6;
  num_stages_ = 6;

  coeff_a_.resize(num_stages_);
  coeff_c_.resize(num_stages_);

  // coefficients
  coeff_gamma_ = 0.3341423670680504;
  coeff_a_[1] = {2.0};
  coeff_a_[2] = {1.751493065942685, -0.1454290536332865};
  coeff_a_[3] = {-1.847093912231436, -2.513756792158473, 1.874707432337999};
  coeff_a_[4] = {10.59634783677141, 1.974951525952609, -1.905211286263863,
                 -3.575118228830491};
  coeff_a_[5] = {2.417642067883312, 0.3050984437044573, -0.2346208879122501,
                 -0.1327038464607418, 0.03912922779645768};
  coeff_c_[1] = {-17.45029492512995};
  coeff_c_[2] = {-12.02359936227844, 1.315910110742745};
  coeff_c_[3] = {23.11230597159272, 12.97893129565445, -8.445374594562038};
  coeff_c_[4] = {-3.147228891330713, -1.761332622909965, 6.115295934038585,
                 14.99319950457112};
  coeff_c_[5] = {-20.15840911262880, -1.603923799800133, 1.155870096920252,
                 0.6304639815292044, -0.1602510215637174};
  coeff_m_ = {33.99347452674165,  -20.91829882847333, -13.75688477471081,
              -11.13925929930077, 2.873406527609468,  38.76609945620840};

  coeff_r_ = {0.334142367068050, -1.614202631302161, -1.718073012358579,
              0.762494753419930, 1.242086134696671,  -1.620894493799871};
  coeff_alpha_ = {0.0, 0.668284734136100, 0.82, 0.219636250757926,
                  0.9, 0.665857632931949};
}
TimeschemeROW6A::~TimeschemeROW6A(){};
void TimeschemeROW6A::BuildData(void) {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeROW6A::BuildData()"));
  TimeschemeRosenbrock::BuildData();
}
void TimeschemeROW6A::Marching(void) {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeROW6A::Marching()"));
  TimeschemeRosenbrock::Marching();
}
}  // namespace deneb