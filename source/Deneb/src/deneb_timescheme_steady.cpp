#include "deneb_timescheme_steady.h"

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
// ------------------------ Steady ----------------------- //
TimeschemeSteady::TimeschemeSteady() : Timescheme(true) {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeSteady"));
  MASTER_MESSAGE("Implicit = " + std::string(is_implicit_ ? "true" : "false") +
                 "\n");

  const auto& config = AVOCADO_CONFIG;
  const int num_steps =
      std::stoi(config->GetConfigValue(CFL_INCREASE_INTERVAL));
  const double amount = std::stod(config->GetConfigValue(CFL_INCREASE_AMOUNT));
  increase_amount_ = std::pow(amount, 1.0 / static_cast<double>(num_steps));
  MASTER_MESSAGE("CFL number increases " + std::to_string(amount) +
                 " times at every " + std::to_string(num_steps) + " steps.\n");
  MASTER_MESSAGE(
      "\t CFL multiply factor = " + std::to_string(increase_amount_) + "\n");

  steady_stop_residual_ =
      std::stod(config->GetConfigValue(STEADY_CONVERGENCE_TOL));
  MASTER_MESSAGE("Steady stop condition: residual < " +
                 std::to_string(steady_stop_residual_) + "\n");

  if (timestep_control_ == TimestepControl::DT) {
    ERROR_MESSAGE(
        "Wrong timestep control: dt=const\n\tPlease use the CFL control.\n");
  }
  std::string text = config->GetConfigValue(STOP_CONDITION);
  if (!text.compare("Time")) {
    ERROR_MESSAGE(
        "Wrong stop condition. Please use the iteration condition.\n");
  }
  text = config->GetConfigValue(POST_CONTROL);
  if (!text.compare("Time")) {
    ERROR_MESSAGE("Wrong post control. Please use the iteration control.\n");
  }
  text = config->GetConfigValue(SAVE_CONTROL);
  if (!text.compare("Time")) {
    ERROR_MESSAGE("Wrong save control. Please use the iteration control.\n");
  }
};
TimeschemeSteady::~TimeschemeSteady() {
  VecDestroy(&rhs_);
  VecDestroy(&delta_);
};

void TimeschemeSteady::BuildData(void) {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeSteady::BuildData()"));
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

  InitSolution();
}
void TimeschemeSteady::Marching(void) {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeSteady::Marching()"));
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
      &solution_[0], iteration_);
  DENEB_CONTOUR->CellSolution(
      dir + RETURN_POST_DIR + "Iter" + std::to_string(iteration_) + ".plt",
      &solution_[0], iteration_);

  std::string filename =
      dir + RETURN_POST_DIR + "history" + std::to_string(iteration_) + ".dat";
  std::ofstream history;
  if (MYRANK == MASTER_NODE) {
    history.open(filename);
    history << "variables = \"iteration\", \"CFL\", \"ComT\", \"ComT_iter\", "
               "\"GMRES subiter\", \"Residual\"\n";
    history.close();
  }

  bool is_stop, is_post, is_save;
  const double* solution = &solution_[0];
  double *rhs_ptr, *delta_ptr;
  double delta_norm, delta_norm_prev = 0.0;

  const int& num_states = DENEB_EQUATION->GetNumStates();
  const int& num_bases = DENEB_DATA->GetNumBases();
  const int& num_cells = DENEB_DATA->GetNumCells();
  const auto& mat_index = DENEB_DATA->GetMatIndex();
  const int sb = num_states * num_bases;
  std::vector<double> block(sb * sb);

  while (iteration_ < max_iteration_) {
    START_TIMER();
    // Computing time step
    unsigned __int64 time_step = 0;
    is_stop = stop_.CheckTimeFinish(current_time_, time_step);
    is_post = post_.CheckTimeEvent(current_time_, time_step);
    is_save = save_.CheckTimeEvent(current_time_, time_step);

    ComputeLocalTimestep(solution, local_timestep_);

    // Updating solution
    DENEB_EQUATION->PreProcess(solution);
    DENEB_ARTIFICIAL_VISCOSITY->ComputeArtificialViscosity(solution);
    VecGetArray(rhs_, &rhs_ptr);
    DENEB_EQUATION->ComputeRHS(solution, rhs_ptr, 0.0);
    cblas_dscal(length_, -1.0, rhs_ptr, 1);
    VecRestoreArray(rhs_, &rhs_ptr);
    DENEB_EQUATION->ComputeSystemMatrix(solution, sysmat_, 0.0);
    for (int icell = 0; icell < num_cells; icell++) {
      memset(&block[0], 0, sb * sb * sizeof(double));
      const double dt_factor = 1.0 / local_timestep_[icell];

      for (int i = 0; i < sb * sb; i += (sb + 1)) block[i] = dt_factor;

      MatSetValuesBlocked(sysmat_, 1, &mat_index[icell], 1, &mat_index[icell],
                          &block[0], ADD_VALUES);
    }
    MatAssemblyBegin(sysmat_, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(sysmat_, MAT_FINAL_ASSEMBLY);

    const int sub_iteration = solver_.Solve(sysmat_, rhs_, delta_);

    // Check divergence
    VecNorm(delta_, NORM_2, &delta_norm);
    if ((!std::isnormal(delta_norm)) || (sub_iteration == 0)) {
      MASTER_MESSAGE("GMRES diverged!\n");
      is_stop = true;
    } else {
      VecGetArray(delta_, &delta_ptr);
      cblas_daxpy(length_, 1.0, delta_ptr, 1, &solution_[0], 1);
      VecRestoreArray(delta_, &delta_ptr);
      DENEB_LIMITER->Limiting(&solution_[0]);     
      DENEB_PRESSUREFIX->Execute(&solution_[0]);  

      // Updating time and iteration
      current_time_ += time_step;
      iteration_++;
    }

    // Measuring computing cost
    const double CFL = timestep_control_value_;
    const double cost = STOP_TIMER();
    computing_cost_ += cost;
    
    if (delta_norm < steady_stop_residual_) is_stop = true;
    if (delta_norm < delta_norm_prev && sub_iteration < 600) {
      timestep_control_value_ *= increase_amount_;
    } else {
      if (delta_norm > delta_norm_prev*1.01)
      timestep_control_value_ /= increase_amount_;
    }
    delta_norm_prev = delta_norm;

    // Printing message
    std::stringstream ss;
    ss << "Iter=" << std::setw(2) << iteration_;
    ss << " |  CFL=" << std::scientific << std::setprecision(3) << CFL;
    ss << " |  ComT=" << std::fixed << std::setprecision(2) << computing_cost_;
    ss << " |  ComT/iter=" << std::scientific << std::setprecision(4) << cost;
    ss << " |  GMRES subiter=" << sub_iteration;
    ss << " |  delta norm=" << std::scientific << std::setprecision(3)
       << delta_norm;
    MASTER_MESSAGE(ss.str() + "\n");

    if (MYRANK == MASTER_NODE) {
      history.open(filename, std::ios::app);
      history << std::scientific << std::setprecision(5);
      history << iteration_ << "\t";
      history << CFL << "\t";
      history << computing_cost_ << "\t";
      history << cost << "\t";
      history << sub_iteration << "\t";
      history << delta_norm << "\n";
      history.close();
    }

    // Checking interruption
    is_stop = is_stop || stop_.CheckIterationFinish(iteration_);
    is_post = is_post || post_.CheckIterationEvent(iteration_);
    is_save = is_save || save_.CheckIterationEvent(iteration_);
    
    if (is_post) {
      DENEB_CONTOUR->FaceSolution(dir + RETURN_POST_DIR + "Face/Iter" +
                                      std::to_string(iteration_) + ".plt",
                                  &solution_[0], iteration_);
      DENEB_CONTOUR->CellSolution(
          dir + RETURN_POST_DIR + "Iter" + std::to_string(iteration_) + ".plt",
          &solution_[0], iteration_);
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
      &solution_[0], iteration_);
  DENEB_CONTOUR->CellSolution(
      dir + RETURN_POST_DIR + "Iter" + std::to_string(iteration_) + ".plt",
      &solution_[0], iteration_);
}
}  // namespace deneb