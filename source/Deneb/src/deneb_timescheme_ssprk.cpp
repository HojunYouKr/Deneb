#include "deneb_timescheme_ssprk.h"

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
// ------------------------ SSPRK ------------------------ //
TimeschemeSSPRK::TimeschemeSSPRK() : Timescheme(false){};
TimeschemeSSPRK::~TimeschemeSSPRK(){};
void TimeschemeSSPRK::BuildData(void) {
  const int& num_cells = DENEB_DATA->GetNumCells();
  const int& num_states = DENEB_EQUATION->GetNumStates();
  const int& num_bases = DENEB_DATA->GetNumBases();
  length_ = num_cells * num_states * num_bases;

  rhs_.resize(length_);
  u_.resize(length_);
  local_timestep_.resize(num_cells);
  computing_cost_ = 0.0;

  InitSolution();
}
void TimeschemeSSPRK::Marching(void) {
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
  double* solution = &solution_[0];

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

    // Updating solution
    DENEB_EQUATION->PreProcess(solution);
    UpdateSolution(solution, t, dt);

    // Updating time and iteration
    current_time_ += time_step;
    iteration_++;

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
  {
    SaveLoad::SaveData data;
    data.iteration_ = iteration_;
    data.strandid_ = DENEB_CONTOUR->GetStrandID();
    data.time_ = current_time_;
    DENEB_SAVELOAD->Save(
        dir + RETURN_SAVE_DIR + "Iter" + std::to_string(iteration_) + ".SAVE",
        &solution_[0], data);
  }
}

// ----------------------- SSPRK11 ----------------------- //
TimeschemeSSPRK11::TimeschemeSSPRK11() : TimeschemeSSPRK() {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeSSPRK11"));
  MASTER_MESSAGE("Implicit = " + std::string(is_implicit_ ? "true" : "false") +
                 "\n");
  c_.resize(1);
  c_[0] = 1.0;
};
TimeschemeSSPRK11::~TimeschemeSSPRK11(){};
void TimeschemeSSPRK11::BuildData(void) {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeSSPRK11::BuildData()"));
  TimeschemeSSPRK::BuildData();
}
void TimeschemeSSPRK11::Marching(void) {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeSSPRK11::Marching()"));
  TimeschemeSSPRK::Marching();
}
void TimeschemeSSPRK11::UpdateSolution(double* solution, const double t,
                                       const double dt) {
  // SSPRK11 first stage
  DENEB_ARTIFICIAL_VISCOSITY->ComputeArtificialViscosity(solution); 
  DENEB_EQUATION->ComputeRHS(solution, &rhs_[0], t + c_[0] * dt);
  for (int i = 0; i < length_; i++) solution[i] = solution[i] - dt * rhs_[i];
  DENEB_LIMITER->Limiting(solution);    
  DENEB_PRESSUREFIX->Execute(solution); 
}

// ----------------------- SSPRK33 ----------------------- //
TimeschemeSSPRK33::TimeschemeSSPRK33() : TimeschemeSSPRK() {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeSSPRK33"));
  MASTER_MESSAGE("Implicit = " + std::string(is_implicit_ ? "true" : "false") +
                 "\n");
  c_.resize(3);
  c_[0] = 0.0;
  c_[1] = 1.0;
  c_[2] = 0.5;
};
TimeschemeSSPRK33::~TimeschemeSSPRK33(){};
void TimeschemeSSPRK33::BuildData(void) {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeSSPRK33::BuildData()"));
  TimeschemeSSPRK::BuildData();
}
void TimeschemeSSPRK33::Marching(void) {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeSSPRK33::Marching()"));
  TimeschemeSSPRK::Marching();
}
void TimeschemeSSPRK33::UpdateSolution(double* solution, const double t,
                                       const double dt) {
  // SSPRK33 first stage
  DENEB_ARTIFICIAL_VISCOSITY->ComputeArtificialViscosity(solution);
  DENEB_EQUATION->ComputeRHS(solution, &rhs_[0], t + c_[0] * dt);
  ptr_ = solution;
  solution = &u_[0];
  for (int i = 0; i < length_; i++) solution[i] = ptr_[i] - dt * rhs_[i];
  DENEB_LIMITER->Limiting(solution);
  DENEB_PRESSUREFIX->Execute(solution); 

  // SSPRK33 second stage
  DENEB_ARTIFICIAL_VISCOSITY->ComputeArtificialViscosity(solution);
  DENEB_EQUATION->ComputeRHS(solution, &rhs_[0], t + c_[1] * dt);
  for (int i = 0; i < length_; i++)
    solution[i] = 0.25 * (3.0 * ptr_[i] + solution[i] - dt * rhs_[i]);
  DENEB_LIMITER->Limiting(solution);
  DENEB_PRESSUREFIX->Execute(solution); 

  // SSPRK33 third stage
  DENEB_ARTIFICIAL_VISCOSITY->ComputeArtificialViscosity(solution);
  DENEB_EQUATION->ComputeRHS(solution, &rhs_[0], t + c_[2] * dt);
  for (int i = 0; i < length_; i++)
    ptr_[i] = (ptr_[i] + 2.0 * solution[i] - 2.0 * dt * rhs_[i]) * c13_;
  solution = ptr_;
  DENEB_LIMITER->Limiting(solution);
  DENEB_PRESSUREFIX->Execute(solution);
}

// ----------------------- SSPRK54 ----------------------- //
TimeschemeSSPRK54::TimeschemeSSPRK54() : TimeschemeSSPRK() {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeSSPRK54"));
  MASTER_MESSAGE("Implicit = " + std::string(is_implicit_ ? "true" : "false") +
                 "\n");
  c_.resize(5);
  c_[0] = 0.0;
  c_[1] = 3.917522265718900e-01;
  c_[2] = 5.860796893115398e-01;
  c_[3] = 4.745423631214000e-01;
  c_[4] = 9.350106309676530e-01;
};
TimeschemeSSPRK54::~TimeschemeSSPRK54(){};
void TimeschemeSSPRK54::BuildData(void) {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeSSPRK54::BuildData()"));
  TimeschemeSSPRK::BuildData();
  rhs2_.resize(length_);
  u2_.resize(length_);
}
void TimeschemeSSPRK54::Marching(void) {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeSSPRK54::Marching()"));
  TimeschemeSSPRK::Marching();
}
void TimeschemeSSPRK54::UpdateSolution(double* solution, const double t,
                                       const double dt) {
  // SSPRK54 first stage
  DENEB_ARTIFICIAL_VISCOSITY->ComputeArtificialViscosity(solution);
  DENEB_EQUATION->ComputeRHS(solution, &rhs_[0], t + c_[0] * dt);
  ptr_ = solution;
  solution = &u_[0];
  for (int i = 0; i < length_; i++)
    solution[i] = ptr_[i] - 0.391752226571890 * dt * rhs_[i];
  DENEB_LIMITER->Limiting(solution);
  DENEB_PRESSUREFIX->Execute(solution);

  // SSPRK54 second stage
  DENEB_ARTIFICIAL_VISCOSITY->ComputeArtificialViscosity(solution);
  DENEB_EQUATION->ComputeRHS(solution, &rhs_[0], t + c_[1] * dt);
  for (int i = 0; i < length_; i++)
    solution[i] = 0.444370493651235 * ptr_[i] +
                  0.555629506348765 * solution[i] -
                  0.368410593050371 * dt * rhs_[i];
  DENEB_LIMITER->Limiting(solution);
  DENEB_PRESSUREFIX->Execute(solution);

  // SSPRK54 third stage
  DENEB_ARTIFICIAL_VISCOSITY->ComputeArtificialViscosity(solution);
  DENEB_EQUATION->ComputeRHS(solution, &rhs_[0], t + c_[2] * dt);
  ptr2_ = solution;
  solution = &u2_[0];
  for (int i = 0; i < length_; i++)
    solution[i] = 0.620101851488403 * ptr_[i] + 0.379898148511597 * ptr2_[i] -
                  0.251891774271694 * dt * rhs_[i];
  DENEB_LIMITER->Limiting(solution);
  DENEB_PRESSUREFIX->Execute(solution);

  // SSPRK54 fourth stage
  DENEB_ARTIFICIAL_VISCOSITY->ComputeArtificialViscosity(solution);
  DENEB_EQUATION->ComputeRHS(solution, &rhs2_[0], t + c_[3] * dt);
  solution = ptr_;
  ptr_ = &u2_[0];
  for (int i = 0; i < length_; i++)
    solution[i] = 0.178079954393132 * solution[i] +
                  0.821920045606868 * ptr_[i] -
                  0.544974750228521 * dt * rhs2_[i];
  DENEB_LIMITER->Limiting(solution);
  DENEB_PRESSUREFIX->Execute(solution);

  // SSPRK54 fifth stage
  DENEB_ARTIFICIAL_VISCOSITY->ComputeArtificialViscosity(solution);
  DENEB_EQUATION->ComputeRHS(solution, &rhs_[0], t + c_[4] * dt);
  for (int i = 0; i < length_; i++)
    solution[i] = 0.517231671970585 * ptr2_[i] + 0.096059710526147 * ptr_[i] -
                  0.063692468666290 * dt * rhs2_[i] +
                  0.386708617503269 * solution[i] -
                  0.226007483236906 * dt * rhs_[i];
  DENEB_LIMITER->Limiting(solution);
  DENEB_PRESSUREFIX->Execute(solution);
}
}  // namespace deneb