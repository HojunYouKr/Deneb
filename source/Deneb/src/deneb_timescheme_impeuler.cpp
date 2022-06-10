#include "deneb_timescheme_impeuler.h"

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
// ----------------------- ImpEuler ---------------------- //
TimeschemeImpEuler::TimeschemeImpEuler() : Timescheme(true) {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeImpEuler"));
  MASTER_MESSAGE("Implicit = " + std::string(is_implicit_ ? "true" : "false") +
                 "\n");
};
TimeschemeImpEuler::~TimeschemeImpEuler() {
  VecDestroy(&rhs_);
  VecDestroy(&delta_);
};

void TimeschemeImpEuler::BuildData(void) {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeImpEuler::BuildData()"));
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
void TimeschemeImpEuler::Marching(void) {
  MASTER_MESSAGE(avocado::GetTitle("TimeschemeImpEuler::Marching()"));
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
  double *rhs_ptr, *delta_ptr;
  double error_norm;

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
    DENEB_ARTIFICIAL_VISCOSITY->ComputeArtificialViscosity(solution); 
    VecGetArray(rhs_, &rhs_ptr);
    DENEB_EQUATION->ComputeRHS(solution, rhs_ptr, t);
    cblas_dscal(length_, -1.0, rhs_ptr, 1);
    VecRestoreArray(rhs_, &rhs_ptr);
    DENEB_EQUATION->ComputeSystemMatrix(solution, sysmat_, t);
    MatShift(sysmat_, 1.0 / dt);

    const int sub_iteration = solver_.Solve(sysmat_, rhs_, delta_);

    // Check divergence
    VecNorm(delta_, NORM_2, &error_norm);
    if ((!std::isnormal(error_norm)) || (sub_iteration == 0)) {
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
    const double cost = STOP_TIMER();
    computing_cost_ += cost;

    // Printing message
    std::stringstream ss;
    ss << "Iter=" << std::setw(2) << iteration_;
    ss << " |  PhyT=" << std::scientific << std::setprecision(6)
       << GetCurrentTime();
    ss << " |  ComT=" << std::fixed << std::setprecision(2) << computing_cost_;
    ss << " |  ComT/iter=" << std::scientific << std::setprecision(4) << cost;
    ss << " |  GMRES subiter=" << sub_iteration;
    ss << " |  delta norm=" << std::scientific << std::setprecision(3)
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
}  // namespace deneb