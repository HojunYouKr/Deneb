#include "deneb_system_matrix.h"

#include "avocado.h"
#include "deneb_config_macro.h"
#include "deneb_data.h"
#include "deneb_equation.h"

namespace deneb {
void InitializeSystemMatrix(Mat& sysmat) {
  const int& num_states = DENEB_EQUATION->GetNumStates();
  const int& num_bases = DENEB_DATA->GetNumBases();
  const int SB = num_states * num_bases;

  const int& num_cells = DENEB_DATA->GetNumCells();
  MatCreate(MPI_COMM_WORLD, &sysmat);
  MatSetSizes(sysmat, num_cells * SB, num_cells * SB, PETSC_DECIDE,
              PETSC_DECIDE);
  MatSetBlockSize(sysmat, SB);
  MatSetFromOptions(sysmat);

  // Matrix preallocation
  std::vector<int> d_nnz_block(num_cells, SB);
  std::vector<int> o_nnz_block(num_cells, 0);

  const int& num_faces = DENEB_DATA->GetNumFaces();
  const int& num_inner_faces = DENEB_DATA->GetNumInnerFaces();
  const auto& face_owner_cell = DENEB_DATA->GetFaceOwnerCell();
  const auto& face_neighbor_cell = DENEB_DATA->GetFaceNeighborCell();
  for (int iface = 0; iface < num_inner_faces; iface++) {
    d_nnz_block[face_owner_cell[iface]] += SB;
    d_nnz_block[face_neighbor_cell[iface]] += SB;
  }
  for (int iface = num_inner_faces; iface < num_faces; iface++)
    o_nnz_block[face_owner_cell[iface]] += SB;

  std::vector<PetscInt> d_nnz(num_cells * SB, 0);
  std::vector<PetscInt> o_nnz(num_cells * SB, 0);
  for (int icell = 0; icell < num_cells; icell++) {
    for (int icont = 0; icont < SB; icont++)
      d_nnz[icell * SB + icont] += d_nnz_block[icell];
    for (int icont = 0; icont < SB; icont++)
      o_nnz[icell * SB + icont] += o_nnz_block[icell];
  }
  d_nnz_block.clear();
  o_nnz_block.clear();

  if (NDOMAIN == 1)
    MatSeqAIJSetPreallocation(sysmat, NULL, &d_nnz[0]);
  else
    MatMPIAIJSetPreallocation(sysmat, NULL, &d_nnz[0], NULL, &o_nnz[0]);
  d_nnz.clear();
  o_nnz.clear();

  // Matrix check
  std::vector<PetscReal> initial_values(SB * SB, 0.0);
  const auto& mat_index = DENEB_DATA->GetMatIndex();
  for (int icell = 0; icell < num_cells; icell++)
    MatSetValuesBlocked(sysmat, 1, &mat_index[icell], 1, &mat_index[icell],
                        &initial_values[0], INSERT_VALUES);
  for (int iface = 0; iface < num_inner_faces; iface++) {
    const int& owner_cell = face_owner_cell[iface];
    const int& neighbor_cell = face_neighbor_cell[iface];
    MatSetValuesBlocked(sysmat, 1, &mat_index[owner_cell], 1,
                        &mat_index[owner_cell], &initial_values[0],
                        INSERT_VALUES);
    MatSetValuesBlocked(sysmat, 1, &mat_index[owner_cell], 1,
                        &mat_index[neighbor_cell], &initial_values[0],
                        INSERT_VALUES);
    MatSetValuesBlocked(sysmat, 1, &mat_index[neighbor_cell], 1,
                        &mat_index[owner_cell], &initial_values[0],
                        INSERT_VALUES);
    MatSetValuesBlocked(sysmat, 1, &mat_index[neighbor_cell], 1,
                        &mat_index[neighbor_cell], &initial_values[0],
                        INSERT_VALUES);
  }
  for (int iface = num_inner_faces; iface < num_faces; iface++) {
    const int& owner_cell = face_owner_cell[iface];
    const int& neighbor_cell = face_neighbor_cell[iface];
    MatSetValuesBlocked(sysmat, 1, &mat_index[owner_cell], 1,
                        &mat_index[owner_cell], &initial_values[0],
                        INSERT_VALUES);
    MatSetValuesBlocked(sysmat, 1, &mat_index[owner_cell], 1,
                        &mat_index[neighbor_cell + num_cells],
                        &initial_values[0], INSERT_VALUES);
  }
  MatAssemblyBegin(sysmat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(sysmat, MAT_FINAL_ASSEMBLY);
}

GMRES::GMRES() { KSPCreate(MPI_COMM_WORLD, &krylov_solver_); };
GMRES::~GMRES() { KSPDestroy(&krylov_solver_); };

void GMRES::Initialize(const Mat& mat) {
  const auto& config = AVOCADO_CONFIG;
  const int num_search_direction =
      std::stoi(config->GetConfigValue(GMRES_SEARCH_DIRECTION));
  const double relative_error_tol =
      std::stod(config->GetConfigValue(GMRES_ERROR_TOL));
  const int maximum_iteration =
      std::stoi(config->GetConfigValue(GMRES_MAX_ITERATION));
  MASTER_MESSAGE("Number of GMRES Krylov vectors = " +
                 std::to_string(num_search_direction) + "\n");
  MASTER_MESSAGE("Relative error tolerance of GMRES = " +
                 std::to_string(relative_error_tol) + "\n");
  MASTER_MESSAGE("Maximum iterations of GMRES = " +
                 std::to_string(maximum_iteration) + "\n");

  // Krylov solver
  const double absolute_error_tol = 1.0e-50;
  const double diverge_tol = 1.0e5;
  const PetscBool nonzero_guess = PETSC_TRUE;

  KSPSetFromOptions(krylov_solver_);
  KSPSetType(krylov_solver_, KSPGMRES);
  KSPSetInitialGuessNonzero(krylov_solver_, nonzero_guess);
  KSPGMRESSetRestart(krylov_solver_, num_search_direction);
  KSPGMRESSetOrthogonalization(krylov_solver_,
                               KSPGMRESModifiedGramSchmidtOrthogonalization);
  KSPSetTolerances(krylov_solver_, relative_error_tol, absolute_error_tol,
                   diverge_tol, maximum_iteration);

  // Preconditioner
  const PCSide preconditioning_side = PC_RIGHT;
  const PCType preconditioner_type = PCBJACOBI;

  KSPGetPC(krylov_solver_, &preconditioner_);
  KSPSetPCSide(krylov_solver_, preconditioning_side);
  PCSetType(preconditioner_, preconditioner_type);
  KSPSetPC(krylov_solver_, preconditioner_);

  KSPSetOperators(krylov_solver_, mat, mat);
}
int GMRES::Solve(const Mat& A, const Vec& b, Vec& solution) {
  KSPSetOperators(krylov_solver_, A, A);
  KSPSolve(krylov_solver_, b, solution);
  return GetIterationNumber();
}
}  // namespace deneb