#pragma once

#include <petsc.h>

namespace deneb {
void InitializeSystemMatrix(Mat& sysmat);

class GMRES {
 private:
  KSP krylov_solver_;
  PC preconditioner_;

 public:
  GMRES();
  ~GMRES();

  inline int GetIterationNumber() const {
    int iteration = 0;
    KSPGetIterationNumber(krylov_solver_, &iteration);
    return iteration;
  }

  void Initialize(const Mat& mat);
  int Solve(const Mat& A, const Vec& b, Vec& solution);
};
}  // namespace deneb