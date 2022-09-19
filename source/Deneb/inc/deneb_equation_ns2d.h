#pragma once

#include <cmath>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "deneb_equation.h"

// 2-D compressible Navier-Stokes equations

namespace deneb {
class ProblemNS2D;
class BoundaryNS2D;

// ------------------------------- Constants ------------------------------- //
class ConstantsNS2D {
 protected:
  static constexpr int D_ = 2;
  static constexpr int S_ = 4;
  static constexpr int DS_ = D_ * S_;
  static constexpr int SS_ = S_ * S_;
  static constexpr int DSS_ = DS_ * S_;
  static constexpr int DDSS_ = D_ * DSS_;

  const double r_;
  const double r_inv_;
  const double rm1_;
  const double rm1_inv_;

  const double c23_ = 2.0 / 3.0;
  const double pi_ = 4.0 * std::atan(1.0);

  double Re_;
  double Ma_;
  double Pr_;
  double alpha_;
  double beta_;

  int max_num_points_;
  int max_num_cell_points_;
  int max_num_face_points_;
  int max_num_bdry_points_;

 public:
  ConstantsNS2D();
  virtual ~ConstantsNS2D(){};

  inline int GetMaxNumBdryPoints() const { return max_num_bdry_points_; };

 protected:
  double ComputePressure(const double* sol) const;
  double ComputeTotalEnergy(const double* pri) const;
};

// ------------------------------- Equation -------------------------------- //
class EquationNS2D : public Equation, public ConstantsNS2D {
 private:
  std::shared_ptr<ProblemNS2D> problem_;
  void (EquationNS2D::*compute_numflux_)(const int num_points,
                                         std::vector<double>& flux,
                                         FACE_INPUTS);
  void (EquationNS2D::*compute_numflux_jacobi_)(const int num_points,
                                                FACE_JACOBI_OUTPUTS,
                                                FACE_INPUTS);

  std::unordered_map<int, std::shared_ptr<BoundaryNS2D>> boundary_registry_;
  std::vector<std::shared_ptr<BoundaryNS2D>> boundaries_;

 public:
  EquationNS2D();
  virtual ~EquationNS2D();

  virtual void RegistBoundary(const std::vector<int>& bdry_tag);
  virtual void BuildData(void);
  virtual void PreProcess(const double* solution) { return; };
  virtual void ComputeRHS(const double* solution, double* rhs, const double t);
  virtual bool ComputeRHSdt(const double* solution, double* rhs_dt,
                            const double t) {
    return false;
  };
  virtual void ComputeSystemMatrix(const double* solution, Mat& sysmat,
                                   const double t);
  virtual void GetCellPostSolution(const int icell, const int num_points,
                                   const std::vector<double>& solution,
                                   const std::vector<double>& solution_grad,
                                   std::vector<double>& post_solution);
  virtual void GetFacePostSolution(const int num_points,
                                   const std::vector<double>& solution,
                                   const std::vector<double>& solution_grad,
                                   const std::vector<double>& normal,
                                   std::vector<double>& post_solution);
  virtual void ComputeInitialSolution(double* solution, const double t);
  virtual void ComputeLocalTimestep(
      const double* solution, std::vector<double>& local_timestep);
  virtual void ComputePressureCoefficient(
      std::vector<double>& pressure_coefficient, const int num_points,
      const std::vector<double>& solution){};
  virtual void ComputeViscousStress(std::vector<double>& viscous_stress,
                                    const int num_points,
                                    const std::vector<double>& solution,
                                    const std::vector<double>& solution_grad){};
  virtual bool IsContact(const int& icell,
                         const std::vector<int>& neighbor_cells,
                         const double* solution,
                         const double* total_solution) const;
  virtual double ComputeMaxCharacteristicSpeed(
      const double* input_solution) const;
  virtual const std::vector<double>& ComputePressureFixValues(
      const double* input_solution);

  void ComputeComFlux(const int num_points, std::vector<double>& flux,
                      const int icell,
                      const std::vector<double>& owner_u,
                      const std::vector<double>& owner_div_u);
  void ComputeComFluxJacobi(const int num_points,
                            std::vector<double>& flux_jacobi,
                            std::vector<double>& flux_grad_jacobi,
                            const int icell, const std::vector<double>& owner_u,
                            const std::vector<double>& owner_div_u);
  inline void ComputeNumFlux(const int num_points,
                             std::vector<double>& flux,
                             FACE_INPUTS) {
    (this->*compute_numflux_)(num_points, flux, owner_cell, neighbor_cell,
                              owner_u, owner_div_u, neighbor_u, neighbor_div_u,
                              normal);
  };
  inline void ComputeNumFluxJacobi(const int num_points,
                                   FACE_JACOBI_OUTPUTS,
                                   FACE_INPUTS) {
    (this->*compute_numflux_jacobi_)(
        num_points, flux_owner_jacobi, flux_neighbor_jacobi,
        flux_owner_grad_jacobi, flux_neighbor_grad_jacobi, owner_cell,
        neighbor_cell, owner_u, owner_div_u, neighbor_u, neighbor_div_u,
        normal);
  };

  DEFINE_FLUX(LLF);
  DEFINE_FLUX(Roe);
};

// ------------------------------- Boundary -------------------------------- //
class BoundaryNS2D : public ConstantsNS2D {
 public:
  static std::shared_ptr<BoundaryNS2D> GetBoundary(const std::string& type,
                                                   const int bdry_tag,
                                                   EquationNS2D* equation);

 protected:
  int bdry_tag_;
  EquationNS2D* equation_;

 public:
  BoundaryNS2D(const int bdry_tag, EquationNS2D* equation)
      : bdry_tag_(bdry_tag), equation_(equation){};
  virtual ~BoundaryNS2D(){};

  virtual void ComputeBdrySolution(
      const int num_points, std::vector<double>& bdry_u,
      std::vector<double>& bdry_div_u, const std::vector<double>& owner_u,
      const std::vector<double>& owner_div_u, const std::vector<double>& normal,
      const std::vector<double>& coords, const double& time) = 0;
  virtual void ComputeBdryFlux(const int num_points, std::vector<double>& flux,
                               FACE_INPUTS, const std::vector<double>& coords,
                               const double& time) = 0;
  virtual void ComputeBdrySolutionJacobi(const int num_points,
                                         double* bdry_u_jacobi,
                                         const std::vector<double>& owner_u,
                                         const std::vector<double>& owner_div_u,
                                         const std::vector<double>& normal,
                                         const std::vector<double>& coords,
                                         const double& time) = 0;
  virtual void ComputeBdryFluxJacobi(const int num_points,
                                     std::vector<double>& flux_jacobi,
                                     std::vector<double>& flux_grad_jacobi,
                                     FACE_INPUTS,
                                     const std::vector<double>& coords,
                                     const double& time) = 0;
};
// Boundary = AdiabaticWall
// Dependency: -
class AdiabaticWallNS2D : public BoundaryNS2D {
 public:
  AdiabaticWallNS2D(const int bdry_tag, EquationNS2D* equation);
  virtual ~AdiabaticWallNS2D() {}

  BOUNDARY_METHODS;
};
// Boundary = Riemann
// Dependency: Ma, AOA
class RiemannNS2D : public BoundaryNS2D {
 private:
  double Ma_;
  double AoA_;  // degree
  double uf_;
  double vf_;

 public:
  RiemannNS2D(const int bdry_tag, EquationNS2D* equation);
  virtual ~RiemannNS2D() {}

  BOUNDARY_METHODS;
};
// Boundary = Constant
// Dependency: BdryInput(num)
// BdryInput(num) = rho, rhoU, rhoV, rhoE
class ConstantBdryNS2D : public BoundaryNS2D {
 private:
  std::vector<double> values_;

 public:
  ConstantBdryNS2D(const int bdry_tag, EquationNS2D* equation);
  virtual ~ConstantBdryNS2D() {}

  BOUNDARY_METHODS;
};

// -------------------------------- Problem -------------------------------- //
class ProblemNS2D : public ConstantsNS2D {
 public:
  static std::shared_ptr<ProblemNS2D> GetProblem(const std::string& name);

  ProblemNS2D() : ConstantsNS2D(){};
  virtual ~ProblemNS2D(){};

  virtual void Problem(const int num_points, std::vector<double>& solutions,
                       const std::vector<double>& coord,
                       const double time = 0.0) const = 0;
};
// Problem = FreeStream
// ProblemInput = Ma, AOA
class FreeStreamNS2D : public ProblemNS2D {
 private:
  double Ma_;
  double AoA_;  // degree
  double uf_;
  double vf_;

 public:
  FreeStreamNS2D();
  virtual ~FreeStreamNS2D(){};

  virtual void Problem(const int num_points, std::vector<double>& solutions,
                       const std::vector<double>& coord,
                       const double time = 0.0) const;
};
// Problem = ShockTube
// ProblemInput = x-split,
//                left_rho, left_rhoU, left_rhoE,
//                right_rho, right_rhoU, right_rhoE
class ShockTubeNS2D : public ProblemNS2D {
 private:
  double split_;
  std::vector<double> left_values_;
  std::vector<double> right_values_;

 public:
  ShockTubeNS2D();
  virtual ~ShockTubeNS2D(){};

  virtual void Problem(const int num_points, std::vector<double>& solutions,
                       const std::vector<double>& coord,
                       const double time = 0.0) const;
};
// Problem = DoubleSine
// ProblemInput = -
class DoubleSineNS2D : public ProblemNS2D {
 private:
  std::vector<double> velocity_;
  std::vector<double> wave_number_;

 public:
  DoubleSineNS2D();
  virtual ~DoubleSineNS2D(){};

  virtual void Problem(const int num_points, std::vector<double>& solutions,
                       const std::vector<double>& coord,
                       const double time = 0.0) const;
};
}  // namespace deneb