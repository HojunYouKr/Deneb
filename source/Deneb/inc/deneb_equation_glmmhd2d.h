#pragma once

#include <cmath>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "deneb_equation.h"

// 2-D ideal magnetohydrodynamic equations with a non-conservative
// magnetic-divergence (eight-wave) term and hyperbolic ansatz of GLM approach

namespace deneb {
class ProblemGLMMHD2D;
class BoundaryGLMMHD2D;

// ------------------------------- Constants ------------------------------- //
class ConstantsGLMMHD2D {
 protected:
  static constexpr int D_ = 2;
  static constexpr int S_ = 7;
  static constexpr int DD_ = D_ * D_;
  static constexpr int DS_ = D_ * S_;
  static constexpr int SS_ = S_ * S_;
  static constexpr int DSS_ = D_ * SS_;
  static constexpr int DDSS_ = D_ * DSS_;

  const double r_;
  const double r_inv_;
  const double rm1_;
  const double rm1_inv_;

  const double pi_ = 4.0 * std::atan(1.0);

  int max_num_points_;
  int max_num_cell_points_;
  int max_num_face_points_;
  int max_num_bdry_points_;

 public:
  ConstantsGLMMHD2D();
  virtual ~ConstantsGLMMHD2D(){};

  inline int GetMaxNumBdryPoints() const { return max_num_bdry_points_; };

 protected:
  // compute p
  double ComputePressure(const double* sol) const;
  // compute dE
  double ComputeTotalEnergy(const double* sol) const;
};

// ------------------------------- Equation -------------------------------- //
class EquationGLMMHD2D : public Equation, public ConstantsGLMMHD2D {
 private:
  std::shared_ptr<ProblemGLMMHD2D> problem_;
  void (EquationGLMMHD2D::*compute_numflux_)(const int num_points,
                                             std::vector<double>& flux,
                                             FACE_INPUTS);
  void (EquationGLMMHD2D::*compute_numflux_jacobi_)(const int num_points,
                                                    FACE_JACOBI_OUTPUTS,
                                                    FACE_INPUTS);

  std::unordered_map<int, std::shared_ptr<BoundaryGLMMHD2D>> boundary_registry_;
  std::vector<std::shared_ptr<BoundaryGLMMHD2D>> boundaries_;

  double ch_;

 public:
  EquationGLMMHD2D();
  virtual ~EquationGLMMHD2D();

  virtual void RegistBoundary(const std::vector<int>& bdry_tag);
  virtual void BuildData(void);
  virtual void PreProcess(const double* solution);
  virtual void ComputeRHS(const double* solution, double* rhs, const double t);
  virtual bool ComputeRHSdt(const double* solution, double* rhs_dt,
                            const double t) {
    return false;
  };
  virtual void ComputeSystemMatrix(const double* solution, Mat& sysmat,
                                   const double t);
  virtual void GetCellPostSolution(const int num_points,
                                   const std::vector<double>& solution,
                                   const std::vector<double>& solution_grad,
                                   std::vector<double>& post_solution);
  virtual void GetFacePostSolution(const int num_points,
                                   const std::vector<double>& solution,
                                   const std::vector<double>& solution_grad,
                                   const std::vector<double>& normal,
                                   std::vector<double>& post_solution);
  virtual void ComputeInitialSolution(double* solution, const double t);
  virtual void ComputeLocalTimestep(const double* solution,
                                    std::vector<double>& local_timestep);
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

   void ComputeComFlux(const int num_points,
                      std::vector<double>& flux, 
                      const int icell, const std::vector<double>& owner_u,
                      const std::vector<double>& owner_div_u);
  void ComputeComFluxJacobi(const int num_points,
                            std::vector<double>& flux_jacobi,
                            std::vector<double>& flux_grad_jacobi,
                            const int icell, const std::vector<double>& owner_u,
                            const std::vector<double>& owner_div_u);
  // source : ps
  void ComputeSource(const int num_points, std::vector<double>& source,
                     const int icell, const std::vector<double>& owner_u,
                     const std::vector<double>& owner_div_u);
  inline double GetCh(void) const { return ch_; };
  inline void ComputeNumFlux(const int num_points, std::vector<double>& flux,
                             FACE_INPUTS) {
    (this->*compute_numflux_)(num_points, flux, owner_cell, neighbor_cell,
                              owner_u, owner_div_u,
                              neighbor_u, neighbor_div_u, normal);
  };
  inline void ComputeNumFluxJacobi(const int num_points, FACE_JACOBI_OUTPUTS,
                                   FACE_INPUTS) {
    (this->*compute_numflux_jacobi_)(
        num_points, flux_owner_jacobi, flux_neighbor_jacobi,
        flux_owner_grad_jacobi, flux_neighbor_grad_jacobi, owner_cell,
        neighbor_cell, 
        owner_u, owner_div_u, neighbor_u, neighbor_div_u, normal);
  };

  DEFINE_FLUX(LLF);
};

// ------------------------------- Boundary -------------------------------- //
class BoundaryGLMMHD2D : public ConstantsGLMMHD2D {
 public:
  static std::shared_ptr<BoundaryGLMMHD2D> GetBoundary(
      const std::string& type, const int bdry_tag, EquationGLMMHD2D* equation);

 protected:
  int bdry_tag_;
  EquationGLMMHD2D* equation_;

 public:
  BoundaryGLMMHD2D(const int bdry_tag, EquationGLMMHD2D* equation)
      : bdry_tag_(bdry_tag), equation_(equation){};
  virtual ~BoundaryGLMMHD2D(){};

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
// Boundary =  Extrapolate
// Dependency: -
class ExtrapolateGLMMHD2D : public BoundaryGLMMHD2D {
 public:
  ExtrapolateGLMMHD2D(const int bdry_tag, EquationGLMMHD2D* equation);
  virtual ~ExtrapolateGLMMHD2D() {}

  BOUNDARY_METHODS;
};
// Boundary = Constant
// Dependency: BdryInput(num)
// BdryInput(num) = rho, rhoU, rhoV, rhotE, B1, B2
class ConstantBdryGLMMHD2D : public BoundaryGLMMHD2D {
 private:
  std::vector<double> values_;

 public:
  ConstantBdryGLMMHD2D(const int bdry_tag, EquationGLMMHD2D* equation);
  virtual ~ConstantBdryGLMMHD2D() {}

  BOUNDARY_METHODS;
};

// -------------------------------- Problem -------------------------------- //
class ProblemGLMMHD2D : public ConstantsGLMMHD2D {
 public:
  static std::shared_ptr<ProblemGLMMHD2D> GetProblem(const std::string& name);

  ProblemGLMMHD2D() : ConstantsGLMMHD2D(){};
  virtual ~ProblemGLMMHD2D(){};

  virtual void Problem(const int num_points, std::vector<double>& solutions,
                       const std::vector<double>& coord,
                       const double time = 0.0) const = 0;
};
// Problem = DoubleSine
// ProblemInput = -
class DoubleSineGLMMHD2D : public ProblemGLMMHD2D {
 private:
  std::vector<double> velocity_;
  std::vector<double> wave_number_;

 public:
  DoubleSineGLMMHD2D();
  virtual ~DoubleSineGLMMHD2D(){};

  virtual void Problem(const int num_points, std::vector<double>& solutions,
                       const std::vector<double>& coord,
                       const double time = 0.0) const;
};
// Problem = ShockTube
// ProblemInput = x-split,
//                left_rho, left_rhoU, left_rhotE, left_B1, left_B2,
//                right_rho, right_rhoU, right_rhotE, right_B1, right_B2
class ShockTubeGLMMHD2D : public ProblemGLMMHD2D {
 private:
  double split_;
  std::vector<double> left_values_;
  std::vector<double> right_values_;

 public:
  ShockTubeGLMMHD2D();
  virtual ~ShockTubeGLMMHD2D(){};

  virtual void Problem(const int num_points, std::vector<double>& solutions,
                       const std::vector<double>& coord,
                       const double time = 0.0) const;
};
// Problem = OrszagTangVortex
// ProblemInput = -
class OrszagTangVortexGLMMHD2D : public ProblemGLMMHD2D {
 public:
  OrszagTangVortexGLMMHD2D();
  virtual ~OrszagTangVortexGLMMHD2D(){};

  virtual void Problem(const int num_points, std::vector<double>& solutions,
                       const std::vector<double>& coord,
                       const double time = 0.0) const;
};
// Problem = MHDRotor
// ProblemInput = -
class MHDRotorGLMMHD2D : public ProblemGLMMHD2D {
 public:
  MHDRotorGLMMHD2D();
  virtual ~MHDRotorGLMMHD2D(){};

  virtual void Problem(const int num_points, std::vector<double>& solutions,
                       const std::vector<double>& coord,
                       const double time = 0.0) const;
};
}  // namespace deneb