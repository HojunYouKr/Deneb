#pragma once

#include <cmath>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "avocado.h"
#include "deneb_equation.h"

// 2-D compressible Navier-Stokes equations at thermochemical equilibrium state

namespace deneb {
class ProblemEquilibriumNS2D;
class BoundaryEquilibriumNS2D;

// ------------------------------- Constants ------------------------------- //
class ConstantsEquilibriumNS2D {
 protected:
  static constexpr int D_ = 2;
  static constexpr int S_ = 4;
  static constexpr int DS_ = D_ * S_;
  static constexpr int SS_ = S_ * S_;
  static constexpr int DSS_ = DS_ * S_;
  static constexpr int DDSS_ = D_ * DSS_;

  const double L_ref_;
  const double rho_ref_;
  const double T_ref_;
  const double a_ref_;
  const double V_ref_;
  const double mu_ref_;
  const double k_ref_;

  const double e_ref_;
  const double p_ref_;
  const double Re_;
  const double Ma_;
  // alpha = Ma / Re
  double alpha_;
  // beta = k*T / (L*rho*a**3)
  double beta_;

  const double c23_ = 2.0 / 3.0;
  const double c43_ = 4.0 / 3.0;
  const double pi_ = 4.0 * std::atan(1.0);

  int max_num_points_;
  int max_num_cell_points_;
  int max_num_face_points_;
  int max_num_bdry_points_;

 public:
  ConstantsEquilibriumNS2D();
  virtual ~ConstantsEquilibriumNS2D(){};

  inline int GetMaxNumBdryPoints() const { return max_num_bdry_points_; };

 protected:
  inline double IDEA_Wrapper(double (*IDEA)(double, double),
                             double (*IDEA_Grad)(double*, double, double),
                             double (*IDEA_Hess)(double*, double*, double,
                                                 double),
                             const double x1, const double x2) const {
    return IDEA(x1, x2);
  }
  inline void IDEA_Grad_Wrapper(double (*IDEA)(double, double),
                                double (*IDEA_Grad)(double*, double, double),
                                double (*IDEA_Hess)(double*, double*, double,
                                                    double),
                                double* grad, const double x1,
                                const double x2) const {
    IDEA_Grad(grad, x1, x2);
  }
  template <int dim>
  Dual<dim> IDEA_Wrapper(double (*IDEA)(double, double),
                         double (*IDEA_Grad)(double*, double, double),
                         double (*IDEA_Hess)(double*, double*, double, double),
                         const Dual<dim>& x1, const Dual<dim>& x2) const {
    double grad[2];
    Dual<dim> var(IDEA_Grad(grad, x1.f, x2.f));
    for (int i = 0; i < dim; i++)
      var.df[i] = grad[0] * x1.df[i] + grad[1] * x2.df[i];
    return var;
  }
  template <int dim>
  void IDEA_Grad_Wrapper(double (*IDEA)(double, double),
                         double (*IDEA_Grad)(double*, double, double),
                         double (*IDEA_Hess)(double*, double*, double, double),
                         Dual<dim>* grad, const Dual<dim>& x1,
                         const Dual<dim>& x2) const {
    double grad_dble[2], hess[3];
    IDEA_Hess(hess, grad_dble, x1.f, x2.f);

    grad[0].f = grad_dble[0];
    grad[1].f = grad_dble[1];
    for (int i = 0; i < dim; i++) {
      grad[0].df[i] = hess[0] * x1.df[i] + hess[1] * x2.df[i];
      grad[1].df[i] = hess[1] * x1.df[i] + hess[2] * x2.df[i];
    }
  }
  double ComputeInternalEnergy(const double* sol) const;
  double ComputePressure(const double* sol) const;
};

// ------------------------------- Equation -------------------------------- //
class EquationEquilibriumNS2D : public Equation,
                                public ConstantsEquilibriumNS2D {
 private:
  std::shared_ptr<ProblemEquilibriumNS2D> problem_;
  void (EquationEquilibriumNS2D::*compute_numflux_)(const int num_points,
                                                    std::vector<double>& flux,
                                                    FACE_INPUTS);
  void (EquationEquilibriumNS2D::*compute_numflux_jacobi_)(const int num_points,
                                                           FACE_JACOBI_OUTPUTS,
                                                           FACE_INPUTS);

  std::unordered_map<int, std::shared_ptr<BoundaryEquilibriumNS2D>>
      boundary_registry_;
  std::vector<std::shared_ptr<BoundaryEquilibriumNS2D>> boundaries_;

 public:
  EquationEquilibriumNS2D();
  virtual ~EquationEquilibriumNS2D();

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
};

// ------------------------------- Boundary -------------------------------- //
class BoundaryEquilibriumNS2D : public ConstantsEquilibriumNS2D {
 public:
  static std::shared_ptr<BoundaryEquilibriumNS2D> GetBoundary(
      const std::string& type, const int bdry_tag,
      EquationEquilibriumNS2D* equation);

 protected:
  int bdry_tag_;
  EquationEquilibriumNS2D* equation_;

 public:
  BoundaryEquilibriumNS2D(const int bdry_tag, EquationEquilibriumNS2D* equation)
      : bdry_tag_(bdry_tag), equation_(equation){};
  virtual ~BoundaryEquilibriumNS2D(){};

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
class AdiabaticWallEquilibriumNS2D : public BoundaryEquilibriumNS2D {
 public:
  AdiabaticWallEquilibriumNS2D(const int bdry_tag,
                               EquationEquilibriumNS2D* equation);
  virtual ~AdiabaticWallEquilibriumNS2D() {}

  BOUNDARY_METHODS;
};
// Boundary = IsothermalWall
// Dependency: BdryInput(num)
// BdryInput(num) = Twall
class IsothermalWallEquilibriumNS2D : public BoundaryEquilibriumNS2D {
 private:
  double Twall_;  // K

 public:
  IsothermalWallEquilibriumNS2D(const int bdry_tag,
                                EquationEquilibriumNS2D* equation);
  virtual ~IsothermalWallEquilibriumNS2D() {}

  BOUNDARY_METHODS;
};
// Boundary = SupersonicInflow
// Dependency: BdryInput(num)
// BdryInput(num) = Ma, AOA, Pinflow, Tinflow
class SupersonicInflowEquilibriumNS2D : public BoundaryEquilibriumNS2D {
 private:
  double Ma_;
  double AoA_;  // degree
  double p_inflow_; // Pa
  double T_inflow_; // K

  double d_;
  double u_;
  double v_;
  double p_;
  double du_;
  double dv_;
  double dE_;

 public:
  SupersonicInflowEquilibriumNS2D(const int bdry_tag,
                                  EquationEquilibriumNS2D* equation);
  virtual ~SupersonicInflowEquilibriumNS2D() {}

  BOUNDARY_METHODS;
};
// Boundary = BackPressure
// for subsonic and supersonic
// Dependency: BdryInput(num)
// BdryInput(num) = back pressure
class BackPressureEquilibriumNS2D : public BoundaryEquilibriumNS2D {
 private:
  double p_back_; // Pa

 public:
  BackPressureEquilibriumNS2D(const int bdry_tag,
                              EquationEquilibriumNS2D* equation);
  virtual ~BackPressureEquilibriumNS2D() {}

  BOUNDARY_METHODS;
};

// -------------------------------- Problem -------------------------------- //
class ProblemEquilibriumNS2D : public ConstantsEquilibriumNS2D {
 public:
  static std::shared_ptr<ProblemEquilibriumNS2D> GetProblem(
      const std::string& name);

  ProblemEquilibriumNS2D() : ConstantsEquilibriumNS2D(){};
  virtual ~ProblemEquilibriumNS2D(){};

  virtual void Problem(const int num_points, std::vector<double>& solutions,
                       const std::vector<double>& coord,
                       const double time = 0.0) const = 0;
};
// Problem = Constant
// ProblemInput = value 1, ... , value 4
class ConstantEquilibriumNS2D : public ProblemEquilibriumNS2D {
 private:
  std::vector<double> values_;

 public:
  ConstantEquilibriumNS2D();
  virtual ~ConstantEquilibriumNS2D(){};

  virtual void Problem(const int num_points, std::vector<double>& solutions,
                       const std::vector<double>& coord,
                       const double time = 0.0) const;
};
// Problem = FreeStream
// ProblemInput = Ma, AOA, Pinf, Tinf
class FreeStreamEquilibriumNS2D : public ProblemEquilibriumNS2D {
 private:
  double Ma_;
  double AoA_;       // degree
  double p_inf_;  // Pa
  double T_inf_;  // K

  double d_;
  double du_;
  double dv_;
  double dE_;

 public:
  FreeStreamEquilibriumNS2D();
  virtual ~FreeStreamEquilibriumNS2D(){};

  virtual void Problem(const int num_points, std::vector<double>& solutions,
                       const std::vector<double>& coord,
                       const double time = 0.0) const;
};
// Problem = DoubleSine
// ProblemInput = -
class DoubleSineEquilibriumNS2D : public ProblemEquilibriumNS2D {
 private:
  std::vector<double> velocity_;
  std::vector<double> wave_number_;

 public:
  DoubleSineEquilibriumNS2D();
  virtual ~DoubleSineEquilibriumNS2D(){};

  virtual void Problem(const int num_points, std::vector<double>& solutions,
                       const std::vector<double>& coord,
                       const double time = 0.0) const;
};
}  // namespace deneb