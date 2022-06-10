#pragma once

#include <memory>
#include <string>
#include <vector>

#include <petsc.h>

#define DENEB_EQUATION_NAME equation_global_ptr
#define DENEB_EQUATION deneb::DENEB_EQUATION_NAME
#define DENEB_EQUATION_INITIALIZE(name) \
  DENEB_EQUATION = deneb::Equation::GetEquation(name)
#define DENEB_EQUATION_FINALIZE() DENEB_EQUATION.reset()

#define FACE_INPUTS                              \
  const int owner_cell, const int neighbor_cell, \
      const std::vector<double>&owner_u,         \
      const std::vector<double>&owner_div_u,     \
      const std::vector<double>&neighbor_u,      \
      const std::vector<double>&neighbor_div_u,  \
      const std::vector<double>&normal
#define FACE_JACOBI_OUTPUTS                       \
  std::vector<double>&flux_owner_jacobi,          \
      std::vector<double>&flux_neighbor_jacobi,   \
      std::vector<double>&flux_owner_grad_jacobi, \
      std::vector<double>&flux_neighbor_grad_jacobi
#define ASSIGN_FLUX(equation, fluxtype)                        \
  compute_numflux_ = &##equation## ::ComputeNumFlux##fluxtype; \
  compute_numflux_jacobi_ = &##equation## ::ComputeNumFluxJacobi##fluxtype
#define DEFINE_FLUX(fluxtype)                                        \
  virtual void ComputeNumFlux##fluxtype##(                           \
      const int num_points, std::vector<double>& flux, FACE_INPUTS); \
  virtual void ComputeNumFluxJacobi##fluxtype##(                     \
      const int num_points, FACE_JACOBI_OUTPUTS, FACE_INPUTS)
#define BOUNDARY_METHODS                                                    \
  virtual void ComputeBdrySolution(                                         \
      const int num_points, std::vector<double>& bdry_u,                    \
      std::vector<double>& bdry_div_u, const std::vector<double>& owner_u,  \
      const std::vector<double>& owner_div_u,                               \
      const std::vector<double>& normal, const std::vector<double>& coords, \
      const double& time);                                                  \
  virtual void ComputeBdryFlux(                                             \
      const int num_points, std::vector<double>& flux, FACE_INPUTS,         \
      const std::vector<double>& coords, const double& time);               \
  virtual void ComputeBdrySolutionJacobi(                                   \
      const int num_points, double* bdry_u_jacobi,                          \
      const std::vector<double>& owner_u,                                   \
      const std::vector<double>& owner_div_u,                               \
      const std::vector<double>& normal, const std::vector<double>& coords, \
      const double& time);                                                  \
  virtual void ComputeBdryFluxJacobi(                                       \
      const int num_points, std::vector<double>& flux_jacobi,               \
      std::vector<double>& flux_grad_jacobi, FACE_INPUTS,                   \
      const std::vector<double>& coords, const double& time)

namespace avocado {
class Communicate;
}
namespace deneb {
class Equation {
 public:
  static std::shared_ptr<Equation> GetEquation(const std::string& name);

 protected:
  int dimension_;
  int num_states_;
  bool source_term_;

  std::vector<double> dt_auxiliary_;
  std::vector<double> auxiliary_solution_;
  std::vector<double> pressure_fix_values_;
  std::shared_ptr<avocado::Communicate> communicate_;
  std::vector<double> outer_solution_;

  std::vector<std::string> cell_variable_names_;
  std::vector<std::string> face_variable_names_;

 public:
  Equation(const int dimension, const int num_states, const bool source_term);
  virtual ~Equation(){};

  inline const int& GetDimension(void) const { return dimension_; };
  inline const int& GetNumStates(void) const { return num_states_; };
  inline bool GetSourceTerm(void) const { return source_term_; };
  inline const std::vector<std::string>& GetCellVariableNames(void) const {
    return cell_variable_names_;
  };
  inline const std::vector<std::string>& GetFaceVariableNames(void) const {
    return face_variable_names_;
  };

  virtual void RegistBoundary(const std::vector<int>& bdry_tag) = 0;
  virtual void BuildData(void) = 0;
  virtual void PreProcess(const double* solution) = 0;
  virtual void ComputeRHS(const double* solution, double* rhs,
                          const double t) = 0;
  virtual bool ComputeRHSdt(const double* solution, double* rhs_dt,
                            const double t) = 0;
  virtual void ComputeSystemMatrix(const double* solution, Mat& sysmat,
                                   const double t) = 0;
  virtual void GetCellPostSolution(const int num_points,
                                   const std::vector<double>& solution,
                                   const std::vector<double>& solution_grad,
                                   std::vector<double>& post_solution) = 0;
  virtual void GetFacePostSolution(const int num_points,
                                   const std::vector<double>& solution,
                                   const std::vector<double>& solution_grad,
                                   const std::vector<double>& normal,
                                   std::vector<double>& post_solution) = 0;
  virtual void ComputeInitialSolution(double* solution, const double t) = 0;
  virtual void ComputeLocalTimestep(const double* solution,
                                    std::vector<double>& local_timestep) = 0;
  virtual void ComputePressureCoefficient(
      std::vector<double>& pressure_coefficient, const int num_points,
      const std::vector<double>& solution) = 0;
  virtual void ComputeViscousStress(
      std::vector<double>& viscous_stress, const int num_points,
      const std::vector<double>& solution,
      const std::vector<double>& solution_grad) = 0;
  virtual bool IsContact(const int& icell,
                         const std::vector<int>& neighbor_cells,
                         const double* solution,
                         const double* total_solution) const = 0;
  virtual double ComputeMaxCharacteristicSpeed(
      const double* input_solution) const = 0;
  virtual const std::vector<double>& ComputePressureFixValues(
      const double* input_solution) = 0;
};
extern std::shared_ptr<Equation> DENEB_EQUATION_NAME;
}  // namespace deneb