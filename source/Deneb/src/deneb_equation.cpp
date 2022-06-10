#include "deneb_equation.h"

#include "avocado.h"
#include "deneb_data.h"
#include "deneb_equation_equilibriumns2d.h"
#include "deneb_equation_euler2d.h"
#include "deneb_equation_glmmhd2d.h"
#include "deneb_equation_ns2d.h"
#include "deneb_equation_ns3d.h"

namespace deneb {
std::shared_ptr<Equation> DENEB_EQUATION_NAME = nullptr;
std::shared_ptr<Equation> Equation::GetEquation(const std::string& name) {
  if (!name.compare("NS2D"))
    return std::make_shared<EquationNS2D>();
  else if (!name.compare("EquilibriumNS2D"))
    return std::make_shared<EquationEquilibriumNS2D>();
  else if (!name.compare("Euler2D"))
    return std::make_shared<EquationEuler2D>();
  else if (!name.compare("GLMMHD2D"))
    return std::make_shared<EquationGLMMHD2D>();
  else if (!name.compare("NS3D"))
    return std::make_shared<EquationNS3D>();
  ERROR_MESSAGE("Wrong equation (no-exist):" + name + "\n");
  return nullptr;
}

Equation::Equation(const int dimension, const int num_states,
                   const bool source_term)
    : dimension_(dimension),
      num_states_(num_states),
      source_term_(source_term) {}
}  // namespace deneb