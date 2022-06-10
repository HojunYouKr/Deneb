#include "deneb_artificial_viscosity.h"

#include <algorithm>
#include <cstring>

#include "avocado.h"
#include "deneb_config_macro.h"
#include "deneb_data.h"
#include "deneb_equation.h"
#include "deneb_quadrature.h"

namespace deneb {
std::shared_ptr<ArtificialViscosity> DENEB_ARTIFICIAL_VISCOSITY_NAME = nullptr;
std::shared_ptr<ArtificialViscosity>
ArtificialViscosity::GetArtificialViscosity(const std::string& name) {
  if (!name.compare("None"))
    return std::make_shared<NoArtificialViscosity>();
  else if (!name.compare("LaplacianP0"))
    return std::make_shared<LaplacianP0>();
  ERROR_MESSAGE("Wrong artificial viscosity (no-exist):" + name + "\n");
  return nullptr;
}

NoArtificialViscosity::NoArtificialViscosity() {
  MASTER_MESSAGE(avocado::GetTitle("NoArtificialViscosity"));
}
void NoArtificialViscosity::BuildData(void) {
  MASTER_MESSAGE(avocado::GetTitle("NoArtificialViscosity::BuildData()"));
}

// -------------------- LaplacianP0 ------------------ //
LaplacianP0::LaplacianP0() { MASTER_MESSAGE(avocado::GetTitle("LaplacianP0")); }
void LaplacianP0::BuildData(void) {
  MASTER_MESSAGE(avocado::GetTitle("LaplacianP0::BuildData()"));

  target_state_ = 0;
  const std::vector<std::string>& variable_names =
      DENEB_EQUATION->GetCellVariableNames();
  MASTER_MESSAGE("Target state: " + variable_names[target_state_] + "\n");

  auto& config = AVOCADO_CONFIG;
  const std::string& equation = config->GetConfigValue(EQUATION);
  if (!equation.compare("Euler2D"))
    ERROR_MESSAGE("Artificial viscosity is not compatible with " + equation +
                  ", use NS equation with Re<0\n");
  else if (!equation.compare("ScalarAdvection2D"))
    ERROR_MESSAGE("Artificial viscosity is not compatible with " + equation +
                  ", use NS equation with Re<0\n");

  Peclet_ = std::stod(config->GetConfigValue(PECLET));
  kappa_ = std::stod(config->GetConfigValue(KAPPA));

  const int& num_outer_cells = DENEB_DATA->GetNumOuterCells();
  artificial_viscosity_.resize(num_outer_cells + 1, 0.0);

  const int& order = DENEB_DATA->GetOrder();
  if (order == 0) ERROR_MESSAGE("P0 doesn't require artificial viscosity");
  S0_ = -3.0 * std::log10(order);

  const int& dimension = DENEB_EQUATION->GetDimension();
  if (dimension == 2)
    num_bases_m1_ = order * (order + 1) / 2;
  else if (dimension == 3)
    num_bases_m1_ = order * (order + 1) * (order + 2) / 6;
  else
    ERROR_MESSAGE("Dimension error\n");

  std::vector<double> GL_points;
  std::vector<double> GL_weights;
  quadrature::Line_Poly1D(order, GL_points, GL_weights);
  for (int ipoint = 0; ipoint < GL_points.size(); ipoint++)
    GL_points[ipoint] = 0.5 * (GL_points[ipoint] + 1.0);
  dLmax_ = 0.0;
  for (int ipoint = 0; ipoint < GL_points.size() - 1; ipoint++) {
    const double dist = GL_points[ipoint] - GL_points[ipoint + 1];
    if (dLmax_ < dist) dLmax_ = dist;
  }

  const int& num_states = DENEB_EQUATION->GetNumStates();
  const int& num_bases = DENEB_DATA->GetNumBases();
  const int sb = num_states * num_bases;
  const int& num_cells = DENEB_DATA->GetNumCells();
  communicate_ = std::make_shared<avocado::Communicate>(
      1, DENEB_DATA->GetOuterSendCellList(),
      DENEB_DATA->GetOuterRecvCellList());
}
void LaplacianP0::ComputeArtificialViscosity(const double* solution) {
  static const int& order = DENEB_DATA->GetOrder();
  if (order == 0) return;
  static const int& num_cells = DENEB_DATA->GetNumCells();
  static const int& num_outer_cells = DENEB_DATA->GetNumOuterCells();
  static const int& num_states = DENEB_EQUATION->GetNumStates();
  static const int& num_bases = DENEB_DATA->GetNumBases();
  static const int sb = num_states * num_bases;
  static const std::vector<double>& cell_volumes = DENEB_DATA->GetCellVolumes();
  static const std::vector<std::vector<double>>& cell_basis_value =
      DENEB_DATA->GetCellBasisValue();

  for (int icell = 0; icell < num_cells; icell++) {
    const double E0 = MaxArtificialViscosity(
        &solution[icell * sb], cell_volumes[icell], cell_basis_value[icell][0]);
    const double Se = SmoothnessIndicator(&solution[icell * sb]);

    if (Se <= S0_ - kappa_)
      artificial_viscosity_[icell] = 0.0;
    else if (Se <= S0_ + kappa_)
      artificial_viscosity_[icell] =
          0.5 * E0 * (1.0 + std::sin(M_PI * (Se - S0_) / (2.0 * kappa_)));
    else
      artificial_viscosity_[icell] = E0;
  }
  communicate_->CommunicateBegin(&artificial_viscosity_[0]);
  communicate_->CommunicateEnd(&artificial_viscosity_[num_cells]);
}
double LaplacianP0::SmoothnessIndicator(const double* solution) {
  static const int& num_bases = DENEB_DATA->GetNumBases();
  const double Pn_value =
      avocado::VecInnerProd(num_bases, &solution[target_state_ * num_bases],
                            &solution[target_state_ * num_bases]);
  const double Pn_minus_1_value =
      avocado::VecInnerProd(num_bases_m1_, &solution[target_state_ * num_bases],
                            &solution[target_state_ * num_bases]);

  return std::log10((Pn_value - Pn_minus_1_value) / Pn_value);
}

double LaplacianP0::MaxArtificialViscosity(const double* solution,
                                           const double cell_volumes,
                                           const double cell_basis_value) {
  static const int& dimension = DENEB_EQUATION->GetDimension();
  static const int& num_states = DENEB_EQUATION->GetNumStates();
  static const int& num_bases = DENEB_DATA->GetNumBases();
  static const int& num_cells = DENEB_DATA->GetNumCells();
  static std::vector<double> input_solutions(num_states);

  for (int istate = 0; istate < num_states; istate++)
    input_solutions[istate] = solution[istate * num_bases] * cell_basis_value;

  const double max_speed =
      DENEB_EQUATION->ComputeMaxCharacteristicSpeed(&input_solutions[0]);
  const double length_scale = std::pow(cell_volumes, 1.0 / dimension);
  return max_speed * length_scale * (2.0 - dLmax_) / Peclet_;
}
}  // namespace deneb