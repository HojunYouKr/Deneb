#include "deneb_equation_equilibriumns2d.h"

#include <algorithm>
#include <cstring>
#include <string>
#include <unordered_set>

#include "idea.h"

#include "avocado.h"
#include "deneb_artificial_viscosity.h"
#include "deneb_config_macro.h"
#include "deneb_data.h"
#include "deneb_timescheme.h"

#define IDEA_ROUTINES(input) \
  IDEA_##input, IDEA_##input##_Grad, IDEA_##input##_Hess

#define GET_NORMAL_PD(data)       \
  ind = D_ * ipoint;              \
  const double& nx = data[ind++]; \
  const double& ny = data[ind]

#define GET_SOLUTION_PS(tag, data)     \
  ind = S_ * ipoint;                   \
  const double& d##tag = data[ind++];  \
  const double& du##tag = data[ind++]; \
  const double& dv##tag = data[ind++]; \
  const double& dE##tag = data[ind]

#define GET_SOLUTION_GRAD_PDS(tag, data) \
  ind = DS_ * ipoint;                    \
  const double& dx##tag = data[ind++];   \
  const double& dux##tag = data[ind++];  \
  const double& dvx##tag = data[ind++];  \
  const double& dEx##tag = data[ind++];  \
  const double& dy##tag = data[ind++];   \
  const double& duy##tag = data[ind++];  \
  const double& dvy##tag = data[ind++];  \
  const double& dEy##tag = data[ind]

#define COMPUTE_FLUX_BASICS(tag, datatype)                                   \
  const auto d_inv##tag = 1.0 / d##tag;                                      \
  const auto u##tag = du##tag * d_inv##tag;                                  \
  const auto v##tag = dv##tag * d_inv##tag;                                  \
  const auto ux##tag = (dux##tag - dx##tag * u##tag) * d_inv##tag;           \
  const auto uy##tag = (duy##tag - dy##tag * u##tag) * d_inv##tag;           \
  const auto vx##tag = (dvx##tag - dx##tag * v##tag) * d_inv##tag;           \
  const auto vy##tag = (dvy##tag - dy##tag * v##tag) * d_inv##tag;           \
  const auto e##tag =                                                        \
      (dE##tag - 0.5 * (du##tag * u##tag + dv##tag * v##tag)) * d_inv##tag;  \
  const auto p##tag =                                                        \
      IDEA_Wrapper(IDEA_ROUTINES(DEP), d##tag * rho_ref_, e##tag * e_ref_) / \
      p_ref_;                                                                \
  const auto mu##tag =                                                       \
      IDEA_Wrapper(IDEA_ROUTINES(DEV), d##tag * rho_ref_, e##tag * e_ref_) / \
      mu_ref_;                                                               \
  const auto k##tag =                                                        \
      IDEA_Wrapper(IDEA_ROUTINES(DEK), d##tag * rho_ref_, e##tag * e_ref_) / \
      k_ref_;                                                                \
  const auto txx##tag = c23_ * mu##tag * (2.0 * ux##tag - vy##tag);          \
  const auto txy##tag = mu##tag * (uy##tag + vx##tag);                       \
  const auto tyy##tag = c23_ * mu##tag * (2.0 * vy##tag - ux##tag);          \
  const auto ex##tag =                                                       \
      (dEx##tag - dx##tag * dE##tag * d_inv##tag) * d_inv##tag -             \
      u##tag * ux##tag - v##tag * vx##tag;                                   \
  const auto ey##tag =                                                       \
      (dEy##tag - dy##tag * dE##tag * d_inv##tag) * d_inv##tag -             \
      u##tag * uy##tag - v##tag * vy##tag;                                   \
  datatype T_grad##tag[2];                                                   \
  IDEA_Grad_Wrapper(IDEA_ROUTINES(DET), T_grad##tag, d##tag* rho_ref_,       \
                    e##tag* e_ref_);                                         \
  T_grad##tag[0] = T_grad##tag[0] * (rho_ref_ / T_ref_);                     \
  T_grad##tag[1] = T_grad##tag[1] * (e_ref_ / T_ref_);                       \
  const auto Qx##tag =                                                       \
      -k##tag * (T_grad##tag[0] * dx##tag + T_grad##tag[1] * ex##tag);       \
  const auto Qy##tag =                                                       \
      -k##tag * (T_grad##tag[0] * dy##tag + T_grad##tag[1] * ey##tag)

#define COMPUTE_VOLUME_FLUX_PDS(flux, icell, ipoint)                          \
  ind = DS_ * ipoint;                                                         \
  AVcoeff =                                                                   \
      DENEB_ARTIFICIAL_VISCOSITY->GetArtificialViscosityValue(icell, ipoint); \
  flux[ind++] = du - AVcoeff * dx;                                            \
  flux[ind++] = du * u + p - alpha_ * txx - AVcoeff * dux;                    \
  flux[ind++] = du * v - alpha_ * txy - AVcoeff * dvx;                        \
  flux[ind++] = (dE + p) * u - alpha_ * (txx * u + txy * v) + beta_ * Qx -    \
                AVcoeff * dEx;                                                \
  flux[ind++] = dv - AVcoeff * dy;                                            \
  flux[ind++] = dv * u - alpha_ * txy - AVcoeff * duy;                        \
  flux[ind++] = dv * v + p - alpha_ * tyy - AVcoeff * dvy;                    \
  flux[ind] =                                                                 \
      (dE + p) * v - alpha_ * (txy * u + tyy * v) + beta_ * Qy - AVcoeff * dEy

#define COMPUTE_NUMERICAL_FLUX_PDS(flux, owner_cell, neighbor_cell, ipoint)    \
  ind = DS_ * ipoint;                                                          \
  AVcoeff_o = DENEB_ARTIFICIAL_VISCOSITY->GetArtificialViscosityValue(         \
      owner_cell, ipoint);                                                     \
  AVcoeff_n = DENEB_ARTIFICIAL_VISCOSITY->GetArtificialViscosityValue(         \
      neighbor_cell, ipoint);                                                  \
  flux[ind++] = 0.5 * (du_o + du_n - nx * diff1 -                              \
                       (AVcoeff_o * dx_o + AVcoeff_n * dx_n));                 \
  flux[ind++] = 0.5 * (du_o * u_o + p_o + du_n * u_n + p_n - nx * diff2 -      \
                       alpha_ * (txx_o + txx_n) -                              \
                       (AVcoeff_o * dux_o + AVcoeff_n * dux_n));               \
  flux[ind++] =                                                                \
      0.5 * (du_o * v_o + du_n * v_n - nx * diff3 - alpha_ * (txy_o + txy_n) - \
             (AVcoeff_o * dvx_o + AVcoeff_n * dvx_n));                         \
  flux[ind++] =                                                                \
      0.5 *                                                                    \
      ((dE_o + p_o) * u_o + (dE_n + p_n) * u_n - nx * diff4 -                  \
       alpha_ * (txx_o * u_o + txy_o * v_o + txx_n * u_n + txy_n * v_n) +      \
       beta_ * (Qx_o + Qx_n) - (AVcoeff_o * dEx_o + AVcoeff_n * dEx_n));       \
  flux[ind++] = 0.5 * (dv_o + dv_n - ny * diff1 -                              \
                       (AVcoeff_o * dy_o + AVcoeff_n * dy_n));                 \
  flux[ind++] =                                                                \
      0.5 * (dv_o * u_o + dv_n * u_n - ny * diff2 - alpha_ * (txy_o + txy_n) - \
             (AVcoeff_o * duy_o + AVcoeff_n * duy_n));                         \
  flux[ind++] = 0.5 * (dv_o * v_o + p_o + dv_n * v_n + p_n - ny * diff3 -      \
                       alpha_ * (tyy_o + tyy_n) -                              \
                       (AVcoeff_o * dvy_o + AVcoeff_n * dvy_n));               \
  flux[ind] =                                                                  \
      0.5 *                                                                    \
      ((dE_o + p_o) * v_o + (dE_n + p_n) * v_n - ny * diff4 -                  \
       alpha_ * (txy_o * u_o + tyy_o * v_o + txy_n * u_n + tyy_n * v_n) +      \
       beta_ * (Qy_o + Qy_n) - (AVcoeff_o * dEy_o + AVcoeff_n * dEy_n))

namespace deneb {
// ------------------------------- Constants -------------------------------
ConstantsEquilibriumNS2D::ConstantsEquilibriumNS2D()
    : L_ref_(std::stod(AVOCADO_CONFIG->GetConfigValue(REFERENCE_LENGTH))),
      rho_ref_(std::stod(AVOCADO_CONFIG->GetConfigValue(REFERENCE_DENSITY))),
      T_ref_(std::stod(AVOCADO_CONFIG->GetConfigValue(REFERENCE_TEMPERATURE))),
      a_ref_(std::stod(AVOCADO_CONFIG->GetConfigValue(REFERENCE_SOUNDSPEED))),
      V_ref_(std::stod(AVOCADO_CONFIG->GetConfigValue(REFERENCE_VELOCITY))),
      mu_ref_(std::stod(AVOCADO_CONFIG->GetConfigValue(REFERENCE_VISCOSITY))),
      k_ref_(std::stod(AVOCADO_CONFIG->GetConfigValue(REFERENCE_CONDUCTIVITY))),
      e_ref_(a_ref_ * a_ref_),
      p_ref_(rho_ref_ * a_ref_ * a_ref_),
      Re_(rho_ref_ * V_ref_ * L_ref_ / mu_ref_),
      Ma_(V_ref_ / a_ref_),
      alpha_(Ma_ / Re_),
      beta_(k_ref_ * T_ref_ / L_ref_ / rho_ref_ / a_ref_ / a_ref_ / a_ref_) {}
double ConstantsEquilibriumNS2D::ComputeInternalEnergy(
    const double* sol) const {
  return (sol[S_ - 1] -
          0.5 / sol[0] * avocado::VecInnerProd(D_, &sol[1], &sol[1])) /
         sol[0];
}
double ConstantsEquilibriumNS2D::ComputePressure(const double* sol) const {
  const double& D = sol[0];
  const double E = ComputeInternalEnergy(sol);
  return IDEA_DEP(D * rho_ref_, E * a_ref_) / p_ref_;
}

// ------------------------------- Equation --------------------------------
EquationEquilibriumNS2D::EquationEquilibriumNS2D()
    : ConstantsEquilibriumNS2D(), Equation(D_, S_, false) {
  MASTER_MESSAGE(avocado::GetTitle("EquationEquilibriumNS2D"));
  MASTER_MESSAGE("Dimension = " + std::to_string(D_) + "\n");
  MASTER_MESSAGE("Number of state variables = " + std::to_string(S_) + "\n");
  MASTER_MESSAGE(
      "Source term = " + std::string(source_term_ ? "true" : "false") + "\n");

  MASTER_MESSAGE("References\n");
  MASTER_MESSAGE("\tLength = " + std::to_string(L_ref_) + " m\n");
  MASTER_MESSAGE("\tDensity = " + std::to_string(rho_ref_) + " kg/m3\n");
  MASTER_MESSAGE("\tTemperature = " + std::to_string(T_ref_) + " K\n");
  MASTER_MESSAGE("\tSpeed of sound = " + std::to_string(a_ref_) + " m/s\n");
  MASTER_MESSAGE("\tVelocity = " + std::to_string(V_ref_) + " m/s\n");
  MASTER_MESSAGE("\tViscosity = " + std::to_string(mu_ref_) + " Pa*s\n");
  MASTER_MESSAGE("\tConductivity = " + std::to_string(k_ref_) + " W/m/K\n");

  MASTER_MESSAGE("Dimensionless Parameters\n");
  MASTER_MESSAGE("\tRe = " + std::to_string(Re_) + "\n");
  MASTER_MESSAGE("\tMa = " + std::to_string(Ma_) + "\n");
  MASTER_MESSAGE("\talpha = " + std::to_string(alpha_) + "\n");
  MASTER_MESSAGE("\tbeta = " + std::to_string(beta_) + "\n");

  const auto& config = AVOCADO_CONFIG;
  const std::string dir = config->GetConfigValue(IDEA_RESOURCE_DIRECTORY);
  const std::string errlev = "ALL=1, DEG=0.5";
  const char* initlog = IDEA_Init(dir.c_str(), errlev.c_str());
  MASTER_MESSAGE(std::string(initlog) + "\n");

  problem_ =
      ProblemEquilibriumNS2D::GetProblem(config->GetConfigValue(PROBLEM));
  const std::string& numflux = config->GetConfigValue(CONVECTIVE_FLUX);
  if (!numflux.compare("LLF")) {
    ASSIGN_FLUX(EquationEquilibriumNS2D, LLF);
  } else
    ERROR_MESSAGE("Wrong numerical flux (no-exist):" + numflux + "\n");
  MASTER_MESSAGE("Problem: " + config->GetConfigValue(PROBLEM) + "\n");
  MASTER_MESSAGE("Convective flux: " + numflux + "\n");
}
EquationEquilibriumNS2D::~EquationEquilibriumNS2D() {
  problem_.reset();
  boundaries_.clear();
  boundary_registry_.clear();
  const char* finlog = IDEA_Finalize();
  MASTER_MESSAGE(std::string(finlog) + "\n");
}

void EquationEquilibriumNS2D::RegistBoundary(const std::vector<int>& bdry_tag) {
  auto& config = AVOCADO_CONFIG;
  std::vector<int> all_bdry_tag;
  {
    std::unordered_set<int> temp(bdry_tag.begin(), bdry_tag.end());
    std::vector<int> bdry_tag_new(temp.begin(), temp.end());
    temp.clear();

    std::vector<std::vector<int>> send_data(NDOMAIN, bdry_tag_new);
    bdry_tag_new.clear();
    std::vector<std::vector<int>> recv_data;
    AVOCADO_MPI->CommunicateData(send_data, recv_data);
    send_data.clear();

    for (auto&& data : recv_data) temp.insert(data.begin(), data.end());
    recv_data.clear();
    all_bdry_tag = std::vector<int>(temp.begin(), temp.end());
  }

  for (auto&& tag : all_bdry_tag) {
    const std::string& bdry_type = config->GetConfigValue(BDRY_TYPE(tag));
    if (boundary_registry_.find(tag) == boundary_registry_.end())
      boundary_registry_[tag] =
          BoundaryEquilibriumNS2D::GetBoundary(bdry_type, tag, this);
  }
}
void EquationEquilibriumNS2D::BuildData(void) {
  MASTER_MESSAGE(avocado::GetTitle("EquationEquilibriumNS2D::BuildData()"));

  const int& num_bdries = DENEB_DATA->GetNumBdries();
  boundaries_.resize(num_bdries, nullptr);
  {
    const std::vector<int>& bdry_tag = DENEB_DATA->GetBdryTag();
    RegistBoundary(bdry_tag);
    for (int ibdry = 0; ibdry < num_bdries; ibdry++)
      boundaries_[ibdry] = boundary_registry_[bdry_tag[ibdry]];
  }

  const int& num_bases = DENEB_DATA->GetNumBases();

  max_num_points_ = 0;
  max_num_cell_points_ = 0;
  max_num_face_points_ = 0;
  max_num_bdry_points_ = 0;
  const std::vector<int>& num_cell_points = DENEB_DATA->GetNumCellPoints();
  const std::vector<int>& num_face_points = DENEB_DATA->GetNumFacePoints();
  const std::vector<int>& num_bdry_points = DENEB_DATA->GetNumBdryPoints();
  if (num_cell_points.size() != 0)
    max_num_cell_points_ =
        *std::max_element(num_cell_points.begin(), num_cell_points.end());
  if (num_face_points.size() != 0)
    max_num_face_points_ =
        *std::max_element(num_face_points.begin(), num_face_points.end());
  if (num_bdry_points.size() != 0)
    max_num_bdry_points_ =
        *std::max_element(num_bdry_points.begin(), num_bdry_points.end());
  max_num_face_points_ = std::max(max_num_bdry_points_, max_num_face_points_);
  max_num_points_ = std::max(max_num_cell_points_, max_num_face_points_);

  const int& order = DENEB_DATA->GetOrder();
  const int& num_cells = DENEB_DATA->GetNumCells();
  static const auto& cell_volumes = DENEB_DATA->GetCellVolumes();
  static const auto& cell_proj_volumes = DENEB_DATA->GetCellProjVolumes();
  dt_auxiliary_.resize(num_cells);
  for (int icell = 0; icell < num_cells; icell++)
    dt_auxiliary_[icell] =
        D_ * static_cast<double>(2 * order + 1) / cell_volumes[icell] *
        avocado::VecInnerProd(D_, &cell_proj_volumes[icell * D_],
                              &cell_proj_volumes[icell * D_]);
  auxiliary_solution_.resize(num_cells * D_ * S_ * num_bases);
  pressure_fix_values_.resize(3);

  const int& num_outer_cells = DENEB_DATA->GetNumOuterCells();
  outer_solution_.resize(
      std::max((num_outer_cells - num_cells) * S_ * num_bases, 1));
  communicate_ = std::make_shared<avocado::Communicate>(
      S_ * num_bases, DENEB_DATA->GetOuterSendCellList(),
      DENEB_DATA->GetOuterRecvCellList());

  // dimensional
  cell_variable_names_ = {"rho",
                          "rhoU",
                          "rhoV",
                          "rhoE",
                          "SpecificInternalEnergy",
                          "Pressure",
                          "Temperature",
                          "Viscosity",
                          "Conductivity",
                          "SoundSpeed",
                          "SpecificEntropy",
                          "SpecificEnthalpy",
                          "SpecificHeatRatio",
                          "MachNumber"};
  face_variable_names_ = cell_variable_names_;

  SYNCRO();
  {
    std::string message = "Contour variables (cell): ";
    for (auto&& name : cell_variable_names_) message += (name + ", ");
    message.pop_back();
    message.pop_back();
    MASTER_MESSAGE(message + "\n");
    message = "Contour variables (face): ";
    for (auto&& name : face_variable_names_) message += (name + ", ");
    message.pop_back();
    message.pop_back();
    MASTER_MESSAGE(message + "\n");
  }
}
void EquationEquilibriumNS2D::GetCellPostSolution(
    const int icell, const int num_points, const std::vector<double>& solution,
    const std::vector<double>& solution_grad,
    std::vector<double>& post_solution) {
  int ind = 0;
  for (int ipoint = 0; ipoint < num_points; ipoint++) {
    const double d = solution[ipoint * S_] * rho_ref_;
    const double u = solution[ipoint * S_ + 1] / solution[ipoint * S_] * a_ref_;
    const double v = solution[ipoint * S_ + 2] / solution[ipoint * S_] * a_ref_;
    const double e = ComputeInternalEnergy(&solution[ipoint * S_]) * e_ref_;
    const double p = IDEA_DEP(d, e);
    const double T = IDEA_DET(d, e);
    const double mu = IDEA_DEV(d, e);
    const double k = IDEA_DEK(d, e);
    const double a = IDEA_DEA(d, e);
    const double s = IDEA_DES(d, e);
    const double h = IDEA_DEH(d, e);
    const double r = IDEA_DEG(d, e);
    const double V = std::sqrt(u * u + v * v);

    post_solution[ind++] = d;
    post_solution[ind++] = d * u;
    post_solution[ind++] = d * v;
    post_solution[ind++] = solution[ipoint * S_ + 3] * rho_ref_ * e_ref_;
    post_solution[ind++] = e;
    post_solution[ind++] = p;
    post_solution[ind++] = T;
    post_solution[ind++] = mu;
    post_solution[ind++] = k;
    post_solution[ind++] = a;
    post_solution[ind++] = s;
    post_solution[ind++] = h;
    post_solution[ind++] = r;
    post_solution[ind++] = V / a;
  }
}
void EquationEquilibriumNS2D::GetFacePostSolution(
    const int num_points, const std::vector<double>& solution,
    const std::vector<double>& solution_grad, const std::vector<double>& normal,
    std::vector<double>& post_solution) {
  int ind = 0;
  for (int ipoint = 0; ipoint < num_points; ipoint++) {
    const double d = solution[ipoint * S_] * rho_ref_;
    const double u = solution[ipoint * S_ + 1] / solution[ipoint * S_] * a_ref_;
    const double v = solution[ipoint * S_ + 2] / solution[ipoint * S_] * a_ref_;
    const double e = ComputeInternalEnergy(&solution[ipoint * S_]) * e_ref_;
    const double p = IDEA_DEP(d, e);
    const double T = IDEA_DET(d, e);
    const double mu = IDEA_DEV(d, e);
    const double k = IDEA_DEK(d, e);
    const double a = IDEA_DEA(d, e);
    const double s = IDEA_DES(d, e);
    const double h = IDEA_DEH(d, e);
    const double r = IDEA_DEG(d, e);
    const double V = std::sqrt(u * u + v * v);

    post_solution[ind++] = d;
    post_solution[ind++] = d * u;
    post_solution[ind++] = d * v;
    post_solution[ind++] = solution[ipoint * S_ + 3] * rho_ref_ * e_ref_;
    post_solution[ind++] = e;
    post_solution[ind++] = p;
    post_solution[ind++] = T;
    post_solution[ind++] = mu;
    post_solution[ind++] = k;
    post_solution[ind++] = a;
    post_solution[ind++] = s;
    post_solution[ind++] = h;
    post_solution[ind++] = r;
    post_solution[ind++] = V / a;
  }
}
void EquationEquilibriumNS2D::ComputeInitialSolution(double* solution,
                                                     const double t) {
  const int& order = DENEB_DATA->GetOrder();
  const int& num_bases = DENEB_DATA->GetNumBases();
  const int& num_cells = DENEB_DATA->GetNumCells();
  static const int sb = S_ * num_bases;

  memset(solution, 0, num_cells * sb * sizeof(double));

  std::vector<double> quad_points;
  std::vector<double> quad_weights;
  std::vector<double> initial_values;
  std::vector<double> basis_values;

  for (int icell = 0; icell < num_cells; icell++) {
    DENEB_DATA->GetCellQuadrature(icell, 2 * order, quad_points, quad_weights);
    DENEB_DATA->GetCellBasisValues(icell, quad_points, basis_values);

    const int num_points = static_cast<int>(quad_weights.size());
    initial_values.resize(num_points * S_);
    problem_->Problem(num_points, initial_values, quad_points, t);

    // solution reconstruction
    for (int istate = 0; istate < S_; istate++)
      for (int ibasis = 0; ibasis < num_bases; ibasis++)
        for (int ipoint = 0; ipoint < num_points; ipoint++)
          solution[icell * sb + istate * num_bases + ibasis] +=
              initial_values[ipoint * S_ + istate] *
              basis_values[ipoint * num_bases + ibasis] * quad_weights[ipoint];
  }
}
void EquationEquilibriumNS2D::ComputeLocalTimestep(
    const double* solution, std::vector<double>& local_timestep) {
  static const int& order = DENEB_DATA->GetOrder();
  static const int& num_bases = DENEB_DATA->GetNumBases();
  static const int sb = S_ * num_bases;

  static const int& num_cells = DENEB_DATA->GetNumCells();
  static const auto& cell_volumes = DENEB_DATA->GetCellVolumes();
  static const auto& cell_proj_volumes = DENEB_DATA->GetCellProjVolumes();
  static const auto& cell_basis_value = DENEB_DATA->GetCellBasisValue();
  static const double visradii = alpha_ / mu_ref_;
  double Tgrad[2];
  for (int icell = 0; icell < num_cells; icell++) {
    const double& Vx = cell_proj_volumes[icell * D_];
    const double& Vy = cell_proj_volumes[icell * D_ + 1];
    const double d = solution[icell * sb] * cell_basis_value[icell][0];
    const double d_inv = 1.0 / d;
    const double du =
        solution[icell * sb + 1 * num_bases] * cell_basis_value[icell][0];
    const double dv =
        solution[icell * sb + 2 * num_bases] * cell_basis_value[icell][0];
    const double dE =
        solution[icell * sb + 3 * num_bases] * cell_basis_value[icell][0];
    const double e = (dE - 0.5 * (du * du + dv * dv) * d_inv) * d_inv;
    const double dd = d * rho_ref_;
    const double ee = e * e_ref_;
    const double k = IDEA_DEK(dd, ee);
    const double mu = IDEA_DEV(dd, ee);
    const double a = IDEA_DEA(dd, ee) / a_ref_;
    IDEA_DET_Grad(Tgrad, dd, ee);
    const double var = d_inv * std::max(k * Tgrad[1], c43_ * mu);

    const double AVcoeff =
        DENEB_ARTIFICIAL_VISCOSITY->GetArtificialViscosityValue(icell, 0);

    local_timestep[icell] =
        cell_volumes[icell] /
        ((std::abs(du) * d_inv + a) * Vx + (std::abs(dv) * d_inv + a) * Vy +
         (var * visradii + AVcoeff) * dt_auxiliary_[icell]);
  }
  avocado::VecScale(num_cells, 1.0 / static_cast<double>(2 * order + 1),
                    &local_timestep[0]);
}
bool EquationEquilibriumNS2D::IsContact(const int& icell,
                                        const std::vector<int>& neighbor_cells,
                                        const double* solution,
                                        const double* total_solution) const {
  static const int& num_bases = DENEB_DATA->GetNumBases();
  static const int sb = S_ * num_bases;

  static const int& num_cells = DENEB_DATA->GetNumCells();
  static const auto& volume = DENEB_DATA->GetCellVolumes();

  std::vector<double> sol(S_);
  cblas_daxpby(S_, 1.0 / std::sqrt(volume[icell]), &solution[icell * sb],
               num_bases, 0.0, &sol[0], 1);

  const double d = sol[0] * rho_ref_;
  const double e = ComputeInternalEnergy(&sol[0]) * e_ref_;
  const double owner_p = IDEA_DEP(d, e) / p_ref_;

  for (auto&& cell : neighbor_cells) {
    const double* target = (cell < num_cells)
                               ? &solution[cell * sb]
                               : &total_solution[(cell - num_cells) * sb];
    cblas_daxpby(S_, 1.0 / std::sqrt(volume[cell]), target, num_bases, 0.0,
                 &sol[0], 1);

    const double d = sol[0] * rho_ref_;
    const double e = ComputeInternalEnergy(&sol[0]) * e_ref_;
    const double neighbor_p = IDEA_DEP(d, e) / p_ref_;

    if (std::abs(owner_p - neighbor_p) > 1.0E-2 * owner_p) return false;
  }
  return true;
}
double EquationEquilibriumNS2D::ComputeMaxCharacteristicSpeed(
    const double* input_solution) const {
  const double d_inv = 1.0 / input_solution[0];
  const double p = ComputePressure(&input_solution[0]);
  const double e = ComputeInternalEnergy(&input_solution[0]);
  const double a = IDEA_Wrapper(IDEA_ROUTINES(DEA),
                                input_solution[0] * rho_ref_, e * e_ref_) /
                   a_ref_;
  const double V = std::sqrt(avocado::VecInnerProd(D_, &input_solution[1],
                                                   &input_solution[1])) *
                   d_inv;
  return V + a;
}
const std::vector<double>& EquationEquilibriumNS2D::ComputePressureFixValues(
    const double* input_solution) {
  pressure_fix_values_[0] = input_solution[0];
  pressure_fix_values_[1] = ComputePressure(&input_solution[0]);
  pressure_fix_values_[2] = ComputeInternalEnergy(&input_solution[0]);
  return pressure_fix_values_;
}
void EquationEquilibriumNS2D::ComputeRHS(const double* solution, double* rhs,
                                         const double t) {
  static const int& num_bases = DENEB_DATA->GetNumBases();
  static const int& num_cells = DENEB_DATA->GetNumCells();
  static const int sb = S_ * num_bases;
  static const int dsb = DS_ * num_bases;

  static std::vector<double> owner_solution(S_ * max_num_points_);
  static std::vector<double> owner_solution_grad(DS_ * max_num_points_);
  static std::vector<double> flux(DS_ * max_num_points_);
  static std::vector<double> neighbor_solution(S_ * max_num_face_points_);
  static std::vector<double> neighbor_solution_grad(DS_ * max_num_face_points_);
  static std::vector<double> solution_difference(S_ * max_num_face_points_);
  static std::vector<double> local_auxiliary(dsb);

  communicate_->CommunicateBegin(solution);

  memset(&rhs[0], 0, num_cells * sb * sizeof(double));
  memset(&auxiliary_solution_[0], 0, num_cells * dsb * sizeof(double));

  // inner face sweep
  static const int& num_inner_faces = DENEB_DATA->GetNumInnerFaces();
  static const auto& num_face_points = DENEB_DATA->GetNumFacePoints();
  static const auto& face_owner_cell = DENEB_DATA->GetFaceOwnerCell();
  static const auto& face_neighbor_cell = DENEB_DATA->GetFaceNeighborCell();
  static const auto& face_normals = DENEB_DATA->GetFaceNormals();
  static const auto& face_owner_basis_value =
      DENEB_DATA->GetFaceOwnerBasisValue();
  static const auto& face_neighbor_basis_value =
      DENEB_DATA->GetFaceNeighborBasisValue();
  static const auto& face_owner_basis_grad_value =
      DENEB_DATA->GetFaceOwnerBasisGradValue();
  static const auto& face_neighbor_basis_grad_value =
      DENEB_DATA->GetFaceNeighborBasisGradValue();
  static const auto& face_owner_coefficients =
      DENEB_DATA->GetFaceOwnerCoefficients();
  static const auto& face_neighbor_coefficients =
      DENEB_DATA->GetFaceNeighborCoefficients();
  for (int iface = 0; iface < num_inner_faces; iface++) {
    const int& num_points = num_face_points[iface];
    const int& owner_cell = face_owner_cell[iface];
    const int& neighbor_cell = face_neighbor_cell[iface];

    avocado::Kernel0::f4(&solution[owner_cell * sb],
                         &face_owner_basis_value[iface][0], &owner_solution[0],
                         S_, num_bases, num_points, 1.0, 0.0);
    avocado::Kernel0::f4(
        &solution[neighbor_cell * sb], &face_neighbor_basis_value[iface][0],
        &neighbor_solution[0], S_, num_bases, num_points, 1.0, 0.0);

    vdSub(num_points * S_, &neighbor_solution[0], &owner_solution[0],
          &solution_difference[0]);

    avocado::Kernel1::f59(&face_owner_coefficients[iface][0],
                          &solution_difference[0], &local_auxiliary[0], D_,
                          num_bases, num_points, S_, 0.5, 0.0);
    avocado::Kernel1::f42(&face_owner_basis_grad_value[iface][0],
                          &solution[owner_cell * sb], &owner_solution_grad[0],
                          D_, num_points, num_bases, S_, 1.0, 0.0);
    avocado::Kernel1::f18(
        &local_auxiliary[0], &face_owner_basis_value[iface][0],
        &owner_solution_grad[0], D_, S_, num_bases, num_points, 1.0, 1.0);
    cblas_daxpy(dsb, 1.0, &local_auxiliary[0], 1,
                &auxiliary_solution_[owner_cell * dsb], 1);

    avocado::Kernel1::f59(&face_neighbor_coefficients[iface][0],
                          &solution_difference[0], &local_auxiliary[0], D_,
                          num_bases, num_points, S_, 0.5, 0.0);
    avocado::Kernel1::f42(&face_neighbor_basis_grad_value[iface][0],
                          &solution[neighbor_cell * sb],
                          &neighbor_solution_grad[0], D_, num_points, num_bases,
                          S_, 1.0, 0.0);
    avocado::Kernel1::f18(
        &local_auxiliary[0], &face_neighbor_basis_value[iface][0],
        &neighbor_solution_grad[0], D_, S_, num_bases, num_points, 1.0, 1.0);
    cblas_daxpy(dsb, 1.0, &local_auxiliary[0], 1,
                &auxiliary_solution_[neighbor_cell * dsb], 1);

    ComputeNumFlux(num_points, flux, owner_cell, neighbor_cell, owner_solution,
                   owner_solution_grad, neighbor_solution,
                   neighbor_solution_grad, face_normals[iface]);

    avocado::Kernel2::f67(&flux[0], &face_owner_coefficients[iface][0],
                          &rhs[owner_cell * sb], S_, D_, num_points, num_bases,
                          1.0, 1.0);
    avocado::Kernel2::f67(&flux[0], &face_neighbor_coefficients[iface][0],
                          &rhs[neighbor_cell * sb], S_, D_, num_points,
                          num_bases, -1.0, 1.0);
  }

  // bdry sweep
  static const int& num_bdries = DENEB_DATA->GetNumBdries();
  static const auto& num_bdry_points = DENEB_DATA->GetNumBdryPoints();
  static const auto& bdry_owner_cell = DENEB_DATA->GetBdryOwnerCell();
  static const auto& bdry_normals = DENEB_DATA->GetBdryNormals();
  static const auto& bdry_points = DENEB_DATA->GetBdryPoints();
  static const auto& bdry_owner_basis_value =
      DENEB_DATA->GetBdryOwnerBasisValue();
  static const auto& bdry_owner_basis_grad_value =
      DENEB_DATA->GetBdryOwnerBasisGradValue();
  static const auto& bdry_owner_coefficients =
      DENEB_DATA->GetBdryOwnerCoefficients();
  for (int ibdry = 0; ibdry < num_bdries; ibdry++) {
    const int& num_points = num_bdry_points[ibdry];
    const int& owner_cell = bdry_owner_cell[ibdry];

    avocado::Kernel0::f4(&solution[owner_cell * sb],
                         &bdry_owner_basis_value[ibdry][0], &owner_solution[0],
                         S_, num_bases, num_points, 1.0, 0.0);
    boundaries_[ibdry]->ComputeBdrySolution(
        num_points, neighbor_solution, neighbor_solution_grad, owner_solution,
        owner_solution_grad, bdry_normals[ibdry], bdry_points[ibdry], t);

    vdSub(num_points * S_, &neighbor_solution[0], &owner_solution[0],
          &solution_difference[0]);

    avocado::Kernel1::f59(&bdry_owner_coefficients[ibdry][0],
                          &solution_difference[0], &local_auxiliary[0], D_,
                          num_bases, num_points, S_, 1.0, 0.0);
    avocado::Kernel1::f42(&bdry_owner_basis_grad_value[ibdry][0],
                          &solution[owner_cell * sb], &owner_solution_grad[0],
                          D_, num_points, num_bases, S_, 1.0, 0.0);
    avocado::Kernel1::f18(
        &local_auxiliary[0], &bdry_owner_basis_value[ibdry][0],
        &owner_solution_grad[0], D_, S_, num_bases, num_points, 1.0, 1.0);
    cblas_daxpy(dsb, 1.0, &local_auxiliary[0], 1,
                &auxiliary_solution_[owner_cell * dsb], 1);

    boundaries_[ibdry]->ComputeBdryFlux(
        num_points, flux, owner_cell, -1, owner_solution, owner_solution_grad,
        neighbor_solution, neighbor_solution_grad, bdry_normals[ibdry],
        bdry_points[ibdry], t);

    avocado::Kernel2::f67(&flux[0], &bdry_owner_coefficients[ibdry][0],
                          &rhs[owner_cell * sb], S_, D_, num_points, num_bases,
                          1.0, 1.0);
  }

  communicate_->CommunicateEnd(&outer_solution_[0]);

  // outer face sweep
  static const int& num_faces = DENEB_DATA->GetNumFaces();
  for (int iface = num_inner_faces; iface < num_faces; iface++) {
    const int& num_points = num_face_points[iface];
    const int& owner_cell = face_owner_cell[iface];
    const int& neighbor_cell = face_neighbor_cell[iface];

    avocado::Kernel0::f4(&solution[owner_cell * sb],
                         &face_owner_basis_value[iface][0], &owner_solution[0],
                         S_, num_bases, num_points, 1.0, 0.0);
    avocado::Kernel0::f4(&outer_solution_[neighbor_cell * sb],
                         &face_neighbor_basis_value[iface][0],
                         &neighbor_solution[0], S_, num_bases, num_points, 1.0,
                         0.0);

    vdSub(num_points * S_, &neighbor_solution[0], &owner_solution[0],
          &solution_difference[0]);

    avocado::Kernel1::f59(&face_owner_coefficients[iface][0],
                          &solution_difference[0], &local_auxiliary[0], D_,
                          num_bases, num_points, S_, 0.5, 0.0);
    avocado::Kernel1::f42(&face_owner_basis_grad_value[iface][0],
                          &solution[owner_cell * sb], &owner_solution_grad[0],
                          D_, num_points, num_bases, S_, 1.0, 0.0);
    avocado::Kernel1::f18(
        &local_auxiliary[0], &face_owner_basis_value[iface][0],
        &owner_solution_grad[0], D_, S_, num_bases, num_points, 1.0, 1.0);
    cblas_daxpy(dsb, 1.0, &local_auxiliary[0], 1,
                &auxiliary_solution_[owner_cell * dsb], 1);

    avocado::Kernel1::f59(&face_neighbor_coefficients[iface][0],
                          &solution_difference[0], &local_auxiliary[0], D_,
                          num_bases, num_points, S_, 0.5, 0.0);
    avocado::Kernel1::f42(&face_neighbor_basis_grad_value[iface][0],
                          &outer_solution_[neighbor_cell * sb],
                          &neighbor_solution_grad[0], D_, num_points, num_bases,
                          S_, 1.0, 0.0);
    avocado::Kernel1::f18(
        &local_auxiliary[0], &face_neighbor_basis_value[iface][0],
        &neighbor_solution_grad[0], D_, S_, num_bases, num_points, 1.0, 1.0);

    ComputeNumFlux(num_points, flux, owner_cell, neighbor_cell + num_cells,
                   owner_solution, owner_solution_grad, neighbor_solution,
                   neighbor_solution_grad, face_normals[iface]);

    avocado::Kernel2::f67(&flux[0], &face_owner_coefficients[iface][0],
                          &rhs[owner_cell * sb], S_, D_, num_points, num_bases,
                          1.0, 1.0);
  }

  // cell sweep
  static const auto& num_cell_points = DENEB_DATA->GetNumCellPoints();
  static const auto& cell_basis_value = DENEB_DATA->GetCellBasisValue();
  static const auto& cell_basis_grad_value =
      DENEB_DATA->GetCellBasisGradValue();
  static const auto& cell_coefficients = DENEB_DATA->GetCellCoefficients();
  for (int icell = 0; icell < num_cells; icell++) {
    const int& num_points = num_cell_points[icell];

    avocado::Kernel0::f4(&solution[icell * sb], &cell_basis_value[icell][0],
                         &owner_solution[0], S_, num_bases, num_points, 1.0,
                         0.0);
    avocado::Kernel1::f42(&cell_basis_grad_value[icell][0],
                          &solution[icell * sb], &owner_solution_grad[0], D_,
                          num_points, num_bases, S_, 1.0, 0.0);
    avocado::Kernel1::f18(&auxiliary_solution_[icell * dsb],
                          &cell_basis_value[icell][0], &owner_solution_grad[0],
                          D_, S_, num_bases, num_points, 1.0, 1.0);

    ComputeComFlux(num_points, flux, icell, owner_solution,
                   owner_solution_grad);

    avocado::Kernel2::f67(&flux[0], &cell_coefficients[icell][0],
                          &rhs[icell * sb], S_, D_, num_points, num_bases, -1.0,
                          1.0);
  }
}
void EquationEquilibriumNS2D::ComputeSystemMatrix(const double* solution,
                                                  Mat& sysmat, const double t) {
  static const int& num_bases = DENEB_DATA->GetNumBases();
  static const int& num_cells = DENEB_DATA->GetNumCells();
  static const int& num_bdries = DENEB_DATA->GetNumBdries();

  static const int db = D_ * num_bases;
  static const int sb = S_ * num_bases;
  static const int dsb = DS_ * num_bases;
  static const int dbb = db * num_bases;
  static const int ssb = SS_ * num_bases;
  static const int ssp = SS_ * max_num_bdry_points_;
  static const int dssb = DSS_ * num_bases;
  static const int ssbb = ssb * num_bases;

  static const std::vector<double> identity(S_, 1.0);

  static std::vector<double> owner_solution(S_ * max_num_points_);
  static std::vector<double> owner_solution_grad(DS_ * max_num_points_);
  static std::vector<double> flux_owner_jacobi(DSS_ * max_num_points_);
  static std::vector<double> flux_owner_grad_jacobi(DDSS_ * max_num_points_);
  static std::vector<double> owner_a(db * max_num_points_);
  static std::vector<double> owner_b(db * max_num_points_);
  static std::vector<double> largeB(dssb * max_num_points_);
  static std::vector<double> flux_derivative(dssb * max_num_points_);

  static std::vector<double> neighbor_solution(S_ * max_num_face_points_);
  static std::vector<double> neighbor_solution_grad(DS_ * max_num_face_points_);
  static std::vector<double> flux_neighbor_jacobi(DSS_ * max_num_face_points_);
  static std::vector<double> flux_neighbor_grad_jacobi(DDSS_ *
                                                       max_num_face_points_);
  static std::vector<double> solution_difference(S_ * max_num_face_points_);
  static std::vector<double> alpha(D_ * max_num_points_ * max_num_face_points_);
  static std::vector<double> largeA(std::max(ssb * max_num_face_points_, dbb));
  static std::vector<double> neighbor_a(db * max_num_face_points_);
  static std::vector<double> neighbor_b(db * max_num_face_points_);
  static std::vector<double> local_auxiliary(dsb);
  static std::vector<double> block(ssbb);
  static std::vector<double> bdry_solution_jacobi(num_bdries * ssp);

  communicate_->CommunicateBegin(solution);

  MatZeroEntries(sysmat);
  memset(&auxiliary_solution_[0], 0, num_cells * dsb * sizeof(double));

  // inner face sweep
  static const auto& mat_index = DENEB_DATA->GetMatIndex();
  static const int& num_inner_faces = DENEB_DATA->GetNumInnerFaces();
  static const auto& num_face_points = DENEB_DATA->GetNumFacePoints();
  static const auto& face_owner_cell = DENEB_DATA->GetFaceOwnerCell();
  static const auto& face_neighbor_cell = DENEB_DATA->GetFaceNeighborCell();
  static const auto& face_normals = DENEB_DATA->GetFaceNormals();
  static const auto& face_owner_basis_value =
      DENEB_DATA->GetFaceOwnerBasisValue();
  static const auto& face_neighbor_basis_value =
      DENEB_DATA->GetFaceNeighborBasisValue();
  static const auto& face_owner_basis_grad_value =
      DENEB_DATA->GetFaceOwnerBasisGradValue();
  static const auto& face_neighbor_basis_grad_value =
      DENEB_DATA->GetFaceNeighborBasisGradValue();
  static const auto& face_owner_coefficients =
      DENEB_DATA->GetFaceOwnerCoefficients();
  static const auto& face_neighbor_coefficients =
      DENEB_DATA->GetFaceNeighborCoefficients();
  for (int iface = 0; iface < num_inner_faces; iface++) {
    const int& num_points = num_face_points[iface];
    const int& owner_cell = face_owner_cell[iface];
    const int& neighbor_cell = face_neighbor_cell[iface];

    avocado::Kernel0::f4(&solution[owner_cell * sb],
                         &face_owner_basis_value[iface][0], &owner_solution[0],
                         S_, num_bases, num_points, 1.0, 0.0);
    avocado::Kernel0::f4(
        &solution[neighbor_cell * sb], &face_neighbor_basis_value[iface][0],
        &neighbor_solution[0], S_, num_bases, num_points, 1.0, 0.0);

    vdSub(num_points * S_, &neighbor_solution[0], &owner_solution[0],
          &solution_difference[0]);

    avocado::Kernel1::f59(&face_owner_coefficients[iface][0],
                          &solution_difference[0], &local_auxiliary[0], D_,
                          num_bases, num_points, S_, 0.5, 0.0);
    avocado::Kernel1::f42(&face_owner_basis_grad_value[iface][0],
                          &solution[owner_cell * sb], &owner_solution_grad[0],
                          D_, num_points, num_bases, S_, 1.0, 0.0);
    avocado::Kernel1::f18(
        &local_auxiliary[0], &face_owner_basis_value[iface][0],
        &owner_solution_grad[0], D_, S_, num_bases, num_points, 1.0, 1.0);
    cblas_daxpy(dsb, 1.0, &local_auxiliary[0], 1,
                &auxiliary_solution_[owner_cell * dsb], 1);

    avocado::Kernel1::f59(&face_neighbor_coefficients[iface][0],
                          &solution_difference[0], &local_auxiliary[0], D_,
                          num_bases, num_points, S_, 0.5, 0.0);
    avocado::Kernel1::f42(&face_neighbor_basis_grad_value[iface][0],
                          &solution[neighbor_cell * sb],
                          &neighbor_solution_grad[0], D_, num_points, num_bases,
                          S_, 1.0, 0.0);
    avocado::Kernel1::f18(
        &local_auxiliary[0], &face_neighbor_basis_value[iface][0],
        &neighbor_solution_grad[0], D_, S_, num_bases, num_points, 1.0, 1.0);
    cblas_daxpy(dsb, 1.0, &local_auxiliary[0], 1,
                &auxiliary_solution_[neighbor_cell * dsb], 1);

    ComputeNumFluxJacobi(num_points, flux_owner_jacobi, flux_neighbor_jacobi,
                         flux_owner_grad_jacobi, flux_neighbor_grad_jacobi,
                         owner_cell, neighbor_cell, owner_solution,
                         owner_solution_grad, neighbor_solution,
                         neighbor_solution_grad, face_normals[iface]);

    avocado::Kernel1::f47(&face_owner_coefficients[iface][0],
                          &face_owner_basis_value[iface][0], &alpha[0], D_,
                          num_points, num_bases, num_points, 1.0, 0.0);
    cblas_dcopy(db * num_points, &face_owner_basis_grad_value[iface][0], 1,
                &owner_a[0], 1);
    avocado::Kernel1::f50(&alpha[0], &face_owner_basis_value[iface][0],
                          &owner_a[0], D_, num_points, num_points, num_bases,
                          -0.5, 1.0);
    avocado::Kernel1::f50(&alpha[0], &face_neighbor_basis_value[iface][0],
                          &neighbor_b[0], D_, num_points, num_points, num_bases,
                          0.5, 0.0);

    avocado::Kernel1::f47(&face_neighbor_coefficients[iface][0],
                          &face_neighbor_basis_value[iface][0], &alpha[0], D_,
                          num_points, num_bases, num_points, 1.0, 0.0);
    cblas_dcopy(db * num_points, &face_neighbor_basis_grad_value[iface][0], 1,
                &neighbor_a[0], 1);
    avocado::Kernel1::f50(&alpha[0], &face_neighbor_basis_value[iface][0],
                          &neighbor_a[0], D_, num_points, num_points, num_bases,
                          0.5, 1.0);
    avocado::Kernel1::f50(&alpha[0], &face_owner_basis_value[iface][0],
                          &owner_b[0], D_, num_points, num_points, num_bases,
                          -0.5, 0.0);

    for (int ipoint = 0; ipoint < num_points; ipoint++) {
      gemmAB(1.0, &flux_owner_grad_jacobi[ipoint * DDSS_],
             &owner_a[ipoint * db], 0.0, &flux_derivative[ipoint * dssb], DSS_,
             D_, num_bases);
      gemmAB(1.0, &flux_neighbor_grad_jacobi[ipoint * DDSS_],
             &owner_b[ipoint * db], 1.0, &flux_derivative[ipoint * dssb], DSS_,
             D_, num_bases);
      cblas_dger(CBLAS_LAYOUT::CblasRowMajor, DSS_, num_bases, 1.0,
                 &flux_owner_jacobi[ipoint * DSS_], 1,
                 &face_owner_basis_value[iface][ipoint * num_bases], 1,
                 &flux_derivative[ipoint * dssb], num_bases);
    }
    avocado::Kernel1::f59(&flux_derivative[0],
                          &face_owner_coefficients[iface][0], &block[0], S_, sb,
                          num_points * D_, num_bases, 1.0, 0.0);
    MatSetValuesBlocked(sysmat, 1, &mat_index[owner_cell], 1,
                        &mat_index[owner_cell], &block[0], ADD_VALUES);
    avocado::Kernel1::f59(&flux_derivative[0],
                          &face_neighbor_coefficients[iface][0], &block[0], S_,
                          sb, num_points * D_, num_bases, -1.0, 0.0);
    MatSetValuesBlocked(sysmat, 1, &mat_index[neighbor_cell], 1,
                        &mat_index[owner_cell], &block[0], ADD_VALUES);

    for (int ipoint = 0; ipoint < num_points; ipoint++) {
      gemmAB(1.0, &flux_owner_grad_jacobi[ipoint * DDSS_],
             &neighbor_b[ipoint * db], 0.0, &flux_derivative[ipoint * dssb],
             DSS_, D_, num_bases);
      gemmAB(1.0, &flux_neighbor_grad_jacobi[ipoint * DDSS_],
             &neighbor_a[ipoint * db], 1.0, &flux_derivative[ipoint * dssb],
             DSS_, D_, num_bases);
      cblas_dger(CBLAS_LAYOUT::CblasRowMajor, DSS_, num_bases, 1.0,
                 &flux_neighbor_jacobi[ipoint * DSS_], 1,
                 &face_neighbor_basis_value[iface][ipoint * num_bases], 1,
                 &flux_derivative[ipoint * dssb], num_bases);
    }
    avocado::Kernel1::f59(&flux_derivative[0],
                          &face_owner_coefficients[iface][0], &block[0], S_, sb,
                          num_points * D_, num_bases, 1.0, 0.0);
    MatSetValuesBlocked(sysmat, 1, &mat_index[owner_cell], 1,
                        &mat_index[neighbor_cell], &block[0], ADD_VALUES);
    avocado::Kernel1::f59(&flux_derivative[0],
                          &face_neighbor_coefficients[iface][0], &block[0], S_,
                          sb, num_points * D_, num_bases, -1.0, 0.0);
    MatSetValuesBlocked(sysmat, 1, &mat_index[neighbor_cell], 1,
                        &mat_index[neighbor_cell], &block[0], ADD_VALUES);
  }

  // bdry sweep
  static const auto& num_bdry_points = DENEB_DATA->GetNumBdryPoints();
  static const auto& bdry_owner_cell = DENEB_DATA->GetBdryOwnerCell();
  static const auto& bdry_normals = DENEB_DATA->GetBdryNormals();
  static const auto& bdry_points = DENEB_DATA->GetBdryPoints();
  static const auto& bdry_owner_basis_value =
      DENEB_DATA->GetBdryOwnerBasisValue();
  static const auto& bdry_owner_basis_grad_value =
      DENEB_DATA->GetBdryOwnerBasisGradValue();
  static const auto& bdry_owner_coefficients =
      DENEB_DATA->GetBdryOwnerCoefficients();
  for (int ibdry = 0; ibdry < num_bdries; ibdry++) {
    const int& num_points = num_bdry_points[ibdry];
    const int& owner_cell = bdry_owner_cell[ibdry];

    avocado::Kernel0::f4(&solution[owner_cell * sb],
                         &bdry_owner_basis_value[ibdry][0], &owner_solution[0],
                         S_, num_bases, num_points, 1.0, 0.0);
    boundaries_[ibdry]->ComputeBdrySolution(
        num_points, neighbor_solution, neighbor_solution_grad, owner_solution,
        owner_solution_grad, bdry_normals[ibdry], bdry_points[ibdry], t);

    vdSub(num_points * S_, &neighbor_solution[0], &owner_solution[0],
          &solution_difference[0]);

    avocado::Kernel1::f59(&bdry_owner_coefficients[ibdry][0],
                          &solution_difference[0], &local_auxiliary[0], D_,
                          num_bases, num_points, S_, 1.0, 0.0);
    avocado::Kernel1::f42(&bdry_owner_basis_grad_value[ibdry][0],
                          &solution[owner_cell * sb], &owner_solution_grad[0],
                          D_, num_points, num_bases, S_, 1.0, 0.0);
    avocado::Kernel1::f18(
        &local_auxiliary[0], &bdry_owner_basis_value[ibdry][0],
        &owner_solution_grad[0], D_, S_, num_bases, num_points, 1.0, 1.0);
    cblas_daxpy(dsb, 1.0, &local_auxiliary[0], 1,
                &auxiliary_solution_[owner_cell * dsb], 1);

    double* solution_jacobi = &bdry_solution_jacobi[ibdry * ssp];
    boundaries_[ibdry]->ComputeBdrySolutionJacobi(
        num_points, solution_jacobi, owner_solution, owner_solution_grad,
        bdry_normals[ibdry], bdry_points[ibdry], t);
    boundaries_[ibdry]->ComputeBdryFluxJacobi(
        num_points, flux_owner_jacobi, flux_owner_grad_jacobi, owner_cell, -1,
        owner_solution, owner_solution_grad, neighbor_solution,
        neighbor_solution_grad, bdry_normals[ibdry], bdry_points[ibdry], t);

    memset(&largeA[0], 0, num_points * ssb * sizeof(double));
    for (int ipoint = 0; ipoint < num_points; ipoint++) {
      cblas_daxpy(S_, -1.0, &identity[0], 1, &solution_jacobi[ipoint * SS_],
                  S_ + 1);
      cblas_dger(CBLAS_LAYOUT::CblasRowMajor, SS_, num_bases, 1.0,
                 &solution_jacobi[ipoint * SS_], 1,
                 &bdry_owner_basis_value[ibdry][ipoint * num_bases], 1,
                 &largeA[ipoint * ssb], num_bases);
    }
    avocado::Kernel1::f47(&bdry_owner_coefficients[ibdry][0],
                          &bdry_owner_basis_value[ibdry][0], &alpha[0], D_,
                          num_points, num_bases, num_points, 1.0, 0.0);
    avocado::Kernel1::f21(&alpha[0], &largeA[0], &largeB[0], num_points, D_,
                          num_points, ssb, 1.0, 0.0);

    for (int ipoint = 0; ipoint < num_points; ipoint++) {
      avocado::Kernel2::f13(
          &flux_owner_grad_jacobi[ipoint * DDSS_], &largeB[ipoint * dssb],
          &flux_derivative[ipoint * dssb], DS_, D_, S_, sb, 1.0, 0.0);
      gemmAB(1.0, &flux_owner_grad_jacobi[ipoint * DDSS_],
             &bdry_owner_basis_grad_value[ibdry][ipoint * db], 1.0,
             &flux_derivative[ipoint * dssb], DSS_, D_, num_bases);
      cblas_dger(CBLAS_LAYOUT::CblasRowMajor, DSS_, num_bases, 1.0,
                 &flux_owner_jacobi[ipoint * DSS_], 1,
                 &bdry_owner_basis_value[ibdry][ipoint * num_bases], 1,
                 &flux_derivative[ipoint * dssb], num_bases);
    }
    avocado::Kernel1::f59(&flux_derivative[0],
                          &bdry_owner_coefficients[ibdry][0], &block[0], S_, sb,
                          num_points * D_, num_bases, 1.0, 0.0);
    MatSetValuesBlocked(sysmat, 1, &mat_index[owner_cell], 1,
                        &mat_index[owner_cell], &block[0], ADD_VALUES);
  }

  communicate_->CommunicateEnd(&outer_solution_[0]);

  // outer face sweep
  static const int& num_faces = DENEB_DATA->GetNumFaces();
  for (int iface = num_inner_faces; iface < num_faces; iface++) {
    const int& num_points = num_face_points[iface];
    const int& owner_cell = face_owner_cell[iface];
    const int& neighbor_cell = face_neighbor_cell[iface];

    avocado::Kernel0::f4(&solution[owner_cell * sb],
                         &face_owner_basis_value[iface][0], &owner_solution[0],
                         S_, num_bases, num_points, 1.0, 0.0);
    avocado::Kernel0::f4(&outer_solution_[neighbor_cell * sb],
                         &face_neighbor_basis_value[iface][0],
                         &neighbor_solution[0], S_, num_bases, num_points, 1.0,
                         0.0);

    vdSub(num_points * S_, &neighbor_solution[0], &owner_solution[0],
          &solution_difference[0]);

    avocado::Kernel1::f59(&face_owner_coefficients[iface][0],
                          &solution_difference[0], &local_auxiliary[0], D_,
                          num_bases, num_points, S_, 0.5, 0.0);
    avocado::Kernel1::f42(&face_owner_basis_grad_value[iface][0],
                          &solution[owner_cell * sb], &owner_solution_grad[0],
                          D_, num_points, num_bases, S_, 1.0, 0.0);
    avocado::Kernel1::f18(
        &local_auxiliary[0], &face_owner_basis_value[iface][0],
        &owner_solution_grad[0], D_, S_, num_bases, num_points, 1.0, 1.0);
    cblas_daxpy(dsb, 1.0, &local_auxiliary[0], 1,
                &auxiliary_solution_[owner_cell * dsb], 1);

    avocado::Kernel1::f59(&face_neighbor_coefficients[iface][0],
                          &solution_difference[0], &local_auxiliary[0], D_,
                          num_bases, num_points, S_, 0.5, 0.0);
    avocado::Kernel1::f42(&face_neighbor_basis_grad_value[iface][0],
                          &outer_solution_[neighbor_cell * sb],
                          &neighbor_solution_grad[0], D_, num_points, num_bases,
                          S_, 1.0, 0.0);
    avocado::Kernel1::f18(
        &local_auxiliary[0], &face_neighbor_basis_value[iface][0],
        &neighbor_solution_grad[0], D_, S_, num_bases, num_points, 1.0, 1.0);

    ComputeNumFluxJacobi(num_points, flux_owner_jacobi, flux_neighbor_jacobi,
                         flux_owner_grad_jacobi, flux_neighbor_grad_jacobi,
                         owner_cell, neighbor_cell + num_cells, owner_solution,
                         owner_solution_grad, neighbor_solution,
                         neighbor_solution_grad, face_normals[iface]);

    avocado::Kernel1::f47(&face_owner_coefficients[iface][0],
                          &face_owner_basis_value[iface][0], &alpha[0], D_,
                          num_points, num_bases, num_points, 1.0, 0.0);
    cblas_dcopy(db * num_points, &face_owner_basis_grad_value[iface][0], 1,
                &owner_a[0], 1);
    avocado::Kernel1::f50(&alpha[0], &face_owner_basis_value[iface][0],
                          &owner_a[0], D_, num_points, num_points, num_bases,
                          -0.5, 1.0);
    avocado::Kernel1::f50(&alpha[0], &face_neighbor_basis_value[iface][0],
                          &neighbor_b[0], D_, num_points, num_points, num_bases,
                          0.5, 0.0);

    avocado::Kernel1::f47(&face_neighbor_coefficients[iface][0],
                          &face_neighbor_basis_value[iface][0], &alpha[0], D_,
                          num_points, num_bases, num_points, 1.0, 0.0);
    cblas_dcopy(db * num_points, &face_neighbor_basis_grad_value[iface][0], 1,
                &neighbor_a[0], 1);
    avocado::Kernel1::f50(&alpha[0], &face_neighbor_basis_value[iface][0],
                          &neighbor_a[0], D_, num_points, num_points, num_bases,
                          0.5, 1.0);
    avocado::Kernel1::f50(&alpha[0], &face_owner_basis_value[iface][0],
                          &owner_b[0], D_, num_points, num_points, num_bases,
                          -0.5, 0.0);

    for (int ipoint = 0; ipoint < num_points; ipoint++) {
      gemmAB(1.0, &flux_owner_grad_jacobi[ipoint * DDSS_],
             &owner_a[ipoint * db], 0.0, &flux_derivative[ipoint * dssb], DSS_,
             D_, num_bases);
      gemmAB(1.0, &flux_neighbor_grad_jacobi[ipoint * DDSS_],
             &owner_b[ipoint * db], 1.0, &flux_derivative[ipoint * dssb], DSS_,
             D_, num_bases);
      cblas_dger(CBLAS_LAYOUT::CblasRowMajor, DSS_, num_bases, 1.0,
                 &flux_owner_jacobi[ipoint * DSS_], 1,
                 &face_owner_basis_value[iface][ipoint * num_bases], 1,
                 &flux_derivative[ipoint * dssb], num_bases);
    }
    avocado::Kernel1::f59(&flux_derivative[0],
                          &face_owner_coefficients[iface][0], &block[0], S_, sb,
                          num_points * D_, num_bases, 1.0, 0.0);
    MatSetValuesBlocked(sysmat, 1, &mat_index[owner_cell], 1,
                        &mat_index[owner_cell], &block[0], ADD_VALUES);

    for (int ipoint = 0; ipoint < num_points; ipoint++) {
      gemmAB(1.0, &flux_owner_grad_jacobi[ipoint * DDSS_],
             &neighbor_b[ipoint * db], 0.0, &flux_derivative[ipoint * dssb],
             DSS_, D_, num_bases);
      gemmAB(1.0, &flux_neighbor_grad_jacobi[ipoint * DDSS_],
             &neighbor_a[ipoint * db], 1.0, &flux_derivative[ipoint * dssb],
             DSS_, D_, num_bases);
      cblas_dger(CBLAS_LAYOUT::CblasRowMajor, DSS_, num_bases, 1.0,
                 &flux_neighbor_jacobi[ipoint * DSS_], 1,
                 &face_neighbor_basis_value[iface][ipoint * num_bases], 1,
                 &flux_derivative[ipoint * dssb], num_bases);
    }
    avocado::Kernel1::f59(&flux_derivative[0],
                          &face_owner_coefficients[iface][0], &block[0], S_, sb,
                          num_points * D_, num_bases, 1.0, 0.0);
    MatSetValuesBlocked(sysmat, 1, &mat_index[owner_cell], 1,
                        &mat_index[neighbor_cell + num_cells], &block[0],
                        ADD_VALUES);
  }

  // cell sweep
  static const auto& num_cell_points = DENEB_DATA->GetNumCellPoints();
  static const auto& cell_basis_value = DENEB_DATA->GetCellBasisValue();
  static const auto& cell_basis_grad_value =
      DENEB_DATA->GetCellBasisGradValue();
  static const auto& cell_coefficients = DENEB_DATA->GetCellCoefficients();
  static const auto& subface_ptr = DENEB_DATA->GetSubfacePtr();
  static const auto& subface_sign = DENEB_DATA->GetSubfaceSign();
  static const auto& subface_neighbor_cells =
      DENEB_DATA->GetSubfaceNeighborCell();
  static const auto& subface_num_points = DENEB_DATA->GetSubfaceNumPoints();
  static const auto& subface_owner_coefficients =
      DENEB_DATA->GetSubfaceOwnerCoefficients();
  static const auto& subface_owner_basis_value =
      DENEB_DATA->GetSubfaceOwnerBasisValue();
  static const auto& subface_neighbor_basis_value =
      DENEB_DATA->GetSubfaceNeighborBasisValue();
  static const auto& subbdry_ptr = DENEB_DATA->GetSubbdryPtr();
  static const auto& subbdry_ind = DENEB_DATA->GetSubbdryInd();
  for (int icell = 0; icell < num_cells; icell++) {
    const int& num_points = num_cell_points[icell];

    avocado::Kernel0::f4(&solution[icell * sb], &cell_basis_value[icell][0],
                         &owner_solution[0], S_, num_bases, num_points, 1.0,
                         0.0);
    avocado::Kernel1::f42(&cell_basis_grad_value[icell][0],
                          &solution[icell * sb], &owner_solution_grad[0], D_,
                          num_points, num_bases, S_, 1.0, 0.0);
    avocado::Kernel1::f18(&auxiliary_solution_[icell * dsb],
                          &cell_basis_value[icell][0], &owner_solution_grad[0],
                          D_, S_, num_bases, num_points, 1.0, 1.0);

    ComputeComFluxJacobi(num_points, flux_owner_jacobi, flux_owner_grad_jacobi,
                         icell, owner_solution, owner_solution_grad);

    memset(&largeA[0], 0, dbb * sizeof(double));
    for (int index = subface_ptr[icell]; index < subface_ptr[icell + 1];
         index++) {
      const int& num_subface_points = subface_num_points[index];
      const int& neighbor_cell = subface_neighbor_cells[index];
      const double& sign = subface_sign[index];

      avocado::Kernel1::f61(
          subface_owner_coefficients[index], subface_owner_basis_value[index],
          &largeA[0], D_, num_bases, num_subface_points, num_bases, -sign, 1.0);
      avocado::Kernel1::f47(
          subface_owner_coefficients[index], &cell_basis_value[icell][0],
          &alpha[0], D_, num_subface_points, num_bases, num_points, sign, 0.0);
      avocado::Kernel1::f50(&alpha[0], subface_neighbor_basis_value[index],
                            &owner_b[0], D_, num_points, num_subface_points,
                            num_bases, 0.5, 0.0);

      for (int ipoint = 0; ipoint < num_points; ipoint++)
        gemmAB(1.0, &flux_owner_grad_jacobi[ipoint * DDSS_],
               &owner_b[ipoint * db], 0.0, &flux_derivative[ipoint * dssb],
               DSS_, D_, num_bases);
      avocado::Kernel1::f59(&flux_derivative[0], &cell_coefficients[icell][0],
                            &block[0], S_, sb, num_points * D_, num_bases, -1.0,
                            0.0);
      MatSetValuesBlocked(sysmat, 1, &mat_index[icell], 1,
                          &mat_index[neighbor_cell], &block[0], ADD_VALUES);
    }
    cblas_dcopy(db * num_points, &cell_basis_grad_value[icell][0], 1,
                &owner_a[0], 1);
    avocado::Kernel1::f46(&largeA[0], &cell_basis_value[icell][0], &owner_a[0],
                          D_, num_bases, num_bases, num_points, 0.5, 1.0);

    memset(&largeB[0], 0, num_points * dssb * sizeof(double));
    for (int index = subbdry_ptr[icell]; index < subbdry_ptr[icell + 1];
         index++) {
      const int& ibdry = subbdry_ind[index];
      const int& num_subbdry_points = num_bdry_points[ibdry];
      const double* solution_jacobi = &bdry_solution_jacobi[ibdry * ssp];

      avocado::Kernel1::f47(
          &bdry_owner_coefficients[ibdry][0], &cell_basis_value[icell][0],
          &alpha[0], D_, num_subbdry_points, num_bases, num_points, 1.0, 0.0);
      memset(&largeA[0], 0, num_subbdry_points * ssb * sizeof(double));
      for (int ipoint = 0; ipoint < num_subbdry_points; ipoint++)
        cblas_dger(CBLAS_LAYOUT::CblasRowMajor, SS_, num_bases, 1.0,
                   &solution_jacobi[ipoint * SS_], 1,
                   &bdry_owner_basis_value[ibdry][ipoint * num_bases], 1,
                   &largeA[ipoint * ssb], num_bases);
      avocado::Kernel1::f21(&alpha[0], &largeA[0], &largeB[0], num_points, D_,
                            num_subbdry_points, ssb, 1.0, 1.0);
    }

    for (int ipoint = 0; ipoint < num_points; ipoint++) {
      avocado::Kernel2::f13(
          &flux_owner_grad_jacobi[ipoint * DDSS_], &largeB[ipoint * dssb],
          &flux_derivative[ipoint * dssb], DS_, D_, S_, sb, 1.0, 0.0);
      gemmAB(1.0, &flux_owner_grad_jacobi[ipoint * DDSS_],
             &owner_a[ipoint * db], 1.0, &flux_derivative[ipoint * dssb], DSS_,
             D_, num_bases);
      cblas_dger(CBLAS_LAYOUT::CblasRowMajor, DSS_, num_bases, 1.0,
                 &flux_owner_jacobi[ipoint * DSS_], 1,
                 &cell_basis_value[icell][ipoint * num_bases], 1,
                 &flux_derivative[ipoint * dssb], num_bases);
    }
    avocado::Kernel1::f59(&flux_derivative[0], &cell_coefficients[icell][0],
                          &block[0], S_, sb, num_points * D_, num_bases, -1.0,
                          0.0);
    MatSetValuesBlocked(sysmat, 1, &mat_index[icell], 1, &mat_index[icell],
                        &block[0], ADD_VALUES);
  }

  MatAssemblyBegin(sysmat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(sysmat, MAT_FINAL_ASSEMBLY);
}
void EquationEquilibriumNS2D::ComputeComFlux(
    const int num_points, std::vector<double>& flux, const int icell,
    const std::vector<double>& owner_u,
    const std::vector<double>& owner_div_u) {
  int ind = 0;
  double AVcoeff = 0.0;
  for (int ipoint = 0; ipoint < num_points; ipoint++) {
    GET_SOLUTION_PS(, owner_u);

    GET_SOLUTION_GRAD_PDS(, owner_div_u);

    COMPUTE_FLUX_BASICS(, double);

    COMPUTE_VOLUME_FLUX_PDS(flux, icell, ipoint);
  }
}
void EquationEquilibriumNS2D::ComputeComFluxJacobi(
    const int num_points, std::vector<double>& flux_jacobi,
    std::vector<double>& flux_grad_jacobi, const int icell,
    const std::vector<double>& owner_u,
    const std::vector<double>& owner_div_u) {
  // flux_jacobi(ds1s2) = F(ds1) over U(s2)
  // flux_grad_jacobi(d1s1 s2d2) = F(d1s1) over gradU(d2s2)
  int ind = 0;
  double AVcoeff = 0.0;
  static std::vector<Dual<S_>> flux1(DS_);
  static std::vector<Dual<DS_>> flux2(DS_);
  for (int ipoint = 0; ipoint < num_points; ipoint++) {
    {
      // ps
      ind = S_ * ipoint;
      const Dual<S_> d(owner_u[ind++], 0);
      const Dual<S_> du(owner_u[ind++], 1);
      const Dual<S_> dv(owner_u[ind++], 2);
      const Dual<S_> dE(owner_u[ind], 3);

      GET_SOLUTION_GRAD_PDS(, owner_div_u);

      COMPUTE_FLUX_BASICS(, Dual<S_>);

      COMPUTE_VOLUME_FLUX_PDS(flux1, icell, 0);

      // pdss
      ind = DSS_ * ipoint;
      for (int ds = 0; ds < DS_; ds++)
        for (int istate = 0; istate < S_; istate++)
          flux_jacobi[ind++] = flux1[ds].df[istate];
    }
    {
      GET_SOLUTION_PS(, owner_u);

      // pds over psd
      ind = DS_ * ipoint;
      const Dual<DS_> dx(owner_div_u[ind++], 0);
      const Dual<DS_> dux(owner_div_u[ind++], 2);
      const Dual<DS_> dvx(owner_div_u[ind++], 4);
      const Dual<DS_> dEx(owner_div_u[ind++], 6);
      const Dual<DS_> dy(owner_div_u[ind++], 1);
      const Dual<DS_> duy(owner_div_u[ind++], 3);
      const Dual<DS_> dvy(owner_div_u[ind++], 5);
      const Dual<DS_> dEy(owner_div_u[ind], 7);

      COMPUTE_FLUX_BASICS(, double);

      COMPUTE_VOLUME_FLUX_PDS(flux2, icell, 0);

      // pdssd
      ind = DDSS_ * ipoint;
      for (int ds = 0; ds < DS_; ds++)
        for (int sd = 0; sd < DS_; sd++)
          flux_grad_jacobi[ind++] = flux2[ds].df[sd];
    }
  }
}
void EquationEquilibriumNS2D::ComputeNumFluxLLF(const int num_points,
                                                std::vector<double>& flux,
                                                FACE_INPUTS) {
  int ind = 0;
  double AVcoeff_o = 0.0;
  double AVcoeff_n = 0.0;
  for (int ipoint = 0; ipoint < num_points; ipoint++) {
    GET_NORMAL_PD(normal);

    GET_SOLUTION_PS(_o, owner_u);
    GET_SOLUTION_PS(_n, neighbor_u);

    GET_SOLUTION_GRAD_PDS(_o, owner_div_u);
    GET_SOLUTION_GRAD_PDS(_n, neighbor_div_u);

    COMPUTE_FLUX_BASICS(_o, double);
    const double V_o = u_o * nx + v_o * ny;
    const double a_o =
        IDEA_Wrapper(IDEA_ROUTINES(DEA), d_o * rho_ref_, e_o * e_ref_) / a_ref_;

    COMPUTE_FLUX_BASICS(_n, double);
    const double V_n = u_n * nx + v_n * ny;
    const double a_n =
        IDEA_Wrapper(IDEA_ROUTINES(DEA), d_n * rho_ref_, e_n * e_ref_) / a_ref_;
    const double r_max = std::max(std::abs(V_o) + a_o, std::abs(V_n) + a_n);

    const auto diff1 = r_max * (d_n - d_o);
    const auto diff2 = r_max * (du_n - du_o);
    const auto diff3 = r_max * (dv_n - dv_o);
    const auto diff4 = r_max * (dE_n - dE_o);

    COMPUTE_NUMERICAL_FLUX_PDS(flux, owner_cell, neighbor_cell, ipoint);
  }
}
void EquationEquilibriumNS2D::ComputeNumFluxJacobiLLF(const int num_points,
                                                      FACE_JACOBI_OUTPUTS,
                                                      FACE_INPUTS) {
  // flux_jacobi(ds1s2) = F(ds1) over U(s2)
  // flux_grad_jacobi(d1s1 s2d2) = F(d1s1) over gradU(d2s2)
  int ind = 0;
  double AVcoeff_o = 0.0;
  double AVcoeff_n = 0.0;
  static std::vector<Dual<S_>> flux1(DS_);
  static std::vector<Dual<DS_>> flux2(DS_);
  for (int ipoint = 0; ipoint < num_points; ipoint++) {
    GET_NORMAL_PD(normal);
    {
      GET_SOLUTION_PS(_n, neighbor_u);

      GET_SOLUTION_GRAD_PDS(_n, neighbor_div_u);

      COMPUTE_FLUX_BASICS(_n, double);
      const double V_n = u_n * nx + v_n * ny;
      const double a_n =
          IDEA_Wrapper(IDEA_ROUTINES(DEA), d_n * rho_ref_, e_n * e_ref_) /
          a_ref_;
      {
        // ps
        ind = S_ * ipoint;
        const Dual<S_> d_o(owner_u[ind++], 0);
        const Dual<S_> du_o(owner_u[ind++], 1);
        const Dual<S_> dv_o(owner_u[ind++], 2);
        const Dual<S_> dE_o(owner_u[ind], 3);

        GET_SOLUTION_GRAD_PDS(_o, owner_div_u);

        COMPUTE_FLUX_BASICS(_o, Dual<S_>);
        const auto V_o = u_o * nx + v_o * ny;
        const auto a_o =
            IDEA_Wrapper(IDEA_ROUTINES(DEA), d_o * rho_ref_, e_o * e_ref_) /
            a_ref_;
        const auto r_max = std::max(std::abs(V_o) + a_o, std::abs(V_n) + a_n);

        const auto diff1 = r_max * (d_n - d_o);
        const auto diff2 = r_max * (du_n - du_o);
        const auto diff3 = r_max * (dv_n - dv_o);
        const auto diff4 = r_max * (dE_n - dE_o);

        COMPUTE_NUMERICAL_FLUX_PDS(flux1, owner_cell, neighbor_cell, 0);

        // pdss
        ind = DSS_ * ipoint;
        for (int ds = 0; ds < DS_; ds++)
          for (int istate = 0; istate < S_; istate++)
            flux_owner_jacobi[ind++] = flux1[ds].df[istate];
      }
      {
        GET_SOLUTION_PS(_o, owner_u);

        // pds over psd
        ind = DS_ * ipoint;
        const Dual<DS_> dx_o(owner_div_u[ind++], 0);
        const Dual<DS_> dux_o(owner_div_u[ind++], 2);
        const Dual<DS_> dvx_o(owner_div_u[ind++], 4);
        const Dual<DS_> dEx_o(owner_div_u[ind++], 6);
        const Dual<DS_> dy_o(owner_div_u[ind++], 1);
        const Dual<DS_> duy_o(owner_div_u[ind++], 3);
        const Dual<DS_> dvy_o(owner_div_u[ind++], 5);
        const Dual<DS_> dEy_o(owner_div_u[ind], 7);

        COMPUTE_FLUX_BASICS(_o, double);
        const auto V_o = u_o * nx + v_o * ny;
        const auto a_o =
            IDEA_Wrapper(IDEA_ROUTINES(DEA), d_o * rho_ref_, e_o * e_ref_) /
            a_ref_;
        const auto r_max = std::max(std::abs(V_o) + a_o, std::abs(V_n) + a_n);

        const auto diff1 = r_max * (d_n - d_o);
        const auto diff2 = r_max * (du_n - du_o);
        const auto diff3 = r_max * (dv_n - dv_o);
        const auto diff4 = r_max * (dE_n - dE_o);

        COMPUTE_NUMERICAL_FLUX_PDS(flux2, owner_cell, neighbor_cell, 0);

        // pdssd
        ind = DDSS_ * ipoint;
        for (int ds = 0; ds < DS_; ds++)
          for (int sd = 0; sd < DS_; sd++)
            flux_owner_grad_jacobi[ind++] = flux2[ds].df[sd];
      }
    }
    {
      GET_SOLUTION_PS(_o, owner_u);

      GET_SOLUTION_GRAD_PDS(_o, owner_div_u);

      COMPUTE_FLUX_BASICS(_o, double);
      const double V_o = u_o * nx + v_o * ny;
      const double a_o =
          IDEA_Wrapper(IDEA_ROUTINES(DEA), d_o * rho_ref_, e_o * e_ref_) /
          a_ref_;
      {
        // ps
        ind = S_ * ipoint;
        const Dual<S_> d_n(neighbor_u[ind++], 0);
        const Dual<S_> du_n(neighbor_u[ind++], 1);
        const Dual<S_> dv_n(neighbor_u[ind++], 2);
        const Dual<S_> dE_n(neighbor_u[ind], 3);

        GET_SOLUTION_GRAD_PDS(_n, neighbor_div_u);

        COMPUTE_FLUX_BASICS(_n, Dual<S_>);
        const auto V_n = u_n * nx + v_n * ny;
        const auto a_n =
            IDEA_Wrapper(IDEA_ROUTINES(DEA), d_n * rho_ref_, e_n * e_ref_) /
            a_ref_;
        const auto r_max = std::max(std::abs(V_o) + a_o, std::abs(V_n) + a_n);

        const auto diff1 = r_max * (d_n - d_o);
        const auto diff2 = r_max * (du_n - du_o);
        const auto diff3 = r_max * (dv_n - dv_o);
        const auto diff4 = r_max * (dE_n - dE_o);

        COMPUTE_NUMERICAL_FLUX_PDS(flux1, owner_cell, neighbor_cell, 0);

        // pdss
        ind = DSS_ * ipoint;
        for (int ds = 0; ds < DS_; ds++)
          for (int istate = 0; istate < S_; istate++)
            flux_neighbor_jacobi[ind++] = flux1[ds].df[istate];
      }
      {
        GET_SOLUTION_PS(_n, neighbor_u);

        // pds over psd
        ind = DS_ * ipoint;
        const Dual<DS_> dx_n(neighbor_div_u[ind++], 0);
        const Dual<DS_> dux_n(neighbor_div_u[ind++], 2);
        const Dual<DS_> dvx_n(neighbor_div_u[ind++], 4);
        const Dual<DS_> dEx_n(neighbor_div_u[ind++], 6);
        const Dual<DS_> dy_n(neighbor_div_u[ind++], 1);
        const Dual<DS_> duy_n(neighbor_div_u[ind++], 3);
        const Dual<DS_> dvy_n(neighbor_div_u[ind++], 5);
        const Dual<DS_> dEy_n(neighbor_div_u[ind], 7);

        COMPUTE_FLUX_BASICS(_n, double);
        const auto V_n = u_n * nx + v_n * ny;
        const auto a_n =
            IDEA_Wrapper(IDEA_ROUTINES(DEA), d_n * rho_ref_, e_n * e_ref_) /
            a_ref_;
        const auto r_max = std::max(std::abs(V_o) + a_o, std::abs(V_n) + a_n);

        const auto diff1 = r_max * (d_n - d_o);
        const auto diff2 = r_max * (du_n - du_o);
        const auto diff3 = r_max * (dv_n - dv_o);
        const auto diff4 = r_max * (dE_n - dE_o);

        COMPUTE_NUMERICAL_FLUX_PDS(flux2, owner_cell, neighbor_cell, 0);

        // pdssd
        ind = DDSS_ * ipoint;
        for (int ds = 0; ds < DS_; ds++)
          for (int sd = 0; sd < DS_; sd++)
            flux_neighbor_grad_jacobi[ind++] = flux2[ds].df[sd];
      }
    }
  }
}

// ------------------------------- Boundary -------------------------------- //
std::shared_ptr<BoundaryEquilibriumNS2D> BoundaryEquilibriumNS2D::GetBoundary(
    const std::string& type, const int bdry_tag,
    EquationEquilibriumNS2D* equation) {
  if (!type.compare("AdiabaticWall"))
    return std::make_shared<AdiabaticWallEquilibriumNS2D>(bdry_tag, equation);
  else if (!type.compare("IsothermalWall"))
    return std::make_shared<IsothermalWallEquilibriumNS2D>(bdry_tag, equation);
  else if (!type.compare("SupersonicInflow"))
    return std::make_shared<SupersonicInflowEquilibriumNS2D>(bdry_tag,
                                                             equation);
  else if (!type.compare("BackPressure"))
    return std::make_shared<BackPressureEquilibriumNS2D>(bdry_tag, equation);
  ERROR_MESSAGE("Wrong boundary condition (no-exist):" + type + "\n");
  return nullptr;
}
// Boundary = AdiabaticWall
// Dependency: -
AdiabaticWallEquilibriumNS2D::AdiabaticWallEquilibriumNS2D(
    const int bdry_tag, EquationEquilibriumNS2D* equation)
    : BoundaryEquilibriumNS2D(bdry_tag, equation) {
  MASTER_MESSAGE("AdiabaticWall (tag=" + std::to_string(bdry_tag) + ")\n");
}
void AdiabaticWallEquilibriumNS2D::ComputeBdrySolution(
    const int num_points, std::vector<double>& bdry_u,
    std::vector<double>& bdry_div_u, const std::vector<double>& owner_u,
    const std::vector<double>& owner_div_u, const std::vector<double>& normal,
    const std::vector<double>& coords, const double& time) {
  int ind = 0;
  for (int ipoint = 0; ipoint < num_points; ipoint++) {
    GET_SOLUTION_PS(, owner_u);

    // ps
    ind = S_ * ipoint;
    bdry_u[ind++] = d;
    bdry_u[ind++] = 0.0;
    bdry_u[ind++] = 0.0;
    bdry_u[ind] = dE - 0.5 * (du * du + dv * dv) / d;
  }
  // Don't touch bdry_div_u.
}
void AdiabaticWallEquilibriumNS2D::ComputeBdryFlux(
    const int num_points, std::vector<double>& flux, FACE_INPUTS,
    const std::vector<double>& coords, const double& time) {
  int ind = 0;
  double AVcoeff = 0.0;
  for (int ipoint = 0; ipoint < num_points; ipoint++) {
    GET_SOLUTION_PS(, owner_u);

    GET_SOLUTION_GRAD_PDS(, owner_div_u);

    COMPUTE_FLUX_BASICS(, double);

    // pds
    AVcoeff = DENEB_ARTIFICIAL_VISCOSITY->GetArtificialViscosityValue(
        owner_cell, ipoint);
    flux[ind++] = -AVcoeff * dx;
    flux[ind++] = p - alpha_ * txx - AVcoeff * dux;
    flux[ind++] = -alpha_ * txy - AVcoeff * dvx;
    flux[ind++] = -AVcoeff * dEx;

    flux[ind++] = -AVcoeff * dy;
    flux[ind++] = -alpha_ * txy - AVcoeff * duy;
    flux[ind++] = p - alpha_ * tyy - AVcoeff * dvy;
    flux[ind] = -AVcoeff * dEy;
  }
}
void AdiabaticWallEquilibriumNS2D::ComputeBdrySolutionJacobi(
    const int num_points, double* bdry_u_jacobi,
    const std::vector<double>& owner_u, const std::vector<double>& owner_div_u,
    const std::vector<double>& normal, const std::vector<double>& coords,
    const double& time) {
  int ind = 0;
  memset(bdry_u_jacobi, 0, num_points * SS_ * sizeof(double));
  for (int ipoint = 0; ipoint < num_points; ipoint++) {
    GET_SOLUTION_PS(, owner_u);

    const double d_inv = 1.0 / d;
    const double u = du * d_inv;
    const double v = dv * d_inv;

    ind = SS_ * ipoint;
    bdry_u_jacobi[ind] = 1.0;
    ind += (S_ * (S_ - 1));
    bdry_u_jacobi[ind++] = 0.5 * (u * u + v * v);
    bdry_u_jacobi[ind++] = -u;
    bdry_u_jacobi[ind++] = -v;
    bdry_u_jacobi[ind] = 1.0;
  }
}
void AdiabaticWallEquilibriumNS2D::ComputeBdryFluxJacobi(
    const int num_points, std::vector<double>& flux_jacobi,
    std::vector<double>& flux_grad_jacobi, FACE_INPUTS,
    const std::vector<double>& coords, const double& time) {
  // flux_jacobi(ds1s2) = F(ds1) over U(s2)
  // flux_grad_jacobi(d1s1 s2d2) = F(d1s1) over gradU(d2s2)
  int ind = 0;
  double AVcoeff = 0.0;
  std::vector<Dual<S_>> flux1(DS_);
  std::vector<Dual<DS_>> flux2(DS_);
  for (int ipoint = 0; ipoint < num_points; ipoint++) {
    AVcoeff = DENEB_ARTIFICIAL_VISCOSITY->GetArtificialViscosityValue(
        owner_cell, ipoint);
    {
      // ps
      ind = S_ * ipoint;
      const Dual<S_> d(owner_u[ind++], 0);
      const Dual<S_> du(owner_u[ind++], 1);
      const Dual<S_> dv(owner_u[ind++], 2);
      const Dual<S_> dE(owner_u[ind], 3);

      GET_SOLUTION_GRAD_PDS(, owner_div_u);

      COMPUTE_FLUX_BASICS(, Dual<S_>);

      // pds
      flux1[0] = -AVcoeff * dx;
      flux1[1] = p - alpha_ * txx - AVcoeff * dux;
      flux1[2] = -alpha_ * txy - AVcoeff * dvx;
      flux1[3] = -AVcoeff * dEx;
      flux1[4] = -AVcoeff * dy;
      flux1[5] = -alpha_ * txy - AVcoeff * duy;
      flux1[6] = p - alpha_ * tyy - AVcoeff * dvy;
      flux1[7] = -AVcoeff * dEy;

      // pdss
      ind = DSS_ * ipoint;
      for (int ds = 0; ds < DS_; ds++)
        for (int istate = 0; istate < S_; istate++)
          flux_jacobi[ind++] = flux1[ds].df[istate];
    }
    {
      GET_SOLUTION_PS(, owner_u);

      // pds over psd
      ind = DS_ * ipoint;
      const Dual<DS_> dx(owner_div_u[ind++], 0);
      const Dual<DS_> dux(owner_div_u[ind++], 2);
      const Dual<DS_> dvx(owner_div_u[ind++], 4);
      const Dual<DS_> dEx(owner_div_u[ind++], 6);
      const Dual<DS_> dy(owner_div_u[ind++], 1);
      const Dual<DS_> duy(owner_div_u[ind++], 3);
      const Dual<DS_> dvy(owner_div_u[ind++], 5);
      const Dual<DS_> dEy(owner_div_u[ind], 7);

      COMPUTE_FLUX_BASICS(, double);

      // pds
      flux2[0] = -AVcoeff * dx;
      flux2[1] = p - alpha_ * txx - AVcoeff * dux;
      flux2[2] = -alpha_ * txy - AVcoeff * dvx;
      flux2[3] = -AVcoeff * dEx;
      flux2[4] = -AVcoeff * dy;
      flux2[5] = -alpha_ * txy - AVcoeff * duy;
      flux2[6] = p - alpha_ * tyy - AVcoeff * dvy;
      flux2[7] = -AVcoeff * dEy;

      // pdssd
      ind = DDSS_ * ipoint;
      for (int ds = 0; ds < DS_; ds++)
        for (int sd = 0; sd < DS_; sd++)
          flux_grad_jacobi[ind++] = flux2[ds].df[sd];
    }
  }
}
// Boundary = IsothermalWall
// Dependency: BdryInput(num)
// BdryInput(num) = Twall
IsothermalWallEquilibriumNS2D::IsothermalWallEquilibriumNS2D(
    const int bdry_tag, EquationEquilibriumNS2D* equation)
    : BoundaryEquilibriumNS2D(bdry_tag, equation) {
  auto& config = AVOCADO_CONFIG;
  Twall_ = std::stod(config->GetConfigValue(BDRY_INPUT_I(bdry_tag, 0)));

  MASTER_MESSAGE("IsothermalWall (tag=" + std::to_string(bdry_tag) +
                 ")\n\tWall temperature = " + std::to_string(Twall_) + " K\n");
}
void IsothermalWallEquilibriumNS2D::ComputeBdrySolution(
    const int num_points, std::vector<double>& bdry_u,
    std::vector<double>& bdry_div_u, const std::vector<double>& owner_u,
    const std::vector<double>& owner_div_u, const std::vector<double>& normal,
    const std::vector<double>& coords, const double& time) {
  int ind = 0;
  for (int ipoint = 0; ipoint < num_points; ipoint++) {
    GET_SOLUTION_PS(, owner_u);

    const double e = (dE - 0.5 * (du * du + dv * dv) / d) / d;
    const double dim_p = IDEA_DEP(d * rho_ref_, e * e_ref_);
    const double ew = IDEA_PTE(dim_p, Twall_) / e_ref_;
    const double dw = IDEA_PTD(dim_p, Twall_) / rho_ref_;

    // ps
    ind = S_ * ipoint;
    bdry_u[ind++] = dw;
    bdry_u[ind++] = 0.0;
    bdry_u[ind++] = 0.0;
    bdry_u[ind] = dw * ew;
  }
  // Don't touch bdry_div_u.
}
void IsothermalWallEquilibriumNS2D::ComputeBdryFlux(
    const int num_points, std::vector<double>& flux, FACE_INPUTS,
    const std::vector<double>& coords, const double& time) {
  int ind = 0;
  double AVcoeff = 0.0;
  for (int ipoint = 0; ipoint < num_points; ipoint++) {
    GET_SOLUTION_PS(, owner_u);

    GET_SOLUTION_GRAD_PDS(, owner_div_u);

    const auto d_inv = 1.0 / d;
    const auto u = du * d_inv;
    const auto v = dv * d_inv;
    const auto ux = (dux - dx * u) * d_inv;
    const auto uy = (duy - dy * u) * d_inv;
    const auto vx = (dvx - dx * v) * d_inv;
    const auto vy = (dvy - dy * v) * d_inv;
    const auto e = (dE - 0.5 * (du * u + dv * v)) * d_inv;
    const auto dim_p = IDEA_DEP(d * rho_ref_, e * e_ref_);
    const double p = dim_p / p_ref_;
    const double muw = IDEA_PTV(dim_p, Twall_) / mu_ref_;
    const double kw = IDEA_PTK(dim_p, Twall_) / k_ref_;
    const auto txx = c23_ * muw * (2.0 * ux - vy);
    const auto txy = muw * (uy + vx);
    const auto tyy = c23_ * muw * (2.0 * vy - ux);
    const auto ex = (dEx - dx * dE * d_inv) * d_inv - u * ux - v * vx;
    const auto ey = (dEy - dy * dE * d_inv) * d_inv - u * uy - v * vy;
    double T_grad[2];
    IDEA_DET_Grad(T_grad, d * rho_ref_, e * e_ref_);
    T_grad[0] = T_grad[0] * (rho_ref_ / T_ref_);
    T_grad[1] = T_grad[1] * (e_ref_ / T_ref_);
    const auto Qx = -kw * (T_grad[0] * dx + T_grad[1] * ex);
    const auto Qy = -kw * (T_grad[0] * dy + T_grad[1] * ey);

    // pds
    ind = DS_ * ipoint;
    AVcoeff = DENEB_ARTIFICIAL_VISCOSITY->GetArtificialViscosityValue(
        owner_cell, ipoint);
    flux[ind++] = -AVcoeff * dx;
    flux[ind++] = p - alpha_ * txx - AVcoeff * dux;
    flux[ind++] = -alpha_ * txy - AVcoeff * dvx;
    flux[ind++] = beta_ * Qx - AVcoeff * dEx;

    flux[ind++] = -AVcoeff * dy;
    flux[ind++] = -alpha_ * txy - AVcoeff * duy;
    flux[ind++] = p - alpha_ * tyy - AVcoeff * dvy;
    flux[ind] = beta_ * Qy - AVcoeff * dEy;
  }
}
void IsothermalWallEquilibriumNS2D::ComputeBdrySolutionJacobi(
    const int num_points, double* bdry_u_jacobi,
    const std::vector<double>& owner_u, const std::vector<double>& owner_div_u,
    const std::vector<double>& normal, const std::vector<double>& coords,
    const double& time) {
  ERROR_MESSAGE("Not supported!\n");
}
void IsothermalWallEquilibriumNS2D::ComputeBdryFluxJacobi(
    const int num_points, std::vector<double>& flux_jacobi,
    std::vector<double>& flux_grad_jacobi, FACE_INPUTS,
    const std::vector<double>& coords, const double& time) {
  ERROR_MESSAGE("Not supported!\n");
}
// Boundary = SupersonicInflow
// Dependency: BdryInput(num)
// BdryInput(num) = Ma, AOA, Pinflow, Tinflow
SupersonicInflowEquilibriumNS2D::SupersonicInflowEquilibriumNS2D(
    const int bdry_tag, EquationEquilibriumNS2D* equation)
    : BoundaryEquilibriumNS2D(bdry_tag, equation) {
  auto& config = AVOCADO_CONFIG;
  Ma_ = std::stod(config->GetConfigValue(BDRY_INPUT_I(bdry_tag, 0)));
  AoA_ = std::stod(config->GetConfigValue(BDRY_INPUT_I(bdry_tag, 1)));
  p_inflow_ = std::stod(config->GetConfigValue(BDRY_INPUT_I(bdry_tag, 2)));
  T_inflow_ = std::stod(config->GetConfigValue(BDRY_INPUT_I(bdry_tag, 3)));

  u_ = Ma_ * std::cos(AoA_ * pi_ / 180.0);
  v_ = Ma_ * std::sin(AoA_ * pi_ / 180.0);
  p_ = p_inflow_ / p_ref_;
  const double e = IDEA_PTE(p_inflow_, T_inflow_) / e_ref_;
  d_ = IDEA_PTD(p_inflow_, T_inflow_) / rho_ref_;
  du_ = d_ * u_;
  dv_ = d_ * v_;
  dE_ = d_ * (e + 0.5 * (u_ * u_ + v_ * v_));

  MASTER_MESSAGE(
      "Supersonic (tag=" + std::to_string(bdry_tag) +
      ")\n\tMa = " + std::to_string(Ma_) + "\n\tAOA = " + std::to_string(AoA_) +
      "\n\tInflow pressure = " + std::to_string(p_inflow_) +
      "\n\tInflow temperature = " + std::to_string(T_inflow_) + "\n");
}
void SupersonicInflowEquilibriumNS2D::ComputeBdrySolution(
    const int num_points, std::vector<double>& bdry_u,
    std::vector<double>& bdry_div_u, const std::vector<double>& owner_u,
    const std::vector<double>& owner_div_u, const std::vector<double>& normal,
    const std::vector<double>& coords, const double& time) {
  int ind = 0;
  for (int ipoint = 0; ipoint < num_points; ipoint++) {
    // ps
    ind = S_ * ipoint;
    bdry_u[ind++] = d_;
    bdry_u[ind++] = du_;
    bdry_u[ind++] = dv_;
    bdry_u[ind] = dE_;
  }
  // Don't touch bdry_div_u.
}
void SupersonicInflowEquilibriumNS2D::ComputeBdryFlux(
    const int num_points, std::vector<double>& flux, FACE_INPUTS,
    const std::vector<double>& coords, const double& time) {
  int ind = 0;
  for (int ipoint = 0; ipoint < num_points; ipoint++) {
    // pds
    ind = DS_ * ipoint;
    flux[ind++] = du_;
    flux[ind++] = du_ * u_ + p_;
    flux[ind++] = du_ * v_;
    flux[ind++] = (dE_ + p_) * u_;

    flux[ind++] = dv_;
    flux[ind++] = dv_ * u_;
    flux[ind++] = dv_ * v_ + p_;
    flux[ind] = (dE_ + p_) * v_;
  }
}
void SupersonicInflowEquilibriumNS2D::ComputeBdrySolutionJacobi(
    const int num_points, double* bdry_u_jacobi,
    const std::vector<double>& owner_u, const std::vector<double>& owner_div_u,
    const std::vector<double>& normal, const std::vector<double>& coords,
    const double& time) {
  ERROR_MESSAGE("Not supported!\n");
}
void SupersonicInflowEquilibriumNS2D::ComputeBdryFluxJacobi(
    const int num_points, std::vector<double>& flux_jacobi,
    std::vector<double>& flux_grad_jacobi, FACE_INPUTS,
    const std::vector<double>& coords, const double& time) {
  ERROR_MESSAGE("Not supported!\n");
}
// Boundary = BackPressure
// for subsonic and supersonic
// Dependency: BdryInput(num)
// BdryInput(num) = back pressure
BackPressureEquilibriumNS2D::BackPressureEquilibriumNS2D(
    const int bdry_tag, EquationEquilibriumNS2D* equation)
    : BoundaryEquilibriumNS2D(bdry_tag, equation) {
  auto& config = AVOCADO_CONFIG;
  p_back_ = std::stod(config->GetConfigValue(BDRY_INPUT_I(bdry_tag, 0)));

  MASTER_MESSAGE("BackPressure (tag=" + std::to_string(bdry_tag) +
                 ")\n\tBack pressure = " + std::to_string(p_back_) + " Pa\n");
}
void BackPressureEquilibriumNS2D::ComputeBdrySolution(
    const int num_points, std::vector<double>& bdry_u,
    std::vector<double>& bdry_div_u, const std::vector<double>& owner_u,
    const std::vector<double>& owner_div_u, const std::vector<double>& normal,
    const std::vector<double>& coords, const double& time) {
  int ind = 0;
  for (int ipoint = 0; ipoint < num_points; ipoint++) {
    GET_NORMAL_PD(normal);

    GET_SOLUTION_PS(, owner_u);

    const double d_inv = 1.0 / d;
    const double u = du * d_inv;
    const double v = dv * d_inv;
    const double absV = u * u + v * v;
    const double e = dE * d_inv - 0.5 * absV;
    const double a = IDEA_DEA(d * rho_ref_, e * e_ref_) / a_ref_;
    const double V = u * nx + v * ny;
    const double local_M = V / a;

    // ps
    ind = S_ * ipoint;
    if (local_M >= 1.0) {  // supersonic
      bdry_u[ind++] = d;
      bdry_u[ind++] = du;
      bdry_u[ind++] = dv;
      bdry_u[ind] = dE;
    } else {  // subsonic
      const double dim_p = IDEA_DEP(d * rho_ref_, e * e_ref_);
      const double dim_T = IDEA_DPT(d * rho_ref_, dim_p);
      const double d_back = IDEA_PTD(p_back_, dim_T) / rho_ref_;
      const double e_back = IDEA_DPE(d_back * rho_ref_, p_back_) / e_ref_;

      bdry_u[ind++] = d_back;
      bdry_u[ind++] = d_back * u;
      bdry_u[ind++] = d_back * v;
      bdry_u[ind] = d_back * (e_back + 0.5 * absV);
    }
  }
  memset(&bdry_div_u[0], 0, DS_ * num_points * sizeof(double));
}
void BackPressureEquilibriumNS2D::ComputeBdryFlux(
    const int num_points, std::vector<double>& flux, FACE_INPUTS,
    const std::vector<double>& coords, const double& time) {
  equation_->ComputeNumFlux(num_points, flux, owner_cell, neighbor_cell,
                            owner_u, owner_div_u, neighbor_u, neighbor_div_u,
                            normal);
}
void BackPressureEquilibriumNS2D::ComputeBdrySolutionJacobi(
    const int num_points, double* bdry_u_jacobi,
    const std::vector<double>& owner_u, const std::vector<double>& owner_div_u,
    const std::vector<double>& normal, const std::vector<double>& coords,
    const double& time) {
  ERROR_MESSAGE("Not supported!\n");
}
void BackPressureEquilibriumNS2D::ComputeBdryFluxJacobi(
    const int num_points, std::vector<double>& flux_jacobi,
    std::vector<double>& flux_grad_jacobi, FACE_INPUTS,
    const std::vector<double>& coords, const double& time) {
  ERROR_MESSAGE("Not supported!\n");
}
// -------------------------------- Problem -------------------------------- //
std::shared_ptr<ProblemEquilibriumNS2D> ProblemEquilibriumNS2D::GetProblem(
    const std::string& name) {
  if (!name.compare("Constant"))
    return std::make_shared<ConstantEquilibriumNS2D>();
  else if (!name.compare("FreeStream"))
    return std::make_shared<FreeStreamEquilibriumNS2D>();
  else if (!name.compare("DoubleSine")) return std::make_shared<
        DoubleSineEquilibriumNS2D>();
  ERROR_MESSAGE("Wrong problem (no-exist):" + name + "\n");
  return nullptr;
}
// Problem = Constant
// ProblemInput = value 1, ... , value 4
ConstantEquilibriumNS2D::ConstantEquilibriumNS2D() {
  MASTER_MESSAGE("Constant problem\n");

  auto& config = AVOCADO_CONFIG;
  std::vector<double> values(S_);
  for (int i = 0; i < S_; i++)
    values[i] = std::stod(config->GetConfigValue(PROBLEM_INPUT_I(i)));

  std::stringstream str;
  str << "\tInput:\n";
  for (int i = 0; i < S_; i++)
    str << "\t\tvalue[" << i << "] = " << values[i] << "\n";
  MASTER_MESSAGE(str.str());
  values_ = move(values);
}
void ConstantEquilibriumNS2D::Problem(const int num_points,
                                      std::vector<double>& solutions,
                                      const std::vector<double>& coord,
                                      const double time) const {
  int ind = 0;
  for (int ipoint = 0; ipoint < num_points; ipoint++)
    for (int istate = 0; istate < S_; istate++)
      solutions[ind++] = values_[istate];
}
// Problem = FreeStream
// ProblemInput = Ma, AOA, Pinf, Tinf
FreeStreamEquilibriumNS2D::FreeStreamEquilibriumNS2D(){
  MASTER_MESSAGE("FreeStream problem\n");

  auto& config = AVOCADO_CONFIG;
  Ma_ = std::stod(config->GetConfigValue(PROBLEM_INPUT_I(0)));
  AoA_ = std::stod(config->GetConfigValue(PROBLEM_INPUT_I(1)));
  p_inf_ = std::stod(config->GetConfigValue(PROBLEM_INPUT_I(2)));
  T_inf_ = std::stod(config->GetConfigValue(PROBLEM_INPUT_I(3)));

  std::stringstream str;
  str << "\tInput (fixed):\n";
  str << "\t\tMach number = " << Ma_ << "\n";
  str << "\t\tIncident angle = " << AoA_ << "\n";
  str << "\t\tFreestream pressure = " << p_inf_ << "\n";
  str << "\t\tFreestream temperature = " << T_inf_ << "\n";
  MASTER_MESSAGE(str.str());

  const double u = Ma_ * std::cos(AoA_ * pi_ / 180.0);
  const double v = Ma_ * std::sin(AoA_ * pi_ / 180.0);
  const double p = p_inf_ / p_ref_;
  const double e = IDEA_PTE(p_inf_, T_inf_) / e_ref_;
  d_ = IDEA_PTD(p_inf_, T_inf_) / rho_ref_;
  du_ = d_ * u;
  dv_ = d_ * v;
  dE_ = d_ * (e + 0.5 * (u * u + v * v));
}
void FreeStreamEquilibriumNS2D::Problem(const int num_points,
                                        std::vector<double>& solutions,
                                        const std::vector<double>& coord,
                                        const double time) const {
  int ind = 0;
  for (int ipoint = 0; ipoint < num_points; ipoint++) {
    solutions[ind++] = d_;
    solutions[ind++] = du_;
    solutions[ind++] = dv_;
    solutions[ind++] = dE_;
  }
}
// Problem = DoubleSine
// ProblemInput = -
DoubleSineEquilibriumNS2D::DoubleSineEquilibriumNS2D()
    : velocity_({0.5, 0.5}), wave_number_({1.0, 1.0}) {
  MASTER_MESSAGE("DoubleSine problem\n");

  std::stringstream str;
  str << "\tInput (fixed):\n";
  str << "\t\tx-velocity = " << velocity_[0] << "\n";
  str << "\t\ty-velocity = " << velocity_[1] << "\n";
  str << "\t\tx-wavenumber = " << wave_number_[0] << "\n";
  str << "\t\ty-wavenumber = " << wave_number_[1] << "\n";
  MASTER_MESSAGE(str.str());
}
void DoubleSineEquilibriumNS2D::Problem(const int num_points,
                                        std::vector<double>& solutions,
                                        const std::vector<double>& coord,
                                        const double time) const {
  int ind = 0;
  const double& u = velocity_[0];
  const double& v = velocity_[1];
  const double p = 1.0;  // pressure
  for (int ipoint = 0; ipoint < num_points; ipoint++) {
    const double& x = coord[ipoint * D_];
    const double& y = coord[ipoint * D_ + 1];
    const double rho =
        1.0 + 0.2 * std::sin(2.0 * pi_ * wave_number_[0] * (x - u * time)) *
                  std::sin(2.0 * pi_ * wave_number_[1] * (y - v * time));
    const double e = IDEA_DPE(rho * rho_ref_, p * p_ref_) / e_ref_;

    solutions[ind++] = rho;
    solutions[ind++] = rho * u;
    solutions[ind++] = rho * v;
    solutions[ind++] = rho * (e + 0.5 * (u * u + v * v));
  }
}
}  // namespace deneb