#include "deneb_equation_glmmhd2d.h"

#include <algorithm>
#include <cstring>
#include <string>
#include <unordered_set>

#include "avocado.h"
#include "deneb_config_macro.h"
#include "deneb_data.h"
#include "deneb_timescheme.h"

#define GET_NORMAL_PD(data)       \
  ind = D_ * ipoint;              \
  const double& nx = data[ind++]; \
  const double& ny = data[ind]

#define GET_SOLUTION_PS(tag, data)      \
  ind = S_ * ipoint;                    \
  const double& d##tag = data[ind++];   \
  const double& du##tag = data[ind++];  \
  const double& dv##tag = data[ind++];  \
  const double& dtE##tag = data[ind++]; \
  const double& B1##tag = data[ind++];  \
  const double& B2##tag = data[ind++];  \
  const double& psi##tag = data[ind]

#define COMPUTE_FLUX_BASICS(tag)                                      \
  const auto d_inv##tag = 1.0 / d##tag;                               \
  const auto u##tag = du##tag * d_inv##tag;                           \
  const auto v##tag = dv##tag * d_inv##tag;                           \
  const auto BB##tag = B1##tag * B1##tag + B2##tag * B2##tag;         \
  const auto VB##tag = u##tag * B1##tag + v##tag * B2##tag;           \
  const auto p##tag =                                                 \
      rm1_ * (dtE##tag - 0.5 * (du##tag * u##tag + dv##tag * v##tag + \
                                BB##tag + psi##tag * psi##tag));      \
  const auto tp##tag = p##tag + 0.5 * BB##tag

#define COMPUTE_VOLUME_FLUX_PDS(flux, ipoint)                                \
  ind = DS_ * ipoint;                                                        \
  flux[ind++] = du;                                                          \
  flux[ind++] = du * u - B1 * B1 + tp;                                       \
  flux[ind++] = du * v - B1 * B2;                                            \
  flux[ind++] = (dtE + tp - 0.5 * psi * psi) * u - B1 * VB + ch_ * psi * B1; \
  flux[ind++] = ch_ * psi;                                                   \
  flux[ind++] = u * B2 - v * B1;                                             \
  flux[ind++] = ch_ * B1;                                                    \
  flux[ind++] = dv;                                                          \
  flux[ind++] = dv * u - B2 * B1;                                            \
  flux[ind++] = dv * v - B2 * B2 + tp;                                       \
  flux[ind++] = (dtE + tp - 0.5 * psi * psi) * v - B2 * VB + ch_ * psi * B2; \
  flux[ind++] = v * B1 - u * B2;                                             \
  flux[ind++] = ch_ * psi;                                                   \
  flux[ind] = ch_ * B2

#define COMPUTE_VOLUME_SOURCE_PS(source, ipoint)                \
  ind = S_ * ipoint;                                            \
  source[ind++] = 0.0;                                          \
  source[ind++] = -divB * B1;                                   \
  source[ind++] = -divB * B2;                                   \
  source[ind++] = -divB * VB - divpsiV * psi - 2.0 * psi * psi; \
  source[ind++] = -divB * u;                                    \
  source[ind++] = -divB * v;                                    \
  source[ind] = -divpsiV - 2.0 * psi

#define COMPUTE_NUMERICAL_FLUX_PDS(flux, ipoint)                               \
  ind = DS_ * ipoint;                                                          \
  flux[ind++] = 0.5 * (du_o + du_n - nx * diff1);                              \
  flux[ind++] = 0.5 * (du_o * u_o - B1_o * B1_o + du_n * u_n - B1_n * B1_n -   \
                       nx * diff2 + tp_o + tp_n);                              \
  flux[ind++] = 0.5 * (du_o * v_o - B1_o * B2_o + du_n * v_n - B1_n * B2_n -   \
                       nx * diff3);                                            \
  flux[ind++] =                                                                \
      0.5 * ((dtE_o + tp_o - 0.5 * psi_o * psi_o) * u_o - B1_o * VB_o +        \
             ch_ * psi_o * B1_o + (dtE_n + tp_n - 0.5 * psi_n * psi_n) * u_n - \
             B1_n * VB_n + ch_ * psi_n * B1_n - nx * diff4);                   \
  flux[ind++] = 0.5 * (ch_ * psi_o + ch_ * psi_n - nx * diff5);                \
  flux[ind++] =                                                                \
      0.5 * (u_o * B2_o - v_o * B1_o + u_n * B2_n - v_n * B1_n - nx * diff6);  \
  flux[ind++] = 0.5 * (ch_ * B1_o + ch_ * B1_n - nx * diff7);                  \
  flux[ind++] = 0.5 * (dv_o + dv_n - ny * diff1);                              \
  flux[ind++] = 0.5 * (dv_o * u_o - B2_o * B1_o + dv_n * u_n - B2_n * B1_n -   \
                       ny * diff2);                                            \
  flux[ind++] = 0.5 * (dv_o * v_o - B2_o * B2_o + dv_n * v_n - B2_n * B2_n -   \
                       ny * diff3 + tp_o + tp_n);                              \
  flux[ind++] =                                                                \
      0.5 * ((dtE_o + tp_o - 0.5 * psi_o * psi_o) * v_o - B2_o * VB_o +        \
             ch_ * psi_o * B2_o + (dtE_n + tp_n - 0.5 * psi_n * psi_n) * v_n - \
             B2_n * VB_n + ch_ * psi_n * B2_n - ny * diff4);                   \
  flux[ind++] =                                                                \
      0.5 * (v_o * B1_o - u_o * B2_o + v_n * B1_n - u_n * B2_n - ny * diff5);  \
  flux[ind++] = 0.5 * (ch_ * psi_o + ch_ * psi_n - ny * diff6);                \
  flux[ind] = 0.5 * (ch_ * B2_o + ch_ * B2_n - ny * diff7)

namespace deneb {
// ------------------------------- Constants -------------------------------
// //
ConstantsGLMMHD2D::ConstantsGLMMHD2D()
    : r_(std::stod(AVOCADO_CONFIG->GetConfigValue(HEAT_CAPACITY_RATIO))),
      r_inv_(1.0 / r_),
      rm1_(r_ - 1.0),
      rm1_inv_(1.0 / rm1_) {}
// compute p
double ConstantsGLMMHD2D::ComputePressure(const double* sol) const {
  const double BB = avocado::VecInnerProd(D_, &sol[4], &sol[4]);
  return rm1_ * (sol[3] - 0.5 * BB - 0.5 * sol[6] * sol[6] -
                 0.5 / sol[0] * avocado::VecInnerProd(D_, &sol[1], &sol[1]));
}
// compute dE
double ConstantsGLMMHD2D::ComputeTotalEnergy(const double* sol) const {
  const double BB = avocado::VecInnerProd(D_, &sol[4], &sol[4]);
  return sol[3] - 0.5 * BB - 0.5 * sol[6] * sol[6];
}
// ------------------------------- Equation ---------------------------- //
EquationGLMMHD2D::EquationGLMMHD2D()
    : ConstantsGLMMHD2D(), Equation(D_, S_, true) {
  MASTER_MESSAGE(avocado::GetTitle("EquationGLMMHD2D"));
  MASTER_MESSAGE("Dimension = " + std::to_string(D_) + "\n");
  MASTER_MESSAGE("Number of state variables = " + std::to_string(S_) + "\n");
  MASTER_MESSAGE(
      "Source term = " + std::string(source_term_ ? "true" : "false") + "\n");

  auto& config = AVOCADO_CONFIG;
  problem_ = ProblemGLMMHD2D::GetProblem(config->GetConfigValue(PROBLEM));
  const std::string& numflux = config->GetConfigValue(CONVECTIVE_FLUX);
  if (!numflux.compare("LLF")) {
    ASSIGN_FLUX(EquationGLMMHD2D, LLF);
  } else
    ERROR_MESSAGE("Wrong numerical flux (no-exist):" + numflux + "\n");
  MASTER_MESSAGE("Problem: " + config->GetConfigValue(PROBLEM) + "\n");
  MASTER_MESSAGE("Convective flux: " + numflux + "\n");

  ch_ = 0.0;
}
EquationGLMMHD2D::~EquationGLMMHD2D() {
  problem_.reset();
  boundaries_.clear();
  boundary_registry_.clear();
}

void EquationGLMMHD2D::RegistBoundary(const std::vector<int>& bdry_tag) {
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
          BoundaryGLMMHD2D::GetBoundary(bdry_type, tag, this);
  }
}
void EquationGLMMHD2D::BuildData(void) {
  MASTER_MESSAGE(avocado::GetTitle("EquationGLMMHD2D::BuildData()"));

  const int& num_bdries = DENEB_DATA->GetNumBdries();
  boundaries_.resize(num_bdries, nullptr);
  {
    const std::vector<int>& bdry_tag = DENEB_DATA->GetBdryTag();
    RegistBoundary(bdry_tag);
    for (int ibdry = 0; ibdry < num_bdries; ibdry++)
      boundaries_[ibdry] = boundary_registry_[bdry_tag[ibdry]];
  }

  const int& num_bases = DENEB_DATA->GetNumBases();

  int max_num_points = 0;
  int max_num_cell_points = 0;
  int max_num_face_points = 0;
  int max_num_bdry_points = 0;
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

  if (DENEB_TIMESCHEME->IsImplicit()) {
    ERROR_MESSAGE("Not supported!\n");
  }

  const int& num_cells = DENEB_DATA->GetNumCells();
  auxiliary_solution_.resize(num_cells * D_ * S_ * num_bases);

  const int& num_outer_cells = DENEB_DATA->GetNumOuterCells();
  outer_solution_.resize(
      std::max((num_outer_cells - num_cells) * S_ * num_bases, 1));
  communicate_ = std::make_shared<avocado::Communicate>(
      S_ * num_bases, DENEB_DATA->GetOuterSendCellList(),
      DENEB_DATA->GetOuterRecvCellList());
  pressure_fix_values_.resize(2);

  cell_variable_names_ = {"rho", "rhoU", "rhoV", "rhoTE", "B1",
                          "B2",  "psi",  "rhoE", "P",     "divB"};
  face_variable_names_ = {"rho", "rhoU", "rhoV", "rhoTE", "B1",
                          "B2",  "psi",  "rhoE", "P",     "divB"};

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
void EquationGLMMHD2D::GetCellPostSolution(
    const int icell, const int num_points, const std::vector<double>& solution,
    const std::vector<double>& solution_grad,
    std::vector<double>& post_solution) {
  int ind = 0;
  for (int ipoint = 0; ipoint < num_points; ipoint++) {
    for (int istate = 0; istate < S_; istate++)
      post_solution[ind++] = solution[ipoint * S_ + istate];
    post_solution[ind++] = ComputeTotalEnergy(&solution[ipoint * S_]);
    post_solution[ind++] = ComputePressure(&solution[ipoint * S_]);
    post_solution[ind++] =
        solution_grad[ipoint * DS_ + 4] + solution_grad[ipoint * DS_ + 12];
  }
}
void EquationGLMMHD2D::GetFacePostSolution(
    const int num_points, const std::vector<double>& solution,
    const std::vector<double>& solution_grad, const std::vector<double>& normal,
    std::vector<double>& post_solution) {
  int ind = 0;
  for (int ipoint = 0; ipoint < num_points; ipoint++) {
    for (int istate = 0; istate < S_; istate++)
      post_solution[ind++] = solution[ipoint * S_ + istate];
    post_solution[ind++] = ComputeTotalEnergy(&solution[ipoint * S_]);
    post_solution[ind++] = ComputePressure(&solution[ipoint * S_]);
    post_solution[ind++] =
        solution_grad[ipoint * DS_ + 4] + solution_grad[ipoint * DS_ + 12];
  }
}
void EquationGLMMHD2D::ComputeInitialSolution(double* solution,
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
void EquationGLMMHD2D::ComputeLocalTimestep(
    const double* solution, std::vector<double>& local_timestep) {
  static const int& order = DENEB_DATA->GetOrder();
  static const int& num_bases = DENEB_DATA->GetNumBases();
  static const int sb = S_ * num_bases;
  static std::vector<double> owner_solution(S_);

  static const int& num_cells = DENEB_DATA->GetNumCells();
  static const auto& cell_volumes = DENEB_DATA->GetCellVolumes();
  static const auto& cell_proj_volumes = DENEB_DATA->GetCellProjVolumes();
  static const auto& cell_basis_value = DENEB_DATA->GetCellBasisValue();

  int ind = 0;
  const int ipoint = 0;
  for (int icell = 0; icell < num_cells; icell++) {
    const double& Vx = cell_proj_volumes[icell * D_];
    const double& Vy = cell_proj_volumes[icell * D_ + 1];

    cblas_dcopy(S_, &solution[icell * sb], num_bases, &owner_solution[0], 1);
    cblas_dscal(S_, cell_basis_value[icell][0], &owner_solution[0], 1);
    GET_SOLUTION_PS(, owner_solution);
    COMPUTE_FLUX_BASICS();

    const double alpha = 0.5 * d_inv * (r_ * p + BB);
    const double beta_x = r_ * p * std::pow(d_inv * B1, 2);
    const double beta_y = r_ * p * std::pow(d_inv * B2, 2);
    const double r_max_x =
        std::abs(u) + std::sqrt(alpha + std::sqrt(alpha * alpha - beta_x));
    const double r_max_y =
        std::abs(v) + std::sqrt(alpha + std::sqrt(alpha * alpha - beta_y));

    local_timestep[icell] =
        cell_volumes[icell] / (r_max_x * Vx + r_max_y * Vy);
  }
  avocado::VecScale(num_cells, 1.0 / static_cast<double>(2 * order + 1),
                    &local_timestep[0]);
}
bool EquationGLMMHD2D::IsContact(const int& icell,
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
  const double owner_p = ComputePressure(&sol[0]);

  for (auto&& cell : neighbor_cells) {
    const double* target = (cell < num_cells)
                               ? &solution[cell * sb]
                               : &total_solution[(cell - num_cells) * sb];
    cblas_daxpby(S_, 1.0 / std::sqrt(volume[cell]), target, num_bases, 0.0,
                 &sol[0], 1);
    const double neighbor_p = ComputePressure(&sol[0]);

    if (std::abs(owner_p - neighbor_p) > 1.0E-2 * owner_p) return false;
  }
  return true;
}
double EquationGLMMHD2D::ComputeMaxCharacteristicSpeed(
    const double* input_solution) const {
  ERROR_MESSAGE("Not Supported!");
  return 0.0;
}
const std::vector<double>& EquationGLMMHD2D::ComputePressureFixValues(
    const double* input_solution)
{
  pressure_fix_values_[0] = input_solution[0];
  pressure_fix_values_[1] = ComputePressure(&input_solution[0]);
  return pressure_fix_values_;
}
void EquationGLMMHD2D::PreProcess(const double* solution) {
  static const int& num_bases = DENEB_DATA->GetNumBases();
  static const int sb = S_ * num_bases;
  static std::vector<double> owner_solution(S_);

  static const int& num_cells = DENEB_DATA->GetNumCells();
  static const auto& cell_basis_value = DENEB_DATA->GetCellBasisValue();

  // compute ch_
  int ind = 0;
  const int ipoint = 0;
  double r_max = 0.0;
  double vel_max = 0.0;
  for (int icell = 0; icell < num_cells; icell++) {
    cblas_dcopy(S_, &solution[icell * sb], num_bases, &owner_solution[0], 1);
    cblas_dscal(S_, cell_basis_value[icell][0], &owner_solution[0], 1);
    GET_SOLUTION_PS(, owner_solution);
    COMPUTE_FLUX_BASICS();

    const double alpha = 0.5 * d_inv * (r_ * p + BB);
    const double beta_x = r_ * p * std::pow(d_inv * B1, 2);
    const double beta_y = r_ * p * std::pow(d_inv * B2, 2);
    const double abs_u = std::abs(u);
    const double abs_v = std::abs(v);
    const double r_max_x =
        abs_u + std::sqrt(alpha + std::sqrt(alpha * alpha - beta_x));
    const double r_max_y =
        abs_v + std::sqrt(alpha + std::sqrt(alpha * alpha - beta_y));

    vel_max = std::max(vel_max, std::max(abs_u, abs_v));
    r_max = std::max(r_max, std::max(r_max_x, r_max_y));
  }
  vel_max = AVOCADO_MPI->Reduce(vel_max, avocado::MPI::Op::MAX);
  r_max = AVOCADO_MPI->Reduce(r_max, avocado::MPI::Op::MAX);
  ch_ = r_max - vel_max;
}
void EquationGLMMHD2D::ComputeRHS(const double* solution, double* rhs,
                                  const double t) {
  static const int& num_bases = DENEB_DATA->GetNumBases();
  static const int& num_cells = DENEB_DATA->GetNumCells();
  static const int sb = S_ * num_bases;
  static const int dsb = DS_ * num_bases;

  static std::vector<double> owner_solution(S_ * max_num_points_);
  static std::vector<double> owner_solution_grad(DS_ * max_num_points_);
  static std::vector<double> flux(DS_ * max_num_points_);
  static std::vector<double> neighbor_solution(S_ * max_num_face_points_);
  static std::vector<double> solution_difference(S_ * max_num_face_points_);
  static std::vector<double> neighbor_solution_grad;
  static std::vector<double> source(S_ * max_num_points_);

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
  static const auto& face_owner_coefficients =
      DENEB_DATA->GetFaceOwnerCoefficients();
  static const auto& face_neighbor_coefficients =
      DENEB_DATA->GetFaceNeighborCoefficients();
  for (int iface = 0; iface < num_inner_faces; iface++) {
    const int& num_points = num_face_points[iface];
    const int& owner_cell = face_owner_cell[iface];
    const int& neighbor_cell = face_neighbor_cell[iface];

    avocado::Kernel0::f4(
        &solution[owner_cell * sb], &face_owner_basis_value[iface][0],
        &owner_solution[0], S_, num_bases, num_points, 1.0, 0.0);
    avocado::Kernel0::f4(
        &solution[neighbor_cell * sb], &face_neighbor_basis_value[iface][0],
        &neighbor_solution[0], S_, num_bases, num_points, 1.0, 0.0);

    vdSub(num_points * S_, &neighbor_solution[0], &owner_solution[0],
          &solution_difference[0]);

    avocado::Kernel1::f59(&face_owner_coefficients[iface][0],
                          &solution_difference[0],
                          &auxiliary_solution_[owner_cell * dsb], D_,
                          num_bases, num_points, S_, 0.5, 1.0);

    avocado::Kernel1::f59(&face_neighbor_coefficients[iface][0],
                          &solution_difference[0],
                          &auxiliary_solution_[neighbor_cell * dsb], D_,
                          num_bases, num_points, S_, 0.5, 1.0);

    ComputeNumFlux(num_points, flux, owner_cell, neighbor_cell,
                   owner_solution, owner_solution_grad, neighbor_solution,
                   neighbor_solution_grad, face_normals[iface]);

    avocado::Kernel2::f67(&flux[0], &face_owner_coefficients[iface][0],
                          &rhs[owner_cell * sb], S_, D_, num_points,
                          num_bases, 1.0, 1.0);
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
  static const auto& bdry_owner_coefficients =
      DENEB_DATA->GetBdryOwnerCoefficients();
  for (int ibdry = 0; ibdry < num_bdries; ibdry++) {
    const int& num_points = num_bdry_points[ibdry];
    const int& owner_cell = bdry_owner_cell[ibdry];

    avocado::Kernel0::f4(
        &solution[owner_cell * sb], &bdry_owner_basis_value[ibdry][0],
        &owner_solution[0], S_, num_bases, num_points, 1.0, 0.0);
    boundaries_[ibdry]->ComputeBdrySolution(
        num_points, neighbor_solution, neighbor_solution_grad,
        owner_solution, owner_solution_grad, bdry_normals[ibdry],
        bdry_points[ibdry], t);

    vdSub(num_points * S_, &neighbor_solution[0], &owner_solution[0],
          &solution_difference[0]);

    avocado::Kernel1::f59(&bdry_owner_coefficients[ibdry][0],
                          &solution_difference[0],
                          &auxiliary_solution_[owner_cell * dsb], D_,
                          num_bases, num_points, S_, 1.0, 1.0);

    boundaries_[ibdry]->ComputeBdryFlux(
        num_points, flux, owner_cell, -1, owner_solution, owner_solution_grad,
        neighbor_solution, neighbor_solution_grad, bdry_normals[ibdry],
        bdry_points[ibdry], t);

    avocado::Kernel2::f67(&flux[0], &bdry_owner_coefficients[ibdry][0],
                          &rhs[owner_cell * sb], S_, D_, num_points,
                          num_bases, 1.0, 1.0);
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
                          &solution_difference[0],
                          &auxiliary_solution_[owner_cell * dsb], D_,
                          num_bases, num_points, S_, 0.5, 1.0);

    ComputeNumFlux(num_points, flux, owner_cell, neighbor_cell, 
                   owner_solution, owner_solution_grad, neighbor_solution,
                   neighbor_solution_grad, face_normals[iface]);

    avocado::Kernel2::f67(&flux[0], &face_owner_coefficients[iface][0],
                          &rhs[owner_cell * sb], S_, D_, num_points,
                          num_bases, 1.0, 1.0);
  }

  // cell sweep
  static const auto& num_cell_points = DENEB_DATA->GetNumCellPoints();
  static const auto& cell_basis_value = DENEB_DATA->GetCellBasisValue();
  static const auto& cell_basis_grad_value =
      DENEB_DATA->GetCellBasisGradValue();
  static const auto& cell_coefficients = DENEB_DATA->GetCellCoefficients();
  static const auto& cell_source_coefficients =
      DENEB_DATA->GetCellSourceCoefficients();
  for (int icell = 0; icell < num_cells; icell++) {
    const int& num_points = num_cell_points[icell];

    avocado::Kernel0::f4(&solution[icell * sb], &cell_basis_value[icell][0],
                         &owner_solution[0], S_, num_bases, num_points, 1.0,
                         0.0);
    avocado::Kernel1::f42(&cell_basis_grad_value[icell][0],
                          &solution[icell * sb], &owner_solution_grad[0], D_,
                          num_points, num_bases, S_, 1.0, 0.0);
    avocado::Kernel1::f18(
        &auxiliary_solution_[icell * dsb], &cell_basis_value[icell][0],
        &owner_solution_grad[0], D_, S_, num_bases, num_points, 1.0, 1.0);

    ComputeComFlux(num_points, flux, icell, owner_solution,
                   owner_solution_grad);

    avocado::Kernel2::f67(&flux[0], &cell_coefficients[icell][0],
                          &rhs[icell * sb], S_, D_, num_points, num_bases,
                          -1.0, 1.0);

    ComputeSource(num_points, source, icell, owner_solution, owner_solution_grad);

    gemmATB(-1.0, &source[0], &cell_source_coefficients[icell][0], 1.0,
            &rhs[icell * sb], num_points, S_, num_bases);
  }
}
void EquationGLMMHD2D::ComputeSystemMatrix(const double* solution, Mat& sysmat,
                                           const double t) {
  ERROR_MESSAGE("Not supported!\n");
}
void EquationGLMMHD2D::ComputeComFlux(const int num_points,
                                     std::vector<double>& flux, const int icell,
                                     const std::vector<double>& owner_u,
                                     const std::vector<double>& owner_div_u) {
  int ind = 0;
  for (int ipoint = 0; ipoint < num_points; ipoint++) {
    GET_SOLUTION_PS(, owner_u);

    COMPUTE_FLUX_BASICS();

    COMPUTE_VOLUME_FLUX_PDS(flux, ipoint);
  }
}
void EquationGLMMHD2D::ComputeComFluxJacobi( 
    const int num_points, std::vector<double>& flux_jacobi,
    std::vector<double>& flux_grad_jacobi, const int icell,
    const std::vector<double>& owner_u,
    const std::vector<double>& owner_div_u) {
  ERROR_MESSAGE("Not supported!\n");
}
void EquationGLMMHD2D::ComputeSource(const int num_points,
                                     std::vector<double>& source,
                                     const int cell,
                                     const std::vector<double>& owner_u,
                                     const std::vector<double>& owner_div_u) {
  int ind = 0;
  for (int ipoint = 0; ipoint < num_points; ipoint++) {
    GET_SOLUTION_PS(, owner_u);

    const double& B1x = owner_div_u[DS_ * ipoint + 4];
    const double& B2y = owner_div_u[DS_ * ipoint + 12];
    const double& psix = owner_div_u[DS_ * ipoint + 6];
    const double& psiy = owner_div_u[DS_ * ipoint + 13];

    const auto d_inv = 1.0 / d;
    const auto u = du * d_inv;
    const auto v = dv * d_inv;
    const auto VB = u * B1 + v * B2;
    const auto divB = B1x + B2y;
    const auto divpsiV = psix * u + psiy * v;

    COMPUTE_VOLUME_SOURCE_PS(source, ipoint);
  }
}
void EquationGLMMHD2D::ComputeNumFluxLLF(const int num_points,
                                         std::vector<double>& flux,
                                         FACE_INPUTS) {
  int ind = 0;
  for (int ipoint = 0; ipoint < num_points; ipoint++) {
    GET_NORMAL_PD(normal);

    GET_SOLUTION_PS(_o, owner_u);
    GET_SOLUTION_PS(_n, neighbor_u);

    COMPUTE_FLUX_BASICS(_o);
    const double V_o = u_o * nx + v_o * ny;
    const double min_BB_o = std::min(B1_o * B1_o, B2_o * B2_o);
    const double alpha_o = 0.5 * d_inv_o * (r_ * p_o + BB_o);
    const double beta_o = r_ * p_o * min_BB_o * std::pow(d_inv_o, 2);
    const double r_max_o =
        std::abs(V_o) +
        std::sqrt(alpha_o + std::sqrt(alpha_o * alpha_o - beta_o));

    COMPUTE_FLUX_BASICS(_n);
    const double V_n = u_n * nx + v_n * ny;
    const double min_BB_n = std::min(B1_n * B1_n, B2_n * B2_n);
    const double alpha_n = 0.5 * d_inv_n * (r_ * p_n + BB_n);
    const double beta_n = r_ * p_n * min_BB_n * std::pow(d_inv_n, 2);
    const double r_max_n =
        std::abs(V_n) +
        std::sqrt(alpha_n + std::sqrt(alpha_n * alpha_n - beta_n));
    const double r_max = std::max(r_max_o, r_max_n);

    const auto diff1 = r_max * (d_n - d_o);
    const auto diff2 = r_max * (du_n - du_o);
    const auto diff3 = r_max * (dv_n - dv_o);
    const auto diff4 = r_max * (dtE_n - dtE_o);
    const auto diff5 = r_max * (B1_n - B1_o);
    const auto diff6 = r_max * (B2_n - B2_o);
    const auto diff7 = r_max * (psi_n - psi_o);

    COMPUTE_NUMERICAL_FLUX_PDS(flux, ipoint);
  }
}
void EquationGLMMHD2D::ComputeNumFluxJacobiLLF(const int num_points,
                                               FACE_JACOBI_OUTPUTS,
                                               FACE_INPUTS) {
  ERROR_MESSAGE("Not supported!\n");
}

// ------------------------------- Boundary --------------------------------
// //
std::shared_ptr<BoundaryGLMMHD2D> BoundaryGLMMHD2D::GetBoundary(
    const std::string& type, const int bdry_tag, EquationGLMMHD2D* equation) {
  if (!type.compare("Extrapolate"))
    return std::make_shared<ExtrapolateGLMMHD2D>(bdry_tag, equation);
  else if (!type.compare("Constant"))
    return std::make_shared<ConstantBdryGLMMHD2D>(bdry_tag, equation);
  ERROR_MESSAGE("Wrong boundary condition (no-exist):" + type + "\n");
  return nullptr;
}
// Boundary =  Extrapolate
// Dependency: -
ExtrapolateGLMMHD2D::ExtrapolateGLMMHD2D(const int bdry_tag,
                                         EquationGLMMHD2D* equation)
    : BoundaryGLMMHD2D(bdry_tag, equation) {
  MASTER_MESSAGE("Extrapolate (tag=" + std::to_string(bdry_tag) + ")\n");
}
void ExtrapolateGLMMHD2D::ComputeBdrySolution(
    const int num_points, std::vector<double>& bdry_u,
    std::vector<double>& bdry_div_u, const std::vector<double>& owner_u,
    const std::vector<double>& owner_div_u, const std::vector<double>& normal,
    const std::vector<double>& coords, const double& time) {
  cblas_dcopy(num_points * S_, &owner_u[0], 1, &bdry_u[0], 1);
}
void ExtrapolateGLMMHD2D::ComputeBdryFlux(const int num_points,
                                          std::vector<double>& flux,
                                          FACE_INPUTS,
                                          const std::vector<double>& coords,
                                          const double& time) {
  int ind = 0;
  const double ch_ = equation_->GetCh();
  for (int ipoint = 0; ipoint < num_points; ipoint++) {
    GET_SOLUTION_PS(, owner_u);

    COMPUTE_FLUX_BASICS();

    COMPUTE_VOLUME_FLUX_PDS(flux, ipoint);
  }
}
void ExtrapolateGLMMHD2D::ComputeBdrySolutionJacobi(
    const int num_points, 
    double* bdry_u_jacobi, const std::vector<double>& owner_u,
    const std::vector<double>& owner_div_u, const std::vector<double>& normal,
    const std::vector<double>& coords, const double& time) {
  ERROR_MESSAGE("Not supported!\n");
}
void ExtrapolateGLMMHD2D::ComputeBdryFluxJacobi(
    const int num_points, std::vector<double>& flux_jacobi,
    std::vector<double>& flux_grad_jacobi, FACE_INPUTS,
    const std::vector<double>& coords, const double& time) {
  ERROR_MESSAGE("Not supported!\n");
}
// Boundary = Constant
// Dependency: BdryInput(num)
// BdryInput(num) = rho, rhoU, rhoV, rhotE, B1, B2
ConstantBdryGLMMHD2D::ConstantBdryGLMMHD2D(const int bdry_tag,
                                           EquationGLMMHD2D* equation)
    : BoundaryGLMMHD2D(bdry_tag, equation) {
  MASTER_MESSAGE("Constant boundary (tag=" + std::to_string(bdry_tag) + ")\n");

  auto& config = AVOCADO_CONFIG;
  values_.resize(S_);
  for (int i = 0; i < S_; i++)
    values_[i] = std::stod(config->GetConfigValue(BDRY_INPUT_I(bdry_tag, i)));

  std::vector<std::string> names = {"rho", "rhoU", "rhoV", "rhotE",
                                    "B1",  "B2",   "psi"};
  std::stringstream str;
  str << "\tInput: ";
  for (int i = 0; i < S_; i++)
    str << names[i] << " = " << std::scientific << values_[i] << "\n\t";
  MASTER_MESSAGE(str.str());
}
void ConstantBdryGLMMHD2D::ComputeBdrySolution(
    const int num_points, std::vector<double>& bdry_u,
    std::vector<double>& bdry_div_u, const std::vector<double>& owner_u,
    const std::vector<double>& owner_div_u, const std::vector<double>& normal,
    const std::vector<double>& coords, const double& time) {
  int ind = 0;
  for (int ipoint = 0; ipoint < num_points; ipoint++)
    for (int istate = 0; istate < S_; istate++) bdry_u[ind++] = values_[istate];
}
void ConstantBdryGLMMHD2D::ComputeBdryFlux(const int num_points,
                                           std::vector<double>& flux,
                                           FACE_INPUTS,
                                           const std::vector<double>& coords,
                                           const double& time) {
  equation_->ComputeNumFlux(num_points, flux, owner_cell, neighbor_cell, owner_u, owner_div_u, neighbor_u,
                            neighbor_div_u, normal);
}
void ConstantBdryGLMMHD2D::ComputeBdrySolutionJacobi(
    const int num_points, 
    double* bdry_u_jacobi, const std::vector<double>& owner_u,
    const std::vector<double>& owner_div_u, const std::vector<double>& normal,
    const std::vector<double>& coords, const double& time) {
  ERROR_MESSAGE("Not supported!\n");
}
void ConstantBdryGLMMHD2D::ComputeBdryFluxJacobi(
    const int num_points, std::vector<double>& flux_jacobi,
    std::vector<double>& flux_grad_jacobi, FACE_INPUTS,
    const std::vector<double>& coords, const double& time) {
  ERROR_MESSAGE("Not supported!\n");
}
// -------------------------------- Problem --------------------------------
// //
std::shared_ptr<ProblemGLMMHD2D> ProblemGLMMHD2D::GetProblem(
    const std::string& name) {
  if (!name.compare("DoubleSine"))
    return std::make_shared<DoubleSineGLMMHD2D>();
  else if (!name.compare("ShockTube"))
    return std::make_shared<ShockTubeGLMMHD2D>();
  else if (!name.compare("OrszagTangVortex"))
    return std::make_shared<OrszagTangVortexGLMMHD2D>();
  else if (!name.compare("MHDRotor"))
    return std::make_shared<MHDRotorGLMMHD2D>();
  ERROR_MESSAGE("Wrong problem (no-exist):" + name + "\n");
  return nullptr;
}
// Problem = DoubleSine
// ProblemInput = -
DoubleSineGLMMHD2D::DoubleSineGLMMHD2D()
    : velocity_({0.5, 0.5}), wave_number_({1.0, 1.0}) {
  MASTER_MESSAGE("DoubleSine problem\n");

  std::stringstream str;
  str << "\tInput: ";
  str << "x-velocity = " << velocity_[0] << "\n\t";
  str << "y-velocity = " << velocity_[1] << "\n\t";
  str << "x-wavenumber = " << wave_number_[0] << "\n\t";
  str << "y-wavenumber = " << wave_number_[1] << "\n\t";
  MASTER_MESSAGE(str.str());
}
void DoubleSineGLMMHD2D::Problem(const int num_points,
                                 std::vector<double>& solutions,
                                 const std::vector<double>& coord,
                                 const double time) const {
  const double p = 1.0;
  const double B1 = 1.0;
  const double B2 = 0.3;
  const double psi = 0.0;
  const double BB = B1 * B1 + B2 * B2;
  for (int ipoint = 0; ipoint < num_points; ipoint++) {
    const double& x = coord[ipoint * D_];
    const double& y = coord[ipoint * D_ + 1];

    const double& u = velocity_[0];
    const double& v = velocity_[1];
    const double d =
        1.0 +
        0.2 *
            std::sin(2.0 * pi_ * wave_number_[0] * (x - velocity_[0] * time)) *
            std::sin(2.0 * pi_ * wave_number_[1] * (y - velocity_[1] * time));

    double* sol = &solutions[ipoint * S_];
    sol[0] = d;
    sol[1] = d * u;
    sol[2] = d * v;
    sol[3] =
        rm1_inv_ * p + 0.5 * d * (u * u + v * v) + 0.5 * BB + 0.5 * psi * psi;
    sol[4] = B1;
    sol[5] = B2;
    sol[6] = psi;
  }
}
// Problem = ShockTube
// ProblemInput = x-split,
//                left_rho, left_rhoU, left_rhotE, left_B1, left_B2,
//                right_rho, right_rhoU, right_rhotE, right_B1, right_B2
ShockTubeGLMMHD2D::ShockTubeGLMMHD2D() {
  MASTER_MESSAGE("ShockTube problem\n");
  auto& config = AVOCADO_CONFIG;

  int ind = 0;
  split_ = std::stod(config->GetConfigValue(PROBLEM_INPUT_I(ind++)));
  left_values_.resize(S_, 0.0);
  left_values_[0] = std::stod(config->GetConfigValue(PROBLEM_INPUT_I(ind++)));
  left_values_[1] = std::stod(config->GetConfigValue(PROBLEM_INPUT_I(ind++)));
  left_values_[3] = std::stod(config->GetConfigValue(PROBLEM_INPUT_I(ind++)));
  left_values_[4] = std::stod(config->GetConfigValue(PROBLEM_INPUT_I(ind++)));
  left_values_[5] = std::stod(config->GetConfigValue(PROBLEM_INPUT_I(ind++)));

  right_values_.resize(S_, 0.0);
  right_values_[0] = std::stod(config->GetConfigValue(PROBLEM_INPUT_I(ind++)));
  right_values_[1] = std::stod(config->GetConfigValue(PROBLEM_INPUT_I(ind++)));
  right_values_[3] = std::stod(config->GetConfigValue(PROBLEM_INPUT_I(ind++)));
  right_values_[4] = std::stod(config->GetConfigValue(PROBLEM_INPUT_I(ind++)));
  right_values_[5] = std::stod(config->GetConfigValue(PROBLEM_INPUT_I(ind++)));

  std::vector<std::string> names = {"rho", "rhoU", "rhoV", "rhotE",
                                    "B1",  "B2",   "psi"};
  std::stringstream str;
  str << "\tInput: split = " << split_ << "\n\t";
  for (int i = 0; i < S_; i++)
    str << "left_" << names[i] << " = " << left_values_[i] << "\t"
        << "right_" << names[i] << " = " << right_values_[i] << "\n\t";
  MASTER_MESSAGE(str.str());
}
void ShockTubeGLMMHD2D::Problem(const int num_points,
                                std::vector<double>& solutions,
                                const std::vector<double>& coord,
                                const double time) const {
  int ind = 0;
  for (int ipoint = 0; ipoint < num_points; ipoint++) {
    const double& x = coord[ipoint * D_];

    if (x < split_)
      for (int i = 0; i < S_; i++) solutions[ind++] = left_values_[i];
    else
      for (int i = 0; i < S_; i++) solutions[ind++] = right_values_[i];
  }
}
// Problem = OrszagTangVortex
// ProblemInput = -
OrszagTangVortexGLMMHD2D::OrszagTangVortexGLMMHD2D() {
  MASTER_MESSAGE("Orszag-Tang GLMMHD Vortex problem\n");
}
void OrszagTangVortexGLMMHD2D::Problem(const int num_points,
                                       std::vector<double>& solutions,
                                       const std::vector<double>& coord,
                                       const double time) const {
  const double d = 1.0;
  const double p = r_inv_;
  const double psi = 0.0;
  for (int ipoint = 0; ipoint < num_points; ipoint++) {
    const double& x = coord[ipoint * D_];
    const double& y = coord[ipoint * D_ + 1];

    const double u = -std::sin(2.0 * pi_ * y);
    const double v = std::sin(2.0 * pi_ * x);

    const double B1 = -r_inv_ * std::sin(2.0 * pi_ * y);
    const double B2 = r_inv_ * std::sin(4.0 * pi_ * x);
    const double BB = B1 * B1 + B2 * B2;

    double* sol = &solutions[ipoint * S_];
    sol[0] = d;
    sol[1] = d * u;
    sol[2] = d * v;
    sol[3] =
        rm1_inv_ * p + 0.5 * d * (u * u + v * v) + 0.5 * BB + 0.5 * psi * psi;
    sol[4] = B1;
    sol[5] = B2;
    sol[6] = psi;
  }
}
// Problem = MHDRotor
// ProblemInput = -
MHDRotorGLMMHD2D::MHDRotorGLMMHD2D() { MASTER_MESSAGE("MHD Rotor problem\n"); }
void MHDRotorGLMMHD2D::Problem(const int num_points,
                               std::vector<double>& solutions,
                               const std::vector<double>& coord,
                               const double time) const {
  const double xc = 0.5;
  const double yc = 0.5;
  const double r0 = 0.1;
  const double r1 = 0.115;

  const double p = 0.5;
  const double B1 = 2.5 / std::sqrt(4.0 * pi_);
  const double B2 = 0.0;
  const double BB = B1 * B1 + B2 * B2;

  double d, u, v;
  for (int ipoint = 0; ipoint < num_points; ipoint++) {
    const double& x = coord[ipoint * D_];
    const double& y = coord[ipoint * D_ + 1];

    const double dx = x - xc;
    const double dy = y - yc;
    const double rr = dx * dx + dy * dy;

    if (rr <= r0 * r0) {
      d = 10.0;
      u = -dy / r0;
      v = dx / r0;
    } else if (rr <= r1 * r1) {
      const double r = std::sqrt(rr);
      const double f = (r1 - r) / (r1 - r0);
      d = 1.0 + 9.0 * f;
      u = -f * dy / r;
      v = f * dx / r;
    } else {
      d = 1.0;
      u = 0.0;
      v = 0.0;
    }

    double* sol = &solutions[ipoint * S_];
    sol[0] = d;
    sol[1] = d * u;
    sol[2] = d * v;
    sol[3] = rm1_inv_ * p + 0.5 * d * (u * u + v * v) + 0.5 * BB;
    sol[4] = B1;
    sol[5] = B2;
    sol[6] = 0.0;
  }
}
}  // namespace deneb