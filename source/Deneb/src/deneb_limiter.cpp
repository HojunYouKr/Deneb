#include "deneb_limiter.h"

#include <algorithm>
#include <cstring>
#include <unordered_set>

#include "avocado.h"
#include "deneb_data.h"
#include "deneb_equation.h"

namespace deneb {
std::shared_ptr<Limiter> DENEB_LIMITER_NAME = nullptr;
std::shared_ptr<Limiter> Limiter::GetLimiter(const std::string& name) {
  if (!name.compare("None"))
    return std::make_shared<NoLimiter>();
  else if (!name.compare("hMLP"))
    return std::make_shared<hMLP>();
  else if (!name.compare("hMLP_BD"))
    return std::make_shared<hMLP_BD>();
  ERROR_MESSAGE("Wrong limiter (no-exist):" + name + "\n");
  return nullptr;
}

NoLimiter::NoLimiter() { MASTER_MESSAGE(avocado::GetTitle("NoLimiter")); }
void NoLimiter::BuildData(void) {
  MASTER_MESSAGE(avocado::GetTitle("NoLimiter::BuildData()"));
}

// --------------------- hMLP --------------------- //
hMLP::hMLP() : eps_(1.0e-14), deact_eps_(1.0e-3) {
  MASTER_MESSAGE(avocado::GetTitle("hMLP"));
}
void hMLP::BuildData(void) {
  MASTER_MESSAGE(avocado::GetTitle("hMLP::BuildData()"));

  ConstructNodeCells();

  target_state_ = 0;
  const std::vector<std::string>& variable_names =
      DENEB_EQUATION->GetCellVariableNames();
  MASTER_MESSAGE("Target state: " + variable_names[target_state_] + "\n");

  const int& order = DENEB_DATA->GetOrder();
  num_bases_list_.resize(order + 1);
  for (int i = 0; i < order + 1; i++)
    num_bases_list_[i] = DENEB_DATA->GetNumBases(i);

  const int& num_nodes = DENEB_DATA->GetNumNodes();
  const std::vector<int>& cell_node_ind = DENEB_DATA->GetCellNodeInd();
  const std::vector<int>& cell_node_ptr = DENEB_DATA->GetCellNodePtr();
  nodetypes_.resize(num_nodes, NODETYPE::NORMAL);
  const int& num_bdries = DENEB_DATA->GetNumBdries();
  const auto& cell_element = DENEB_DATA->GetCellElement();
  const std::vector<int>& bdry_owner_cell = DENEB_DATA->GetBdryOwnerCell();
  const std::vector<int>& bdry_owner_type = DENEB_DATA->GetBdryOwnerType();
  for (int ibdry = 0; ibdry < num_bdries; ibdry++) {
    const int& owner_cell = bdry_owner_cell[ibdry];
    const int& owner_type = bdry_owner_type[ibdry];

    std::vector<int> bdry_nodes =
        cell_element[owner_cell]->GetFacetypeNodes(1)[owner_type];
    const int& start_index = cell_node_ptr[owner_cell];
    for (auto&& inode : bdry_nodes) inode = cell_node_ind[start_index + inode];
    for (auto&& inode : bdry_nodes) nodetypes_[inode] = NODETYPE::BOUNDARY;
  }

  const int& num_cells = DENEB_DATA->GetNumCells();
  const int& num_total_cells = DENEB_DATA->GetNumTotalCells();
  const int& num_states = DENEB_EQUATION->GetNumStates();
  const int& num_bases = DENEB_DATA->GetNumBases();
  const int sb = num_states * num_bases;
  foreign_solution_.resize(std::max((num_total_cells - num_cells) * sb, 1));
  communicate_ = std::make_shared<avocado::Communicate>(
      sb, DENEB_DATA->GetTotalSendCellList(),
      DENEB_DATA->GetTotalRecvCellList());

  const int& dimension = DENEB_EQUATION->GetDimension();
  const std::vector<double>& node_coords = DENEB_DATA->GetNodeCoords();
  cell_vertex_basis_value_.resize(num_cells);
  for (int icell = 0; icell < num_cells; icell++) {
    const int num_vertex = cell_node_ptr[icell + 1] - cell_node_ptr[icell];

    std::vector<double> coords;
    for (int v = 0; v < num_vertex; v++) {
      const int& node = cell_node_ind[cell_node_ptr[icell] + v];
      const double* coord = &node_coords[node * dimension];
      for (int idim = 0; idim < dimension; idim++)
        coords.push_back(coord[idim]);
    }
    DENEB_DATA->GetCellBasisValues(icell, coords,
                                   cell_vertex_basis_value_[icell]);
  }

  foreign_cell_basis_value_.resize(std::max(num_total_cells - num_cells, 1));
  {
    std::shared_ptr<avocado::Communicate> communicate =
        std::make_shared<avocado::Communicate>(
            1, DENEB_DATA->GetTotalSendCellList(),
            DENEB_DATA->GetTotalRecvCellList());
    std::vector<double> cell_basis_value(num_cells);
    for (int icell = 0; icell < num_cells; icell++)
      cell_basis_value[icell] = cell_vertex_basis_value_[icell][0];
    communicate->CommunicateBegin(&cell_basis_value[0]);
    communicate->CommunicateEnd(&foreign_cell_basis_value_[0]);
  }

  cell_average_.resize(num_total_cells * num_states);
  vertex_min_.resize(num_states * num_nodes);
  vertex_max_.resize(num_states * num_nodes);
}
void hMLP::ConstructNodeCells(void) {
  static const std::vector<int>& node_global_index =
      DENEB_DATA->GetNodeGlobalIndex();
  std::unordered_map<int, int> node_mapping;  // global to local index
  {
    int ind = 0;
    for (auto&& global_index : node_global_index)
      node_mapping[global_index] = ind++;
  }

  const std::unordered_map<int, std::vector<int>>&
      periodic_matching_global_node_index =
          DENEB_DATA->GetPeriodicMatchingGlobalNodeIndex();
  const auto& node_grps = periodic_matching_global_node_index;

  static const int& num_nodes = DENEB_DATA->GetNumNodes();
  static const int& num_total_cells = DENEB_DATA->GetNumTotalCells();
  const std::vector<int>& cell_node_ptr = DENEB_DATA->GetCellNodePtr();
  const std::vector<int>& cell_node_ind = DENEB_DATA->GetCellNodeInd();
  std::vector<std::vector<int>> node_cells;
  {
    std::unordered_map<int, std::unordered_set<int>> node_cells_temp;
    for (int i = 0; i < num_total_cells; i++)
      for (int ptr = cell_node_ptr[i]; ptr < cell_node_ptr[i + 1]; ptr++)
        node_cells_temp[cell_node_ind[ptr]].insert(i);

    for (auto&& iterator : node_grps) {
      const int& mynode = node_mapping.at(iterator.first);

      std::unordered_set<int> new_cells;
      for (auto&& node : iterator.second) {
        const auto& cells = node_cells_temp.find(node_mapping.at(node));
        if (cells == node_cells_temp.end()) continue;
        new_cells.insert(cells->second.begin(), cells->second.end());
      }
      node_cells_temp[mynode].insert(new_cells.begin(), new_cells.end());
    }

    node_cells.resize(num_nodes);
    for (auto&& iterator : node_cells_temp) {
      if (iterator.first >= num_nodes) continue;
      node_cells[iterator.first] =
          std::vector<int>(iterator.second.begin(), iterator.second.end());
    }
  }

  std::vector<std::vector<int>> node_vertices(num_nodes);
  for (int i = 0; i < num_nodes; i++) {
    const auto& matching_nodes = node_grps.find(node_global_index[i]);
    for (auto&& cell : node_cells[i]) {
      for (int ptr = cell_node_ptr[cell]; ptr < cell_node_ptr[cell + 1];
           ptr++) {
        if (cell_node_ind[ptr] == i) {
          node_vertices[i].push_back(ptr - cell_node_ptr[cell]);
          break;
        }
        if (matching_nodes == node_grps.end()) continue;
        const int& node = node_global_index[cell_node_ind[ptr]];
        if (std::find(matching_nodes->second.begin(),
                      matching_nodes->second.end(),
                      node) != matching_nodes->second.end()) {
          node_vertices[i].push_back(ptr - cell_node_ptr[cell]);
          break;
        }
      }
    }
  }

  node_cells_ = std::move(node_cells);
  node_vertices_ = std::move(node_vertices);
}
void hMLP::Limiting(double* solution) {
  static const int& order = DENEB_DATA->GetOrder();
  if (order == 0) return;

  ComputeVertexMinMax(solution);

  static const int& num_states = DENEB_EQUATION->GetNumStates();
  static const int& num_bases = DENEB_DATA->GetNumBases();
  static const int& num_cells = DENEB_DATA->GetNumCells();

  const int sb = num_states * num_bases;
  for (int icell = 0; icell < num_cells; icell++) {
    bool is_trouble = true;
    int temp_order = order;

    while (true) {
      if (temp_order == 0) {
        is_trouble = false;
        break;
      }

      is_trouble = Indicator(solution, icell, temp_order);

      if (is_trouble == true) {
        if (temp_order == 1)
          break;
        else if (temp_order == 2) {
          temp_order = 1;
          Projection(solution, icell, 1);
          break;
        } else {
          is_trouble = false;
          temp_order = temp_order - 1;
          Projection(solution, icell, temp_order);
        }
      } else
        break;
    }
    if (is_trouble == true) {
      MLPu1(solution, icell);
    }
  }

  return;
}

void hMLP::ComputeVertexMinMax(double* solution) {
  static const int& num_cells = DENEB_DATA->GetNumCells();
  static const int& num_states = DENEB_EQUATION->GetNumStates();
  static const int& num_bases = DENEB_DATA->GetNumBases();
  const int sb = num_states * num_bases;

  communicate_->CommunicateBegin(solution);
  for (int icell = 0; icell < num_cells; icell++)
    for (int istate = 0; istate < num_states; istate++)
      cell_average_[icell * num_states + istate] =
          solution[icell * sb + istate * num_bases] *
          cell_vertex_basis_value_[icell][0];

  static const int& num_total_cells = DENEB_DATA->GetNumTotalCells();
  communicate_->CommunicateEnd(&foreign_solution_[0]);
  for (int icell = num_cells; icell < num_total_cells; icell++)
    for (int istate = 0; istate < num_states; istate++)
      cell_average_[icell * num_states + istate] =
          foreign_solution_[(icell - num_cells) * sb + istate * num_bases] *
          foreign_cell_basis_value_[icell - num_cells];

  static const int& num_nodes = DENEB_DATA->GetNumNodes();
  for (int inode = 0; inode < num_nodes; inode++) {
    for (int istate = 0; istate < num_states; istate++) {
      std::vector<double> avgs;
      for (auto&& icell : node_cells_[inode])
        avgs.push_back(cell_average_[icell * num_states + istate]);

      vertex_min_[istate * num_nodes + inode] =
          *std::min_element(avgs.begin(), avgs.end());
      vertex_max_[istate * num_nodes + inode] =
          *std::max_element(avgs.begin(), avgs.end());
    }
  }
}
bool hMLP::Indicator(const double* solution_ptr, const int icell,
                     const int order) {
  static const int& num_states = DENEB_EQUATION->GetNumStates();
  static const int& num_basis = DENEB_DATA->GetNumBases();
  static const int& num_cells = DENEB_DATA->GetNumCells();
  static const int sb = num_states * num_basis;
  static const std::vector<int>& cell_node_ptr = DENEB_DATA->GetCellNodePtr();
  static const std::vector<int>& cell_node_ind = DENEB_DATA->GetCellNodeInd();

  const int num_vertex = cell_node_ptr[icell + 1] - cell_node_ptr[icell];
  const double average_value =
      solution_ptr[icell * sb + target_state_ * num_basis] *
      cell_vertex_basis_value_[icell][0];
  for (int ivertex = 0; ivertex < num_vertex; ivertex++) {
    const int& node = cell_node_ind[cell_node_ptr[icell] + ivertex];
    if (nodetypes_[node] == NODETYPE::BOUNDARY) continue;

    const double vertex_value = cblas_ddot(
        num_basis, &solution_ptr[icell * sb + target_state_ * num_basis], 1,
        &cell_vertex_basis_value_[icell][ivertex * num_basis], 1);
    const double P1var =
        cblas_ddot(num_bases_list_[1],
                   &solution_ptr[icell * sb + target_state_ * num_basis], 1,
                   &cell_vertex_basis_value_[icell][ivertex * num_basis], 1);

    bool flag = MLPcondition(solution_ptr, icell, ivertex, target_state_, P1var,
                             vertex_value, average_value);
    if (flag == true) {
      if (order > 1) {
        flag = SmoothDetect(solution_ptr, icell, ivertex, target_state_, P1var,
                            vertex_value, average_value);
      }
    }
    if (flag == true) {
      return true;
    }
  }
  return false;
}

bool hMLP::MLPcondition(const double* solution_ptr, const int icell,
                        const int ivertex, const int istate, const double P1var,
                        const double vertex_value, const double average_value) {
  static const int& num_nodes = DENEB_DATA->GetNumNodes();
  static const std::vector<int>& cell_node_ptr = DENEB_DATA->GetCellNodePtr();
  static const std::vector<int>& cell_node_ind = DENEB_DATA->GetCellNodeInd();
  static const std::vector<double>& cell_volumes = DENEB_DATA->GetCellVolumes();

  if (Deactivation(solution_ptr, cell_volumes[icell], vertex_value,
                   average_value))
    return false;

  if (P1var > vertex_max_[istate * num_nodes +
                          cell_node_ind[cell_node_ptr[icell] + ivertex]])
    return true;
  else if (P1var < vertex_min_[istate * num_nodes +
                               cell_node_ind[cell_node_ptr[icell] + ivertex]])
    return true;
  else
    return false;
}
bool hMLP::Deactivation(const double* solution_ptr, const double volume,
                        const double vertex_value, const double average_value) {
  if (std::abs(vertex_value - average_value) <=
      std::max(deact_eps_ * std::abs(average_value), volume))
    return true;
  return false;
}
bool hMLP::SmoothDetect(const double* solution_ptr, const int icell,
                        const int ivertex, const int istate, const double P1var,
                        const double vertex_value, const double average_value) {
  static const int num_nodes = DENEB_DATA->GetNumNodes();
  static const std::vector<int>& cell_node_ptr = DENEB_DATA->GetCellNodePtr();
  static const std::vector<int>& cell_node_ind = DENEB_DATA->GetCellNodeInd();

  const double high_term = vertex_value - P1var;
  const double P1term = P1var - average_value;

  if ((P1term > 0.0) && (high_term < eps_) &&
      (vertex_value > vertex_min_[istate * num_nodes +
                                  cell_node_ind[cell_node_ptr[icell] +
                                                ivertex]]))  // c1 condition
    return false;
  else if ((P1term < 0.0) && (high_term > -eps_) &&
           (vertex_value <
            vertex_max_[istate * num_nodes +
                        cell_node_ind[cell_node_ptr[icell] +
                                      ivertex]]))  // c2 condition
    return false;

  return true;
}
void hMLP::Projection(double* solution_ptr, const int icell, const int iorder) {
  static const int& num_states = DENEB_EQUATION->GetNumStates();
  static const int& num_bases = DENEB_DATA->GetNumBases();
  static const int sb = num_states * num_bases;

  for (int istate = 0; istate < num_states; istate++)
    for (int ibasis = num_bases_list_[iorder]; ibasis < num_bases; ibasis++)
      solution_ptr[icell * sb + istate * num_bases + ibasis] = 0.0;
}
void hMLP::MLPu1(double* solution_ptr, const int icell) {
  static const int& num_states = DENEB_EQUATION->GetNumStates();
  static const int& num_bases = DENEB_DATA->GetNumBases();
  static const int& num_nodes = DENEB_DATA->GetNumNodes();
  static const int sb = num_states * num_bases;
  static const std::vector<int>& cell_node_ptr = DENEB_DATA->GetCellNodePtr();
  static const std::vector<int>& cell_node_ind = DENEB_DATA->GetCellNodeInd();

  const int num_vertex = cell_node_ptr[icell + 1] - cell_node_ptr[icell];

  std::vector<double> p1_value(num_states * num_vertex, 0.0);
  for (int istate = 0; istate < num_states; istate++)
    for (int ivertex = 0; ivertex < num_vertex; ivertex++)
      p1_value[istate * num_vertex + ivertex] = cblas_ddot(
          num_bases_list_[1], &solution_ptr[icell * sb + istate * num_bases], 1,
          &cell_vertex_basis_value_[icell][ivertex * num_bases], 1);

  for (int istate = 0; istate < num_states; istate++) {
    double temp, phi = 1.0;
    for (int ivertex = 0; ivertex < num_vertex; ivertex++) {
      const int& node = cell_node_ind[cell_node_ptr[icell] + ivertex];
      if (nodetypes_[node] == NODETYPE::BOUNDARY) continue;

      temp = solution_ptr[icell * sb + istate * num_bases] *
                 cell_vertex_basis_value_[icell][0] -
             p1_value[istate * num_vertex + ivertex];
      if (p1_value[istate * num_vertex + ivertex] <
          vertex_min_[istate * num_nodes +
                      cell_node_ind[cell_node_ptr[icell] + ivertex]]) {
        if (std::abs(temp) < 1.0e-14)
          temp = 0.0;
        else
          temp = (solution_ptr[icell * sb + istate * num_bases] *
                      cell_vertex_basis_value_[icell][0] -
                  vertex_min_[istate * num_nodes +
                              cell_node_ind[cell_node_ptr[icell] + ivertex]]) /
                 temp;
        if (temp < phi) phi = temp;
      } else if (vertex_max_[istate * num_nodes +
                             cell_node_ind[cell_node_ptr[icell] + ivertex]] <
                 p1_value[istate * num_vertex + ivertex]) {
        if (std::abs(temp) < 1.0e-14)
          temp = 0.0;
        else
          temp = (solution_ptr[icell * sb + istate * num_bases] *
                      cell_vertex_basis_value_[icell][0] -
                  vertex_max_[istate * num_nodes +
                              cell_node_ind[cell_node_ptr[icell] + ivertex]]) /
                 temp;
        if (temp < phi) phi = temp;
      }
    }
    cblas_dscal(num_bases_list_[1] - 1, phi,
                &solution_ptr[icell * sb + istate * num_bases + 1], 1);
  }
}

hMLP_BD::hMLP_BD() : eps_(1.0e-14), deact_eps_(1.0e-3) {
  MASTER_MESSAGE(avocado::GetTitle("hMLP_BD"));
}
void hMLP_BD::BuildData(void) {
  SYNCRO();
  START_TIMER_TAG("BuildData");
  MASTER_MESSAGE(avocado::GetTitle("hMLP_BD::BuildData()"));

  START_TIMER();
  MASTER_MESSAGE("Constructing cell to cells data... ");
  ConstructCellCells();
  SYNCRO();
  MASTER_MESSAGE("Complete. (Time: " + std::to_string(STOP_TIMER()) + "s)\n");

  START_TIMER();
  MASTER_MESSAGE("Constructing node to cells data... ");
  ConstructNodeCells();
  SYNCRO();
  MASTER_MESSAGE("Complete. (Time: " + std::to_string(STOP_TIMER()) + "s)\n");

  START_TIMER();
  MASTER_MESSAGE("Constructing cell to faces data... ");
  {
    const int& num_cells = DENEB_DATA->GetNumCells();
    const int& num_faces = DENEB_DATA->GetNumFaces();
    const int& num_inner_faces = DENEB_DATA->GetNumInnerFaces();
    const auto& face_owner_cell = DENEB_DATA->GetFaceOwnerCell();
    const auto& face_neighbor_cell = DENEB_DATA->GetFaceNeighborCell();
    cell_faces_.resize(num_cells);
    for (int iface = 0; iface < num_faces; iface++)
      cell_faces_[face_owner_cell[iface]].push_back(iface);
    for (int iface = 0; iface < num_inner_faces; iface++)
      cell_faces_[face_neighbor_cell[iface]].push_back(iface);
  }
  SYNCRO();
  MASTER_MESSAGE("Complete. (Time: " + std::to_string(STOP_TIMER()) + "s)\n");

  target_state_ = 0;
  const auto& variable_names = DENEB_EQUATION->GetCellVariableNames();
  MASTER_MESSAGE("Target state: " + variable_names[target_state_] + "\n");

  const int& order = DENEB_DATA->GetOrder();
  num_bases_list_.resize(order + 1);
  for (int i = 0; i < order + 1; i++)
    num_bases_list_[i] = DENEB_DATA->GetNumBases(i);

  const int& num_nodes = DENEB_DATA->GetNumNodes();
  const std::vector<int>& cell_node_ind = DENEB_DATA->GetCellNodeInd();
  const std::vector<int>& cell_node_ptr = DENEB_DATA->GetCellNodePtr();
  nodetypes_.resize(num_nodes, NODETYPE::NORMAL);
  const int& num_bdries = DENEB_DATA->GetNumBdries();
  const auto& cell_element = DENEB_DATA->GetCellElement();
  const std::vector<int>& bdry_owner_cell = DENEB_DATA->GetBdryOwnerCell();
  const std::vector<int>& bdry_owner_type = DENEB_DATA->GetBdryOwnerType();
  for (int ibdry = 0; ibdry < num_bdries; ibdry++) {
    const int& owner_cell = bdry_owner_cell[ibdry];
    const int& owner_type = bdry_owner_type[ibdry];

    std::vector<int> bdry_nodes =
        cell_element[owner_cell]->GetFacetypeNodes(1)[owner_type];
    const int& start_index = cell_node_ptr[owner_cell];
    for (auto&& inode : bdry_nodes) inode = cell_node_ind[start_index + inode];
    for (auto&& inode : bdry_nodes) nodetypes_[inode] = NODETYPE::BOUNDARY;
  }

  const int& num_cells = DENEB_DATA->GetNumCells();
  const int& num_total_cells = DENEB_DATA->GetNumTotalCells();
  const int& num_states = DENEB_EQUATION->GetNumStates();
  const int& num_bases = DENEB_DATA->GetNumBases();
  const int sb = num_states * num_bases;
  foreign_solution_.resize(std::max((num_total_cells - num_cells) * sb, 1));
  communicate_ = std::make_shared<avocado::Communicate>(
      sb, DENEB_DATA->GetTotalSendCellList(),
      DENEB_DATA->GetTotalRecvCellList());

  const int& dimension = DENEB_EQUATION->GetDimension();
  const std::vector<double>& node_coords = DENEB_DATA->GetNodeCoords();
  cell_vertex_basis_value_.resize(num_cells);
  for (int icell = 0; icell < num_cells; icell++) {
    const int num_vertex = cell_node_ptr[icell + 1] - cell_node_ptr[icell];

    std::vector<double> coords;
    for (int v = 0; v < num_vertex; v++) {
      const int& node = cell_node_ind[cell_node_ptr[icell] + v];
      const double* coord = &node_coords[node * dimension];
      for (int idim = 0; idim < dimension; idim++)
        coords.push_back(coord[idim]);
    }
    DENEB_DATA->GetCellBasisValues(icell, coords,
                                   cell_vertex_basis_value_[icell]);
  }

  cell_num_for_type_one_.resize(num_cells);
  cell_num_for_type_two_.resize(num_cells);
  const auto& cell_elemtype = DENEB_DATA->GetCellElemtype();
  for (int icell = 0; icell < num_cells; icell++) {
    if (cell_elemtype[icell] == ElemType::TRIS) {
      cell_num_for_type_one_[icell] = 2;
      cell_num_for_type_two_[icell] = 1;
    } else if (cell_elemtype[icell] == ElemType::QUAD) {
      cell_num_for_type_one_[icell] = 2;
      cell_num_for_type_two_[icell] = 1;
    } else if (cell_elemtype[icell] == ElemType::TETS) {
      cell_num_for_type_one_[icell] = 2;
      cell_num_for_type_two_[icell] = 1;
    } else if (cell_elemtype[icell] == ElemType::HEXA) {
      cell_num_for_type_one_[icell] = 2;
      cell_num_for_type_two_[icell] = 1;
    } else if (cell_elemtype[icell] == ElemType::PRIS) {
      cell_num_for_type_one_[icell] = 2;
      cell_num_for_type_two_[icell] = 1;
    } else if (cell_elemtype[icell] == ElemType::PYRA) {
      cell_num_for_type_one_[icell] = 2;
      cell_num_for_type_two_[icell] = 1;
    } else
      ERROR_MESSAGE("(hMLPBD) ElemType error\n");
  }

  const int& num_faces = DENEB_DATA->GetNumFaces();
  const int db = dimension * num_bases;
  const std::vector<int>& num_face_points = DENEB_DATA->GetNumFacePoints();
  const std::vector<int>& face_owner_cell = DENEB_DATA->GetFaceOwnerCell();
  const auto& cell_volumes = DENEB_DATA->GetCellVolumes();
  const auto& face_normals = DENEB_DATA->GetFaceNormals();
  const auto& face_owner_coefficients =
      DENEB_DATA->GetFaceOwnerCoefficients();  // check here
  std::vector<std::vector<double>> face_coefficients(num_faces);
  face_areas_.resize(num_faces);
  face_characteristic_length_.resize(num_faces);
  for (int iface = 0; iface < num_faces; iface++) {
    const int& num_points = num_face_points[iface];
    const int& owner_cell = face_owner_cell[iface];

    face_coefficients[iface].resize(num_points);
    for (int ipoint = 0; ipoint < num_points; ipoint++)
      face_coefficients[iface][ipoint] =
          cblas_ddot(dimension, &face_owner_coefficients[iface][ipoint * db],
                     num_bases, &face_normals[iface][ipoint * dimension], 1);
    cblas_dscal(num_points, std::sqrt(cell_volumes[owner_cell]),
                &face_coefficients[iface][0], 1);

    double area = 0.0;
    for (int ipoint = 0; ipoint < num_points; ipoint++)
      area += face_coefficients[iface][ipoint];

    double length = area;
    if (dimension == 3) length = std::sqrt(area);

    face_areas_[iface] = area;
    face_characteristic_length_[iface] = length;
  }
  face_coefficients_ = std::move(face_coefficients);

  START_TIMER();
  MASTER_MESSAGE("Constructing simplex decomposition data... ");
  ConstructSimplex();
  SYNCRO();
  MASTER_MESSAGE("Complete. (Time: " + std::to_string(STOP_TIMER()) + "s)\n");

  SYNCRO();
  STOP_TIMER_TAG("BuildData");
  MASTER_MESSAGE("hMLP_BD::BuildData() completes. (Time: " +
                 std::to_string(GET_RECORD_TAG("BuildData")) + "s)\n");
}
void hMLP_BD::ConstructSimplex(void) {
  const int& dimension = DENEB_EQUATION->GetDimension();
  ElemType simplex_type;
  if (dimension == 2)
    simplex_type = ElemType::TRIS;
  else if (dimension == 3)
    simplex_type = ElemType::TETS;
  const Element* simplex_element = DENEB_ELEMENT->GetElement(simplex_type);
  const int num_basis_p1 = dimension + 1;

  const int& num_total_cells = DENEB_DATA->GetNumTotalCells();
  simplex_average_coefficients_.resize(num_total_cells);
  simplex_p1_coefficients_.resize(num_total_cells);
  simplex_volume_.resize(num_total_cells);

  const auto& order = DENEB_DATA->GetOrder();
  const auto& num_bases = DENEB_DATA->GetNumBases();
  const auto& node_coords = DENEB_DATA->GetNodeCoords();
  const auto& cell_element = DENEB_DATA->GetCellElement();
  const auto& cell_elem_type = DENEB_DATA->GetCellElemtype();
  const auto& cell_basis = DENEB_DATA->GetCellBasis();
  const auto& cell_volumes = DENEB_DATA->GetCellVolumes();
  const std::vector<int>& cell_node_ind = DENEB_DATA->GetCellNodeInd();
  const std::vector<int>& cell_node_ptr = DENEB_DATA->GetCellNodePtr();
  const std::vector<int>& cell_subnode_ind = DENEB_DATA->GetCellSubnodeInd();
  const std::vector<int>& cell_subnode_ptr = DENEB_DATA->GetCellSubnodePtr();
  for (int icell = 0; icell < num_total_cells; icell++) {
    int num_vertex = cell_element[icell]->GetNumNodes(1);
    if (static_cast<ElemType>(cell_elem_type[icell]) == ElemType::PYRA)
      num_vertex -= 1;

    std::vector<double> vertex_coords(num_vertex * dimension);
    avocado::VecCopy(num_vertex, &cell_node_ind[cell_node_ptr[icell]],
                     &node_coords[0], dimension, &vertex_coords[0], dimension);
    std::vector<double> cell_vertex_basis_value(num_vertex * num_bases);
    cell_basis[icell]->GetBasis(num_vertex, &vertex_coords[0],
                                &cell_vertex_basis_value[0]);

    simplex_average_coefficients_[icell].resize(num_vertex * num_bases, 0.0);
    simplex_p1_coefficients_[icell].resize(num_vertex * num_bases, 0.0);
    simplex_volume_[icell].resize(num_vertex, 0.0);

    if (static_cast<ElemType>(cell_elem_type[icell]) == ElemType::TRIS ||
        static_cast<ElemType>(cell_elem_type[icell]) == ElemType::TETS) {
      for (int ivertex = 0; ivertex < num_vertex; ivertex++)
        simplex_volume_[icell][ivertex] = cell_volumes[icell];
      for (int ivertex = 0; ivertex < num_vertex; ivertex++) {
        simplex_average_coefficients_[icell][ivertex * num_bases] =
            cell_vertex_basis_value[ivertex * num_bases];
        for (int ibasis = 0; ibasis < num_basis_p1; ibasis++)
          simplex_p1_coefficients_[icell][ivertex * num_bases + ibasis] =
              cell_vertex_basis_value[ivertex * num_bases + ibasis];
      }

      continue;
    }

    std::vector<std::vector<int>> simplex_nodes_index;
    switch (cell_elem_type[icell]) {
      case ElemType::QUAD:
        simplex_nodes_index = SED_simplex_node_ind_[0];
        break;
      case ElemType::HEXA:
        simplex_nodes_index = SED_simplex_node_ind_[1];
        break;
      case ElemType::PRIS:
        simplex_nodes_index = SED_simplex_node_ind_[2];
        break;
      case ElemType::PYRA:
        simplex_nodes_index = SED_simplex_node_ind_[3];
        break;
    }

    for (int ivertex = 0; ivertex < num_vertex; ivertex++) {
      const int num_nodes = simplex_nodes_index[ivertex].size();

      std::vector<double> cell_coords;
      const int start_index = cell_subnode_ptr[icell];
      cell_coords.reserve(dimension * num_nodes);
      for (int inode = 0; inode < num_nodes; inode++)
        for (int idim = 0; idim < dimension; idim++)
          cell_coords.push_back(
              node_coords[cell_subnode_ind[start_index + simplex_nodes_index
                                                             [ivertex][inode]] *
                              dimension +
                          idim]);

      std::shared_ptr<Jacobian> jacobian = Jacobian::GetJacobian(simplex_type);
      jacobian->SetTopology(1, &cell_coords[0]);

      double volume = 0.0;
      std::vector<double> quad_points, quad_weights;
      {
        std::shared_ptr<Basis> basis = Basis::GetStandardBasis(simplex_type);
        std::vector<double> ref_points;
        basis->GetBasisPolynomial()->GetVolumeQuadrature(
            2 * order, simplex_type, 1, ref_points, quad_weights);
        const int num_points = static_cast<int>(quad_weights.size());
        std::vector<double> jacobian_det(num_points);
        jacobian->CalJacobianDet(num_points, &ref_points[0], &jacobian_det[0]);
        for (int ipoint = 0; ipoint < num_points; ipoint++)
          quad_weights[ipoint] *= std::abs(jacobian_det[ipoint]);
        quad_points.resize(num_points * dimension);
        jacobian->TransformToPhyCoords(num_points, &ref_points[0],
                                       &quad_points[0]);
        for (int ipoint = 0; ipoint < num_points; ipoint++)
          volume += quad_weights[ipoint];
      }
      simplex_volume_[icell][ivertex] = volume;

      // std::vector<double> center(dimension, 0.0);
      // for (int inode = 0; inode < num_nodes; inode++)
      //	for (int idim = 0; idim < dimension; idim++)
      //		center[idim] += cell_coords[inode*dimension + idim];
      // for (auto&& var : center)
      //	var /= double(num_nodes);

      std::shared_ptr<Basis> basis = Basis::GetStandardBasis(simplex_type);
      basis->SetTransform(cell_coords);
      basis->ComputeConnectingMatrix(1, quad_points, quad_weights);

      const int num_points = quad_weights.size();
      std::vector<double> cell_basis_value(num_points * num_bases);
      std::vector<double> simplex_basis_value(num_points * num_basis_p1);
      cell_basis[icell]->GetBasis(num_points, &quad_points[0],
                                  &cell_basis_value[0]);
      basis->GetBasis(num_points, &quad_points[0], &simplex_basis_value[0]);

      std::vector<double> transition_matrix(num_bases * num_basis_p1, 0.0);
      for (int jbasis = 0; jbasis < num_basis_p1; jbasis++)
        for (int ibasis = 0; ibasis < num_bases; ibasis++)
          for (int ipoint = 0; ipoint < quad_weights.size(); ipoint++)
            transition_matrix[jbasis * num_bases + ibasis] +=
                cell_basis_value[ipoint * num_bases + ibasis] *
                simplex_basis_value[ipoint * num_basis_p1 + jbasis] *
                quad_weights[ipoint];

      std::vector<double> com_coord(
          simplex_element->GetNodesCoord(1).begin(),
          simplex_element->GetNodesCoord(1).begin() + dimension);
      std::vector<double> phy_coord(dimension);
      jacobian->TransformToPhyCoords(1, &com_coord[0], &phy_coord[0]);

      std::vector<double> basis_vertex_value(1 * num_basis_p1);
      basis->GetBasis(1, &phy_coord[0], &basis_vertex_value[0]);
      for (int ibasis = 0; ibasis < num_bases; ibasis++)
        simplex_average_coefficients_[icell][ivertex * num_bases + ibasis] =
            transition_matrix[ibasis] * basis_vertex_value[0];
      for (int jbasis = 0; jbasis < num_basis_p1; jbasis++)
        for (int ibasis = 0; ibasis < num_bases; ibasis++)
          simplex_p1_coefficients_[icell][ivertex * num_bases + ibasis] +=
              transition_matrix[jbasis * num_bases + ibasis] *
              basis_vertex_value[jbasis];
    }
  }
}
void hMLP_BD::ConstructCellCells(void) {
  const std::vector<int>& face_owner_cell = DENEB_DATA->GetFaceOwnerCell();
  const std::vector<int>& face_neighbor_cell =
      DENEB_DATA->GetFaceNeighborCell();

  const int& num_cells = DENEB_DATA->GetNumCells();
  const int& num_faces = DENEB_DATA->GetNumFaces();
  const int& num_inner_faces = DENEB_DATA->GetNumInnerFaces();
  std::vector<std::vector<int>> cell_cells(num_cells);
  for (int i = 0; i < num_inner_faces; i++) {
    const int& owner_cell = face_owner_cell[i];
    const int& neighbor_cell = face_neighbor_cell[i];

    cell_cells[owner_cell].push_back(neighbor_cell);
    cell_cells[neighbor_cell].push_back(owner_cell);
  }
  for (int i = num_inner_faces; i < num_faces; i++) {
    const int& owner_cell = face_owner_cell[i];
    const int& neighbor_cell = face_neighbor_cell[i];

    cell_cells[owner_cell].push_back(neighbor_cell + num_cells);
  }

  cell_cells_ = std::move(cell_cells);
}
void hMLP_BD::ConstructNodeCells(void) {
  const std::vector<int>& node_global_index = DENEB_DATA->GetNodeGlobalIndex();
  std::unordered_map<int, int> node_mapping;  // global to local index
  {
    int ind = 0;
    for (auto&& global_index : node_global_index)
      node_mapping[global_index] = ind++;
  }

  const std::unordered_map<int, std::vector<int>>&
      periodic_matching_global_node_index =
          DENEB_DATA->GetPeriodicMatchingGlobalNodeIndex();
  const auto& node_grps = periodic_matching_global_node_index;

  const int& num_nodes = DENEB_DATA->GetNumNodes();
  const int& num_cells = DENEB_DATA->GetNumCells();
  const int& num_total_cells = DENEB_DATA->GetNumTotalCells();
  const std::vector<int>& cell_node_ptr = DENEB_DATA->GetCellNodePtr();
  const std::vector<int>& cell_node_ind = DENEB_DATA->GetCellNodeInd();
  std::vector<std::vector<int>> node_cells;
  {
    std::unordered_map<int, std::unordered_set<int>> node_cells_temp;
    for (int i = 0; i < num_total_cells; i++)
      for (int ptr = cell_node_ptr[i]; ptr < cell_node_ptr[i + 1]; ptr++)
        node_cells_temp[cell_node_ind[ptr]].insert(i);

    for (auto&& iterator : node_grps) {
      const int& mynode = node_mapping.at(iterator.first);

      std::unordered_set<int> new_cells;
      for (auto&& node : iterator.second) {
        const auto& cells = node_cells_temp.find(node_mapping.at(node));
        if (cells == node_cells_temp.end()) continue;
        new_cells.insert(cells->second.begin(), cells->second.end());
      }
      node_cells_temp[mynode].insert(new_cells.begin(), new_cells.end());
    }

    node_cells.resize(num_nodes);
    for (auto&& iterator : node_cells_temp) {
      if (iterator.first >= num_nodes) continue;
      node_cells[iterator.first] =
          std::vector<int>(iterator.second.begin(), iterator.second.end());
    }
  }

  std::vector<std::vector<int>> node_vertices(num_nodes);
  for (int i = 0; i < num_nodes; i++) {
    const auto& matching_nodes = node_grps.find(node_global_index[i]);
    for (auto&& cell : node_cells[i]) {
      for (int ptr = cell_node_ptr[cell]; ptr < cell_node_ptr[cell + 1];
           ptr++) {
        if (cell_node_ind[ptr] == i) {
          node_vertices[i].push_back(ptr - cell_node_ptr[cell]);
          break;
        }
        if (matching_nodes == node_grps.end()) continue;
        const int& node = node_global_index[cell_node_ind[ptr]];
        if (std::find(matching_nodes->second.begin(),
                      matching_nodes->second.end(),
                      node) != matching_nodes->second.end()) {
          node_vertices[i].push_back(ptr - cell_node_ptr[cell]);
          break;
        }
      }
    }
  }

  node_cells_ = std::move(node_cells);
  node_vertices_ = std::move(node_vertices);
}
void hMLP_BD::Limiting(double* solution_ptr) {
  static const int& order = DENEB_DATA->GetOrder();
  if (order == 0) return;

  static const int& num_cells = DENEB_DATA->GetNumCells();

  VertexMinMax(solution_ptr);
  for (int icell = 0; icell < num_cells; icell++) {
    bool is_trouble = true;
    int temp_order = order;

    while (true) {
      if (temp_order == 0) {
        is_trouble = false;
        break;
      }

      is_trouble = Indicator(solution_ptr, icell, temp_order);

      if (is_trouble == true) {
        if (temp_order == 1)
          break;
        else if (temp_order == 2) {
          temp_order = 1;
          Projection(solution_ptr, icell, 1);
          break;
        } else {
          is_trouble = false;
          temp_order = temp_order - 1;
          Projection(solution_ptr, icell, temp_order);
        }
      } else {
        if (BoundaryDetect(icell, cell_num_for_type_one_[icell])) {
          if (DENEB_EQUATION->IsContact(icell, cell_cells_[icell], solution_ptr,
                                        &foreign_solution_[0]) == false) {
            temp_order = 1;
            Projection(solution_ptr, icell, 1);  // type one boundary detector
          }
        }
        break;
      }
    }
    if (temp_order == 1 && is_trouble == true) {
      MLPu1(solution_ptr, icell);
    }
  }

  return;
}
void hMLP_BD::VertexMinMax(const double* solution_ptr) {
  static const int& num_total_cells = DENEB_DATA->GetNumTotalCells();
  static const int& num_cells = DENEB_DATA->GetNumCells();
  static const int& num_bases = DENEB_DATA->GetNumBases();
  static const int& num_nodes = DENEB_DATA->GetNumNodes();
  static const int& num_states = DENEB_EQUATION->GetNumStates();
  static const int sb = num_states * num_bases;

  const std::vector<int>& cell_node_ind = DENEB_DATA->GetCellNodeInd();
  const std::vector<int>& cell_node_ptr = DENEB_DATA->GetCellNodePtr();

  // hMLP_BD_sim
  communicate_->CommunicateBegin(solution_ptr);
  std::vector<std::vector<double>> cell_vertex_averages(num_total_cells);
  for (int icell = 0; icell < num_cells; icell++) {
    const int num_vertex = cell_node_ptr[icell + 1] - cell_node_ptr[icell];
    cell_vertex_averages[icell].resize(num_vertex * num_states, 0.0);
    for (int ivertex = 0; ivertex < num_vertex; ivertex++)
      for (int istate = 0; istate < num_states; istate++)
        for (int ibasis = 0; ibasis < num_bases; ibasis++)
          cell_vertex_averages[icell][ivertex * num_states + istate] +=
              solution_ptr[icell * sb + istate * num_bases + ibasis] *
              simplex_average_coefficients_[icell]
                                           [ivertex * num_bases + ibasis];
  }

  communicate_->CommunicateEnd(&foreign_solution_[0]);
  for (int icell = num_cells; icell < num_total_cells; icell++) {
    const int num_vertex = cell_node_ptr[icell + 1] - cell_node_ptr[icell];
    cell_vertex_averages[icell].resize(num_vertex * num_states, 0.0);
    for (int ivertex = 0; ivertex < num_vertex; ivertex++)
      for (int istate = 0; istate < num_states; istate++)
        for (int ibasis = 0; ibasis < num_bases; ibasis++)
          cell_vertex_averages[icell][ivertex * num_states + istate] +=
              foreign_solution_[(icell - num_cells) * sb + istate * num_bases +
                                ibasis] *
              simplex_average_coefficients_[icell]
                                           [ivertex * num_bases + ibasis];
  }

  std::vector<double> vertex_max(num_states * num_nodes);
  std::vector<double> vertex_min(num_states * num_nodes);
  std::vector<double> vertex_average(num_states * num_nodes);
  const auto& cell_elem_type = DENEB_DATA->GetCellElemtype();
  for (int inode = 0; inode < num_nodes; inode++) {
    std::vector<std::vector<double>> local_cell_averages(num_states);
    std::vector<double> state_sum(num_states, 0.0);
    double area_sum = 0.0;

    for (int icell = 0; icell < node_cells_[inode].size(); icell++) {
      const int cell_index = node_cells_[inode][icell];

      for (int istate = 0; istate < num_states; istate++) {
        switch (static_cast<ElemType>(cell_elem_type[cell_index])) {
          case ElemType::TRIS:
            local_cell_averages[istate].push_back(
                cell_vertex_averages[cell_index][istate]);
            state_sum[istate] += cell_vertex_averages[cell_index][istate] *
                                 simplex_volume_[cell_index][0];
            if (istate == 0) area_sum += simplex_volume_[cell_index][0];
            break;
          case ElemType::QUAD:
            local_cell_averages[istate].push_back(
                cell_vertex_averages[cell_index]
                                    [node_vertices_[inode][icell] * num_states +
                                     istate]);
            state_sum[istate] +=
                cell_vertex_averages[cell_index]
                                    [node_vertices_[inode][icell] * num_states +
                                     istate] *
                simplex_volume_[cell_index][node_vertices_[inode][icell]];
            if (istate == 0)
              area_sum +=
                  simplex_volume_[cell_index][node_vertices_[inode][icell]];
            break;
          case ElemType::TETS:
            local_cell_averages[istate].push_back(
                cell_vertex_averages[cell_index][istate]);
            state_sum[istate] += cell_vertex_averages[cell_index][istate] *
                                 simplex_volume_[cell_index][0];
            if (istate == 0) area_sum += simplex_volume_[cell_index][0];
            break;
          case ElemType::HEXA:
            local_cell_averages[istate].push_back(
                cell_vertex_averages[cell_index]
                                    [node_vertices_[inode][icell] * num_states +
                                     istate]);
            state_sum[istate] +=
                cell_vertex_averages[cell_index]
                                    [node_vertices_[inode][icell] * num_states +
                                     istate] *
                simplex_volume_[cell_index][node_vertices_[inode][icell]];
            if (istate == 0)
              area_sum +=
                  simplex_volume_[cell_index][node_vertices_[inode][icell]];
            break;
          case ElemType::PRIS:
            local_cell_averages[istate].push_back(
                cell_vertex_averages[cell_index]
                                    [node_vertices_[inode][icell] * num_states +
                                     istate]);
            state_sum[istate] +=
                cell_vertex_averages[cell_index]
                                    [node_vertices_[inode][icell] * num_states +
                                     istate] *
                simplex_volume_[cell_index][node_vertices_[inode][icell]];
            if (istate == 0)
              area_sum +=
                  simplex_volume_[cell_index][node_vertices_[inode][icell]];
            break;
          case ElemType::PYRA:
            switch (node_vertices_[inode][icell]) {
              case 0:
              case 1:
              case 2:
              case 3:
                local_cell_averages[istate].push_back(
                    cell_vertex_averages[cell_index]
                                        [node_vertices_[inode][icell] *
                                             num_states +
                                         istate]);
                state_sum[istate] +=
                    cell_vertex_averages[cell_index]
                                        [node_vertices_[inode][icell] *
                                             num_states +
                                         istate] *
                    simplex_volume_[cell_index][node_vertices_[inode][icell]];
                if (istate == 0)
                  area_sum +=
                      simplex_volume_[cell_index][node_vertices_[inode][icell]];
                break;
              case 4:
                local_cell_averages[istate].push_back(
                    cell_vertex_averages[cell_index][istate]);
                local_cell_averages[istate].push_back(
                    cell_vertex_averages[cell_index][2 * num_states + istate]);
                state_sum[istate] += cell_vertex_averages[cell_index][istate] *
                                     simplex_volume_[cell_index][0];
                state_sum[istate] +=
                    cell_vertex_averages[cell_index][2 * num_states + istate] *
                    simplex_volume_[cell_index][2];
                if (istate == 0) {
                  area_sum += simplex_volume_[cell_index][0];
                  area_sum += simplex_volume_[cell_index][2];
                }
                break;
            }
            break;
        }
      }
    }

    for (int istate = 0; istate < num_states; istate++) {
      vertex_max[istate * num_nodes + inode] =
          *max_element(local_cell_averages[istate].begin(),
                       local_cell_averages[istate].end());
      vertex_min[istate * num_nodes + inode] =
          *min_element(local_cell_averages[istate].begin(),
                       local_cell_averages[istate].end());
      vertex_average[istate * num_nodes + inode] = state_sum[istate] / area_sum;
    }
  }

  vertex_min_ = move(vertex_min);
  vertex_max_ = move(vertex_max);
  vertex_average_ = move(vertex_average);

  // face sweep
  static const int& order = DENEB_DATA->GetOrder();
  static const int& num_faces = DENEB_DATA->GetNumFaces();
  static const int& num_inner_faces = DENEB_DATA->GetNumInnerFaces();
  static const std::vector<int>& num_face_points =
      DENEB_DATA->GetNumFacePoints();
  static const std::vector<int>& face_owner_cell =
      DENEB_DATA->GetFaceOwnerCell();
  static const std::vector<int>& face_neighbor_cell =
      DENEB_DATA->GetFaceNeighborCell();
  static const std::vector<std::vector<double>>& face_owner_basis_value =
      DENEB_DATA->GetFaceOwnerBasisValue();
  static const std::vector<std::vector<double>>& face_neighbor_basis_value =
      DENEB_DATA->GetFaceNeighborBasisValue();
  std::vector<std::vector<double>> face_difference(num_states);
  for (int iface = 0; iface < num_inner_faces; iface++) {
    const int num_points = num_face_points[iface];
    const int owner_cell = face_owner_cell[iface];
    const int neighbor_cell = face_neighbor_cell[iface];

    std::vector<double> owner_solution(num_points * num_states, 0.0);
    std::vector<double> neighbor_solution(num_points * num_states, 0.0);
    std::vector<double> state_sum(num_states, 0.0);
    for (int ipoint = 0; ipoint < num_points; ipoint++)
      for (int istate = 0; istate < num_states; istate++) {
        for (int ibasis = 0; ibasis < num_bases; ibasis++) {
          owner_solution[ipoint * num_states + istate] +=
              face_owner_basis_value[iface][ipoint * num_bases + ibasis] *
              solution_ptr[owner_cell * sb + istate * num_bases + ibasis];
          neighbor_solution[ipoint * num_states + istate] +=
              face_neighbor_basis_value[iface][ipoint * num_bases + ibasis] *
              solution_ptr[neighbor_cell * sb + istate * num_bases + ibasis];
        }
        state_sum[istate] +=
            std::abs(owner_solution[ipoint * num_states + istate] -
                     neighbor_solution[ipoint * num_states + istate]) *
            face_coefficients_[iface][ipoint];
      }
    const double denominator = std::pow(face_characteristic_length_[iface],
                                        0.5 * static_cast<double>(order + 1));
    for (int istate = 0; istate < num_states; istate++) {
      state_sum[istate] /= (face_areas_[iface] * denominator);
      face_difference[istate].push_back(state_sum[istate]);
    }
  }
  for (int iface = num_inner_faces; iface < num_faces; iface++) {
    const int num_points = num_face_points[iface];
    const int owner_cell = face_owner_cell[iface];
    const int neighbor_cell = face_neighbor_cell[iface];

    std::vector<double> owner_solution(num_points * num_states, 0.0);
    std::vector<double> neighbor_solution(num_points * num_states, 0.0);
    std::vector<double> state_sum(num_states, 0.0);
    for (int ipoint = 0; ipoint < num_points; ipoint++)
      for (int istate = 0; istate < num_states; istate++) {
        for (int ibasis = 0; ibasis < num_bases; ibasis++) {
          owner_solution[ipoint * num_states + istate] +=
              face_owner_basis_value[iface][ipoint * num_bases + ibasis] *
              solution_ptr[owner_cell * sb + istate * num_bases + ibasis];
          neighbor_solution[ipoint * num_states + istate] +=
              face_neighbor_basis_value[iface][ipoint * num_bases + ibasis] *
              foreign_solution_[neighbor_cell * sb + istate * num_bases +
                                ibasis];
        }
        state_sum[istate] +=
            std::abs(owner_solution[ipoint * num_states + istate] -
                     neighbor_solution[ipoint * num_states + istate]) *
            face_coefficients_[iface][ipoint];
      }
    const double denominator = std::pow(face_characteristic_length_[iface],
                                        0.5 * static_cast<double>(order + 1));
    for (int istate = 0; istate < num_states; istate++) {
      state_sum[istate] /= (face_areas_[iface] * denominator);
      face_difference[istate].push_back(state_sum[istate]);
    }
  }
  face_difference_ = move(face_difference);

  std::vector<int> cell_num_troubled_boundary(num_cells, 0);
  for (int icell = 0; icell < num_cells; icell++) {
    const int num_faces = cell_faces_[icell].size();
    for (int iface = 0; iface < num_faces; iface++) {
      for (int istate = 0; istate < num_states; istate++) {
        if (face_difference_[istate][cell_faces_[icell][iface]] > 1.0) {
          cell_num_troubled_boundary[icell]++;
          break;
        }
      }
    }
  }
  cell_num_troubled_boundary_ = move(cell_num_troubled_boundary);
}
bool hMLP_BD::Indicator(const double* solution_ptr, const int icell,
                        const int order) {
  static const int& num_states = DENEB_EQUATION->GetNumStates();
  static const int& num_basis = DENEB_DATA->GetNumBases();
  static const int& num_cells = DENEB_DATA->GetNumCells();
  static const int sb = num_states * num_basis;
  static const std::vector<int>& cell_node_ptr = DENEB_DATA->GetCellNodePtr();
  static const std::vector<int>& cell_node_ind = DENEB_DATA->GetCellNodeInd();

  const int num_vertex = cell_node_ptr[icell + 1] - cell_node_ptr[icell];
  const double average_value =
      solution_ptr[icell * sb + target_state_ * num_basis] *
      cell_vertex_basis_value_[icell][0];

  for (int ivertex = 0; ivertex < num_vertex; ivertex++) {
    const int& node = cell_node_ind[cell_node_ptr[icell] + ivertex];
    if (nodetypes_[node] == NODETYPE::BOUNDARY) continue;

    const double vertex_value = cblas_ddot(
        num_basis, &solution_ptr[icell * sb + target_state_ * num_basis], 1,
        &cell_vertex_basis_value_[icell][ivertex * num_basis], 1);
    const double P1var =
        cblas_ddot(num_bases_list_[1],
                   &solution_ptr[icell * sb + target_state_ * num_basis], 1,
                   &cell_vertex_basis_value_[icell][ivertex * num_basis], 1);

    bool flag = MLPcondition(solution_ptr, icell, ivertex, target_state_, P1var,
                             vertex_value, average_value);
    if (flag == true) {
      if (order > 1) {
        flag = SmoothDetect(solution_ptr, icell, ivertex, target_state_, P1var,
                            vertex_value, average_value);
        if (flag == false)
          if (BoundaryDetect(icell, cell_num_for_type_two_[icell]))
            flag = true;  // type two boundary detector
      }
    }
    if (flag == true) {
      return true;
    }
  }
  return false;
}
bool hMLP_BD::MLPcondition(const double* solution_ptr, const int icell,
                           const int ivertex, const int istate,
                           const double P1var, const double vertex_value,
                           const double average_value) {
  static const int& num_nodes = DENEB_DATA->GetNumNodes();
  static const std::vector<int>& cell_node_ptr = DENEB_DATA->GetCellNodePtr();
  static const std::vector<int>& cell_node_ind = DENEB_DATA->GetCellNodeInd();
  static const std::vector<double>& cell_volumes = DENEB_DATA->GetCellVolumes();

  if (Deactivation(solution_ptr, cell_volumes[icell], vertex_value,
                   average_value))
    return false;

  if (P1var > vertex_max_[istate * num_nodes +
                          cell_node_ind[cell_node_ptr[icell] + ivertex]])
    return true;
  else if (P1var < vertex_min_[istate * num_nodes +
                               cell_node_ind[cell_node_ptr[icell] + ivertex]])
    return true;
  else
    return false;
}
bool hMLP_BD::Deactivation(const double* solution_ptr, const double volume,
                           const double vertex_value,
                           const double average_value) {
  if (std::abs(vertex_value - average_value) <=
      std::max(deact_eps_ * std::abs(average_value), volume))
    return true;
  return false;
}
bool hMLP_BD::SmoothDetect(const double* solution_ptr, const int icell,
                           const int ivertex, const int istate,
                           const double P1var, const double vertex_value,
                           const double average_value) {
  static const int num_nodes = DENEB_DATA->GetNumNodes();
  static const std::vector<int>& cell_node_ptr = DENEB_DATA->GetCellNodePtr();
  static const std::vector<int>& cell_node_ind = DENEB_DATA->GetCellNodeInd();

  const double high_term = vertex_value - P1var;
  const double P1term = P1var - average_value;

  if ((P1term > 0.0) && (high_term < eps_) &&
      (vertex_value > vertex_min_[istate * num_nodes +
                                  cell_node_ind[cell_node_ptr[icell] +
                                                ivertex]]))  // c1 condition
    return false;
  else if ((P1term < 0.0) && (high_term > -eps_) &&
           (vertex_value <
            vertex_max_[istate * num_nodes +
                        cell_node_ind[cell_node_ptr[icell] +
                                      ivertex]]))  // c2 condition
    return false;

  return true;
}
void hMLP_BD::Projection(double* solution_ptr, const int icell,
                         const int iorder) {
  static const int& num_states = DENEB_EQUATION->GetNumStates();
  static const int& num_bases = DENEB_DATA->GetNumBases();
  static const int sb = num_states * num_bases;

  for (int istate = 0; istate < num_states; istate++)
    for (int ibasis = num_bases_list_[iorder]; ibasis < num_bases; ibasis++)
      solution_ptr[icell * sb + istate * num_bases + ibasis] = 0.0;
}
void hMLP_BD::MLPu1(double* solution_ptr, const int icell) {
  static const int& num_states = DENEB_EQUATION->GetNumStates();
  static const int& num_bases = DENEB_DATA->GetNumBases();
  static const int& num_nodes = DENEB_DATA->GetNumNodes();
  static const int sb = num_states * num_bases;
  static const std::vector<int>& cell_node_ptr = DENEB_DATA->GetCellNodePtr();
  static const std::vector<int>& cell_node_ind = DENEB_DATA->GetCellNodeInd();

  const int num_vertex = cell_node_ptr[icell + 1] - cell_node_ptr[icell];

  std::vector<double> p1_value(num_states * num_vertex, 0.0);
  for (int istate = 0; istate < num_states; istate++)
    for (int ivertex = 0; ivertex < num_vertex; ivertex++)
      p1_value[istate * num_vertex + ivertex] = cblas_ddot(
          num_bases_list_[1], &solution_ptr[icell * sb + istate * num_bases], 1,
          &cell_vertex_basis_value_[icell][ivertex * num_bases], 1);

  for (int istate = 0; istate < num_states; istate++) {
    double temp, phi = 1.0;
    for (int ivertex = 0; ivertex < num_vertex; ivertex++) {
      const int& node = cell_node_ind[cell_node_ptr[icell] + ivertex];
      if (nodetypes_[node] == NODETYPE::BOUNDARY) continue;

      temp = solution_ptr[icell * sb + istate * num_bases] *
                 cell_vertex_basis_value_[icell][0] -
             p1_value[istate * num_vertex + ivertex];
      if (p1_value[istate * num_vertex + ivertex] <
          vertex_min_[istate * num_nodes +
                      cell_node_ind[cell_node_ptr[icell] + ivertex]]) {
        if (std::abs(temp) < 1.0e-14)
          temp = 0.0;
        else
          temp = (solution_ptr[icell * sb + istate * num_bases] *
                      cell_vertex_basis_value_[icell][0] -
                  vertex_min_[istate * num_nodes +
                              cell_node_ind[cell_node_ptr[icell] + ivertex]]) /
                 temp;
        if (temp < phi) phi = temp;
      } else if (vertex_max_[istate * num_nodes +
                             cell_node_ind[cell_node_ptr[icell] + ivertex]] <
                 p1_value[istate * num_vertex + ivertex]) {
        if (std::abs(temp) < 1.0e-14)
          temp = 0.0;
        else
          temp = (solution_ptr[icell * sb + istate * num_bases] *
                      cell_vertex_basis_value_[icell][0] -
                  vertex_max_[istate * num_nodes +
                              cell_node_ind[cell_node_ptr[icell] + ivertex]]) /
                 temp;
        if (temp < phi) phi = temp;
      }
    }
    cblas_dscal(num_bases_list_[1] - 1, phi,
                &solution_ptr[icell * sb + istate * num_bases + 1], 1);
  }
}
bool hMLP_BD::BoundaryDetect(const int icell, const int num_troubled_boundary) {
  if (cell_num_troubled_boundary_[icell] >= num_troubled_boundary) return true;
  return false;
}
}  // namespace deneb