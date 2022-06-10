#include "deneb_grid_builder.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <memory>

#include "deneb_config_macro.h"
#include "deneb_equation.h"
#include "deneb_grid_reader.h"

#include "parmetis.h"

#undef MPI_Allgather
#define MOVE_DATA(data) data##_ = std::move(data)

namespace deneb {
void GridBuilder::BuildGrid(void) {
  SYNCRO();
  START_TIMER_TAG("BuildGrid");
  MASTER_MESSAGE(avocado::GetTitle("GridBuilder::BuildGrid()"));

  const std::string& grid_format =
      AVOCADO_CONFIG->GetConfigValue(GRID_FILE_FORMAT);
  std::string grid_file_path = AVOCADO_CONFIG->GetConfigValue(GRID_FILE_PATH);
  avocado::ToAbsolutePath(grid_file_path);
  std::shared_ptr<GridReader> grid_reader =
      GridReader::GetGridReader(grid_format);

  START_TIMER();
  grid_reader->Open(dimension_, grid_file_path);
  num_global_nodes_ = grid_reader->GetNumGlobalNodes();
  num_global_cells_ = grid_reader->GetNumGlobalCells();
  num_global_bdries_ = grid_reader->GetNumGlobalBdries();
  num_global_peribdries_ = grid_reader->GetNumGlobalPeribdries();
  DENEB_EQUATION->RegistBoundary(grid_reader->GetAllBdryTag());
  if (num_global_cells_ < 5 * NDOMAIN)
    ERROR_MESSAGE(
        "Too many processors for the given number of cells\n\tNumber of "
        "processors: " +
        std::to_string(NDOMAIN) +
        "\n\tNumber of cells: " + std::to_string(num_global_cells_) + "\n");

  SYNCRO();
  START_TIMER();
  MASTER_MESSAGE("Loading cell data in parallel... ");
  std::vector<int> cell_dist(NDOMAIN + 1, 0);
  {
    const int remain = num_global_cells_ % NDOMAIN;
    const int num_local_cells = (num_global_cells_ - remain) / NDOMAIN;
    for (int idomain = 0; idomain < remain; idomain++)
      cell_dist[idomain + 1] = cell_dist[idomain] + (num_local_cells + 1);
    for (int idomain = remain; idomain < NDOMAIN; idomain++)
      cell_dist[idomain + 1] = cell_dist[idomain] + num_local_cells;
  }

  std::vector<ElemType> cell_elemtype;
  std::vector<int> cell_global_index;
  std::vector<int> cell_node_ptr;
  std::vector<int> cell_node_ind;
  grid_reader->ReadCellData(cell_dist[MYRANK], cell_dist[MYRANK + 1],
                            cell_elemtype, cell_global_index, cell_node_ptr,
                            cell_node_ind);
  SYNCRO();
  MASTER_MESSAGE("Complete. (Time: " + std::to_string(STOP_TIMER()) + "s)\n");

  START_TIMER();
  MASTER_MESSAGE("Coloring cell data... ");
  std::vector<int> cell_color;
  SetCellColor(cell_node_ptr, cell_node_ind, cell_color);
  {
    std::vector<std::vector<int>> send_cell_list(NDOMAIN);
    for (size_t i = 0, len = cell_color.size(); i < len; i++)
      send_cell_list[cell_color[i]].push_back(static_cast<int>(i));
    RedistributeCellData(send_cell_list, cell_global_index);
    RedistributeCellData(send_cell_list, cell_elemtype);
    RedistributeCellData(send_cell_list, cell_node_ptr, cell_node_ind);
    send_cell_list.clear();
  }
  cell_color.clear();
  SYNCRO();
  MASTER_MESSAGE("Complete. (Time: " + std::to_string(STOP_TIMER()) + "s)\n");

  START_TIMER();
  MASTER_MESSAGE("Constructing face data... ");
  std::vector<int> face_owner_cell;
  std::vector<int> face_neighbor_cell;
  std::vector<int> face_node_ptr;
  std::vector<int> face_node_ind;
  std::vector<int> outer_face_owner_cell;
  std::vector<int> outer_face_node_ptr;
  std::vector<int> outer_face_node_ind;
  ConstructFaceData(cell_global_index, cell_elemtype, cell_node_ptr,
                    cell_node_ind, face_owner_cell, face_neighbor_cell,
                    face_node_ptr, face_node_ind, outer_face_owner_cell,
                    outer_face_node_ptr, outer_face_node_ind);
  SYNCRO();
  MASTER_MESSAGE("Complete. (Time: " + std::to_string(STOP_TIMER()) + "s)\n");

  START_TIMER();
  MASTER_MESSAGE("Constructing boundary data... ");
  std::vector<int> all_bdry_global_index;
  std::vector<int> all_bdry_tag;
  std::vector<int> all_bdry_node_ptr;
  std::vector<int> all_bdry_node_ind;
  grid_reader->ReadAllBdryData(all_bdry_global_index, all_bdry_tag,
                               all_bdry_node_ptr, all_bdry_node_ind);

  std::vector<int> bdry_global_index;
  std::vector<int> bdry_tag;
  std::vector<int> bdry_owner_cell;
  std::vector<int> bdry_node_ptr;
  std::vector<int> bdry_node_ind;
  std::vector<int> peribdry_global_index;
  std::vector<int> peribdry_tag;
  std::vector<int> peribdry_owner_cell;
  std::vector<int> peribdry_node_ptr;
  std::vector<int> peribdry_node_ind;
  ConstructBdryData(all_bdry_global_index, all_bdry_tag, all_bdry_node_ptr,
                    all_bdry_node_ind, outer_face_owner_cell,
                    outer_face_node_ptr, outer_face_node_ind, bdry_global_index,
                    bdry_tag, bdry_owner_cell, bdry_node_ptr, bdry_node_ind,
                    peribdry_global_index, peribdry_tag, peribdry_owner_cell,
                    peribdry_node_ptr, peribdry_node_ind);
  all_bdry_global_index.clear();
  all_bdry_tag.clear();
  all_bdry_node_ptr.clear();
  all_bdry_node_ind.clear();
  int num_bdries = static_cast<int>(bdry_global_index.size());
  SYNCRO();
  MASTER_MESSAGE("Complete. (Time: " + std::to_string(STOP_TIMER()) + "s)\n");

  START_TIMER();
  MASTER_MESSAGE("Matching periodic boundary... ");
  std::vector<int> all_periodic_nodes;
  std::vector<double> all_periodic_node_coords;
  grid_reader->ReadAllPeriodicData(all_periodic_nodes,
                                   all_periodic_node_coords);

  std::vector<int> peribdry_pair_node_ind;
  PeriodicMatching(all_periodic_nodes, all_periodic_node_coords,
                   peribdry_global_index, peribdry_tag, peribdry_node_ptr,
                   peribdry_node_ind, peribdry_pair_node_ind);
  all_periodic_nodes.clear();
  all_periodic_node_coords.clear();

  std::unordered_map<int, std::vector<int>> periodic_matching_nodes;
  PeriodicMatching(peribdry_node_ind, peribdry_pair_node_ind,
                   periodic_matching_nodes);
  SYNCRO();
  MASTER_MESSAGE("Complete. (Time: " + std::to_string(STOP_TIMER()) + "s)\n");

  START_TIMER();
  MASTER_MESSAGE("Constructing outer grid data... ");
  std::vector<int> outer_cell_global_index;
  std::vector<ElemType> outer_cell_elemtype;
  std::vector<int> outer_cell_node_ptr;
  std::vector<int> outer_cell_node_ind;
  {
    std::unordered_set<int> node_list;
    node_list.insert(outer_face_node_ind.begin(), outer_face_node_ind.end());
    node_list.insert(peribdry_node_ind.begin(), peribdry_node_ind.end());
    for (auto&& nodes : periodic_matching_nodes)
      node_list.insert(nodes.second.begin(), nodes.second.end());
    grid_reader->ReadCellData(node_list, outer_cell_global_index,
                              outer_cell_elemtype, outer_cell_node_ptr,
                              outer_cell_node_ind);
    node_list.clear();
  }

  int num_cells, num_outer_cells, num_total_cells;
  int num_inner_faces, num_faces;
  int num_inner_peribdries, num_peribdries, num_total_peribdries;
  std::vector<int> peribdry_neighbor_cell;
  ConstructOuterGrid(
      num_cells, num_outer_cells, num_total_cells, cell_global_index,
      cell_elemtype, cell_node_ptr, cell_node_ind, outer_cell_global_index,
      outer_cell_elemtype, outer_cell_node_ptr, outer_cell_node_ind,
      num_inner_faces, num_faces, face_owner_cell, face_neighbor_cell,
      face_node_ptr, face_node_ind, outer_face_owner_cell, outer_face_node_ptr,
      outer_face_node_ind, num_inner_peribdries, num_peribdries,
      num_total_peribdries, peribdry_global_index, peribdry_tag,
      peribdry_owner_cell, peribdry_neighbor_cell, peribdry_node_ptr,
      peribdry_node_ind, peribdry_pair_node_ind);
  outer_cell_global_index.clear();
  outer_cell_elemtype.clear();
  outer_cell_node_ptr.clear();
  outer_cell_node_ind.clear();
  outer_face_owner_cell.clear();
  outer_face_node_ptr.clear();
  outer_face_node_ind.clear();
  SYNCRO();
  MASTER_MESSAGE("Complete. (Time: " + std::to_string(STOP_TIMER()) + "s)\n");

  START_TIMER();
  MASTER_MESSAGE("Constructing face type... ");
  std::vector<int> face_owner_type;
  std::vector<int> face_neighbor_type;
  std::vector<int> bdry_owner_type;
  std::vector<int> peribdry_owner_type;
  std::vector<int> peribdry_neighbor_type;
  {
    std::unordered_map<int, int> cell_global_mapping;
    cell_global_mapping.reserve(num_outer_cells);
    for (int icell = 0; icell < num_outer_cells; icell++)
      cell_global_mapping[cell_global_index[icell]] = icell;

    ConstructFacetype(cell_global_mapping, cell_elemtype, cell_node_ptr,
                      cell_node_ind, face_owner_cell, face_node_ptr,
                      face_node_ind, face_owner_type);
    ConstructFacetype(cell_global_mapping, cell_elemtype, cell_node_ptr,
                      cell_node_ind, face_neighbor_cell, face_node_ptr,
                      face_node_ind, face_neighbor_type);
    ConstructFacetype(cell_global_mapping, cell_elemtype, cell_node_ptr,
                      cell_node_ind, bdry_owner_cell, bdry_node_ptr,
                      bdry_node_ind, bdry_owner_type);
    ConstructFacetype(cell_global_mapping, cell_elemtype, cell_node_ptr,
                      cell_node_ind, peribdry_owner_cell, peribdry_node_ptr,
                      peribdry_node_ind, peribdry_owner_type);
    ConstructFacetype(cell_global_mapping, cell_elemtype, cell_node_ptr,
                      cell_node_ind, peribdry_neighbor_cell, peribdry_node_ptr,
                      peribdry_pair_node_ind, peribdry_neighbor_type);

    cell_global_mapping.clear();
  }
  peribdry_pair_node_ind.clear();
  SYNCRO();
  MASTER_MESSAGE("Complete. (Time: " + std::to_string(STOP_TIMER()) + "s)\n");

  START_TIMER();
  MASTER_MESSAGE("Reading cell sub-nodes data... ");
  std::vector<int> cell_order;
  std::vector<int> cell_subnode_ptr;
  std::vector<int> cell_subnode_ind;
  {
    std::unordered_map<int, int> cell_global_mapping;
    for (auto&& global_index : cell_global_index)
      cell_global_mapping[global_index] = -1;
    grid_reader->ReadCellData(cell_global_mapping, cell_order, cell_subnode_ptr,
                              cell_subnode_ind);

    std::vector<int> mapped_from;
    mapped_from.reserve(num_total_cells);
    for (int icell = 0; icell < num_total_cells; icell++)
      mapped_from.push_back(cell_global_mapping[cell_global_index[icell]]);
    cell_global_mapping.clear();
    ReorderData(mapped_from, cell_order, 1);
    ReorderData(mapped_from, cell_subnode_ptr, cell_subnode_ind);
    mapped_from.clear();
  }
  SYNCRO();
  MASTER_MESSAGE("Complete. (Time: " + std::to_string(STOP_TIMER()) + "s)\n");

  START_TIMER();
  MASTER_MESSAGE("Reading node coordinates... ");
  int num_total_nodes;
  int num_nodes;
  std::vector<int> node_global_index;
  std::vector<double> node_coords;
  {
    std::unordered_map<int, int> node_global_mapping;
    for (auto&& inode : cell_subnode_ind) node_global_mapping[inode] = -1;
    grid_reader->ReadNodeData(node_global_mapping, node_coords);

    num_total_nodes = static_cast<int>(node_global_mapping.size());
    node_global_index.resize(num_total_nodes);
    for (auto&& iterator : node_global_mapping)
      node_global_index[iterator.second] = iterator.first;
    node_global_mapping.clear();

    std::unordered_set<int> nodes(
        cell_node_ind.begin(),
        cell_node_ind.begin() + cell_node_ptr[num_cells]);
    std::vector<int> node_mapping(num_total_nodes, -1);
    int mapped_node = 0;
    for (int inode = 0; inode < num_total_nodes; inode++)
      if (nodes.find(node_global_index[inode]) != nodes.end())
        node_mapping[inode] = mapped_node++;
    num_nodes = mapped_node;
    nodes.clear();
    for (int inode = 0; inode < num_total_nodes; inode++)
      if (node_mapping[inode] == -1) node_mapping[inode] = mapped_node++;
    std::vector<int> mapped_from(num_total_nodes);
    for (int inode = 0; inode < num_total_nodes; inode++)
      mapped_from[node_mapping[inode]] = inode;
    node_mapping.clear();
    ReorderData(mapped_from, node_global_index, 1);
    ReorderData(mapped_from, node_coords, dimension_);
    mapped_from.clear();
  }
  SYNCRO();
  MASTER_MESSAGE("Complete. (Time: " + std::to_string(STOP_TIMER()) + "s)\n");

  START_TIMER();
  MASTER_MESSAGE("Constructing communication data... ");
  std::vector<std::vector<int>> outer_send_cell_list;
  std::vector<std::vector<int>> outer_recv_cell_list;
  std::vector<std::vector<int>> total_send_cell_list;
  std::vector<std::vector<int>> total_recv_cell_list;
  {
    std::vector<int> target_data;
    target_data.reserve(num_outer_cells - num_cells);
    for (int icell = num_cells; icell < num_outer_cells; icell++)
      target_data.push_back(cell_global_index[icell]);

    std::vector<int> local_data(cell_global_index.begin(),
                                cell_global_index.begin() + num_cells);
    ConstructCommunication(local_data, target_data, outer_send_cell_list,
                           outer_recv_cell_list);
    target_data.clear();

    target_data.reserve(num_total_cells - num_cells);
    for (int icell = num_cells; icell < num_total_cells; icell++)
      target_data.push_back(cell_global_index[icell]);
    ConstructCommunication(local_data, target_data, total_send_cell_list,
                           total_recv_cell_list);
    local_data.clear();
    target_data.clear();
  }
  SYNCRO();
  MASTER_MESSAGE("Complete. (Time: " + std::to_string(STOP_TIMER()) + "s)\n");

  START_TIMER();
  MASTER_MESSAGE("Remapping node and cell... ");
  {
    std::unordered_map<int, int> node_mapping;
    for (int inode = 0; inode < num_total_nodes; inode++)
      node_mapping[node_global_index[inode]] = inode;

    std::unordered_map<int, int> cell_mapping;
    for (int icell = 0; icell < num_total_cells; icell++)
      cell_mapping[cell_global_index[icell]] = icell;

    RemapData(node_mapping, cell_node_ind);
    RemapData(node_mapping, cell_subnode_ind);
    RemapData(cell_mapping, face_owner_cell);
    RemapData(cell_mapping, face_neighbor_cell);
    RemapData(cell_mapping, bdry_owner_cell);
    RemapData(cell_mapping, peribdry_owner_cell);
    RemapData(cell_mapping, peribdry_neighbor_cell);
    node_mapping.clear();
    cell_mapping.clear();
  }
  SYNCRO();
  MASTER_MESSAGE("Complete. (Time: " + std::to_string(STOP_TIMER()) + "s)\n");

  START_TIMER();
  MASTER_MESSAGE("Pre-processing for system matrix construction... ");
  {
    for (int iface = 0; iface < num_inner_faces; iface++) {
      const int& owner_cell = face_owner_cell[iface];
      const int& neighbor_cell = face_neighbor_cell[iface];
      if (owner_cell > neighbor_cell) {
        std::swap(face_owner_cell[iface], face_neighbor_cell[iface]);
        std::swap(face_owner_type[iface], face_neighbor_type[iface]);
      }
    }
    for (int ibdry = 0; ibdry < num_inner_peribdries; ibdry++) {
      const int& owner_cell = peribdry_owner_cell[ibdry];
      const int& neighbor_cell = peribdry_neighbor_cell[ibdry];
      if (owner_cell > neighbor_cell) {
        std::swap(peribdry_owner_cell[ibdry], peribdry_neighbor_cell[ibdry]);
        std::swap(peribdry_owner_type[ibdry], peribdry_neighbor_type[ibdry]);
      }
    }
    for (int ibdry = num_peribdries; ibdry < num_total_peribdries; ibdry++) {
      const int& owner_cell = peribdry_owner_cell[ibdry];
      const int& neighbor_cell = peribdry_neighbor_cell[ibdry];
      if (owner_cell > neighbor_cell) {
        std::swap(peribdry_owner_cell[ibdry], peribdry_neighbor_cell[ibdry]);
        std::swap(peribdry_owner_type[ibdry], peribdry_neighbor_type[ibdry]);
      }
    }

    std::vector<int> forced_cell_mapping(num_total_cells);
    for (int icell = 0; icell < num_total_cells; icell++)
      forced_cell_mapping[icell] = icell;
    int ind = num_cells;
    for (int idomain = 0; idomain < NDOMAIN; idomain++)
      for (auto&& cell : outer_recv_cell_list[idomain])
        forced_cell_mapping[ind++] = cell + num_cells;

    std::unordered_map<int, int> cell_mapping;
    cell_mapping.reserve(num_total_cells - num_cells);
    for (int icell = num_cells; icell < num_total_cells; icell++)
      cell_mapping[forced_cell_mapping[icell] - num_cells] = icell - num_cells;

    ReorderData(forced_cell_mapping, cell_global_index, 1);
    ReorderData(forced_cell_mapping, cell_elemtype, 1);
    ReorderData(forced_cell_mapping, cell_order, 1);
    ReorderData(forced_cell_mapping, cell_node_ptr, cell_node_ind);
    ReorderData(forced_cell_mapping, cell_subnode_ptr, cell_subnode_ind);
    forced_cell_mapping.clear();

    for (auto&& cells : outer_recv_cell_list) RemapData(cell_mapping, cells);
    for (auto&& cells : total_recv_cell_list) RemapData(cell_mapping, cells);
    for (int iface = num_inner_faces; iface < num_faces; iface++)
      face_neighbor_cell[iface] =
          cell_mapping.at(face_neighbor_cell[iface] - num_cells) + num_cells;
    for (int ibdry = num_inner_peribdries; ibdry < num_peribdries; ibdry++)
      peribdry_neighbor_cell[ibdry] =
          cell_mapping.at(peribdry_neighbor_cell[ibdry] - num_cells) +
          num_cells;
    cell_mapping.clear();

    std::unordered_map<int, int> outer_cell_to_domain;
    outer_cell_to_domain.reserve(num_outer_cells - num_cells);
    for (int idomain = 0; idomain < NDOMAIN; idomain++)
      for (auto&& cell : outer_recv_cell_list[idomain])
        outer_cell_to_domain[cell + num_cells] = idomain;

    std::vector<int> forced_face_mapping(num_faces);
    for (int iface = 0; iface < num_inner_faces; iface++)
      forced_face_mapping[iface] = iface;
    ind = num_inner_faces;
    for (int idomain = 0; idomain < NDOMAIN; idomain++)
      for (int iface = num_inner_faces; iface < num_faces; iface++)
        if (outer_cell_to_domain[face_neighbor_cell[iface]] == idomain)
          forced_face_mapping[ind++] = iface;

    ReorderData(forced_face_mapping, face_owner_cell, 1);
    ReorderData(forced_face_mapping, face_neighbor_cell, 1);
    ReorderData(forced_face_mapping, face_owner_type, 1);
    ReorderData(forced_face_mapping, face_neighbor_type, 1);
    forced_face_mapping.clear();

    std::vector<int> forced_peribdry_mapping(num_total_peribdries);
    for (int ibdry = 0; ibdry < num_total_peribdries; ibdry++)
      forced_peribdry_mapping[ibdry] = ibdry;
    ind = num_inner_peribdries;
    for (int idomain = 0; idomain < NDOMAIN; idomain++)
      for (int ibdry = num_inner_peribdries; ibdry < num_peribdries; ibdry++)
        if (outer_cell_to_domain[peribdry_neighbor_cell[ibdry]] == idomain)
          forced_peribdry_mapping[ind++] = ibdry;
    ReorderData(forced_peribdry_mapping, peribdry_global_index, 1);
    ReorderData(forced_peribdry_mapping, peribdry_tag, 1);
    ReorderData(forced_peribdry_mapping, peribdry_owner_cell, 1);
    ReorderData(forced_peribdry_mapping, peribdry_neighbor_cell, 1);
    ReorderData(forced_peribdry_mapping, peribdry_owner_type, 1);
    ReorderData(forced_peribdry_mapping, peribdry_neighbor_type, 1);
    forced_peribdry_mapping.clear();
    outer_cell_to_domain.clear();
  }
  SYNCRO();
  MASTER_MESSAGE("Complete. (Time: " + std::to_string(STOP_TIMER()) + "s)\n");

  START_TIMER();
  {
    MOVE_DATA(num_total_nodes);
    MOVE_DATA(num_nodes);
    MOVE_DATA(node_global_index);
    MOVE_DATA(node_coords);

    MOVE_DATA(num_total_cells);
    MOVE_DATA(num_outer_cells);
    MOVE_DATA(num_cells);
    MOVE_DATA(cell_global_index);
    MOVE_DATA(cell_elemtype);
    MOVE_DATA(cell_order);
    MOVE_DATA(cell_node_ptr);
    MOVE_DATA(cell_node_ind);
    MOVE_DATA(cell_subnode_ptr);
    MOVE_DATA(cell_subnode_ind);

    MOVE_DATA(num_bdries);
    MOVE_DATA(bdry_global_index);
    MOVE_DATA(bdry_tag);
    MOVE_DATA(bdry_owner_cell);
    MOVE_DATA(bdry_owner_type);

    MOVE_DATA(num_total_peribdries);
    MOVE_DATA(num_peribdries);
    MOVE_DATA(num_inner_peribdries);
    MOVE_DATA(peribdry_global_index);
    MOVE_DATA(peribdry_tag);
    MOVE_DATA(peribdry_owner_cell);
    MOVE_DATA(peribdry_neighbor_cell);
    MOVE_DATA(peribdry_owner_type);
    MOVE_DATA(peribdry_neighbor_type);

    MOVE_DATA(num_faces);
    MOVE_DATA(num_inner_faces);
    MOVE_DATA(face_owner_cell);
    MOVE_DATA(face_neighbor_cell);
    MOVE_DATA(face_owner_type);
    MOVE_DATA(face_neighbor_type);

    periodic_matching_global_node_index_ = std::move(periodic_matching_nodes);

    MOVE_DATA(outer_send_cell_list);
    MOVE_DATA(outer_recv_cell_list);
    MOVE_DATA(total_send_cell_list);
    MOVE_DATA(total_recv_cell_list);
  }
  SYNCRO();
  STOP_TIMER_TAG("BuildGrid");
  MASTER_MESSAGE("GridBuilder::BuildGrid() completes. (Time: " +
                 std::to_string(GET_RECORD_TAG("BuildGrid")) + "s)\n");

  grid_reader.reset();
}

bool GridBuilder::IsParent(const int* cell_node_index, const int cell_num_nodes,
                           const int* face_node_index,
                           const int face_num_nodes) const {
  for (int i = 0; i < face_num_nodes; i++)
    if (std::find(cell_node_index, cell_node_index + cell_num_nodes,
                  face_node_index[i]) == cell_node_index + cell_num_nodes)
      return false;
  return true;
}

bool GridBuilder::IsSame(const int* index1, const int num_index1,
                         const int* index2, const int num_index2) const {
  if (num_index1 != num_index2) return false;
  for (int i = 0; i < num_index1; i++)
    if (std::find(index1, index1 + num_index1, index2[i]) ==
        index1 + num_index1)
      return false;
  return true;
}

int GridBuilder::FindFacetype(
    const std::vector<std::vector<int>>& facetype_nodes_index,
    const int* cell_nodes_index, const int* face_nodes_index,
    const int num_face_nodes) const {
  const int num_facetypes = static_cast<int>(facetype_nodes_index.size());
  for (int itype = 0; itype < num_facetypes; itype++) {
    std::vector<int> nodes_index;
    for (auto&& inode : facetype_nodes_index[itype])
      nodes_index.push_back(cell_nodes_index[inode]);
    if (IsSame(&nodes_index[0], static_cast<int>(nodes_index.size()),
               face_nodes_index, num_face_nodes))
      return itype;
  }
  ERROR_MESSAGE("Face type is not found.\n");
  return -1;
}

void GridBuilder::RemapData(const std::unordered_map<int, int>& mapping,
                            std::vector<int>& data) const {
  for (auto&& var : data) {
    auto&& iterator = mapping.find(var);
    if (iterator == mapping.end()) ERROR_MESSAGE("mapping is incomplete.\n");
    var = iterator->second;
  }
}

void GridBuilder::ConstructCommunication(
    const std::vector<int>& local_data, const std::vector<int>& target_data,
    std::vector<std::vector<int>>& send_locations,
    std::vector<std::vector<int>>& recv_locations) const {
  send_locations.clear();
  recv_locations.clear();
  send_locations.resize(NDOMAIN);
  recv_locations.resize(NDOMAIN);
  if (NDOMAIN == 1) return;

  std::vector<std::vector<int>> send_data(NDOMAIN);
  std::vector<std::vector<int>> recv_data;
  for (int idomain = 0; idomain < NDOMAIN; idomain++)
    send_data[idomain] = target_data;
  AVOCADO_MPI->CommunicateData(send_data, recv_data);
  send_data.clear();

  send_data.resize(NDOMAIN);
  std::unordered_map<int, int> data_mapping;
  for (int i = 0, len = static_cast<int>(local_data.size()); i < len; i++)
    data_mapping[local_data[i]] = i;
  for (int idomain = 0; idomain < NDOMAIN; idomain++) {
    for (auto&& data : recv_data[idomain]) {
      auto&& iterator = data_mapping.find(data);
      if (iterator == data_mapping.end()) continue;
      send_data[idomain].push_back(data);
      send_locations[idomain].push_back(iterator->second);
    }
  }
  data_mapping.clear();
  recv_data.clear();

  AVOCADO_MPI->CommunicateData(send_data, recv_data);
  send_data.clear();

  for (int i = 0, len = static_cast<int>(target_data.size()); i < len; i++)
    data_mapping[target_data[i]] = i;
  for (int idomain = 0; idomain < NDOMAIN; idomain++)
    for (auto&& data : recv_data[idomain])
      recv_locations[idomain].push_back(data_mapping[data]);
  data_mapping.clear();
  recv_data.clear();
}

void GridBuilder::SetCellColor(const std::vector<int>& cell_node_ptr,
                               const std::vector<int>& cell_node_ind,
                               std::vector<int>& cell_color) const {
  cell_color.clear();
  const int num_cells = static_cast<int>(cell_node_ptr.size()) - 1;
  if (NDOMAIN != 1) {
    MASTER_MESSAGE("\n\tMETIS integer (idx_t) size: " +
                   std::to_string(static_cast<int>(sizeof(idx_t))) + "\n");

    idx_t wgtflag = 0;
    idx_t numflag = 0;
    idx_t ncon = 1;
    idx_t ncommonnodes = dimension_;
    idx_t npart = NDOMAIN;
    std::vector<double> tpwgts;
    tpwgts.resize(npart);
    for (idx_t ipart = 0; ipart < npart; ++ipart)
      tpwgts[ipart] = 1.0 / double(npart);
    double ubvec = 1.05;

    std::vector<idx_t> options;
    options.resize(METIS_NOPTIONS);
    METIS_SetDefaultOptions(&options[0]);
    options[0] = 0;

    idx_t edgecut;
    std::vector<idx_t> part(num_cells);
    std::vector<idx_t> eptr;
    std::vector<idx_t> eind;
    std::vector<idx_t> elmdist;

    eptr.reserve(cell_node_ptr.size());
    for (auto&& iterator : cell_node_ptr) eptr.push_back(iterator);

    eind.reserve(cell_node_ind.size());
    for (auto&& iterator : cell_node_ind) eind.push_back(iterator);

    elmdist.resize(NDOMAIN + 1, 0);
    std::vector<int> num_cells_gather(NDOMAIN);
    MPI_Allgather(&num_cells, 1, MPI_INT, &num_cells_gather[0], 1, MPI_INT,
                  MPI_COMM_WORLD);
    for (int idomain = 0; idomain < NDOMAIN; idomain++)
      elmdist[idomain + 1] = elmdist[idomain] + num_cells_gather[idomain];

    MPI_Comm comm = MPI_COMM_WORLD;
    int ret = ParMETIS_V3_PartMeshKway(&elmdist[0], &eptr[0], &eind[0], NULL,
                                       &wgtflag, &numflag, &ncon, &ncommonnodes,
                                       &npart, &tpwgts[0], &ubvec, &options[0],
                                       &edgecut, &part[0], &comm);

    switch (ret) {
      case METIS_OK:
        MASTER_MESSAGE("ParMETIS works normally!\n");
        break;
      case METIS_ERROR_INPUT:
        ERROR_MESSAGE("ParMETIS returns ERROR_INPUT!\n");
        break;
      case METIS_ERROR_MEMORY:
        ERROR_MESSAGE("ParMETIS returns ERROR_MEMORY!\n");
        break;
      case METIS_ERROR:
        ERROR_MESSAGE("ParMETIS returns ERROR!\n");
        break;
      default:
        ERROR_MESSAGE("ParMETIS returns unknown! : " + std::to_string(ret) +
                      "\n");
    }

    cell_color.reserve(num_cells);
    for (auto&& it : part) cell_color.push_back(it);

  } else {
    MASTER_MESSAGE("No partition... ");
    cell_color.resize(num_cells, MYRANK);
  }
  SYNCRO();
}

void GridBuilder::ConstructFaceData(
    const std::vector<int>& cell_global_index,
    const std::vector<ElemType>& cell_elemtype,
    const std::vector<int>& cell_node_ptr,
    const std::vector<int>& cell_node_ind, std::vector<int>& face_owner_cell,
    std::vector<int>& face_neighbor_cell, std::vector<int>& face_node_ptr,
    std::vector<int>& face_node_ind, std::vector<int>& outer_face_owner_cell,
    std::vector<int>& outer_face_node_ptr,
    std::vector<int>& outer_face_node_ind) const {
  const int num_cells = static_cast<int>(cell_global_index.size());
  const int num_estimate_faces = 2 * num_cells * dimension_;

  face_owner_cell.clear();
  face_neighbor_cell.clear();
  face_node_ptr.clear();
  face_node_ind.clear();
  face_owner_cell.reserve(num_estimate_faces);
  face_neighbor_cell.reserve(num_estimate_faces);
  face_node_ptr.reserve(num_estimate_faces);
  face_node_ptr.push_back(0);
  face_node_ind.reserve(num_estimate_faces * (dimension_ - 1));
  {
    int num_nodes = 0;
    std::vector<int> index_to_node;
    std::vector<int> cell_node_index;
    {
      std::unordered_map<int, int> node_to_index;
      for (auto&& node : cell_node_ind) {
        auto&& iterator = node_to_index.find(node);
        if (iterator == node_to_index.end()) node_to_index[node] = num_nodes++;
      }

      index_to_node.resize(num_nodes);
      for (auto&& iterator : node_to_index)
        index_to_node[iterator.second] = iterator.first;

      cell_node_index = cell_node_ind;
      for (auto&& node : cell_node_index) node = node_to_index[node];
      node_to_index.clear();
    }

    std::vector<int> first_face(num_nodes, -1);
    std::vector<int> next_face(num_estimate_faces, -1);
    for (int icell = 0; icell < num_cells; icell++) {
      const Element* element = DENEB_ELEMENT->GetElement(cell_elemtype[icell]);
      const std::vector<std::vector<int>>& local_faces =
          element->GetFaceNodes();
      for (auto&& local_face : local_faces) {
        std::vector<int> nodes_index;
        for (auto&& inode : local_face)
          nodes_index.push_back(cell_node_index[cell_node_ptr[icell] + inode]);

        const int minnode =
            *std::min_element(nodes_index.begin(), nodes_index.end());
        int current_face = first_face[minnode];
        int previous_face = -1;
        while (true) {
          if (current_face == -1) {
            face_owner_cell.push_back(cell_global_index[icell]);
            face_neighbor_cell.resize(face_owner_cell.size(), -1);
            for (auto&& inode : nodes_index) face_node_ind.push_back(inode);
            face_node_ptr.push_back(static_cast<int>(face_node_ind.size()));
            if (previous_face == -1)
              first_face[minnode] =
                  static_cast<int>(face_owner_cell.size()) - 1;
            else
              next_face[previous_face] =
                  static_cast<int>(face_owner_cell.size()) - 1;
            break;
          } else if (IsParent(&face_node_ind[face_node_ptr[current_face]],
                              face_node_ptr[current_face + 1] -
                                  face_node_ptr[current_face],
                              &nodes_index[0],
                              static_cast<int>(nodes_index.size()))) {
            face_neighbor_cell[current_face] = cell_global_index[icell];
            break;
          }
          previous_face = current_face;
          current_face = next_face[current_face];
        }
      }
    }
    first_face.clear();
    next_face.clear();
    cell_node_index.clear();

    for (auto&& node : face_node_ind) node = index_to_node[node];
    index_to_node.clear();
  }

  const int num_faces = static_cast<int>(face_owner_cell.size());
  for (int iface = 0; iface < num_faces; iface++) {
    const int& a = face_owner_cell[iface];
    const int& b = face_neighbor_cell[iface];
    if (b == -1) continue;
    if (a > b) std::swap(face_owner_cell[iface], face_neighbor_cell[iface]);
  }

  for (int iface = 0; iface < num_faces; iface++) {
    const std::vector<int> face_nodes(
        face_node_ind.begin() + face_node_ptr[iface],
        face_node_ind.begin() + face_node_ptr[iface + 1]);
    const int mnode = static_cast<int>(
        std::min_element(face_nodes.begin(), face_nodes.end()) -
        face_nodes.begin());
    const int nnode = static_cast<int>(face_nodes.size());

    const int right = (mnode + 1 == nnode) ? 0 : mnode + 1;
    const int left = (mnode == 0) ? nnode - 1 : mnode - 1;
    int ind = 0;
    if (face_nodes[right] <= face_nodes[left])
      for (int jnode = 0; jnode < nnode; jnode++)
        face_node_ind[face_node_ptr[iface] + ind++] =
            face_nodes[(mnode + jnode) % nnode];
    else
      for (int jnode = 0; jnode < nnode; jnode++)
        face_node_ind[face_node_ptr[iface] + ind++] =
            face_nodes[(mnode - jnode + nnode) % nnode];
  }

  outer_face_owner_cell.clear();
  outer_face_node_ptr.clear();
  outer_face_node_ptr.push_back(0);
  outer_face_node_ind.clear();
  {
    std::vector<int> face_node_ind_temp = std::move(face_node_ind);
    std::vector<int> face_node_ptr_temp = std::move(face_node_ptr);
    std::vector<int> face_owner_cell_temp = std::move(face_owner_cell);
    std::vector<int> face_neighbor_cell_temp = std::move(face_neighbor_cell);
    face_node_ptr.push_back(0);
    for (int iface = 0; iface < num_faces; iface++) {
      if (face_neighbor_cell_temp[iface] != -1) {
        face_owner_cell.push_back(face_owner_cell_temp[iface]);
        face_neighbor_cell.push_back(face_neighbor_cell_temp[iface]);
        for (int index = face_node_ptr_temp[iface];
             index < face_node_ptr_temp[iface + 1]; index++)
          face_node_ind.push_back(face_node_ind_temp[index]);
        face_node_ptr.push_back(static_cast<int>(face_node_ind.size()));
      } else {
        outer_face_owner_cell.push_back(face_owner_cell_temp[iface]);
        for (int index = face_node_ptr_temp[iface];
             index < face_node_ptr_temp[iface + 1]; index++)
          outer_face_node_ind.push_back(face_node_ind_temp[index]);
        outer_face_node_ptr.push_back(
            static_cast<int>(outer_face_node_ind.size()));
      }
    }
  }
}

void GridBuilder::ConstructBdryData(
    const std::vector<int>& all_bdry_global_index,
    const std::vector<int>& all_bdry_tag,
    const std::vector<int>& all_bdry_node_ptr,
    const std::vector<int>& all_bdry_node_ind,
    std::vector<int>& outer_face_owner_cell,
    std::vector<int>& outer_face_node_ptr,
    std::vector<int>& outer_face_node_ind, std::vector<int>& bdry_global_index,
    std::vector<int>& bdry_tag, std::vector<int>& bdry_owner_cell,
    std::vector<int>& bdry_node_ptr, std::vector<int>& bdry_node_ind,
    std::vector<int>& peribdry_global_index, std::vector<int>& peribdry_tag,
    std::vector<int>& peribdry_owner_cell, std::vector<int>& peribdry_node_ptr,
    std::vector<int>& peribdry_node_ind) const {
  std::vector<int> target_owner_cell = std::move(outer_face_owner_cell);
  std::vector<int> target_node_ptr = std::move(outer_face_node_ptr);
  std::vector<int> target_node_ind = std::move(outer_face_node_ind);
  const int num_target_bdries = static_cast<int>(target_owner_cell.size());

  std::unordered_set<int> periBC_tag;
  auto& config = AVOCADO_CONFIG;
  for (int idim = 0; idim < 3; idim++) {
    periBC_tag.insert(std::stoi(config->GetConfigValue(PERIODIC_LEFT(idim))));
    periBC_tag.insert(std::stoi(config->GetConfigValue(PERIODIC_RIGHT(idim))));
  }
  periBC_tag.insert(std::stoi(config->GetConfigValue(MATCHING_BC)));

  const int num_all_bdries = static_cast<int>(all_bdry_global_index.size());
  std::unordered_map<int, std::vector<int>> node_to_all_bdries;
  for (int ibdry = 0; ibdry < num_all_bdries; ibdry++)
    for (int index = all_bdry_node_ptr[ibdry];
         index < all_bdry_node_ptr[ibdry + 1]; index++)
      node_to_all_bdries[all_bdry_node_ind[index]].push_back(ibdry);

  outer_face_owner_cell.clear();
  outer_face_node_ptr.clear();
  outer_face_node_ptr.push_back(0);
  outer_face_node_ind.clear();
  bdry_global_index.clear();
  bdry_tag.clear();
  bdry_owner_cell.clear();
  bdry_node_ptr.clear();
  bdry_node_ptr.push_back(0);
  bdry_node_ind.clear();
  peribdry_global_index.clear();
  peribdry_tag.clear();
  peribdry_owner_cell.clear();
  peribdry_node_ptr.clear();
  peribdry_node_ptr.push_back(0);
  peribdry_node_ind.clear();
  std::vector<char> flag(num_all_bdries, false);
  for (int ibdry = 0; ibdry < num_target_bdries; ibdry++) {
    bool ismatch = false;
    for (auto&& jbdry :
         node_to_all_bdries[target_node_ind[target_node_ptr[ibdry]]]) {
      if (flag[jbdry] == false &&
          IsParent(&target_node_ind[target_node_ptr[ibdry]],
                   target_node_ptr[ibdry + 1] - target_node_ptr[ibdry],
                   &all_bdry_node_ind[all_bdry_node_ptr[jbdry]],
                   all_bdry_node_ptr[jbdry + 1] - all_bdry_node_ptr[jbdry])) {
        if (periBC_tag.find(all_bdry_tag[jbdry]) == periBC_tag.end()) {
          bdry_global_index.push_back(all_bdry_global_index[jbdry]);
          bdry_tag.push_back(all_bdry_tag[jbdry]);
          bdry_owner_cell.push_back(target_owner_cell[ibdry]);
          for (int index = target_node_ptr[ibdry];
               index < target_node_ptr[ibdry + 1]; index++)
            bdry_node_ind.push_back(target_node_ind[index]);
          bdry_node_ptr.push_back(static_cast<int>(bdry_node_ind.size()));
        } else {
          peribdry_global_index.push_back(all_bdry_global_index[jbdry]);
          peribdry_tag.push_back(all_bdry_tag[jbdry]);
          peribdry_owner_cell.push_back(target_owner_cell[ibdry]);
          for (int index = target_node_ptr[ibdry];
               index < target_node_ptr[ibdry + 1]; index++)
            peribdry_node_ind.push_back(target_node_ind[index]);
          peribdry_node_ptr.push_back(
              static_cast<int>(peribdry_node_ind.size()));
        }
        flag[jbdry] = true;
        ismatch = true;
        break;
      }
    }
    if (ismatch == false) {
      outer_face_owner_cell.push_back(target_owner_cell[ibdry]);
      for (int index = target_node_ptr[ibdry];
           index < target_node_ptr[ibdry + 1]; index++)
        outer_face_node_ind.push_back(target_node_ind[index]);
      outer_face_node_ptr.push_back(
          static_cast<int>(outer_face_node_ind.size()));
    }
  }
  flag.clear();
  node_to_all_bdries.clear();
  target_owner_cell.clear();
  target_node_ptr.clear();
  target_node_ind.clear();
}

void GridBuilder::PeriodicMatching(
    const std::vector<int>& all_periodic_nodes,
    const std::vector<double>& all_periodic_node_coords,
    const std::vector<int>& peribdry_global_index,
    const std::vector<int>& peribdry_tag,
    const std::vector<int>& peribdry_node_ptr,
    const std::vector<int>& peribdry_node_ind,
    std::vector<int>& peribdry_pair_node_ind) const {
  peribdry_pair_node_ind.clear();
  const int num_peribdries = static_cast<int>(peribdry_global_index.size());
  if (num_peribdries == 0) return;

  std::vector<int> periBC_tag;
  std::vector<double> periBC_offset;
  auto& config = AVOCADO_CONFIG;
  for (int idim = 0; idim < 3; idim++) {
    periBC_tag.push_back(
        std::stoi(config->GetConfigValue(PERIODIC_LEFT(idim))));
    periBC_tag.push_back(
        std::stoi(config->GetConfigValue(PERIODIC_RIGHT(idim))));
    periBC_offset.push_back(
        std::stod(config->GetConfigValue(PERIODIC_OFFSET(idim))));
  }
  periBC_tag.push_back(std::stoi(config->GetConfigValue(MATCHING_BC)));

  std::unordered_map<int, int> node_to_index;
  {
    int index = 0;
    for (auto&& node : all_periodic_nodes) node_to_index[node] = index++;
  }

  std::unordered_map<int, std::unordered_map<int, int>> history;
  peribdry_pair_node_ind.resize(peribdry_node_ptr.back(), -1);
  for (int ibdry = 0; ibdry < num_peribdries; ibdry++) {
    const int tag = static_cast<int>(
        std::find(periBC_tag.begin(), periBC_tag.end(), peribdry_tag[ibdry]) -
        periBC_tag.begin());

    int target_tag = tag;
    std::vector<double> offset(dimension_, 0.0);
    {
      const int sign = (tag % 2 == 0) ? 1 : -1;
      target_tag += sign;
      offset[tag / 2] += sign * periBC_offset[tag / 2];
    }

    for (int ptr = peribdry_node_ptr[ibdry]; ptr < peribdry_node_ptr[ibdry + 1];
         ptr++) {
      const int& node = peribdry_node_ind[ptr];
      const int& index = node_to_index[node];
      std::vector<double> coord(dimension_, 0.0);
      for (int idim = 0; idim < dimension_; idim++)
        coord[idim] =
            all_periodic_node_coords[index * dimension_ + idim] + offset[idim];

      auto&& iter = history.find(node);
      if (iter != history.end()) {
        auto&& jter = iter->second.find(tag);
        if (jter != iter->second.end()) {
          peribdry_pair_node_ind[ptr] = jter->second;
          continue;
        }
      }

      for (auto&& iterator : node_to_index) {
        const int& target_node = iterator.first;
        const int& target_index = iterator.second;
        if (target_node != node &&
            avocado::VecIsSame(
                dimension_, &coord[0],
                &all_periodic_node_coords[target_index * dimension_],
                kPeriodicEps)) {
          peribdry_pair_node_ind[ptr] = target_node;
          history[node][tag] = target_node;
          break;
        }
      }
    }
  }
  history.clear();
  node_to_index.clear();

  for (int ibdry = 0; ibdry < num_peribdries; ibdry++)
    for (int ptr = peribdry_node_ptr[ibdry]; ptr < peribdry_node_ptr[ibdry + 1];
         ptr++)
      if (peribdry_pair_node_ind[ptr] < 0)
        ERROR_MESSAGE(
            "Periodic pair node is not found.\n\tPeribdry global index: " +
            std::to_string(peribdry_global_index[ibdry]) +
            "\n\tPeribdry tag: " + std::to_string(peribdry_tag[ibdry]) + "\n");
}

void GridBuilder::PeriodicMatching(
    const std::vector<int>& peribdry_node_ind,
    const std::vector<int>& peribdry_pair_node_ind,
    std::unordered_map<int, std::vector<int>>& periodic_matching_nodes) const {
  std::unordered_map<int, std::unordered_set<int>> node_grps;
  {
    const int num_periodic_nodes = static_cast<int>(peribdry_node_ind.size());
    for (int i = 0; i < num_periodic_nodes; i++)
      node_grps[peribdry_node_ind[i]].insert(peribdry_pair_node_ind[i]);

    std::vector<int> periodic_pairs;
    periodic_pairs.reserve(num_periodic_nodes * 2);
    for (int i = 0; i < num_periodic_nodes; i++) {
      periodic_pairs.push_back(peribdry_node_ind[i]);
      periodic_pairs.push_back(peribdry_pair_node_ind[i]);
    }

    std::vector<std::vector<int>> send_data(NDOMAIN, periodic_pairs);
    std::vector<std::vector<int>> recv_data;
    AVOCADO_MPI->CommunicateData(send_data, recv_data);
    send_data.clear();

    for (int idim = 0; idim < dimension_; idim++) {
      for (auto&& iterator : node_grps) {
        const int& mynode = iterator.first;

        std::unordered_set<int> new_nodes;
        for (auto&& node : iterator.second) {
          for (int i = 0; i < NDOMAIN; i++) {
            if (recv_data[i].size() == 0) continue;
            const size_t size = recv_data[i].size();
            for (int j = 0; j < size; j += 2) {
              if (recv_data[i][j] == node && recv_data[i][j + 1] != mynode)
                new_nodes.insert(recv_data[i][j + 1]);
              if (recv_data[i][j + 1] == mynode && recv_data[i][j] != mynode)
                new_nodes.insert(recv_data[i][j]);
            }
          }
        }

        iterator.second.insert(new_nodes.begin(), new_nodes.end());
      }
    }
  }

  periodic_matching_nodes.clear();
  for (auto&& grps : node_grps)
    periodic_matching_nodes[grps.first] =
        std::vector<int>(grps.second.begin(), grps.second.end());
}

void GridBuilder::ConstructOuterGrid(
    int& num_cells, int& num_outer_cells, int& num_total_cells,
    std::vector<int>& cell_global_index, std::vector<ElemType>& cell_elemtype,
    std::vector<int>& cell_node_ptr, std::vector<int>& cell_node_ind,
    const std::vector<int>& outer_cell_global_index,
    const std::vector<ElemType>& outer_cell_elemtype,
    const std::vector<int>& outer_cell_node_ptr,
    const std::vector<int>& outer_cell_node_ind, int& num_inner_faces,
    int& num_faces, std::vector<int>& face_owner_cell,
    std::vector<int>& face_neighbor_cell, std::vector<int>& face_node_ptr,
    std::vector<int>& face_node_ind,
    const std::vector<int>& outer_face_owner_cell,
    const std::vector<int>& outer_face_node_ptr,
    const std::vector<int>& outer_face_node_ind, int& num_inner_peribdries,
    int& num_peribdries, int& num_total_peribdries,
    std::vector<int>& peribdry_global_index, std::vector<int>& peribdry_tag,
    std::vector<int>& peribdry_owner_cell,
    std::vector<int>& peribdry_neighbor_cell,
    std::vector<int>& peribdry_node_ptr, std::vector<int>& peribdry_node_ind,
    std::vector<int>& peribdry_pair_node_ind) const {
  std::vector<int> outer_face_neighbor_cell;
  std::vector<int> peribdry_pair_owner_cell;
  {
    num_cells = static_cast<int>(cell_global_index.size());
    const int num_outer_faces = static_cast<int>(outer_face_owner_cell.size());
    const int num_target_cells =
        static_cast<int>(outer_cell_global_index.size());
    num_total_peribdries = static_cast<int>(peribdry_global_index.size());

    outer_face_neighbor_cell.resize(num_outer_faces, -1);
    peribdry_pair_owner_cell.resize(num_total_peribdries, -1);
    enum CellType : char { CELL_LOCAL, CELL_OUTER, CELL_VERTEX };
    std::vector<CellType> cell_types(num_target_cells, CellType::CELL_VERTEX);
    for (int icell = 0; icell < num_target_cells; icell++) {
      CellType& celltype = cell_types[icell];
      for (int jcell = 0; jcell < num_cells; jcell++) {
        if (IsSame(&outer_cell_node_ind[outer_cell_node_ptr[icell]],
                   outer_cell_node_ptr[icell + 1] - outer_cell_node_ptr[icell],
                   &cell_node_ind[cell_node_ptr[jcell]],
                   cell_node_ptr[jcell + 1] - cell_node_ptr[jcell])) {
          celltype = CellType::CELL_LOCAL;
          break;
        }
      }
      if (celltype == CellType::CELL_LOCAL) continue;
      for (int iface = 0; iface < num_outer_faces; iface++) {
        if (IsParent(
                &outer_cell_node_ind[outer_cell_node_ptr[icell]],
                outer_cell_node_ptr[icell + 1] - outer_cell_node_ptr[icell],
                &outer_face_node_ind[outer_face_node_ptr[iface]],
                outer_face_node_ptr[iface + 1] - outer_face_node_ptr[iface])) {
          celltype = CellType::CELL_OUTER;
          outer_face_neighbor_cell[iface] = outer_cell_global_index[icell];
          break;
        }
      }
      if (celltype == CellType::CELL_OUTER) continue;
      for (int ibdry = 0; ibdry < num_total_peribdries; ibdry++) {
        if (IsParent(
                &outer_cell_node_ind[outer_cell_node_ptr[icell]],
                outer_cell_node_ptr[icell + 1] - outer_cell_node_ptr[icell],
                &peribdry_pair_node_ind[peribdry_node_ptr[ibdry]],
                peribdry_node_ptr[ibdry + 1] - peribdry_node_ptr[ibdry])) {
          celltype = CellType::CELL_OUTER;
          peribdry_pair_owner_cell[ibdry] = outer_cell_global_index[icell];
          break;
        }
      }
    }

    for (int icell = 0; icell < num_target_cells; icell++) {
      if (cell_types[icell] != CELL_OUTER) continue;
      cell_global_index.push_back(outer_cell_global_index[icell]);
      cell_elemtype.push_back(outer_cell_elemtype[icell]);
      for (int ptr = outer_cell_node_ptr[icell];
           ptr < outer_cell_node_ptr[icell + 1]; ptr++)
        cell_node_ind.push_back(outer_cell_node_ind[ptr]);
      cell_node_ptr.push_back(static_cast<int>(cell_node_ind.size()));
    }
    num_outer_cells = static_cast<int>(cell_global_index.size());

    for (int icell = 0; icell < num_target_cells; icell++) {
      if (cell_types[icell] != CELL_VERTEX) continue;
      cell_global_index.push_back(outer_cell_global_index[icell]);
      cell_elemtype.push_back(outer_cell_elemtype[icell]);
      for (int ptr = outer_cell_node_ptr[icell];
           ptr < outer_cell_node_ptr[icell + 1]; ptr++)
        cell_node_ind.push_back(outer_cell_node_ind[ptr]);
      cell_node_ptr.push_back(static_cast<int>(cell_node_ind.size()));
    }
    num_total_cells = static_cast<int>(cell_global_index.size());
    cell_types.clear();
  }

  {
    const int num_outer_faces = static_cast<int>(outer_face_owner_cell.size());
    for (int iface = 0; iface < num_outer_faces; iface++) {
      if (outer_face_neighbor_cell[iface] >= 0) continue;
      for (int icell = num_cells; icell < num_outer_cells; icell++) {
        if (IsParent(
                &cell_node_ind[cell_node_ptr[icell]],
                cell_node_ptr[icell + 1] - cell_node_ptr[icell],
                &outer_face_node_ind[outer_face_node_ptr[iface]],
                outer_face_node_ptr[iface + 1] - outer_face_node_ptr[iface])) {
          outer_face_neighbor_cell[iface] = cell_global_index[icell];
          break;
        }
      }
    }

    for (int iface = 0; iface < num_outer_faces; iface++)
      if (outer_face_neighbor_cell[iface] < 0)
        ERROR_MESSAGE("Neighbor cell of outer face is not found.\n");

    num_inner_faces = static_cast<int>(face_owner_cell.size());
    for (int iface = 0; iface < num_outer_faces; iface++) {
      face_owner_cell.push_back(outer_face_owner_cell[iface]);
      face_neighbor_cell.push_back(outer_face_neighbor_cell[iface]);
      for (int ptr = outer_face_node_ptr[iface];
           ptr < outer_face_node_ptr[iface + 1]; ptr++)
        face_node_ind.push_back(outer_face_node_ind[ptr]);
      face_node_ptr.push_back(static_cast<int>(face_node_ind.size()));
    }
    num_faces = static_cast<int>(face_owner_cell.size());
    outer_face_neighbor_cell.clear();
  }

  {
    num_total_peribdries = static_cast<int>(peribdry_global_index.size());
    for (int ibdry = 0; ibdry < num_total_peribdries; ibdry++) {
      if (peribdry_pair_owner_cell[ibdry] >= 0) continue;
      for (int icell = 0; icell < num_outer_cells; icell++) {
        if (IsParent(&cell_node_ind[cell_node_ptr[icell]],
                     cell_node_ptr[icell + 1] - cell_node_ptr[icell],
                     &peribdry_pair_node_ind[peribdry_node_ptr[ibdry]],
                     peribdry_node_ptr[ibdry + 1] - peribdry_node_ptr[ibdry])) {
          peribdry_pair_owner_cell[ibdry] = cell_global_index[icell];
          break;
        }
      }
    }

    for (int ibdry = 0; ibdry < num_total_peribdries; ibdry++)
      if (peribdry_pair_owner_cell[ibdry] < 0)
        ERROR_MESSAGE("Owner cell of paired periodic boundary is not found.\n");

    std::unordered_set<int> local_cell_global_index;
    local_cell_global_index.insert(cell_global_index.begin(),
                                   cell_global_index.begin() + num_cells);

    enum PeribdryType : char { BDRY_INNER, BDRY_OUTER, BDRY_REDUNT };
    std::vector<PeribdryType> peribdry_types(num_total_peribdries,
                                             PeribdryType::BDRY_OUTER);
    for (int ibdry = 0; ibdry < num_total_peribdries; ibdry++) {
      PeribdryType& peribdry_type = peribdry_types[ibdry];
      const int& owner_cell = peribdry_owner_cell[ibdry];
      const int& pair_owner_cell = peribdry_pair_owner_cell[ibdry];
      if (local_cell_global_index.find(owner_cell) !=
              local_cell_global_index.end() &&
          local_cell_global_index.find(pair_owner_cell) !=
              local_cell_global_index.end())
        peribdry_type = (owner_cell < pair_owner_cell)
                            ? PeribdryType::BDRY_INNER
                            : PeribdryType::BDRY_REDUNT;
    }
    local_cell_global_index.clear();

    int ind = 0;
    std::vector<int> mapped_from(num_total_peribdries, -1);
    for (int ibdry = 0; ibdry < num_total_peribdries; ibdry++)
      if (peribdry_types[ibdry] == PeribdryType::BDRY_INNER)
        mapped_from[ind++] = ibdry;
    num_inner_peribdries = ind;

    for (int ibdry = 0; ibdry < num_total_peribdries; ibdry++)
      if (peribdry_types[ibdry] == PeribdryType::BDRY_OUTER)
        mapped_from[ind++] = ibdry;
    num_peribdries = ind;

    for (int ibdry = 0; ibdry < num_total_peribdries; ibdry++)
      if (peribdry_types[ibdry] == PeribdryType::BDRY_REDUNT)
        mapped_from[ind++] = ibdry;
    peribdry_types.clear();

    peribdry_neighbor_cell = std::move(peribdry_pair_owner_cell);
    ReorderData(mapped_from, peribdry_global_index, 1);
    ReorderData(mapped_from, peribdry_tag, 1);
    ReorderData(mapped_from, peribdry_owner_cell, 1);
    ReorderData(mapped_from, peribdry_neighbor_cell, 1);
    std::vector<int> peribdry_pair_node_ptr = peribdry_node_ptr;
    ReorderData(mapped_from, peribdry_node_ptr, peribdry_node_ind);
    ReorderData(mapped_from, peribdry_pair_node_ptr, peribdry_pair_node_ind);
    peribdry_pair_node_ptr.clear();
    mapped_from.clear();
  }
}

void GridBuilder::ConstructFacetype(
    const std::unordered_map<int, int>& cell_global_mapping,
    const std::vector<ElemType>& cell_elemtype,
    const std::vector<int>& cell_node_ptr,
    const std::vector<int>& cell_node_ind,
    const std::vector<int>& face_owner_cell,
    const std::vector<int>& face_node_ptr,
    const std::vector<int>& face_node_ind,
    std::vector<int>& face_owner_type) const {
  face_owner_type.clear();
  const int num_faces = static_cast<int>(face_owner_cell.size());
  face_owner_type.resize(num_faces, -1);

  for (int iface = 0; iface < num_faces; iface++) {
    auto&& iterator = cell_global_mapping.find(face_owner_cell[iface]);
    if (iterator == cell_global_mapping.end())
      ERROR_MESSAGE("The cell containing the target face is not found.\n");
    const int& icell = iterator->second;
    const Element* element = DENEB_ELEMENT->GetElement(cell_elemtype[icell]);
    face_owner_type[iface] = FindFacetype(
        element->GetFacetypeNodes(1), &cell_node_ind[cell_node_ptr[icell]],
        &face_node_ind[face_node_ptr[iface]],
        face_node_ptr[iface + 1] - face_node_ptr[iface]);
  }
}
}  // namespace deneb