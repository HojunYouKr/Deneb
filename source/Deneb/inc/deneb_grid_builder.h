#pragma once

#include <string>
#include <unordered_set>
#include <unordered_map>
#include <vector>

#include "avocado.h"
#include "deneb_element.h"

namespace deneb {
constexpr const double kPeriodicEps = 1.0E-10;
class GridBuilder {
  friend class Data;

 private:
  int dimension_;

  int num_global_nodes_;
  int num_total_nodes_;
  int num_nodes_;
  std::vector<int> node_global_index_;
  std::vector<double> node_coords_;

  int num_global_cells_;
  int num_total_cells_;
  int num_outer_cells_;
  int num_cells_;
  std::vector<int> cell_global_index_;
  std::vector<ElemType> cell_elemtype_;
  std::vector<int> cell_order_;
  std::vector<int> cell_node_ptr_;
  std::vector<int> cell_node_ind_;
  std::vector<int> cell_subnode_ptr_;
  std::vector<int> cell_subnode_ind_;

  int num_global_bdries_;
  int num_bdries_;
  std::vector<int> bdry_global_index_;
  std::vector<int> bdry_tag_;
  std::vector<int> bdry_owner_cell_;
  std::vector<int> bdry_owner_type_;

  int num_global_peribdries_;
  int num_total_peribdries_;
  int num_peribdries_;
  int num_inner_peribdries_;
  std::vector<int> peribdry_global_index_;
  std::vector<int> peribdry_tag_;
  std::vector<int> peribdry_owner_cell_;
  std::vector<int> peribdry_neighbor_cell_;
  std::vector<int> peribdry_owner_type_;
  std::vector<int> peribdry_neighbor_type_;

  int num_faces_;
  int num_inner_faces_;
  std::vector<int> face_owner_cell_;
  std::vector<int> face_neighbor_cell_;
  std::vector<int> face_owner_type_;
  std::vector<int> face_neighbor_type_;

  std::unordered_map<int, std::vector<int>>
      periodic_matching_global_node_index_;

  // solution communicate
  std::vector<std::vector<int>> outer_send_cell_list_;
  std::vector<std::vector<int>> outer_recv_cell_list_;

  // hMLP stencil communicate
  std::vector<std::vector<int>> total_send_cell_list_;
  std::vector<std::vector<int>> total_recv_cell_list_;

 public:
  GridBuilder(const int dimension) : dimension_(dimension){};
  ~GridBuilder(){};

 public:
  void BuildGrid();

 private:
  bool IsParent(const int* cell_node_index, const int cell_num_nodes,
                const int* face_node_index, const int face_num_nodes) const;

  bool IsSame(const int* index1, const int num_index1, const int* index2,
              const int num_index2) const;

  int FindFacetype(const std::vector<std::vector<int>>& facetype_nodes_index,
                   const int* cell_nodes_index, const int* face_nodes_index,
                   const int num_face_nodes) const;

  // data = mapping[data]
  void RemapData(const std::unordered_map<int, int>& mapping,
                 std::vector<int>& data) const;

  // data[i] = data[mapped_from[i]]
  template <typename T>
  void ReorderData(const std::vector<int>& mapped_from, std::vector<T>& data,
                   const int& datasize) const;

  // data[i] = data[mapped_from[i]]
  template <typename T>
  void ReorderData(const std::vector<int>& mapped_from, std::vector<int>& ptr,
                   std::vector<T>& data) const;

  void ConstructCommunication(
      const std::vector<int>& local_data, const std::vector<int>& target_data,
      std::vector<std::vector<int>>& send_locations,
      std::vector<std::vector<int>>& recv_locations) const;

  void SetCellColor(const std::vector<int>& cell_node_ptr,
                    const std::vector<int>& cell_node_ind,
                    std::vector<int>& cell_color) const;

  template <typename T>
  void RedistributeCellData(const std::vector<std::vector<int>>& send_cell_list,
                            std::vector<T>& cell_data) const;
  template <typename T>
  void RedistributeCellData(const std::vector<std::vector<int>>& send_cell_list,
                            std::vector<int>& cell_data_ptr,
                            std::vector<T>& cell_data) const;

  void ConstructFaceData(const std::vector<int>& cell_global_index,
                         const std::vector<ElemType>& cell_elemtype,
                         const std::vector<int>& cell_node_ptr,
                         const std::vector<int>& cell_node_ind,
                         std::vector<int>& face_owner_cell,
                         std::vector<int>& face_neighbor_cell,
                         std::vector<int>& face_node_ptr,
                         std::vector<int>& face_node_ind,
                         std::vector<int>& outer_face_owner_cell,
                         std::vector<int>& outer_face_node_ptr,
                         std::vector<int>& outer_face_node_ind) const;

  void ConstructBdryData(
      const std::vector<int>& all_bdry_global_index,
      const std::vector<int>& all_bdry_tag,
      const std::vector<int>& all_bdry_node_ptr,
      const std::vector<int>& all_bdry_node_ind,
      std::vector<int>& outer_face_owner_cell,
      std::vector<int>& outer_face_node_ptr,
      std::vector<int>& outer_face_node_ind,
      std::vector<int>& bdry_global_index, std::vector<int>& bdry_tag,
      std::vector<int>& bdry_owner_cell, std::vector<int>& bdry_node_ptr,
      std::vector<int>& bdry_node_ind, std::vector<int>& peribdry_global_index,
      std::vector<int>& peribdry_tag, std::vector<int>& peribdry_owner_cell,
      std::vector<int>& peribdry_node_ptr,
      std::vector<int>& peribdry_node_ind) const;

  void PeriodicMatching(const std::vector<int>& all_periodic_nodes,
                        const std::vector<double>& all_periodic_node_coords,
                        const std::vector<int>& peribdry_global_index,
                        const std::vector<int>& peribdry_tag,
                        const std::vector<int>& peribdry_node_ptr,
                        const std::vector<int>& peribdry_node_ind,
                        std::vector<int>& peribdry_pair_node_ind) const;

  void PeriodicMatching(const std::vector<int>& peribdry_node_ind,
                        const std::vector<int>& peribdry_pair_node_ind,
                        std::unordered_map<int, std::vector<int>>&
                            periodic_matching_nodes) const;

  void ConstructOuterGrid(
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
      std::vector<int>& peribdry_pair_node_ind) const;

  void ConstructFacetype(
      const std::unordered_map<int, int>& cell_global_mapping,
      const std::vector<ElemType>& cell_elemtype,
      const std::vector<int>& cell_node_ptr,
      const std::vector<int>& cell_node_ind,
      const std::vector<int>& face_owner_cell,
      const std::vector<int>& face_node_ptr,
      const std::vector<int>& face_node_ind,
      std::vector<int>& face_owner_type) const;
};

template <typename T>
void GridBuilder::ReorderData(const std::vector<int>& mapped_from,
                              std::vector<T>& data, const int& datasize) const {
  const size_t num_mapping = mapped_from.size();
  std::vector<T> data_temp;
  data_temp.reserve(num_mapping * datasize);
  for (size_t imap = 0; imap < num_mapping; imap++)
    for (size_t idata = 0; idata < datasize; idata++)
      data_temp.push_back(data[mapped_from[imap] * datasize + idata]);
  data = std::move(data_temp);
}
template <typename T>
void GridBuilder::ReorderData(const std::vector<int>& mapped_from,
                              std::vector<int>& ptr,
                              std::vector<T>& data) const {
  const size_t num_mapping = mapped_from.size();
  std::vector<int> ptr_temp;
  std::vector<int> data_temp;
  ptr_temp.reserve(num_mapping);
  ptr_temp.push_back(0);
  data_temp.reserve(num_mapping);
  for (size_t imap = 0; imap < num_mapping; imap++) {
    for (int i = ptr[mapped_from[imap]]; i < ptr[mapped_from[imap] + 1]; i++)
      data_temp.push_back(data[i]);
    ptr_temp.push_back(static_cast<int>(data_temp.size()));
  }
  ptr = std::move(ptr_temp);
  data = std::move(data_temp);
}

template <typename T>
void GridBuilder::RedistributeCellData(
    const std::vector<std::vector<int>>& send_cell_list,
    std::vector<T>& cell_data) const {
  std::vector<std::vector<T>> send_data;
  std::vector<std::vector<T>> recv_data;

  send_data.resize(NDOMAIN);
  for (int i = 0; i < NDOMAIN; i++) {
    send_data[i].reserve(send_cell_list[i].size());
    for (auto&& icell : send_cell_list[i])
      send_data[i].push_back(cell_data[icell]);
  }
  cell_data.clear();

  AVOCADO_MPI->CommunicateData(send_data, recv_data);
  send_data.clear();

  size_t size = 0;
  for (int i = 0; i < NDOMAIN; i++) size += recv_data[i].size();
  cell_data.reserve(size);
  for (int i = 0; i < NDOMAIN; i++)
    for (auto& data : recv_data[i]) cell_data.push_back(data);
  recv_data.clear();
}

template <typename T>
void GridBuilder::RedistributeCellData(
    const std::vector<std::vector<int>>& send_cell_list,
    std::vector<int>& cell_data_ptr, std::vector<T>& cell_data) const {
  std::vector<std::vector<T>> send_data;
  std::vector<std::vector<T>> recv_data;

  send_data.resize(NDOMAIN);
  for (int i = 0; i < NDOMAIN; i++) {
    int size = 0;
    for (auto&& icell : send_cell_list[i])
      size += (cell_data_ptr[icell + 1] - cell_data_ptr[icell]);
    send_data[i].reserve(size);
    for (auto&& icell : send_cell_list[i])
      for (int p = cell_data_ptr[icell]; p < cell_data_ptr[icell + 1]; p++)
        send_data[i].push_back(cell_data[p]);
  }
  cell_data.clear();

  AVOCADO_MPI->CommunicateData(send_data, recv_data);
  send_data.clear();

  size_t size = 0;
  for (int i = 0; i < NDOMAIN; i++) size += recv_data[i].size();
  cell_data.reserve(size);
  for (int i = 0; i < NDOMAIN; i++)
    for (auto&& data : recv_data[i]) cell_data.push_back(data);
  recv_data.clear();

  for (size_t i = 1, len = cell_data_ptr.size(); i < len; i++)
    cell_data_ptr[i - 1] = cell_data_ptr[i] - cell_data_ptr[i - 1];
  cell_data_ptr.pop_back();
  RedistributeCellData(send_cell_list, cell_data_ptr);
  for (size_t i = 1, len = cell_data_ptr.size(); i < len; i++)
    cell_data_ptr[i] += cell_data_ptr[i - 1];
  cell_data_ptr.insert(cell_data_ptr.begin(), 0);
}
}  // namespace deneb

// The global mesh is from the gridfile.
// The global mesh is partitioned into non-overlapped local meshes.
// Each processor owns a total mesh that consists of the local and layer meshes.
//		global mesh = sum(local meshes)
//		total mesh = local mesh + layer mesh
//
// Cells in a total mesh consist of local and layer cells.
//		local cells = cells in the local mesh
//		layer cells = cells in the layer mesh = outer cells + vertex
//cells 		total cells = cells in the total mesh = local cells + layer cells
// Layer cells consist of outer and vertex cells.
//		outer cells = cells adjacent to local cells
//		vertex cells = cells sharing vertex with local cells
//		layer cells = outer cells + vertex cells
// Cell data are arranged in order of [local, outer, vertex].
//		num_cells = the number of [local]
//		num_outer_cells = the number of [local, outer]
//		num_total_cells = the number of [local, outer, vertex]
//
// Nodes in a total mesh consist of local and outer nodes.
//		total nodes = nodes in the total mesh
//		local nodes = vertices of local cells
//		outer nodes = total nodes - local nodes
// Node data are arranged in order of [local, outer].
//		num_nodes = the number of [local]
//		num_total_nodes = the number of [local, outer]
//
// Each processor owns boundaries only contained in a local mesh.
//		num_bdries = the number of local boundaries
//
// Each processor owns periodic bounadries only contained in a local mesh.
// The periodic boundaries consist of inner, outer, redundunt boundaries.
//		inner periodic boundaries = boundaries that a neighbor cell is
//                                in the local mesh, and at the same
//                                time, the global index of the neighbor
//                                cell is larger than that of the owner cell
//    outer periodic bounadries = boundaries that a neighbor cell is not
//                                in the local mesh
//    redundunt periodic boundaries = boundaries matching with inner periodic
//                                    boundaries
// Periodic boundary data are arranged in order of [inner, outer, redundunt].
//		num_inner_peribdries = the number of [inner]
//		num_peribdries = the number of [inner, outer]
//		num_total_peribdries = the number of [inner, outer, redundunt]
//
// Faces owned by each processor consist of inner and outer faces.
//		inner faces = faces that both owner and neighbor cells are in
// the  local mesh 		outer faces = faces that a neighbor cells is not
// in the local  mesh
// Face data are arranged in order of [inner, outer].
//		num_inner_faces = the number of [inner]
//		num_faces = the number of [inner, outer]