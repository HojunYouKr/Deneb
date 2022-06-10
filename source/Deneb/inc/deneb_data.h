#pragma once

#include <memory>
#include <unordered_map>
#include <vector>

#include "deneb_element.h"

#define DENEB_DATA_NAME data_global_ptr
#define DENEB_DATA deneb::DENEB_DATA_NAME
#define DENEB_DATA_INITIALIZE() DENEB_DATA = std::make_shared<deneb::Data>()
#define DENEB_DATA_FINALIZE() DENEB_DATA.reset()

#define GETFUNC(name) Get##name##()
#define DEFINE_GETFUNC(type, name) \
  inline type GETFUNC(name) const { return (name); }

namespace deneb {
class Jacobian;
class Basis;

class Data {
 private:
  int order_;
  int volume_flux_order_;
  int surface_flux_order_;
  int num_bases_;

  int dimension_;
  bool has_source_term_;

  // node
  int num_global_nodes_;
  int num_total_nodes_;
  int num_nodes_;
  std::vector<int> node_global_index_;
  std::vector<double> node_coords_;

  // cell
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

  std::vector<const Element*> cell_element_;
  std::vector<std::shared_ptr<Jacobian>> cell_jacobians_;
  std::vector<std::shared_ptr<Basis>> cell_basis_;
  std::vector<int> num_cell_points_;
  std::vector<double> cell_volumes_;
  std::vector<double> cell_proj_volumes_;
  std::vector<std::vector<double>> cell_points_;
  std::vector<std::vector<double>> cell_basis_value_;
  std::vector<std::vector<double>> cell_basis_grad_value_;
  std::vector<std::vector<double>> cell_coefficients_;
  std::vector<std::vector<double>> cell_source_coefficients_;

  // face
  int num_faces_;
  int num_inner_faces_;
  std::vector<int> face_owner_cell_;
  std::vector<int> face_neighbor_cell_;
  std::vector<int> face_owner_type_;
  std::vector<int> face_neighbor_type_;

  std::vector<int> num_face_points_;
  std::vector<std::vector<double>> face_points_;
  std::vector<std::vector<double>> face_normals_;
  std::vector<std::vector<double>> face_owner_basis_value_;
  std::vector<std::vector<double>> face_neighbor_basis_value_;
  std::vector<std::vector<double>> face_owner_basis_grad_value_;
  std::vector<std::vector<double>> face_neighbor_basis_grad_value_;
  std::vector<std::vector<double>> face_owner_coefficients_;
  std::vector<std::vector<double>> face_neighbor_coefficients_;

  // boundary
  int num_global_bdries_;
  int num_bdries_;
  std::vector<int> bdry_global_index_;
  std::vector<int> bdry_tag_;
  std::vector<int> bdry_owner_cell_;
  std::vector<int> bdry_owner_type_;

  std::vector<int> num_bdry_points_;
  std::vector<std::vector<double>> bdry_points_;
  std::vector<std::vector<double>> bdry_normals_;
  std::vector<std::vector<double>> bdry_owner_basis_value_;
  std::vector<std::vector<double>> bdry_owner_basis_grad_value_;
  std::vector<std::vector<double>> bdry_owner_coefficients_;

  // periodic boundary
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

  std::vector<int> num_peribdry_points_;
  std::vector<std::vector<double>> peribdry_points_;
  std::vector<std::vector<double>> peribdry_normals_;
  std::vector<std::vector<double>> peribdry_owner_basis_value_;
  std::vector<std::vector<double>> peribdry_neighbor_basis_value_;
  std::vector<std::vector<double>> peribdry_owner_basis_grad_value_;
  std::vector<std::vector<double>> peribdry_neighbor_basis_grad_value_;
  std::vector<std::vector<double>> peribdry_owner_coefficients_;
  std::vector<std::vector<double>> peribdry_neighbor_coefficients_;

  // periodic matching nodes
  std::unordered_map<int, std::vector<int>>
      periodic_matching_global_node_index_;

  // solution communicate
  std::vector<std::vector<int>> outer_send_cell_list_;
  std::vector<std::vector<int>> outer_recv_cell_list_;

  // hMLP stencil communicate
  std::vector<std::vector<int>> total_send_cell_list_;
  std::vector<std::vector<int>> total_recv_cell_list_;

  // implicit data
  std::vector<int> subface_ptr_;
  std::vector<double> subface_sign_;
  std::vector<int> subface_neighbor_cell_;
  std::vector<int> subface_num_points_;
  std::vector<double*> subface_owner_coefficients_;
  std::vector<double*> subface_owner_basis_value_;
  std::vector<double*> subface_neighbor_basis_value_;
  std::vector<int> subbdry_ptr_;
  std::vector<int> subbdry_ind_;
  std::vector<int> mat_index_;

 public:
  Data();
  ~Data();

 private:
  void BuildCellData(const int icell);
  void BuildFaceData(const int iface);
  void BuildBdryData(const int ibdry);
  void BuildPeribdryData(const int ibdry);
  void CombinePeribdryToFaceData(void);
  void BuildImplicitData(void);
  double ComputeProjVolume(const std::vector<double>& coords,
                           const int idim) const;

  template <typename T>
  void CombineData(std::vector<T>& A, std::vector<T>& B,
                   const std::vector<int>& rules);

 public:
  void BuildData(void);
  inline int GetNumBases(const int& order) const {
    int num_bases = order + 1;
    for (int idim = 1; idim < dimension_; idim++)
      num_bases = num_bases * (order + idim + 1) / (idim + 1);
    return num_bases;
  };
  void GetCellQuadrature(const int icell, const int order,
                         std::vector<double>& quad_points,
                         std::vector<double>& quad_weights) const;
  void GetCellBasisValues(const int icell, const std::vector<double>& coords,
                          std::vector<double>& basis_value) const;
  void GetCellBasisGradValues(const int icell,
                              const std::vector<double>& coords,
                              std::vector<double>& basis_grad_value) const;
  void GetCellRefToPhyCoords(const int icell,
                             const std::vector<double>& ref_coords,
                             std::vector<double>& phy_coords) const;
  void GetFaceCoords(const int icell, const int face_type,
                     std::vector<double>& coords) const;

  inline const std::unordered_map<int, std::vector<int>>&
  GetPeriodicMatchingGlobalNodeIndex(void) const {
    return periodic_matching_global_node_index_;
  };

  // Get functions
  inline const int& GetOrder() const { return order_; };
  inline const int& GetNumBases() const { return num_bases_; };
  inline const std::vector<std::vector<int>>& GetOuterSendCellList() const {
    return outer_send_cell_list_;
  };
  inline const std::vector<std::vector<int>>& GetOuterRecvCellList() const {
    return outer_recv_cell_list_;
  };
  inline const std::vector<std::vector<int>>& GetTotalSendCellList() const {
    return total_send_cell_list_;
  };
  inline const std::vector<std::vector<int>>& GetTotalRecvCellList() const {
    return total_recv_cell_list_;
  };

  // cell
  inline const int& GetNumGlobalNodes() const { return num_global_nodes_; };
  inline const int& GetNumTotalNodes() const { return num_total_nodes_; };
  inline const int& GetNumNodes() const { return num_nodes_; };
  inline const std::vector<int>& GetNodeGlobalIndex() const {
    return node_global_index_;
  };
  inline const std::vector<double>& GetNodeCoords() const {
    return node_coords_;
  };
  inline const int& GetNumGlobalCells() const { return num_global_cells_; };
  inline const int& GetNumTotalCells() const { return num_total_cells_; };
  inline const int& GetNumOuterCells() const { return num_outer_cells_; };
  inline const int& GetNumCells() const { return num_cells_; };
  inline const std::vector<int>& GetCellGlobalIndex() const {
    return cell_global_index_;
  };
  inline const std::vector<ElemType>& GetCellElemtype() const {
    return cell_elemtype_;
  };
  inline const std::vector<int>& GetCellOrder() const { return cell_order_; };
  inline const std::vector<int>& GetCellNodePtr() const {
    return cell_node_ptr_;
  };
  inline const std::vector<int>& GetCellNodeInd() const {
    return cell_node_ind_;
  };
  inline const std::vector<int>& GetCellSubnodePtr() const {
    return cell_subnode_ptr_;
  };
  inline const std::vector<int>& GetCellSubnodeInd() const {
    return cell_subnode_ind_;
  };
  inline const std::vector<const Element*>& GetCellElement() const {
    return cell_element_;
  };
  inline const std::vector<std::shared_ptr<Jacobian>>& GetCellJacobians()
      const {
    return cell_jacobians_;
  };
  inline const std::vector<std::shared_ptr<Basis>>& GetCellBasis() const {
    return cell_basis_;
  };
  inline const std::vector<int>& GetNumCellPoints() const {
    return num_cell_points_;
  };
  inline const std::vector<double>& GetCellVolumes() const {
    return cell_volumes_;
  };
  inline const std::vector<double>& GetCellProjVolumes() const {
    return cell_proj_volumes_;
  };
  inline const std::vector<std::vector<double>>& GetCellPoints() const {
    return cell_points_;
  };
  inline const std::vector<std::vector<double>>& GetCellBasisValue() const {
    return cell_basis_value_;
  };
  inline const std::vector<std::vector<double>>& GetCellBasisGradValue() const {
    return cell_basis_grad_value_;
  };
  inline const std::vector<std::vector<double>>& GetCellCoefficients() const {
    return cell_coefficients_;
  };
  inline const std::vector<std::vector<double>>& GetCellSourceCoefficients()
      const {
    return cell_source_coefficients_;
  };

  // face
  inline const int& GetNumFaces() const { return num_faces_; };
  inline const int& GetNumInnerFaces() const { return num_inner_faces_; };
  inline const std::vector<int>& GetFaceOwnerCell() const {
    return face_owner_cell_;
  };
  inline const std::vector<int>& GetFaceNeighborCell() const {
    return face_neighbor_cell_;
  };
  inline const std::vector<int>& GetFaceOwnerType() const {
    return face_owner_type_;
  };
  inline const std::vector<int>& GetFaceNeighborType() const {
    return face_neighbor_type_;
  };
  inline const std::vector<int>& GetNumFacePoints() const {
    return num_face_points_;
  };
  inline const std::vector<std::vector<double>>& GetFacePoints() const {
    return face_points_;
  };
  inline const std::vector<std::vector<double>>& GetFaceNormals() const {
    return face_normals_;
  };
  inline const std::vector<std::vector<double>>& GetFaceOwnerBasisValue()
      const {
    return face_owner_basis_value_;
  };
  inline const std::vector<std::vector<double>>& GetFaceNeighborBasisValue()
      const {
    return face_neighbor_basis_value_;
  };
  inline const std::vector<std::vector<double>>& GetFaceOwnerBasisGradValue()
      const {
    return face_owner_basis_grad_value_;
  };
  inline const std::vector<std::vector<double>>& GetFaceNeighborBasisGradValue()
      const {
    return face_neighbor_basis_grad_value_;
  };
  inline const std::vector<std::vector<double>>& GetFaceOwnerCoefficients()
      const {
    return face_owner_coefficients_;
  };
  inline const std::vector<std::vector<double>>& GetFaceNeighborCoefficients()
      const {
    return face_neighbor_coefficients_;
  };

  // boundary
  inline const int& GetNumGlobalBdries() const { return num_global_bdries_; };
  inline const int& GetNumBdries() const { return num_bdries_; };
  inline const std::vector<int>& GetBdryGlobalIndex() const {
    return bdry_global_index_;
  };
  inline const std::vector<int>& GetBdryTag() const { return bdry_tag_; };
  inline const std::vector<int>& GetBdryOwnerCell() const {
    return bdry_owner_cell_;
  };
  inline const std::vector<int>& GetBdryOwnerType() const {
    return bdry_owner_type_;
  };
  inline const std::vector<int>& GetNumBdryPoints() const {
    return num_bdry_points_;
  };
  inline const std::vector<std::vector<double>>& GetBdryPoints() const {
    return bdry_points_;
  };
  inline const std::vector<std::vector<double>>& GetBdryNormals() const {
    return bdry_normals_;
  };
  inline const std::vector<std::vector<double>>& GetBdryOwnerBasisValue()
      const {
    return bdry_owner_basis_value_;
  };
  inline const std::vector<std::vector<double>>& GetBdryOwnerBasisGradValue()
      const {
    return bdry_owner_basis_grad_value_;
  };
  inline const std::vector<std::vector<double>>& GetBdryOwnerCoefficients()
      const {
    return bdry_owner_coefficients_;
  };

  // periodic boundary
  inline const int& GetNumGlobalPeribdries() const {
    return num_global_peribdries_;
  };
  inline const int& GetNumTotalPeribdries() const {
    return num_total_peribdries_;
  };
  inline const int& GetNumPeribdries() const { return num_peribdries_; };
  inline const int& GetNumInnerPeribdries() const {
    return num_inner_peribdries_;
  };
  inline const std::vector<int>& GetPeribdryGlobalIndex() const {
    return peribdry_global_index_;
  };
  inline const std::vector<int>& GetPeribdryTag() const {
    return peribdry_tag_;
  };
  inline const std::vector<int>& GetPeribdryOwnerCell() const {
    return peribdry_owner_cell_;
  };
  inline const std::vector<int>& GetPeribdryNeighborCell() const {
    return peribdry_neighbor_cell_;
  };
  inline const std::vector<int>& GetPeribdryOwnerType() const {
    return peribdry_owner_type_;
  };
  inline const std::vector<int>& GetPeribdryNeighborType() const {
    return peribdry_neighbor_type_;
  };
  inline const std::vector<int>& GetNumPeribdryPoints() const {
    return num_peribdry_points_;
  };
  inline const std::vector<std::vector<double>>& GetPeribdryPoints() const {
    return peribdry_points_;
  };
  inline const std::vector<std::vector<double>>& GetPeribdryNormals() const {
    return peribdry_normals_;
  };
  inline const std::vector<std::vector<double>>& GetPeribdryOwnerBasisValue()
      const {
    return peribdry_owner_basis_value_;
  };
  inline const std::vector<std::vector<double>>& GetPeribdryNeighborBasisValue()
      const {
    return peribdry_neighbor_basis_value_;
  };
  inline const std::vector<std::vector<double>>&
  GetPeribdryOwnerBasisGradValue() const {
    return peribdry_owner_basis_grad_value_;
  };
  inline const std::vector<std::vector<double>>&
  GetPeribdryNeighborBasisGradValue() const {
    return peribdry_neighbor_basis_grad_value_;
  };
  inline const std::vector<std::vector<double>>& GetPeribdryOwnerCoefficients()
      const {
    return peribdry_owner_coefficients_;
  };
  inline const std::vector<std::vector<double>>&
  GetPeribdryNeighborCoefficients() const {
    return peribdry_neighbor_coefficients_;
  };

  // implicit data
  inline const std::vector<int>& GetSubfacePtr() const { return subface_ptr_; };
  inline const std::vector<double>& GetSubfaceSign() const {
    return subface_sign_;
  };
  inline const std::vector<int>& GetSubfaceNeighborCell() const {
    return subface_neighbor_cell_;
  };
  inline const std::vector<int>& GetSubfaceNumPoints() const {
    return subface_num_points_;
  };
  inline const std::vector<double*>& GetSubfaceOwnerCoefficients() const {
    return subface_owner_coefficients_;
  };
  inline const std::vector<double*>& GetSubfaceOwnerBasisValue() const {
    return subface_owner_basis_value_;
  };
  inline const std::vector<double*>& GetSubfaceNeighborBasisValue() const {
    return subface_neighbor_basis_value_;
  };
  inline const std::vector<int>& GetSubbdryPtr() const { return subbdry_ptr_; };
  inline const std::vector<int>& GetSubbdryInd() const { return subbdry_ind_; };
  inline const std::vector<int>& GetMatIndex() const { return mat_index_; };
};
extern std::shared_ptr<Data> DENEB_DATA_NAME;

template <typename T>
void Data::CombineData(std::vector<T>& A, std::vector<T>& B,
                       const std::vector<int>& rules) {
  const int num_rules = static_cast<int>(rules.size()) / 2 - 1;
  std::vector<T> combine(rules[num_rules * 2] + rules[num_rules * 2 + 1]);
  int ind = 0;
  for (int i = 0; i < num_rules; i++) {
    for (int j = rules[2 * i], len = rules[2 * i + 2]; j < len; j++)
      combine[ind++] = std::move(A[j]);
    for (int j = rules[2 * i + 1], len = rules[2 * i + 3]; j < len; j++)
      combine[ind++] = std::move(B[j]);
  }
  A = std::move(combine);
}
}  // namespace deneb
