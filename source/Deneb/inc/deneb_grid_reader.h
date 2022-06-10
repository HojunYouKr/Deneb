#pragma once

#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "deneb_element.h"

namespace deneb {
class GridReader {
 public:
  static std::shared_ptr<GridReader> GetGridReader(const std::string& name);

 protected:
  int dimension_;
  std::string grid_file_path_;

  int num_global_nodes_;
  int num_global_cells_;
  int num_global_bdries_;
  int num_global_peribdries_;
  std::vector<int> all_bdry_tag_;

 public:
  GridReader()
      : dimension_(0),
        grid_file_path_("No path input"),
        num_global_nodes_(0),
        num_global_cells_(0),
        num_global_bdries_(0),
        num_global_peribdries_(0){};
  virtual ~GridReader(){};

  inline int GetNumGlobalNodes(void) const { return num_global_nodes_; };
  inline int GetNumGlobalCells(void) const { return num_global_cells_; };
  inline int GetNumGlobalBdries(void) const { return num_global_bdries_; };
  inline int GetNumGlobalPeribdries(void) const {
    return num_global_peribdries_;
  };
  inline const std::vector<int>& GetAllBdryTag(void) const {
    return all_bdry_tag_;
  };

  virtual void Open(const int dimension, const std::string& grid_file_path) = 0;
  virtual void ReadCellData(const int start_index, const int end_index,
                            std::vector<ElemType>& cell_elemtype,
                            std::vector<int>& cell_global_index,
                            std::vector<int>& cell_node_ptr,
                            std::vector<int>& cell_node_ind) const = 0;
  virtual void ReadAllBdryData(std::vector<int>& all_bdry_global_index,
                               std::vector<int>& all_bdry_tag,
                               std::vector<int>& all_bdry_node_ptr,
                               std::vector<int>& all_bdry_node_ind) const = 0;
  virtual void ReadAllPeriodicData(
      std::vector<int>& all_periodic_nodes,
      std::vector<double>& all_periodic_node_coords) const = 0;
  virtual void ReadCellData(const std::unordered_set<int>& node_list,
                            std::vector<int>& cell_global_index,
                            std::vector<ElemType>& cell_elemtype,
                            std::vector<int>& cell_node_ptr,
                            std::vector<int>& cell_node_ind) const = 0;
  virtual void ReadCellData(std::unordered_map<int, int>& cell_mapping,
                            std::vector<int>& cell_order,
                            std::vector<int>& cell_subnode_ptr,
                            std::vector<int>& cell_subnode_ind) const = 0;
  virtual void ReadNodeData(std::unordered_map<int, int>& node_mapping,
                            std::vector<double>& node_coords) const = 0;
};

class GmshGridReader : public GridReader {
 protected:
  enum class GmshElemType : int {
    LINE_P1 = 1,
    LINE_P2 = 8,
    LINE_P3 = 26,
    LINE_P4 = 27,
    LINE_P5 = 28,
    LINE_P6 = 62,
    TRIS_P1 = 2,
    TRIS_P2 = 9,
    TRIS_P3 = 21,
    TRIS_P4 = 23,
    TRIS_P5 = 25,
    QUAD_P1 = 3,
    QUAD_P2 = 10,
    QUAD_P3 = 36,
    QUAD_P4 = 37,
    QUAD_P5 = 38,
    QUAD_P6 = 47,
    TETS_P1 = 4,
    TETS_P2 = 11,
    TETS_P3 = 29,
    TETS_P4 = 30,
    TETS_P5 = 31,
    HEXA_P1 = 5,
    HEXA_P2 = 12,
    HEXA_P3 = 92,
    HEXA_P4 = 93,
    HEXA_P5 = 94,
    PRIS_P1 = 6,
    PRIS_P2 = 13,
    PRIS_P3 = 90,
    PRIS_P4 = 91,
    PRIS_P5 = 106,
    PYRA_P1 = 7,
    PYRA_P2 = 14,
    PYRA_P3 = 118,
    PYRA_P4 = 119
  };

 public:
  GmshGridReader(){};
  virtual ~GmshGridReader(){};

  void Open(const int dimension, const std::string& grid_file_path);
  virtual void ReadCellData(const int start_index, const int end_index,
                            std::vector<ElemType>& cell_elemtype,
                            std::vector<int>& cell_global_index,
                            std::vector<int>& cell_node_ptr,
                            std::vector<int>& cell_node_ind) const;
  virtual void ReadAllBdryData(std::vector<int>& all_bdry_global_index,
                               std::vector<int>& all_bdry_tag,
                               std::vector<int>& all_bdry_node_ptr,
                               std::vector<int>& all_bdry_node_ind) const;
  virtual void ReadAllPeriodicData(
      std::vector<int>& all_periodic_nodes,
      std::vector<double>& all_periodic_node_coords) const;
  virtual void ReadCellData(const std::unordered_set<int>& node_list,
                            std::vector<int>& cell_global_index,
                            std::vector<ElemType>& cell_elemtype,
                            std::vector<int>& cell_node_ptr,
                            std::vector<int>& cell_node_ind) const;
  virtual void ReadCellData(std::unordered_map<int, int>& cell_mapping,
                            std::vector<int>& cell_order,
                            std::vector<int>& cell_subnode_ptr,
                            std::vector<int>& cell_subnode_ind) const;
  virtual void ReadNodeData(std::unordered_map<int, int>& node_mapping,
                            std::vector<double>& node_coords) const;

 protected:
  virtual void GenerateGmshFile(void) const {};
  int GetGmshElemOrder(const GmshElemType gmshtype) const;
  ElemType GetElemTypeFromGmshType(const GmshElemType gmshtype) const;
};

class Gmsh2DTrisGridReader : public GmshGridReader {
 public:
  Gmsh2DTrisGridReader(){};
  virtual ~Gmsh2DTrisGridReader(){};

 protected:
  virtual void GenerateGmshFile(void) const;
};

class Gmsh2DQuadGridReader : public GmshGridReader {
 public:
  Gmsh2DQuadGridReader(){};
  virtual ~Gmsh2DQuadGridReader(){};

 protected:
  virtual void GenerateGmshFile(void) const;
};

class Gmsh3DHexaGridReader : public GmshGridReader {
 public:
  Gmsh3DHexaGridReader(){};
  virtual ~Gmsh3DHexaGridReader(){};

 protected:
  virtual void GenerateGmshFile(void) const;
};
}  // namespace deneb
