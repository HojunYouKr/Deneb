#pragma once

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>

#include "avocado.h"

#define DENEB_ELEMENT_NAME element_global_ptr
#define DENEB_ELEMENT deneb::DENEB_ELEMENT_NAME
#define DENEB_ELEMENT_INITIALIZE(resource_dir)        \
  DENEB_ELEMENT = std::make_shared<deneb::Element>(); \
  DENEB_ELEMENT->LoadElementData(resource_dir)
#define DENEB_ELEMENT_FINALIZE() \
  DENEB_ELEMENT->Clear();        \
  DENEB_ELEMENT.reset()

#define ElemType Element::Type
#define ElemTypeHash Element::ElemtypeHash

namespace deneb {
class Element {
 public:
  enum class Type : int {
    LINE = 0,
    TRIS = 1,
    QUAD = 2,
    TETS = 3,
    HEXA = 4,
    PRIS = 5,
    PYRA = 6
  };

  struct ElemtypeHash {
    size_t operator()(const ElemType& elemtype) const {
      return std::hash<int>()(static_cast<int>(elemtype));
    }
  };

  static const std::string GetElementName(const ElemType elemtype);
  static bool IsFace(const int dimension, const ElemType elemtype);
  static int GetDimension(const ElemType elemtype);

 private:
  std::unordered_map<ElemType, Element*, ElemtypeHash> registry_;

 public:
  Element(){};
  ~Element() { Clear(); };

  const Element* GetElement(const ElemType elemtype) const;
  void GetSubcell(const ElemType elemtype, const int order,
                  std::vector<double>& coords,
                  std::vector<int>& connectivity) const;
  inline void Clear(void) {
    for (auto&& iterator : registry_) delete iterator.second;
    registry_.clear();
  };
  void LoadElementData(const std::string& directory);

 private:
  // input
  ElemType elemtype_;
  int elemorder_support_;

  // hard-code
  int dimension_;
  int num_faces_;
  double volume_;
  std::vector<std::vector<int>> face_nodes_;  // P1 by face
  std::vector<ElemType> face_elemtype_;

  // read from file
  std::vector<int> num_nodes_;                    // by elemorder
  std::vector<std::vector<double>> node_coords_;  // by elemorder

  // compute automatically
  int num_facetypes_;
  std::vector<double> face_areas_;     // by face
  std::vector<double> normals_;        // by face
  std::vector<int> facetype_to_face_;  // by facetype
  std::vector<std::vector<std::vector<int>>>
      facetype_nodes_;  // by elemorder, facetype

 public:
  inline ElemType GetElemType() const { return elemtype_; };
  inline int GetNumFaces() const { return num_faces_; };
  inline double GetVolume() const { return volume_; };
  inline const std::vector<std::vector<int>>& GetFaceNodes() const {
    return face_nodes_;
  };
  inline ElemType GetFaceElemtype(const int facetype) const {
    ASSERT(facetype < num_facetypes_);
    return face_elemtype_[facetype_to_face_[facetype]];
  };

  inline int GetNumNodes(const int elemorder) const {
    ASSERT(elemorder <= elemorder_support_);
    return num_nodes_[elemorder - 1];
  };
  inline const std::vector<double>& GetNodesCoord(const int elemorder) const {
    ASSERT(elemorder <= elemorder_support_);
    return node_coords_[elemorder - 1];
  };

  inline double GetFacetypeArea(const int facetype) const {
    ASSERT(facetype < num_facetypes_);
    return face_areas_[facetype_to_face_[facetype]];
  };
  inline const double* GetFacetypeNormal(const int facetype) const {
    ASSERT(facetype < num_facetypes_);
    return &normals_[facetype_to_face_[facetype] * dimension_];
  };
  inline const std::vector<std::vector<int>>& GetFacetypeNodes(
      const int elemorder) const {
    ASSERT(elemorder <= elemorder_support_);
    return facetype_nodes_[elemorder - 1];
  };

  void TransformToPhyCoords(const int facetype, const int num_nodes,
                            const double* ref_coords, double* phy_coords) const;

 private:
  void ReadFromFile(const std::string& filename_template);
  void ComputeFaceData(void);

  double ComputeFaceArea(const int dimension,
                         const std::vector<double>& face_coords) const;
  void ComputeNormal(const int dimension,
                     const std::vector<double>& face_coords,
                     const std::vector<double>& inside_coord,
                     double* normal) const;
  const std::vector<std::vector<int>> GenerateFacetype(
      const std::vector<int>& node_index) const;
  void TransformToPhyCoords(const ElemType elemtype, const int num_coords,
                            const double* ref_coords, double* phy_coords,
                            const double* phyelem_coords) const;

  const std::vector<int> GetType1Connects(const std::vector<int>& pts) const;
  const std::vector<int> GetType2Connects(const std::vector<int>& pts) const;
  const std::vector<int> GetType3Connects(const std::vector<int>& pts) const;
  const std::vector<int> GetType4Connects(const std::vector<int>& pts) const;

  void GetLineSubcell(const int order, std::vector<double>& coords,
                      std::vector<int>& connectivity) const;
  void GetTrisSubcell(const int order, std::vector<double>& coords,
                      std::vector<int>& connectivity) const;
  void GetQuadSubcell(const int order, std::vector<double>& coords,
                      std::vector<int>& connectivity) const;
  void GetTetsSubcell(const int order, std::vector<double>& coords,
                      std::vector<int>& connectivity) const;
  void GetHexaSubcell(const int order, std::vector<double>& coords,
                      std::vector<int>& connectivity) const;
  void GetPrisSubcell(const int order, std::vector<double>& coords,
                      std::vector<int>& connectivity) const;
  void GetPyraSubcell(const int order, std::vector<double>& coords,
                      std::vector<int>& connectivity) const;
};  // namespace deneb
extern std::shared_ptr<Element> DENEB_ELEMENT_NAME;
}  // namespace deneb