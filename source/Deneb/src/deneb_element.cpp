#include "deneb_element.h"

#include <fstream>

namespace deneb {
std::shared_ptr<Element> DENEB_ELEMENT_NAME = nullptr;
const std::string Element::GetElementName(const ElemType elemtype) {
  switch (elemtype) {
    case ElemType::LINE:
      return "Line";
    case ElemType::TRIS:
      return "Tris";
    case ElemType::QUAD:
      return "Quad";
    case ElemType::TETS:
      return "Tets";
    case ElemType::HEXA:
      return "Hexa";
    case ElemType::PRIS:
      return "Pris";
    case ElemType::PYRA:
      return "Pyra";
  }
  ERROR_MESSAGE("Illegal element type: " +
                std::to_string(static_cast<int>(elemtype)) + "\n");
  return std::string();
}
bool Element::IsFace(const int dimension, const ElemType elemtype) {
  if (dimension == 2 && elemtype == ElemType::LINE)
    return true;
  else if (dimension == 3 &&
           (elemtype == ElemType::TRIS || elemtype == ElemType::QUAD))
    return true;
  return false;
}
int Element::GetDimension(const ElemType elemtype) {
  switch (elemtype) {
    case ElemType::LINE:
      return 1;
    case ElemType::TRIS:
    case ElemType::QUAD:
      return 2;
    case ElemType::TETS:
    case ElemType::HEXA:
    case ElemType::PRIS:
    case ElemType::PYRA:
      return 3;
  }
  ERROR_MESSAGE("Illegal element type: " +
                std::to_string(static_cast<int>(elemtype)) + "\n");
  return -1;
}

void Element::LoadElementData(const std::string& directory) {
  const std::string filepath = directory + "Element/&(Elem)_P&(Order).msh";

  Element* element;
  {
    element = new Element();
    element->elemtype_ = ElemType::LINE;
    element->elemorder_support_ = 6;
    element->dimension_ = 1;
    element->volume_ = 2.0;
    std::string filename = filepath;
    avocado::ReplaceString(filename, "&(Elem)", "LINE");
    element->ReadFromFile(filename);
    registry_[element->elemtype_] = element;
  }
  {
    element = new Element();
    element->elemtype_ = ElemType::TRIS;
    element->elemorder_support_ = 5;
    element->dimension_ = 2;
    element->volume_ = 2.0;
    std::string filename = filepath;
    avocado::ReplaceString(filename, "&(Elem)", "TRIS");
    element->ReadFromFile(filename);

    element->num_faces_ = 3;
    element->face_nodes_ = {{0, 1}, {1, 2}, {2, 0}};
    element->face_elemtype_.resize(3, ElemType::LINE);
    element->ComputeFaceData();
    registry_[element->elemtype_] = element;
  }
  {
    element = new Element();
    element->elemtype_ = ElemType::QUAD;
    element->elemorder_support_ = 6;
    element->dimension_ = 2;
    element->volume_ = 4.0;
    std::string filename = filepath;
    avocado::ReplaceString(filename, "&(Elem)", "QUAD");
    element->ReadFromFile(filename);

    element->num_faces_ = 4;
    element->face_nodes_ = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    element->face_elemtype_.resize(4, ElemType::LINE);
    element->ComputeFaceData();
    registry_[element->elemtype_] = element;
  }
  {
    element = new Element();
    element->elemtype_ = ElemType::TETS;
    element->elemorder_support_ = 5;
    element->dimension_ = 3;
    element->volume_ = 4.0 / 3.0;
    std::string filename = filepath;
    avocado::ReplaceString(filename, "&(Elem)", "TETS");
    element->ReadFromFile(filename);

    element->num_faces_ = 4;
    element->face_nodes_ = {{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}};
    element->face_elemtype_.resize(4, ElemType::TRIS);
    element->ComputeFaceData();
    registry_[element->elemtype_] = element;
  }
  {
    element = new Element();
    element->elemtype_ = ElemType::HEXA;
    element->elemorder_support_ = 5;
    element->dimension_ = 3;
    element->volume_ = 8.0;
    std::string filename = filepath;
    avocado::ReplaceString(filename, "&(Elem)", "HEXA");
    element->ReadFromFile(filename);

    element->num_faces_ = 6;
    element->face_nodes_ = {{0, 1, 2, 3}, {0, 1, 5, 4}, {1, 2, 6, 5},
                            {2, 3, 7, 6}, {3, 0, 4, 7}, {4, 5, 6, 7}};
    element->face_elemtype_.resize(6, ElemType::QUAD);
    element->ComputeFaceData();
    registry_[element->elemtype_] = element;
  }
  {
    element = new Element();
    element->elemtype_ = ElemType::PRIS;
    element->elemorder_support_ = 5;
    element->dimension_ = 3;
    element->volume_ = 4.0;
    std::string filename = filepath;
    avocado::ReplaceString(filename, "&(Elem)", "PRIS");
    element->ReadFromFile(filename);

    element->num_faces_ = 5;
    element->face_nodes_ = {
        {0, 1, 2}, {0, 1, 4, 3}, {0, 2, 5, 3}, {1, 2, 5, 4}, {3, 4, 5}};
    element->face_elemtype_.resize(5, ElemType::QUAD);
    element->face_elemtype_[0] = ElemType::TRIS;
    element->face_elemtype_[4] = ElemType::TRIS;
    element->ComputeFaceData();
    registry_[element->elemtype_] = element;
  }
  {
    element = new Element();
    element->elemtype_ = ElemType::PYRA;
    element->elemorder_support_ = 4;
    element->dimension_ = 3;
    element->volume_ = 8.0 / 3.0;
    std::string filename = filepath;
    avocado::ReplaceString(filename, "&(Elem)", "PYRA");
    element->ReadFromFile(filename);

    element->num_faces_ = 5;
    element->face_nodes_ = {
        {0, 1, 2, 3}, {0, 1, 4}, {0, 3, 4}, {1, 2, 4}, {2, 3, 4}};
    element->face_elemtype_.resize(5, ElemType::TRIS);
    element->face_elemtype_[0] = ElemType::QUAD;
    element->ComputeFaceData();
    registry_[element->elemtype_] = element;
  }
}

const Element* Element::GetElement(const ElemType elemtype) const {
  auto iterator = registry_.find(elemtype);
  if (iterator == registry_.end())
    ERROR_MESSAGE("Illegal element type (no-registered): " +
                  GetElementName(elemtype) + "\n");
  return iterator->second;
}

void Element::GetSubcell(const ElemType elemtype, const int order,
                         std::vector<double>& coords,
                         std::vector<int>& connectivity) const {
  switch (elemtype) {
    case ElemType::LINE:
      GetLineSubcell(order, coords, connectivity);
      break;
    case ElemType::TRIS:
      GetTrisSubcell(order, coords, connectivity);
      break;
    case ElemType::QUAD:
      GetQuadSubcell(order, coords, connectivity);
      break;
    case ElemType::TETS:
      GetTetsSubcell(order, coords, connectivity);
      break;
    case ElemType::HEXA:
      GetHexaSubcell(order, coords, connectivity);
      break;
    case ElemType::PRIS:
      GetPrisSubcell(order, coords, connectivity);
      break;
    case ElemType::PYRA:
      GetPyraSubcell(order, coords, connectivity);
      break;
    default:
      ERROR_MESSAGE("Illegal element type: " +
                    std::to_string(static_cast<int>(elemtype)) + "\n");
  }
}

void Element::TransformToPhyCoords(const int facetype, const int num_nodes,
                                   const double* ref_coords,
                                   double* phy_coords) const {
  const ElemType face_elemtype = GetFaceElemtype(facetype);
  const size_t n = GetFacetypeNodes(1)[facetype].size();
  std::vector<double> face_coord_P1(n * dimension_);
  avocado::VecCopy(static_cast<int>(n), &GetFacetypeNodes(1)[facetype][0],
                   &GetNodesCoord(1)[0], dimension_, &face_coord_P1[0],
                   dimension_);
  TransformToPhyCoords(face_elemtype, num_nodes, &ref_coords[0], &phy_coords[0],
                       &face_coord_P1[0]);
}

void Element::ReadFromFile(const std::string& filename_template) {
  std::vector<int> num_nodes(elemorder_support_);
  std::vector<std::vector<double>> node_coords(elemorder_support_);

  std::ifstream infile;
  for (int iorder = 0; iorder < elemorder_support_; iorder++) {
    std::string filename = filename_template;
    avocado::ReplaceString(filename, "&(Order)", std::to_string(iorder + 1));
    infile.open(filename.c_str());
    if (!infile.is_open())
      ERROR_MESSAGE("Illegal file path (no-exist): " + filename + "\n");

    std::string text;
    std::vector<double> node_coords_temp;
    const int& npt = num_nodes[iorder];
    while (getline(infile, text)) {
      if (text.find("$Nodes", 0) != std::string::npos) {
        infile >> num_nodes[iorder];
        double trash;
        node_coords_temp.resize(npt * dimension_);
        for (int inode = 0; inode < npt; inode++) {
          infile >> trash;
          for (int idim = 0; idim < dimension_; idim++)
            infile >> node_coords_temp[inode * dimension_ + idim];
          for (int idim = dimension_; idim < 3; idim++) infile >> trash;
        }
      }
      if (text.find("$Elements", 0) != std::string::npos) {
        int elemtype, numbering;
        infile >> elemtype;
        node_coords[iorder].resize(npt * dimension_);
        for (int inode = 0; inode < npt; inode++) {
          infile >> numbering;
          for (int idim = 0; idim < dimension_; idim++)
            node_coords[iorder][inode * dimension_ + idim] =
                node_coords_temp[(numbering - 1) * dimension_ + idim];
        }
      }
    }
    infile.close();
  }

  num_nodes_ = std::move(num_nodes);
  node_coords_ = std::move(node_coords);
}

void Element::ComputeFaceData(void) {
  std::vector<double> face_areas(num_faces_);
  std::vector<double> normals(num_faces_ * dimension_);

  std::vector<double> center(dimension_);
  for (int idim = 0; idim < dimension_; idim++)
    center[idim] =
        avocado::VecAverage(num_nodes_[0], &node_coords_[0][idim], dimension_);

  for (int iface = 0; iface < num_faces_; iface++) {
    const size_t num_nodes = face_nodes_[iface].size();
    std::vector<double> face_coords(num_nodes * dimension_);
    avocado::VecCopy(static_cast<int>(num_nodes), &face_nodes_[iface][0],
                     &node_coords_[0][0], dimension_, &face_coords[0],
                     dimension_);
    ComputeNormal(dimension_, face_coords, center,
                  &normals[iface * dimension_]);
    face_areas[iface] = ComputeFaceArea(dimension_, face_coords);
  }

  num_facetypes_ = 0;
  std::vector<int> facetype_to_face;
  std::vector<std::vector<std::vector<int>>> facetype_nodes(elemorder_support_);
  for (int iface = 0; iface < num_faces_; iface++) {
    std::vector<std::vector<int>> nodes_per_facetype =
        GenerateFacetype(face_nodes_[iface]);
    num_facetypes_ += static_cast<int>(nodes_per_facetype.size());
    for (int i = 0; i < nodes_per_facetype.size(); i++) {
      facetype_to_face.push_back(iface);
      facetype_nodes[0].push_back(nodes_per_facetype[i]);
    }
  }

  for (int iorder = 1; iorder < elemorder_support_; iorder++) {
    facetype_nodes[iorder].resize(num_facetypes_);

    for (int itype = 0; itype < num_facetypes_; itype++) {
      const ElemType& face_elemtype = face_elemtype_[facetype_to_face[itype]];
      const Element* face_element = DENEB_ELEMENT->GetElement(face_elemtype);
      const int& num_nodes = face_element->GetNumNodes(iorder + 1);
      const std::vector<double>& nodes_coord =
          face_element->GetNodesCoord(iorder + 1);

      const size_t n = facetype_nodes[0][itype].size();
      std::vector<double> face_coord_P1(n * dimension_);
      avocado::VecCopy(static_cast<int>(n), &facetype_nodes[0][itype][0],
                       &node_coords_[0][0], dimension_, &face_coord_P1[0],
                       dimension_);

      std::vector<double> face_coords(num_nodes * dimension_);
      TransformToPhyCoords(face_elemtype, num_nodes, &nodes_coord[0],
                           &face_coords[0], &face_coord_P1[0]);

      std::vector<int> node_index(num_nodes);
      for (int inode = 0; inode < num_nodes; inode++) {
        for (int jnode = 0; jnode < num_nodes_[iorder]; jnode++) {
          if (avocado::VecIsSame(dimension_, &face_coords[inode * dimension_],
                                 &node_coords_[iorder][jnode * dimension_])) {
            node_index[inode] = jnode;
            break;
          }
        }
      }
      facetype_nodes[iorder][itype] = node_index;
    }
  }

  face_areas_ = std::move(face_areas);
  normals_ = std::move(normals);
  facetype_to_face_ = std::move(facetype_to_face);
  facetype_nodes_ = std::move(facetype_nodes);
}

double Element::ComputeFaceArea(const int dimension,
                                const std::vector<double>& face_coords) const {
  double area = 0.0;
  if (dimension == 1) {
    area = 1.0;
  } else if (dimension == 2) {
    if (face_coords.size() != 2 * dimension)
      ERROR_MESSAGE("Wrong face coordinates\n");
    const std::vector<double> vec = {face_coords[2] - face_coords[0],
                                     face_coords[3] - face_coords[1]};
    area = avocado::VecLength(2, &vec[0]);
  } else if (dimension == 3) {
    if (face_coords.size() < 3 * dimension)
      ERROR_MESSAGE("Wrong face coordinates\n");
    std::vector<double> vec(dimension, 0.0);
    std::vector<double> vec1(dimension, 0.0);
    std::vector<double> vec2(dimension, 0.0);
    for (size_t i = 1, len = face_coords.size() / dimension - 1; i < len; i++) {
      for (int idim = 0; idim < dimension; idim++) {
        vec1[idim] = face_coords[i * dimension + idim] - face_coords[idim];
        vec2[idim] =
            face_coords[(i + 1) * dimension + idim] - face_coords[idim];
      }
      avocado::VecCrossProd(&vec[0], &vec1[0], &vec2[0]);
      area += (0.5 * avocado::VecLength(dimension, &vec[0]));
    }
  } else
    ERROR_MESSAGE("Illegal dimension: " + std::to_string(dimension) + "\n");
  return area;
}

void Element::ComputeNormal(const int dimension,
                            const std::vector<double>& face_coords,
                            const std::vector<double>& inside_coord,
                            double* normal) const {
  std::vector<double> direction(dimension, 0.0);
  for (int idim = 0; idim < dimension; idim++)
    direction[idim] = face_coords[idim] - inside_coord[idim];
  if (dimension == 1) {
    normal[0] = 1.0;
  } else if (dimension == 2) {
    normal[0] = face_coords[3] - face_coords[1];
    normal[1] = face_coords[0] - face_coords[2];
  } else if (dimension == 3) {
    std::vector<double> vec1(dimension, 0.0);
    std::vector<double> vec2(dimension, 0.0);
    for (int idim = 0; idim < dimension; idim++) {
      vec1[idim] = face_coords[2 * dimension + idim] - face_coords[idim];
      vec2[idim] = face_coords[1 * dimension + idim] - face_coords[idim];
    }
    avocado::VecCrossProd(normal, &vec1[0], &vec2[0]);
  } else
    ERROR_MESSAGE("Illegal dimension: " + std::to_string(dimension) + "\n");
  avocado::VecNormal(dimension, normal);
  if (avocado::VecInnerProd(dimension, &direction[0], normal) < 0.0)
    for (int idim = 0; idim < dimension; idim++) normal[idim] = -normal[idim];
}

const std::vector<std::vector<int>> Element::GenerateFacetype(
    const std::vector<int>& face_nodes) const {
  std::vector<std::vector<int>> facetype;
  std::vector<int> facetype_temp = face_nodes;
  // forward direction
  for (size_t i = 0, len = face_nodes.size(); i < len; i++) {
    facetype.push_back(facetype_temp);
    int temp = facetype_temp[0];
    facetype_temp.erase(facetype_temp.begin());
    facetype_temp.push_back(temp);
  }

  if (face_nodes.size() == 2) return facetype;
  std::reverse(facetype_temp.begin(), facetype_temp.end());

  // reverse direction
  for (size_t i = 0, len = face_nodes.size(); i < len; i++) {
    facetype.push_back(facetype_temp);
    int temp = facetype_temp[0];
    facetype_temp.erase(facetype_temp.begin());
    facetype_temp.push_back(temp);
  }
  return facetype;
}
void Element::TransformToPhyCoords(const ElemType elemtype,
                                   const int num_coords,
                                   const double* ref_coords, double* phy_coords,
                                   const double* phyelem_coords) const {
  int num_nodes, ref_dim, phy_dim;
  std::vector<double> rs;
  if (elemtype == ElemType::LINE) {
    num_nodes = 2;
    ref_dim = 1;
    phy_dim = 2;
    rs = {0.5, -0.5, 0.5, 0.5};
  } else if (elemtype == ElemType::TRIS) {
    num_nodes = 3;
    ref_dim = 2;
    phy_dim = 3;
    rs = {0.0, -0.5, -0.5, 0.5, 0.5, 0.0, 0.5, 0.0, 0.5};
  } else if (elemtype == ElemType::QUAD) {
    num_nodes = 4;
    ref_dim = 2;
    phy_dim = 3;
    rs = {0.25, -0.25, -0.25, 0.25, 0.25, 0.25,  -0.25, -0.25,
          0.25, 0.25,  0.25,  0.25, 0.25, -0.25, 0.25,  -0.25};
  } else
    ERROR_MESSAGE("Not supported element type: " + GetElementName(elemtype) +
                  "\n");
  std::vector<double> trans(phy_dim * num_nodes);
  gemmATB(1.0, &phyelem_coords[0], &rs[0], 0.0, &trans[0], num_nodes, phy_dim,
          num_nodes);

  int ind = 0;
  std::vector<double> pt_rs(num_coords * num_nodes);
  if (elemtype == ElemType::LINE) {
    for (int i = 0; i < num_coords; i++) {
      const double& r = ref_coords[i];
      pt_rs[ind++] = 1.0;
      pt_rs[ind++] = r;
    }
  } else if (elemtype == ElemType::TRIS) {
    for (int i = 0; i < num_coords; i++) {
      const double& r = ref_coords[i * ref_dim];
      const double& s = ref_coords[i * ref_dim + 1];
      pt_rs[ind++] = 1.0;
      pt_rs[ind++] = r;
      pt_rs[ind++] = s;
    }
  } else if (elemtype == ElemType::QUAD) {
    for (int i = 0; i < num_coords; i++) {
      const double& r = ref_coords[i * ref_dim];
      const double& s = ref_coords[i * ref_dim + 1];
      pt_rs[ind++] = 1.0;
      pt_rs[ind++] = r;
      pt_rs[ind++] = s;
      pt_rs[ind++] = r * s;
    }
  }
  gemmABT(1.0, &pt_rs[0], &trans[0], 0.0, &phy_coords[0], num_coords, num_nodes,
          phy_dim);
}

const std::vector<int> Element::GetType1Connects(
    const std::vector<int>& pts) const {
  //                   z
  //	                .
  //         	      ,/
  //	             /
  //	           3
  //	         ,/|`\
	//	       ,/  |  `\
	//	     ,/    '.   `\    
	//	   ,/       |     `\
	//   ,/         |       `\
	//	0-----------'.--------2 --> y
  //	 `\.        |        ,/
  //	    `\.     |      ,/
  //	       `\.   '.  ,/
  //	          `\. | /
  //	             `1
  //	               `\.
  //	                  `x
  return std::vector<int>({pts[0], pts[1], pts[2], pts[3]});
}
const std::vector<int> Element::GetType2Connects(
    const std::vector<int>& pts) const {
  //         z
  //	4----------6
  //	|\     ^   |\ 
	//	| \    |   | \
	//	|  \   |   |  \
	//	|   5------+---none
  //	|   |  +-- |-- |->y
  //	0---+---\--3   |
  //	 \  |    \  \  |
  //	  \ |     \  \ |
  //	   \|      x  \|
  //	    1----------2
  return std::vector<int>({pts[0], pts[1], pts[3], pts[4], pts[1],
                           pts[3], pts[5], pts[4], pts[3], pts[4],
                           pts[5], pts[6], pts[1], pts[2], pts[3],
                           pts[5], pts[2], pts[3], pts[5], pts[6]});
}
const std::vector<int> Element::GetType3Connects(
    const std::vector<int>& pts) const {
  //         z
  //	4----------7
  //	|\     ^   |\ 
	//	| \    |   | \
	//	|  \   |   |  \
	//	|   5------+---6
  //	|   |  +-- |-- |->y
  //	0---+---\--3   |
  //	 \  |    \  \  |
  //	  \ |     \  \ |
  //	   \|      x  \|
  //	    1----------2
  return std::vector<int>({pts[0], pts[1], pts[3], pts[4], pts[1], pts[3],
                           pts[5], pts[4], pts[3], pts[4], pts[5], pts[7],
                           pts[1], pts[2], pts[3], pts[5], pts[2], pts[3],
                           pts[5], pts[7], pts[2], pts[5], pts[6], pts[7]});
}
const std::vector<int> Element::GetType4Connects(
    const std::vector<int>& pts) const {
  //         z
  //	3----------5
  //	|\     ^   |\ 
	//	| \    |   | \
	//	|  \   |   |  \
	//	|   4------+---none
  //	|   |  +-- |-- |->y
  //	0---+---\--2   |
  //	 \  |    \  \  |
  //	  \ |     \  \ |
  //	   \|      x  \|
  //	    1----------none
  return std::vector<int>({pts[0], pts[1], pts[2], pts[3], pts[1], pts[2],
                           pts[4], pts[3], pts[2], pts[3], pts[4], pts[5]});
}

void Element::GetLineSubcell(const int order, std::vector<double>& coords,
                             std::vector<int>& connectivity) const {
  std::vector<double> coords_local;
  std::vector<int> connectivity_local;

  coords_local.reserve(order + 1);
  for (int ix = 0; ix <= order; ix++)
    coords_local.push_back(double(2 * ix) / double(order) - 1.0);

  connectivity_local.reserve(2 * order);
  for (int ix = 0; ix < order; ix++) {
    connectivity_local.push_back(ix);
    connectivity_local.push_back(ix + 1);
  }

  coords = std::move(coords_local);
  connectivity = std::move(connectivity_local);
}
void Element::GetTrisSubcell(const int order, std::vector<double>& coords,
                             std::vector<int>& connectivity) const {
  std::vector<double> coords_local;
  std::vector<int> connectivity_local;

  coords_local.reserve((order + 1) * (order + 2));
  for (int iy = 0; iy <= order; iy++) {
    for (int ix = 0; ix <= order - iy; ix++) {
      coords_local.push_back(double(2 * ix) / double(order) - 1.0);
      coords_local.push_back(double(2 * iy) / double(order) - 1.0);
    }
  }

  connectivity_local.reserve(3 * order * order);
  for (int iy = 0; iy < order; iy++) {
    for (int ix = 0; ix < order - iy; ix++) {
      const int var1 = (order + 2) * iy - iy * (iy + 1) / 2 + ix;
      const int var2 = (order + 2) * (iy + 1) - (iy + 1) * (iy + 2) / 2 + ix;
      connectivity_local.push_back(var1);
      connectivity_local.push_back(var1 + 1);
      connectivity_local.push_back(var2);

      if (ix == order - iy - 1) continue;

      connectivity_local.push_back(var1 + 1);
      connectivity_local.push_back(var2 + 1);
      connectivity_local.push_back(var2);
    }
  }

  coords = std::move(coords_local);
  connectivity = std::move(connectivity_local);
}
void Element::GetQuadSubcell(const int order, std::vector<double>& coords,
                             std::vector<int>& connectivity) const {
  std::vector<double> coords_local;
  std::vector<int> connectivity_local;

  coords_local.reserve(2 * (order + 1) * (order + 1));
  for (int iy = 0; iy <= order; iy++) {
    for (int ix = 0; ix <= order; ix++) {
      coords_local.push_back(double(2 * ix) / double(order) - 1.0);
      coords_local.push_back(double(2 * iy) / double(order) - 1.0);
    }
  }

  connectivity_local.reserve(6 * order * order);
  for (int iy = 0; iy < order; iy++) {
    for (int ix = 0; ix < order; ix++) {
      const int var1 = (order + 1) * iy + ix;
      const int var2 = (order + 1) * (iy + 1) + ix;
      connectivity_local.push_back(var1);
      connectivity_local.push_back(var1 + 1);
      connectivity_local.push_back(var2);

      connectivity_local.push_back(var1 + 1);
      connectivity_local.push_back(var2 + 1);
      connectivity_local.push_back(var2);
    }
  }

  coords = std::move(coords_local);
  connectivity = std::move(connectivity_local);
}
void Element::GetTetsSubcell(const int order, std::vector<double>& coords,
                             std::vector<int>& connectivity) const {
  std::vector<double> coords_local;
  std::vector<int> connectivity_local;

  coords_local.reserve((order + 1) * (order + 2) * (order + 3) / 2);
  for (int iz = 0; iz <= order; iz++) {
    for (int iy = 0; iy <= order - iz; iy++) {
      for (int ix = 0; ix <= order - iz - iy; ix++) {
        coords_local.push_back(double(2 * ix) / double(order) - 1.0);
        coords_local.push_back(double(2 * iy) / double(order) - 1.0);
        coords_local.push_back(double(2 * iz) / double(order) - 1.0);
      }
    }
  }

  connectivity_local.reserve(4 * order * order * order);
  for (int iz = 0; iz < order; iz++) {
    for (int iy = 0; iy < order - iz; iy++) {
      for (int ix = 0; ix < order - iz - iy; ix++) {
        const int sum = ix + iy + iz;
        if (sum == order - 1) {
          const int var1 = (order - iz + 2) * iy - iy * (iy + 1) / 2 + ix +
                           iz *
                               (3 * order * order - 3 * order * iz +
                                12 * order + iz * iz - 6 * iz + 11) /
                               6;
          const int var2 = (order - iz + 2) * (iy + 1) -
                           (iy + 1) * (iy + 2) / 2 + ix +
                           iz *
                               (3 * order * order - 3 * order * iz +
                                12 * order + iz * iz - 6 * iz + 11) /
                               6;
          const int var3 = (order - iz + 1) * iy - iy * (iy + 1) / 2 + ix +
                           (iz + 1) *
                               (3 * order * order - 3 * order * iz + 9 * order +
                                iz * iz - 4 * iz + 6) /
                               6;
          const std::vector<int> type1 = {var1, var1 + 1, var2, var3};
          for (auto&& i : GetType1Connects(type1))
            connectivity_local.push_back(i);

        } else if (sum == order - 2) {
          const int var1 = (order - iz + 2) * iy - iy * (iy + 1) / 2 + ix +
                           iz *
                               (3 * order * order - 3 * order * iz +
                                12 * order + iz * iz - 6 * iz + 11) /
                               6;
          const int var2 = (order - iz + 2) * (iy + 1) -
                           (iy + 1) * (iy + 2) / 2 + ix +
                           iz *
                               (3 * order * order - 3 * order * iz +
                                12 * order + iz * iz - 6 * iz + 11) /
                               6;
          const int var3 = (order - iz + 1) * iy - iy * (iy + 1) / 2 + ix +
                           (iz + 1) *
                               (3 * order * order - 3 * order * iz + 9 * order +
                                iz * iz - 4 * iz + 6) /
                               6;
          const int var4 = (order - iz + 1) * (iy + 1) -
                           (iy + 1) * (iy + 2) / 2 + ix +
                           (iz + 1) *
                               (3 * order * order - 3 * order * iz + 9 * order +
                                iz * iz - 4 * iz + 6) /
                               6;
          const std::vector<int> type2 = {var1, var1 + 1, var2 + 1, var2,
                                          var3, var3 + 1, var4};
          for (auto&& i : GetType2Connects(type2))
            connectivity_local.push_back(i);
        } else {
          const int var1 = (order - iz + 2) * iy - iy * (iy + 1) / 2 + ix +
                           iz *
                               (3 * order * order - 3 * order * iz +
                                12 * order + iz * iz - 6 * iz + 11) /
                               6;
          const int var2 = (order - iz + 2) * (iy + 1) -
                           (iy + 1) * (iy + 2) / 2 + ix +
                           iz *
                               (3 * order * order - 3 * order * iz +
                                12 * order + iz * iz - 6 * iz + 11) /
                               6;
          const int var3 = (order - iz + 1) * iy - iy * (iy + 1) / 2 + ix +
                           (iz + 1) *
                               (3 * order * order - 3 * order * iz + 9 * order +
                                iz * iz - 4 * iz + 6) /
                               6;
          const int var4 = (order - iz + 1) * (iy + 1) -
                           (iy + 1) * (iy + 2) / 2 + ix +
                           (iz + 1) *
                               (3 * order * order - 3 * order * iz + 9 * order +
                                iz * iz - 4 * iz + 6) /
                               6;
          const std::vector<int> type3 = {var1, var1 + 1, var2 + 1, var2,
                                          var3, var3 + 1, var4 + 1, var4};
          for (auto&& i : GetType3Connects(type3))
            connectivity_local.push_back(i);
        }
      }
    }
  }

  coords = std::move(coords_local);
  connectivity = std::move(connectivity_local);
}
void Element::GetHexaSubcell(const int order, std::vector<double>& coords,
                             std::vector<int>& connectivity) const {
  std::vector<double> coords_local;
  std::vector<int> connectivity_local;

  coords_local.reserve((order + 1) * (order + 1) * (order + 1) * 3);
  for (int iz = 0; iz <= order; iz++) {
    for (int iy = 0; iy <= order; iy++) {
      for (int ix = 0; ix <= order; ix++) {
        coords_local.push_back(double(2 * ix) / double(order) - 1.0);
        coords_local.push_back(double(2 * iy) / double(order) - 1.0);
        coords_local.push_back(double(2 * iz) / double(order) - 1.0);
      }
    }
  }

  connectivity_local.reserve(24 * order * order * order);
  for (int iz = 0; iz < order; iz++) {
    for (int iy = 0; iy < order; iy++) {
      for (int ix = 0; ix < order; ix++) {
        const int var1 = (order + 1) * (order + 1) * iz + (order + 1) * iy + ix;
        const int var2 =
            (order + 1) * (order + 1) * iz + (order + 1) * (iy + 1) + ix;
        const int var3 =
            (order + 1) * (order + 1) * (iz + 1) + (order + 1) * iy + ix;
        const int var4 =
            (order + 1) * (order + 1) * (iz + 1) + (order + 1) * (iy + 1) + ix;
        const std::vector<int> type1 = {var1, var1 + 1, var2 + 1, var2,
                                        var3, var3 + 1, var4 + 1, var4};
        for (auto&& i : GetType3Connects(type1))
          connectivity_local.push_back(i);
      }
    }
  }

  coords = std::move(coords_local);
  connectivity = std::move(connectivity_local);
}
void Element::GetPrisSubcell(const int order, std::vector<double>& coords,
                             std::vector<int>& connectivity) const {
  std::vector<double> coords_local;
  std::vector<int> connectivity_local;

  coords_local.reserve((order + 1) * (order + 2) * (order + 1) / 2 * 3);
  for (int iz = 0; iz <= order; iz++) {
    for (int iy = 0; iy <= order; iy++) {
      for (int ix = 0; ix <= order - iy; ix++) {
        coords_local.push_back(double(2 * ix) / double(order) - 1.0);
        coords_local.push_back(double(2 * iy) / double(order) - 1.0);
        coords_local.push_back(double(2 * iz) / double(order) - 1.0);
      }
    }
  }

  connectivity_local.reserve(12 * order * order * order);
  for (int iz = 0; iz < order; iz++) {
    for (int iy = 0; iy < order; iy++) {
      for (int ix = 0; ix < order - iy; ix++) {
        const int sum = ix + iy;
        if (sum == order - 1) {
          const int var1 = (order + 1) * (order + 2) / 2 * iz +
                           (2 * order - iy + 3) * iy / 2 + ix;
          const int var2 = (order + 1) * (order + 2) / 2 * iz +
                           (2 * order - iy + 2) * (iy + 1) / 2 + ix;
          const int var3 = (order + 1) * (order + 2) / 2 * (iz + 1) +
                           (2 * order - iy + 3) * iy / 2 + ix;
          const int var4 = (order + 1) * (order + 2) / 2 * (iz + 1) +
                           (2 * order - iy + 2) * (iy + 1) / 2 + ix;
          const std::vector<int> type2 = {var1, var1 + 1, var2,
                                          var3, var3 + 1, var4};
          for (auto&& i : GetType4Connects(type2))
            connectivity_local.push_back(i);
        } else {
          const int var1 = (order + 1) * (order + 2) / 2 * iz +
                           (2 * order - iy + 3) * iy / 2 + ix;
          const int var2 = (order + 1) * (order + 2) / 2 * iz +
                           (2 * order - iy + 2) * (iy + 1) / 2 + ix;
          const int var3 = (order + 1) * (order + 2) / 2 * (iz + 1) +
                           (2 * order - iy + 3) * iy / 2 + ix;
          const int var4 = (order + 1) * (order + 2) / 2 * (iz + 1) +
                           (2 * order - iy + 2) * (iy + 1) / 2 + ix;
          const std::vector<int> type1 = {var1, var1 + 1, var2 + 1, var2,
                                          var3, var3 + 1, var4 + 1, var4};
          for (auto&& i : GetType3Connects(type1))
            connectivity_local.push_back(i);
        }
      }
    }
  }

  coords = std::move(coords_local);
  connectivity = std::move(connectivity_local);
}
void Element::GetPyraSubcell(const int order, std::vector<double>& coords,
                             std::vector<int>& connectivity) const {
  std::vector<double> coords_local;
  std::vector<int> connectivity_local;

  coords_local.reserve((order + 1) * (order + 2) * (2 * order + 3) / 2);
  for (int iz = 0; iz <= order; iz++) {
    for (int iy = 0; iy <= order - iz; iy++) {
      for (int ix = 0; ix <= order - iz - iy; ix++) {
        const double r = double(2 * ix) / double(order) - 1.0;
        const double s = double(2 * iy) / double(order) - 1.0;
        const double t = double(2 * iz) / double(order) - 1.0;
        coords_local.push_back(0.5 + r + 0.5 * t);
        coords_local.push_back(0.5 + s + 0.5 * t);
        coords_local.push_back(t);
      }
    }
  }
  for (int iz = 0; iz < order; iz++) {
    for (int iy = 0; iy < order - iz; iy++) {
      for (int ix = 0; ix < order - iz - iy; ix++) {
        const double r = double(2 * ix) / double(order) - 1.0;
        const double s = double(2 * iy) / double(order) - 1.0;
        const double t = double(2 * iz) / double(order) - 1.0;
        coords_local.push_back(-0.5 - s - 0.5 * t);
        coords_local.push_back(-0.5 - r - 0.5 * t);
        coords_local.push_back(t);
      }
    }
  }

  connectivity_local.reserve(8 * order * order * order);
  for (int iz = 0; iz < order; iz++) {
    for (int iy = 0; iy < order - iz; iy++) {
      for (int ix = 0; ix < order - iz - iy; ix++) {
        const int sum = ix + iy + iz;
        if (sum == order - 1) {
          const int var1 = (order - iz + 2) * iy - iy * (iy + 1) / 2 + ix +
                           iz *
                               (3 * order * order - 3 * order * iz +
                                12 * order + iz * iz - 6 * iz + 11) /
                               6;
          const int var2 = (order - iz + 2) * (iy + 1) -
                           (iy + 1) * (iy + 2) / 2 + ix +
                           iz *
                               (3 * order * order - 3 * order * iz +
                                12 * order + iz * iz - 6 * iz + 11) /
                               6;
          const int var3 = (order - iz + 1) * iy - iy * (iy + 1) / 2 + ix +
                           (iz + 1) *
                               (3 * order * order - 3 * order * iz + 9 * order +
                                iz * iz - 4 * iz + 6) /
                               6;
          const std::vector<int> type1 = {var1, var1 + 1, var2, var3};
          for (auto&& i : GetType1Connects(type1))
            connectivity_local.push_back(i);

        } else if (sum == order - 2) {
          const int var1 = (order - iz + 2) * iy - iy * (iy + 1) / 2 + ix +
                           iz *
                               (3 * order * order - 3 * order * iz +
                                12 * order + iz * iz - 6 * iz + 11) /
                               6;
          const int var2 = (order - iz + 2) * (iy + 1) -
                           (iy + 1) * (iy + 2) / 2 + ix +
                           iz *
                               (3 * order * order - 3 * order * iz +
                                12 * order + iz * iz - 6 * iz + 11) /
                               6;
          const int var3 = (order - iz + 1) * iy - iy * (iy + 1) / 2 + ix +
                           (iz + 1) *
                               (3 * order * order - 3 * order * iz + 9 * order +
                                iz * iz - 4 * iz + 6) /
                               6;
          const int var4 = (order - iz + 1) * (iy + 1) -
                           (iy + 1) * (iy + 2) / 2 + ix +
                           (iz + 1) *
                               (3 * order * order - 3 * order * iz + 9 * order +
                                iz * iz - 4 * iz + 6) /
                               6;
          const std::vector<int> type2 = {var1, var1 + 1, var2 + 1, var2,
                                          var3, var3 + 1, var4};
          for (auto&& i : GetType2Connects(type2))
            connectivity_local.push_back(i);
        } else {
          const int var1 = (order - iz + 2) * iy - iy * (iy + 1) / 2 + ix +
                           iz *
                               (3 * order * order - 3 * order * iz +
                                12 * order + iz * iz - 6 * iz + 11) /
                               6;
          const int var2 = (order - iz + 2) * (iy + 1) -
                           (iy + 1) * (iy + 2) / 2 + ix +
                           iz *
                               (3 * order * order - 3 * order * iz +
                                12 * order + iz * iz - 6 * iz + 11) /
                               6;
          const int var3 = (order - iz + 1) * iy - iy * (iy + 1) / 2 + ix +
                           (iz + 1) *
                               (3 * order * order - 3 * order * iz + 9 * order +
                                iz * iz - 4 * iz + 6) /
                               6;
          const int var4 = (order - iz + 1) * (iy + 1) -
                           (iy + 1) * (iy + 2) / 2 + ix +
                           (iz + 1) *
                               (3 * order * order - 3 * order * iz + 9 * order +
                                iz * iz - 4 * iz + 6) /
                               6;
          const std::vector<int> type3 = {var1, var1 + 1, var2 + 1, var2,
                                          var3, var3 + 1, var4 + 1, var4};
          for (auto&& i : GetType3Connects(type3))
            connectivity_local.push_back(i);
        }
      }
    }
  }
  const int gap = (order + 1) * (order + 2) * (order + 3) / 6;
  const int order_1 = order - 1;
  for (int iz = 0; iz < order; iz++) {
    for (int iy = 0; iy < order - iz; iy++) {
      for (int ix = 0; ix < order - iz - iy; ix++) {
        const int sum = ix + iy + iz;
        if (sum == order - 1) {
          const int var1_1 = gap + (order_1 - iz + 2) * iy - iy * (iy + 1) / 2 +
                             ix +
                             iz *
                                 (3 * order_1 * order_1 - 3 * order_1 * iz +
                                  12 * order_1 + iz * iz - 6 * iz + 11) /
                                 6;
          const int var1 = (order - iz + 2) * iy - iy * (iy + 1) / 2 + ix +
                           iz *
                               (3 * order * order - 3 * order * iz +
                                12 * order + iz * iz - 6 * iz + 11) /
                               6;
          const int var2 = (order - iz + 2) * (iy + 1) -
                           (iy + 1) * (iy + 2) / 2 + ix +
                           iz *
                               (3 * order * order - 3 * order * iz +
                                12 * order + iz * iz - 6 * iz + 11) /
                               6;
          const int var3 = (order - iz + 1) * iy - iy * (iy + 1) / 2 + ix +
                           (iz + 1) *
                               (3 * order * order - 3 * order * iz + 9 * order +
                                iz * iz - 4 * iz + 6) /
                               6;
          const std::vector<int> type1 = {var1_1, var1 + 1, var2, var3};
          for (auto&& i : GetType1Connects(type1))
            connectivity_local.push_back(i);

        } else if (sum == order - 2) {
          const int var1_1 = gap + (order_1 - iz + 2) * iy - iy * (iy + 1) / 2 +
                             ix +
                             iz *
                                 (3 * order_1 * order_1 - 3 * order_1 * iz +
                                  12 * order_1 + iz * iz - 6 * iz + 11) /
                                 6;
          const int var2 = (order - iz + 2) * (iy + 1) -
                           (iy + 1) * (iy + 2) / 2 + ix +
                           iz *
                               (3 * order * order - 3 * order * iz +
                                12 * order + iz * iz - 6 * iz + 11) /
                               6;
          const int var2_1 = gap + (order_1 - iz + 2) * (iy + 1) -
                             (iy + 1) * (iy + 2) / 2 + ix +
                             iz *
                                 (3 * order_1 * order_1 - 3 * order_1 * iz +
                                  12 * order_1 + iz * iz - 6 * iz + 11) /
                                 6;
          const int var3 = (order - iz + 1) * iy - iy * (iy + 1) / 2 + ix +
                           (iz + 1) *
                               (3 * order * order - 3 * order * iz + 9 * order +
                                iz * iz - 4 * iz + 6) /
                               6;
          const int var3_1 = gap + (order_1 - iz + 1) * iy - iy * (iy + 1) / 2 +
                             ix +
                             (iz + 1) *
                                 (3 * order_1 * order_1 - 3 * order_1 * iz +
                                  9 * order_1 + iz * iz - 4 * iz + 6) /
                                 6;
          const int var4 = (order - iz + 1) * (iy + 1) -
                           (iy + 1) * (iy + 2) / 2 + ix +
                           (iz + 1) *
                               (3 * order * order - 3 * order * iz + 9 * order +
                                iz * iz - 4 * iz + 6) /
                               6;
          const std::vector<int> type2 = {var1_1, var1_1 + 1, var2 + 1, var2_1,
                                          var3_1, var3 + 1,   var4};
          for (auto&& i : GetType2Connects(type2))
            connectivity_local.push_back(i);
        } else {
          const int var4 = (order - iz + 1) * (iy + 1) -
                           (iy + 1) * (iy + 2) / 2 + ix +
                           (iz + 1) *
                               (3 * order * order - 3 * order * iz + 9 * order +
                                iz * iz - 4 * iz + 6) /
                               6;
          const int var1_1 = gap + (order_1 - iz + 2) * iy - iy * (iy + 1) / 2 +
                             ix +
                             iz *
                                 (3 * order_1 * order_1 - 3 * order_1 * iz +
                                  12 * order_1 + iz * iz - 6 * iz + 11) /
                                 6;
          const int var2_1 = gap + (order_1 - iz + 2) * (iy + 1) -
                             (iy + 1) * (iy + 2) / 2 + ix +
                             iz *
                                 (3 * order_1 * order_1 - 3 * order_1 * iz +
                                  12 * order_1 + iz * iz - 6 * iz + 11) /
                                 6;
          const int var3_1 = gap + (order_1 - iz + 1) * iy - iy * (iy + 1) / 2 +
                             ix +
                             (iz + 1) *
                                 (3 * order_1 * order_1 - 3 * order_1 * iz +
                                  9 * order_1 + iz * iz - 4 * iz + 6) /
                                 6;
          const int var4_1 = gap + (order_1 - iz + 1) * (iy + 1) -
                             (iy + 1) * (iy + 2) / 2 + ix +
                             (iz + 1) *
                                 (3 * order_1 * order_1 - 3 * order_1 * iz +
                                  9 * order_1 + iz * iz - 4 * iz + 6) /
                                 6;
          int var5 = var4_1 + 1;
          if (sum == order - 3) var5 = var4 + 1;
          const std::vector<int> type3 = {var1_1, var1_1 + 1, var2_1 + 1,
                                          var2_1, var3_1,     var3_1 + 1,
                                          var5,   var4_1};
          for (auto&& i : GetType3Connects(type3))
            connectivity_local.push_back(i);
        }
      }
    }
  }

  coords = std::move(coords_local);
  connectivity = std::move(connectivity_local);
}

size_t Hash_ElemType(const ElemType& elemtype) {
  return std::hash<int>()(static_cast<int>(elemtype));
}
}  // namespace deneb