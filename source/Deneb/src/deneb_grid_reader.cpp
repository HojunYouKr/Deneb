#include "deneb_grid_reader.h"

#include <algorithm>
#include <fstream>
#include <iomanip>

#include "avocado.h"
#include "deneb_config_macro.h"

namespace deneb {
std::shared_ptr<GridReader> GridReader::GetGridReader(const std::string& name) {
  if (!name.compare("Gmsh"))
    return std::make_shared<GmshGridReader>();
  else if (!name.compare("Gmsh2DTris"))
    return std::make_shared<Gmsh2DTrisGridReader>();
  else if (!name.compare("Gmsh2DQuad"))
    return std::make_shared<Gmsh2DQuadGridReader>();
  else if (!name.compare("Gmsh3DHexa"))
    return std::make_shared<Gmsh3DHexaGridReader>();
  ERROR_MESSAGE("Invalid grid-type (no-exist):" + name + "\n");
  return nullptr;
}

void GmshGridReader::Open(const int dimension,
                          const std::string& grid_file_path) {
  dimension_ = dimension;
  grid_file_path_ = grid_file_path;
  num_global_nodes_ = 0;
  num_global_cells_ = 0;
  num_global_bdries_ = 0;
  num_global_peribdries_ = 0;

  GenerateGmshFile();
  SYNCRO();
  std::ifstream infile(grid_file_path_.c_str());
  if (infile.is_open())
    MASTER_MESSAGE("Open Gmsh-type grid file: " + grid_file_path_ + "\n");
  else
    ERROR_MESSAGE("Invalid Gmsh-type grid file (no-exist): " + grid_file_path_ +
                  "\n");

  std::unordered_set<int> periBC_tag;
  auto& config = AVOCADO_CONFIG;
  for (int idim = 0; idim < 3; idim++) {
    periBC_tag.insert(std::stoi(config->GetConfigValue(PERIODIC_LEFT(idim))));
    periBC_tag.insert(std::stoi(config->GetConfigValue(PERIODIC_RIGHT(idim))));
  }
  periBC_tag.insert(std::stoi(config->GetConfigValue(MATCHING_BC)));

  std::string text;
  std::unordered_set<int> physical_bc_tag;
  std::unordered_set<int> periodic_bc_tag;
  std::unordered_map<ElemType, int, ElemTypeHash> cell_elemtypes;
  std::unordered_map<ElemType, int, ElemTypeHash> bdry_elemtypes;
  std::unordered_map<ElemType, int, ElemTypeHash> peribdry_elemtypes;
  std::unordered_map<ElemType, int, ElemTypeHash>* ptr = nullptr;
  while (getline(infile, text)) {
    if (text.find("$Nodes", 0) != std::string::npos) {
      infile >> num_global_nodes_;
      MASTER_MESSAGE("Number of nodes: " + std::to_string(num_global_nodes_) +
                     "\n");
    }
    if (text.find("$Elements", 0) != std::string::npos) {
      int num_total;
      infile >> num_total;
      MASTER_MESSAGE("Number of connectivities: " + std::to_string(num_total) +
                     "\n");

      ElemType elemtype;
      int cell_index, gmshtype, garbage, bc_tag;
      for (int icell = 0; icell < num_total; icell++) {
        infile >> cell_index >> gmshtype >> garbage >> bc_tag;

        elemtype = GetElemTypeFromGmshType(static_cast<GmshElemType>(gmshtype));
        if (Element::IsFace(dimension_, elemtype)) {
          auto&& iterator = periBC_tag.find(bc_tag);
          if (iterator == periBC_tag.end()) {
            num_global_bdries_++;
            physical_bc_tag.insert(bc_tag);
            ptr = &bdry_elemtypes;
          } else {
            num_global_peribdries_++;
            periodic_bc_tag.insert(bc_tag);
            ptr = &peribdry_elemtypes;
          }
        } else {
          num_global_cells_++;
          ptr = &cell_elemtypes;
        }

        auto&& iterator = (*ptr).find(elemtype);
        if (iterator == (*ptr).end())
          (*ptr)[elemtype] = 1;
        else
          iterator->second++;

        getline(infile, text);
      }
    }
  }
  infile.close();

  MASTER_MESSAGE("Number of cells: " + std::to_string(num_global_cells_) +
                 "\n");
  for (auto&& iterator : cell_elemtypes)
    MASTER_MESSAGE("\t" + Element::GetElementName(iterator.first) + ": " +
                   std::to_string(iterator.second) + "\n");

  MASTER_MESSAGE("Number of boundaries: " + std::to_string(num_global_bdries_) +
                 "\n");
  for (auto&& iterator : bdry_elemtypes)
    MASTER_MESSAGE("\t" + Element::GetElementName(iterator.first) + ": " +
                   std::to_string(iterator.second) + "\n");

  MASTER_MESSAGE("Number of periodic boundaries: " +
                 std::to_string(num_global_peribdries_) + "\n");
  for (auto&& iterator : peribdry_elemtypes)
    MASTER_MESSAGE("\t" + Element::GetElementName(iterator.first) + ": " +
                   std::to_string(iterator.second) + "\n");

  {
    std::string message = "Boundary conditions: ";
    for (auto&& bc_tag : physical_bc_tag) {
      const std::string& bc_type = config->GetConfigValue(BDRY_TYPE(bc_tag));
      const std::string& bc_name = config->GetConfigValue(BDRY_NAME(bc_tag));
      message += ("\n\tTag = " + std::to_string(bc_tag) +
                  "\tType = " + bc_type + "\tName = " + bc_name);
    }
    MASTER_MESSAGE(message + "\n");
    for (auto&& tag : physical_bc_tag) all_bdry_tag_.push_back(tag);
  }
}
void GmshGridReader::ReadCellData(const int start_index, const int end_index,
                                  std::vector<ElemType>& cell_elemtype,
                                  std::vector<int>& cell_global_index,
                                  std::vector<int>& cell_node_ptr,
                                  std::vector<int>& cell_node_ind) const {
  cell_elemtype.clear();
  cell_global_index.clear();
  cell_node_ind.clear();
  cell_node_ptr.clear();
  cell_elemtype.reserve(end_index - start_index);
  cell_global_index.reserve(end_index - start_index);
  cell_node_ptr.reserve(end_index - start_index);
  cell_node_ptr.push_back(0);
  cell_node_ind.reserve(dimension_ * (end_index - start_index));

  std::ifstream infile(grid_file_path_.c_str());
  std::string text;
  while (getline(infile, text)) {
    if (text.find("$Elements", 0) != std::string::npos) {
      int num_total;
      infile >> num_total;
      ElemType elemtype;
      int cell_index, gmshtype, garbage, node_index;

      int num_cells = 0;
      for (int icell = 0; icell < num_total; icell++) {
        infile >> cell_index >> gmshtype >> garbage >> garbage >> garbage;
        elemtype = GetElemTypeFromGmshType(static_cast<GmshElemType>(gmshtype));
        if (!Element::IsFace(dimension_, elemtype)) {
          if ((start_index <= num_cells) && (num_cells < end_index)) {
            const int num_nodes =
                DENEB_ELEMENT->GetElement(elemtype)->GetNumNodes(1);
            for (int inode = 0; inode < num_nodes; inode++) {
              infile >> node_index;
              cell_node_ind.push_back(node_index);
            }
            cell_elemtype.push_back(elemtype);
            cell_node_ptr.push_back(static_cast<int>(cell_node_ind.size()));
            cell_global_index.push_back(cell_index);  // note here
          }
          num_cells++;
        }
        if (end_index <= num_cells) break;
        getline(infile, text);
      }
      break;
    }
  }
  infile.close();
}

void GmshGridReader::ReadAllBdryData(
    std::vector<int>& all_bdry_global_index, std::vector<int>& all_bdry_tag,
    std::vector<int>& all_bdry_node_ptr,
    std::vector<int>& all_bdry_node_ind) const {
  const int num_all_bdries = num_global_bdries_ + num_global_peribdries_;
  all_bdry_global_index.clear();
  all_bdry_tag.clear();
  all_bdry_node_ptr.clear();
  all_bdry_node_ind.clear();
  all_bdry_global_index.reserve(num_all_bdries);
  all_bdry_tag.reserve(num_all_bdries);
  all_bdry_node_ptr.reserve(num_all_bdries);
  all_bdry_node_ptr.push_back(0);
  all_bdry_node_ind.reserve(num_all_bdries * dimension_);

  std::ifstream infile(grid_file_path_.c_str());
  std::string text;
  while (getline(infile, text)) {
    if (text.find("$Elements", 0) != std::string::npos) {
      int num_total;
      infile >> num_total;
      ElemType elemtype;
      int bdry_index, gmshtype, tag, garbage, node_index;

      int num_bdries = 0;
      for (int icell = 0; icell < num_total; icell++) {
        infile >> bdry_index >> gmshtype >> garbage >> tag >> garbage;
        elemtype = GetElemTypeFromGmshType(static_cast<GmshElemType>(gmshtype));
        if (Element::IsFace(dimension_, elemtype)) {
          const int num_nodes =
              DENEB_ELEMENT->GetElement(elemtype)->GetNumNodes(1);
          for (int inode = 0; inode < num_nodes; inode++) {
            infile >> node_index;
            all_bdry_node_ind.push_back(node_index);
          }
          all_bdry_node_ptr.push_back(
              static_cast<int>(all_bdry_node_ind.size()));
          all_bdry_global_index.push_back(bdry_index);  // note here
          all_bdry_tag.push_back(tag);
          num_bdries++;
        }
        if (num_bdries == num_all_bdries) break;
        getline(infile, text);
      }
      break;
    }
  }
  infile.close();
}
void GmshGridReader::ReadAllPeriodicData(
    std::vector<int>& all_periodic_nodes,
    std::vector<double>& all_periodic_node_coords) const {
  all_periodic_nodes.clear();
  all_periodic_node_coords.clear();
  if (num_global_peribdries_ == 0) return;

  std::unordered_set<int> periBC_tag;
  auto& config = AVOCADO_CONFIG;
  for (int idim = 0; idim < 3; idim++) {
    periBC_tag.insert(std::stoi(config->GetConfigValue(PERIODIC_LEFT(idim))));
    periBC_tag.insert(std::stoi(config->GetConfigValue(PERIODIC_RIGHT(idim))));
  }

  std::ifstream infile(grid_file_path_.c_str());
  std::string text;
  std::unordered_map<int, int> peribdry_nodes;
  while (getline(infile, text)) {
    if (text.find("$Elements", 0) != std::string::npos) {
      int num_total;
      infile >> num_total;
      ElemType elemtype;
      int bdry_index, gmshtype, tag, garbage, node_index;

      int num_peribdries = 0;
      for (int icell = 0; icell < num_total; icell++) {
        infile >> bdry_index >> gmshtype >> garbage >> tag >> garbage;
        elemtype = GetElemTypeFromGmshType(static_cast<GmshElemType>(gmshtype));
        if (periBC_tag.find(tag) != periBC_tag.end() &&
            Element::IsFace(dimension_, elemtype)) {
          const int num_nodes =
              DENEB_ELEMENT->GetElement(elemtype)->GetNumNodes(1);
          for (int inode = 0; inode < num_nodes; inode++) {
            infile >> node_index;
            peribdry_nodes[node_index] = -1;
          }
          num_peribdries++;
        }
        if (num_peribdries == num_global_peribdries_) break;
        getline(infile, text);
      }
      break;
    }
  }
  infile.close();

  const int num_periodic_nodes = static_cast<int>(peribdry_nodes.size());
  all_periodic_node_coords.resize(num_periodic_nodes * dimension_);

  infile.open(grid_file_path_.c_str());
  while (getline(infile, text)) {
    if (text.find("$Nodes", 0) != std::string::npos) {
      int num_total;
      infile >> num_total;
      int node_index;
      double coord[3];

      int num_nodes = 0;
      for (int inode = 0; inode < num_total; inode++) {
        infile >> node_index >> coord[0] >> coord[1] >> coord[2];
        auto&& iterator = peribdry_nodes.find(node_index);
        if (iterator != peribdry_nodes.end()) {
          iterator->second = num_nodes;
          for (int idim = 0; idim < dimension_; idim++)
            all_periodic_node_coords[num_nodes * dimension_ + idim] =
                coord[idim];
          num_nodes++;
        }
        if (num_nodes == num_periodic_nodes) break;
      }
      break;
    }
  }
  infile.close();

  all_periodic_nodes.resize(num_periodic_nodes);
  for (auto&& iterator : peribdry_nodes)
    all_periodic_nodes[iterator.second] = iterator.first;
  peribdry_nodes.clear();
}
void GmshGridReader::ReadCellData(const std::unordered_set<int>& node_list,
                                  std::vector<int>& cell_global_index,
                                  std::vector<ElemType>& cell_elemtype,
                                  std::vector<int>& cell_node_ptr,
                                  std::vector<int>& cell_node_ind) const {
  cell_global_index.clear();
  cell_elemtype.clear();
  cell_node_ptr.clear();
  cell_node_ptr.push_back(0);
  cell_node_ind.clear();

  std::ifstream infile(grid_file_path_.c_str());
  std::string text;
  while (getline(infile, text)) {
    if (text.find("$Elements", 0) != std::string::npos) {
      int num_total;
      infile >> num_total;
      ElemType elemtype;
      int cell_index, gmshtype, garbage, node_index;

      int num_cells = 0;
      for (int icell = 0; icell < num_total; icell++) {
        infile >> cell_index >> gmshtype >> garbage >> garbage >> garbage;
        elemtype = GetElemTypeFromGmshType(static_cast<GmshElemType>(gmshtype));
        if (!Element::IsFace(dimension_, elemtype)) {
          std::vector<int> local_nodes;
          const int num_nodes =
              DENEB_ELEMENT->GetElement(elemtype)->GetNumNodes(1);
          for (int inode = 0; inode < num_nodes; inode++) {
            infile >> node_index;
            local_nodes.push_back(node_index);
          }
          bool istarget = false;
          for (auto&& node : local_nodes) {
            if (node_list.find(node) != node_list.end()) {
              istarget = true;
              break;
            }
          }
          if (istarget == true) {
            cell_global_index.push_back(cell_index);  // note here
            cell_elemtype.push_back(elemtype);
            for (auto&& node : local_nodes) cell_node_ind.push_back(node);
            cell_node_ptr.push_back(static_cast<int>(cell_node_ind.size()));
          }
          num_cells++;
        }
        if (num_cells == num_global_cells_) break;
        getline(infile, text);
      }
      break;
    }
  }
  infile.close();
}
void GmshGridReader::ReadCellData(std::unordered_map<int, int>& cell_mapping,
                                  std::vector<int>& cell_order,
                                  std::vector<int>& cell_subnode_ptr,
                                  std::vector<int>& cell_subnode_ind) const {
  const int num_target_cell = static_cast<int>(cell_mapping.size());
  cell_order.clear();
  cell_subnode_ptr.clear();
  cell_subnode_ind.clear();
  cell_order.reserve(num_target_cell);
  cell_subnode_ptr.reserve(num_target_cell);
  cell_subnode_ptr.push_back(0);
  cell_subnode_ind.reserve(num_target_cell * (dimension_ + 1));

  std::ifstream infile(grid_file_path_.c_str());
  std::string text;
  while (getline(infile, text)) {
    if (text.find("$Elements", 0) != std::string::npos) {
      int num_total;
      infile >> num_total;
      ElemType elemtype;
      int cell_index, gmshtype, garbage, node_index;

      int num_cells = 0;
      for (int icell = 0; icell < num_total; icell++) {
        infile >> cell_index >> gmshtype >> garbage >> garbage >> garbage;
        auto&& iterator = cell_mapping.find(cell_index);
        if (iterator != cell_mapping.end()) {
          iterator->second = num_cells++;
          elemtype =
              GetElemTypeFromGmshType(static_cast<GmshElemType>(gmshtype));
          const int elemorder =
              GetGmshElemOrder(static_cast<GmshElemType>(gmshtype));
          const int num_nodes =
              DENEB_ELEMENT->GetElement(elemtype)->GetNumNodes(elemorder);
          for (int inode = 0; inode < num_nodes; inode++) {
            infile >> node_index;
            cell_subnode_ind.push_back(node_index);  // note here
          }
          cell_subnode_ptr.push_back(static_cast<int>(cell_subnode_ind.size()));
          cell_order.push_back(elemorder);
        }
        if (num_cells == num_target_cell) break;
        getline(infile, text);
      }
      break;
    }
  }
  infile.close();
}
void GmshGridReader::ReadNodeData(std::unordered_map<int, int>& node_mapping,
                                  std::vector<double>& node_coords) const {
  const int num_target_nodes = static_cast<int>(node_mapping.size());
  node_coords.clear();
  node_coords.reserve(num_target_nodes * dimension_);

  std::ifstream infile(grid_file_path_.c_str());
  std::string text;
  while (getline(infile, text)) {
    if (text.find("$Nodes", 0) != std::string::npos) {
      int num_global_nodes;
      infile >> num_global_nodes;
      int node_index;
      double coord[3];

      int num_nodes = 0;
      for (int inode = 0; inode < num_global_nodes; inode++) {
        infile >> node_index >> coord[0] >> coord[1] >> coord[2];
        auto&& iterator = node_mapping.find(node_index);
        if (iterator != node_mapping.end()) {
          for (int idim = 0; idim < dimension_; idim++)
            node_coords.push_back(coord[idim]);
          iterator->second = num_nodes++;
        }
        if (num_nodes == num_target_nodes) break;
      }
      break;
    }
  }
  infile.close();
}
int GmshGridReader::GetGmshElemOrder(const GmshElemType gmshtype) const {
  switch (gmshtype) {
    case GmshElemType::TRIS_P1:
    case GmshElemType::QUAD_P1:
    case GmshElemType::TETS_P1:
    case GmshElemType::HEXA_P1:
    case GmshElemType::PRIS_P1:
    case GmshElemType::PYRA_P1:
      return 1;
    case GmshElemType::TRIS_P2:
    case GmshElemType::QUAD_P2:
    case GmshElemType::TETS_P2:
    case GmshElemType::HEXA_P2:
    case GmshElemType::PRIS_P2:
    case GmshElemType::PYRA_P2:
      return 2;
    case GmshElemType::TRIS_P3:
    case GmshElemType::QUAD_P3:
    case GmshElemType::TETS_P3:
    case GmshElemType::HEXA_P3:
    case GmshElemType::PRIS_P3:
    case GmshElemType::PYRA_P3:
      return 3;
    case GmshElemType::TRIS_P4:
    case GmshElemType::QUAD_P4:
    case GmshElemType::TETS_P4:
    case GmshElemType::HEXA_P4:
    case GmshElemType::PRIS_P4:
    case GmshElemType::PYRA_P4:
      return 4;
    case GmshElemType::TRIS_P5:
    case GmshElemType::QUAD_P5:
    case GmshElemType::TETS_P5:
    case GmshElemType::HEXA_P5:
    case GmshElemType::PRIS_P5:
      return 5;
    case GmshElemType::QUAD_P6:
      return 6;
  }
  ERROR_MESSAGE("Invalid Gmsh-element type\n");
  return -1;
}

ElemType GmshGridReader::GetElemTypeFromGmshType(
    const GmshElemType gmshtype) const {
  switch (gmshtype) {
    case GmshElemType::LINE_P1:
    case GmshElemType::LINE_P2:
    case GmshElemType::LINE_P3:
    case GmshElemType::LINE_P4:
    case GmshElemType::LINE_P5:
    case GmshElemType::LINE_P6:
      return ElemType::LINE;
    case GmshElemType::TRIS_P1:
    case GmshElemType::TRIS_P2:
    case GmshElemType::TRIS_P3:
    case GmshElemType::TRIS_P4:
    case GmshElemType::TRIS_P5:
      return ElemType::TRIS;
    case GmshElemType::QUAD_P1:
    case GmshElemType::QUAD_P2:
    case GmshElemType::QUAD_P3:
    case GmshElemType::QUAD_P4:
    case GmshElemType::QUAD_P5:
    case GmshElemType::QUAD_P6:
      return ElemType::QUAD;
    case GmshElemType::TETS_P1:
    case GmshElemType::TETS_P2:
    case GmshElemType::TETS_P3:
    case GmshElemType::TETS_P4:
    case GmshElemType::TETS_P5:
      return ElemType::TETS;
    case GmshElemType::HEXA_P1:
    case GmshElemType::HEXA_P2:
    case GmshElemType::HEXA_P3:
    case GmshElemType::HEXA_P4:
    case GmshElemType::HEXA_P5:
      return ElemType::HEXA;
    case GmshElemType::PRIS_P1:
    case GmshElemType::PRIS_P2:
    case GmshElemType::PRIS_P3:
    case GmshElemType::PRIS_P4:
    case GmshElemType::PRIS_P5:
      return ElemType::PRIS;
    case GmshElemType::PYRA_P1:
    case GmshElemType::PYRA_P2:
    case GmshElemType::PYRA_P3:
    case GmshElemType::PYRA_P4:
      return ElemType::PYRA;
  }
  ERROR_MESSAGE("Invalid Gmsh-element type\n");
  return static_cast<ElemType>(0);
}

// Gmsh2DTrisGridReader
void Gmsh2DTrisGridReader::GenerateGmshFile(void) const {
  if (MYRANK != MASTER_NODE) return;

  START_TIMER();
  MASTER_MESSAGE(
      "Generate Gmsh2DTris file"
      "\n\tgrid file path: " +
      grid_file_path_ + "\n");

  std::vector<double> start_coords;
  std::vector<double> end_coords;
  std::vector<int> num_division;
  std::vector<int> bdry_tags;
  const std::vector<std::string> dir = {"X", "Y"};
  auto& config = AVOCADO_CONFIG;
  for (int idim = 0; idim < dimension_; idim++) {
    start_coords.push_back(
        std::stod(config->GetConfigValue(GEN_GMSH_LEFT(dir[idim]))));
    end_coords.push_back(
        std::stod(config->GetConfigValue(GEN_GMSH_RIGHT(dir[idim]))));
    num_division.push_back(
        std::stoi(config->GetConfigValue(GEN_GMSH_DIVISION(dir[idim]))));
    bdry_tags.push_back(
        std::stoi(config->GetConfigValue(GEN_GMSH_LEFT_TAG(dir[idim]))));
    bdry_tags.push_back(
        std::stoi(config->GetConfigValue(GEN_GMSH_RIGHT_TAG(dir[idim]))));
    MASTER_MESSAGE("\t" + dir[idim] + "-direction\n");
    MASTER_MESSAGE("\t\tRange: (" + std::to_string(start_coords[idim]) + ", " +
                   std::to_string(end_coords[idim]) + ")\n");
    MASTER_MESSAGE(
        "\t\tNumber of division: " + std::to_string(num_division[idim]) + "\n");
    MASTER_MESSAGE("\t\tBC tags: (" + std::to_string(bdry_tags[idim * 2 + 0]) +
                   ", " + std::to_string(bdry_tags[idim * 2 + 1]) + ")\n");
  }

  avocado::MakeDirectory(grid_file_path_);
  std::ofstream outfile(grid_file_path_);

  outfile << std::setprecision(16);
  outfile << "$MeshFormat" << std::endl;
  outfile << "2.2 0 8" << std::endl;
  outfile << "$EndMeshFormat" << std::endl;

  outfile << "$Comments" << std::endl;
  outfile << "Gmsh2DTris by " << AVOCADO_MPI->GetCodeName() << std::endl;
  for (int idim = 0; idim < dimension_; idim++) {
    outfile << "\t" << dir[idim] + "-direction\n";
    outfile << "\t\tRange: (" + std::to_string(start_coords[idim]) + ", " +
                   std::to_string(end_coords[idim]) + ")\n";
    outfile << "\t\tNumber of division: " + std::to_string(num_division[idim]) +
                   "\n";
    outfile << "\t\tBC tags: (" + std::to_string(bdry_tags[idim * 2 + 0]) +
                   ", " + std::to_string(bdry_tags[idim * 2 + 1]) + ")\n";
  }
  outfile << "$Comments" << std::endl;

  int num_nodes = (num_division[0] + 1) * (num_division[1] + 1);
  const double dx = (end_coords[0] - start_coords[0]) / double(num_division[0]);
  const double dy = (end_coords[1] - start_coords[1]) / double(num_division[1]);

  outfile << "$Nodes" << std::endl;
  outfile << num_nodes << std::endl;
  int index = 0;
  for (int iy = 0; iy < num_division[1] + 1; iy++)
    for (int ix = 0; ix < num_division[0] + 1; ix++)
      outfile << ++index << " " << start_coords[0] + dx * double(ix) << " "
              << start_coords[1] + dy * double(iy) << " 0" << std::endl;
  outfile << "$EndNodes" << std::endl;
  MASTER_MESSAGE("\tNode data are written.\n");

  int num_elements = 2 * (num_nodes - 1);
  outfile << "$Elements" << std::endl;
  outfile << num_elements << std::endl;
  index = 0;
  //          3*
  //    4  --------- 3
  //    |            |
  // 4* |            |  2*
  //    1  --------- 2
  //          1*
  for (int ix = 0; ix < num_division[0]; ix++)  // 1 - 2
    outfile << ++index << " 1 2 " << bdry_tags[2] << " 1 " << ix + 1 << " "
            << ix + 2 << std::endl;
  for (int iy = 0; iy < num_division[1]; iy++)  // 2 - 3
    outfile << ++index << " 1 2 " << bdry_tags[1] << " 2 "
            << (num_division[0] + 1) * (iy + 1) << " "
            << (num_division[0] + 1) * (iy + 2) << std::endl;
  for (int ix = 0; ix < num_division[0]; ix++)  // 4 - 3
    outfile << ++index << " 1 2 " << bdry_tags[3] << " 3 "
            << ix + (num_division[0] + 1) * num_division[1] + 1 << " "
            << ix + (num_division[0] + 1) * num_division[1] + 2 << std::endl;
  for (int iy = 0; iy < num_division[1]; iy++)  // 1 - 4
    outfile << ++index << " 1 2 " << bdry_tags[0] << " 4 "
            << (num_division[0] + 1) * iy + 1 << " "
            << (num_division[0] + 1) * (iy + 1) + 1 << std::endl;

  for (int iy = 0; iy < num_division[1]; iy++) {
    for (int ix = 0; ix < num_division[0]; ix++) {
      outfile << ++index << " 2 2 1 0 " << ix + iy * (num_division[0] + 1) + 1
              << " " << ix + iy * (num_division[0] + 1) + 2 << " "
              << ix + (iy + 1) * (num_division[0] + 1) + 2 << std::endl;
      outfile << ++index << " 2 2 1 0 " << ix + iy * (num_division[0] + 1) + 1
              << " " << ix + (iy + 1) * (num_division[0] + 1) + 2 << " "
              << ix + (iy + 1) * (num_division[0] + 1) + 1 << std::endl;
    }
  }
  outfile << "$EndElements" << std::endl;
  MASTER_MESSAGE("\tElement data are written.\n");
  MASTER_MESSAGE("Gmsh2DTris file is written.(Time: " +
                 std::to_string(STOP_TIMER()) + "s)\n");

  outfile.close();
}

// Gmsh2DQuadGridReader
void Gmsh2DQuadGridReader::GenerateGmshFile(void) const {
  if (MYRANK != MASTER_NODE) return;

  START_TIMER();
  MASTER_MESSAGE("Gmsh2DQuad file is being generated.\n");
  MASTER_MESSAGE("\tgrid file path: " + grid_file_path_ + "\n");

  std::vector<double> start_coords;
  std::vector<double> end_coords;
  std::vector<int> num_division;
  std::vector<int> bdry_tags;
  const std::vector<std::string> dir = {"X", "Y"};
  auto& config = AVOCADO_CONFIG;
  for (int idim = 0; idim < dimension_; idim++) {
    start_coords.push_back(
        std::stod(config->GetConfigValue(GEN_GMSH_LEFT(dir[idim]))));
    end_coords.push_back(
        std::stod(config->GetConfigValue(GEN_GMSH_RIGHT(dir[idim]))));
    num_division.push_back(
        std::stoi(config->GetConfigValue(GEN_GMSH_DIVISION(dir[idim]))));
    bdry_tags.push_back(
        std::stoi(config->GetConfigValue(GEN_GMSH_LEFT_TAG(dir[idim]))));
    bdry_tags.push_back(
        std::stoi(config->GetConfigValue(GEN_GMSH_RIGHT_TAG(dir[idim]))));
    MASTER_MESSAGE("\t" + dir[idim] + "-direction\n");
    MASTER_MESSAGE("\t\tRange: (" + std::to_string(start_coords[idim]) + ", " +
                   std::to_string(end_coords[idim]) + ")\n");
    MASTER_MESSAGE(
        "\t\tNumber of division: " + std::to_string(num_division[idim]) + "\n");
    MASTER_MESSAGE("\t\tBC tags: (" + std::to_string(bdry_tags[idim * 2 + 0]) +
                   ", " + std::to_string(bdry_tags[idim * 2 + 1]) + ")\n");
  }

  avocado::MakeDirectory(grid_file_path_);
  std::ofstream outfile(grid_file_path_);

  outfile << std::setprecision(16);
  outfile << "$MeshFormat" << std::endl;
  outfile << "2.2 0 8" << std::endl;
  outfile << "$EndMeshFormat" << std::endl;

  outfile << "$Comments" << std::endl;
  outfile << "Gmsh2DQuad by " << AVOCADO_MPI->GetCodeName() << std::endl;
  for (int idim = 0; idim < dimension_; idim++) {
    outfile << "\t" << dir[idim] + "-direction\n";
    outfile << "\t\tRange: (" + std::to_string(start_coords[idim]) + ", " +
                   std::to_string(end_coords[idim]) + ")\n";
    outfile << "\t\tNumber of division: " + std::to_string(num_division[idim]) +
                   "\n";
    outfile << "\t\tBC tags: (" + std::to_string(bdry_tags[idim * 2 + 0]) +
                   ", " + std::to_string(bdry_tags[idim * 2 + 1]) + ")\n";
  }
  outfile << "$Comments" << std::endl;

  int num_nodes = (num_division[0] + 1) * (num_division[1] + 1);
  const double dx = (end_coords[0] - start_coords[0]) / double(num_division[0]);
  const double dy = (end_coords[1] - start_coords[1]) / double(num_division[1]);

  outfile << "$Nodes" << std::endl;
  outfile << num_nodes << std::endl;
  int index = 0;
  for (int iy = 0; iy < num_division[1] + 1; iy++)
    for (int ix = 0; ix < num_division[0] + 1; ix++)
      outfile << ++index << " " << start_coords[0] + dx * double(ix) << " "
              << start_coords[1] + dy * double(iy) << " 0" << std::endl;
  outfile << "$EndNodes" << std::endl;
  MASTER_MESSAGE("\tNode data are written.\n");

  int num_elements = (num_division[0] + 2) * (num_division[1] + 2) - 4;
  outfile << "$Elements" << std::endl;
  outfile << num_elements << std::endl;
  index = 0;
  //          3*
  //    4  --------- 3
  //    |            |
  // 4* |            |  2*
  //    1  --------- 2
  //          1*
  for (int ix = 0; ix < num_division[0]; ix++)  // 1 - 2
    outfile << ++index << " 1 2 " << bdry_tags[2] << " 1 " << ix + 1 << " "
            << ix + 2 << std::endl;
  for (int iy = 0; iy < num_division[1]; iy++)  // 2 - 3
    outfile << ++index << " 1 2 " << bdry_tags[1] << " 2 "
            << (num_division[0] + 1) * (iy + 1) << " "
            << (num_division[0] + 1) * (iy + 2) << std::endl;
  for (int ix = 0; ix < num_division[0]; ix++)  // 4 - 3
    outfile << ++index << " 1 2 " << bdry_tags[3] << " 3 "
            << ix + (num_division[0] + 1) * num_division[1] + 1 << " "
            << ix + (num_division[0] + 1) * num_division[1] + 2 << std::endl;
  for (int iy = 0; iy < num_division[1]; iy++)  // 1 - 4
    outfile << ++index << " 1 2 " << bdry_tags[0] << " 4 "
            << (num_division[0] + 1) * iy + 1 << " "
            << (num_division[0] + 1) * (iy + 1) + 1 << std::endl;

  for (int iy = 0; iy < num_division[1]; iy++)
    for (int ix = 0; ix < num_division[0]; ix++)
      outfile << ++index << " 3 2 1 0 " << ix + iy * (num_division[0] + 1) + 1
              << " " << ix + iy * (num_division[0] + 1) + 2 << " "
              << ix + (iy + 1) * (num_division[0] + 1) + 2 << " "
              << ix + (iy + 1) * (num_division[0] + 1) + 1 << std::endl;

  outfile << "$EndElements" << std::endl;
  MASTER_MESSAGE("\tElement data are written.\n");
  MASTER_MESSAGE("Gmsh2DQuad file is written.(Time: " +
                 std::to_string(STOP_TIMER()) + "s)\n");

  outfile.close();
}

// Gmsh3DHexaGridReader
void Gmsh3DHexaGridReader::GenerateGmshFile(void) const {
  if (MYRANK != MASTER_NODE) return;

  START_TIMER();
  MASTER_MESSAGE("Gmsh3DHexa file is being generated.\n");
  MASTER_MESSAGE("\tgrid file path: " + grid_file_path_ + "\n");

  std::vector<double> start_coords;
  std::vector<double> end_coords;
  std::vector<int> num_division;
  std::vector<int> bdry_tags;
  const std::vector<std::string> dir = {"X", "Y", "Z"};
  auto& config = AVOCADO_CONFIG;
  for (int idim = 0; idim < dimension_; idim++) {
    start_coords.push_back(
        std::stod(config->GetConfigValue(GEN_GMSH_LEFT(dir[idim]))));
    end_coords.push_back(
        std::stod(config->GetConfigValue(GEN_GMSH_RIGHT(dir[idim]))));
    num_division.push_back(
        std::stoi(config->GetConfigValue(GEN_GMSH_DIVISION(dir[idim]))));
    bdry_tags.push_back(
        std::stoi(config->GetConfigValue(GEN_GMSH_LEFT_TAG(dir[idim]))));
    bdry_tags.push_back(
        std::stoi(config->GetConfigValue(GEN_GMSH_RIGHT_TAG(dir[idim]))));
    MASTER_MESSAGE("\t" + dir[idim] + "-direction\n");
    MASTER_MESSAGE("\t\tRange: (" + std::to_string(start_coords[idim]) + ", " +
                   std::to_string(end_coords[idim]) + ")\n");
    MASTER_MESSAGE(
        "\t\tNumber of division: " + std::to_string(num_division[idim]) + "\n");
    MASTER_MESSAGE("\t\tBC tags: (" + std::to_string(bdry_tags[idim * 2 + 0]) +
                   ", " + std::to_string(bdry_tags[idim * 2 + 1]) + ")\n");
  }

  avocado::MakeDirectory(grid_file_path_);
  std::ofstream outfile(grid_file_path_);

  outfile << std::setprecision(16);
  outfile << "$MeshFormat" << std::endl;
  outfile << "2.2 0 8" << std::endl;
  outfile << "$EndMeshFormat" << std::endl;

  outfile << "$Comments" << std::endl;
  outfile << "Gmsh3DHexa by " << AVOCADO_MPI->GetCodeName() << std::endl;
  for (int idim = 0; idim < dimension_; idim++) {
    outfile << "\t" << dir[idim] + "-direction\n";
    outfile << "\t\tRange: (" + std::to_string(start_coords[idim]) + ", " +
                   std::to_string(end_coords[idim]) + ")\n";
    outfile << "\t\tNumber of division: " + std::to_string(num_division[idim]) +
                   "\n";
    outfile << "\t\tBC tags: (" + std::to_string(bdry_tags[idim * 2 + 0]) +
                   ", " + std::to_string(bdry_tags[idim * 2 + 1]) + ")\n";
  }
  outfile << "$Comments" << std::endl;

  int num_nodes =
      (num_division[0] + 1) * (num_division[1] + 1) * (num_division[2] + 1);
  const double dx = (end_coords[0] - start_coords[0]) / double(num_division[0]);
  const double dy = (end_coords[1] - start_coords[1]) / double(num_division[1]);
  const double dz = (end_coords[2] - start_coords[2]) / double(num_division[2]);

  outfile << "$Nodes" << std::endl;
  outfile << num_nodes << std::endl;
  int index = 0;
  for (int iz = 0; iz < num_division[2] + 1; iz++)
    for (int iy = 0; iy < num_division[1] + 1; iy++)
      for (int ix = 0; ix < num_division[0] + 1; ix++)
        outfile << ++index << " " << start_coords[0] + dx * double(ix) << " "
                << start_coords[1] + dy * double(iy) << " "
                << start_coords[2] + dz * double(iz) << std::endl;
  outfile << "$EndNodes" << std::endl;
  MASTER_MESSAGE("\tNode data are written.\n");

  int num_elements = num_division[0] * num_division[1] * num_division[2] +
                     2 * (num_division[0] * num_division[1] +
                          num_division[1] * num_division[2] +
                          num_division[2] * num_division[0]);
  outfile << "$Elements" << std::endl;
  outfile << num_elements << std::endl;
  index = 0;
  //		    2*
  //          z
  //	 3----------2
  //	 |\     ^   |\ 
	//	 | \ 3* |   | \
	//	 |  \   |   |  \
	//	 |   7------+---6
  //5* |   |  +---+---|->y 6*
  //	 0---+---\--1   |
  //	 \   |    \  \  |
  //	  \  |     \  \ |
  //	   \ |  1*  x  \|
  //	     4----------5
  //               4*
  for (int iy = 0; iy < num_division[1]; iy++)  // 0 - 1 - 4 - 5
    for (int ix = 0; ix < num_division[0]; ix++)
      outfile << ++index << " 3 2 " << bdry_tags[4] << " 1 "
              << ix + iy * (num_division[0] + 1) + 1 << " "
              << ix + (iy + 1) * (num_division[0] + 1) + 1 << " "
              << ix + 2 + (iy + 1) * (num_division[0] + 1) << " "
              << ix + 2 + iy * (num_division[0] + 1) << std::endl;

  for (int iy = 0; iy < num_division[1]; iy++)  // 3 - 2 - 7 - 6
    for (int ix = 0; ix < num_division[0]; ix++)
      outfile
          << ++index << " 3 2 " << bdry_tags[5] << " 2 "
          << ix + 1 + iy * (num_division[0] + 1) +
                 num_division[2] * (num_division[0] + 1) * (num_division[1] + 1)
          << " "
          << ix + 1 + (iy + 1) * (num_division[0] + 1) +
                 num_division[2] * (num_division[0] + 1) * (num_division[1] + 1)
          << " "
          << ix + 2 + (iy + 1) * (num_division[0] + 1) +
                 num_division[2] * (num_division[0] + 1) * (num_division[1] + 1)
          << " "
          << ix + 2 + iy * (num_division[0] + 1) +
                 num_division[2] * (num_division[0] + 1) * (num_division[1] + 1)
          << std::endl;

  for (int iz = 0; iz < num_division[2]; iz++)  // 1 - 2 - 3 - 0
    for (int iy = 0; iy < num_division[1]; iy++)
      outfile << ++index << " 3 2 " << bdry_tags[0] << " 3 "
              << iy * (num_division[0] + 1) +
                     iz * (num_division[0] + 1) * (num_division[1] + 1) + 1
              << " "
              << (iy + 1) * (num_division[0] + 1) +
                     iz * (num_division[0] + 1) * (num_division[1] + 1) + 1
              << " "
              << (iy + 1) * (num_division[0] + 1) +
                     (iz + 1) * (num_division[0] + 1) * (num_division[1] + 1) +
                     1
              << " "
              << iy * (num_division[0] + 1) +
                     (iz + 1) * (num_division[0] + 1) * (num_division[1] + 1) +
                     1
              << std::endl;

  for (int iz = 0; iz < num_division[2]; iz++)  // 4 - 5 - 6 - 7
    for (int iy = 0; iy < num_division[1]; iy++)
      outfile << ++index << " 3 2 " << bdry_tags[1] << " 4 "
              << num_division[0] + 1 + iy * (num_division[0] + 1) +
                     iz * (num_division[0] + 1) * (num_division[1] + 1)
              << " "
              << num_division[0] + 1 + (iy + 1) * (num_division[0] + 1) +
                     iz * (num_division[0] + 1) * (num_division[1] + 1)
              << " "
              << num_division[0] + 1 + (iy + 1) * (num_division[0] + 1) +
                     (iz + 1) * (num_division[0] + 1) * (num_division[1] + 1)
              << " "
              << num_division[0] + iy * (num_division[0] + 1) +
                     (iz + 1) * (num_division[0] + 1) * (num_division[1] + 1) +
                     1
              << std::endl;

  for (int iz = 0; iz < num_division[2]; iz++)  // 0 - 3 - 4 - 7
    for (int ix = 0; ix < num_division[0]; ix++)
      outfile << ++index << " 3 2 " << bdry_tags[2] << " 5 "
              << ix + 1 + iz * (num_division[0] + 1) * (num_division[1] + 1)
              << " "
              << ix + 2 + iz * (num_division[0] + 1) * (num_division[1] + 1)
              << " "
              << ix + 2 +
                     (iz + 1) * (num_division[0] + 1) * (num_division[1] + 1)
              << " "
              << ix + (iz + 1) * (num_division[0] + 1) * (num_division[1] + 1) +
                     1
              << std::endl;

  for (int iz = 0; iz < num_division[2]; iz++)  // 1 - 2 - 5 - 6
    for (int ix = 0; ix < num_division[0]; ix++)
      outfile << ++index << " 3 2 " << bdry_tags[3] << " 6 "
              << ix + 1 + num_division[1] * (num_division[0] + 1) +
                     iz * (num_division[0] + 1) * (num_division[1] + 1)
              << " "
              << ix + 2 + num_division[1] * (num_division[0] + 1) +
                     iz * (num_division[0] + 1) * (num_division[1] + 1)
              << " "
              << ix + 2 + num_division[1] * (num_division[0] + 1) +
                     (iz + 1) * (num_division[0] + 1) * (num_division[1] + 1)
              << " "
              << ix + num_division[1] * (num_division[0] + 1) +
                     (iz + 1) * (num_division[0] + 1) * (num_division[1] + 1) +
                     1
              << std::endl;

  for (int iz = 0; iz < num_division[2]; iz++) {
    for (int iy = 0; iy < num_division[1]; iy++) {
      for (int ix = 0; ix < num_division[0]; ix++) {
        std::vector<int> nodes_index;
        nodes_index.push_back(ix + 1 + iy * (num_division[0] + 1) +
                              iz * (num_division[0] + 1) *
                                  (num_division[1] + 1));
        nodes_index.push_back(ix + 1 + (iy + 1) * (num_division[0] + 1) +
                              iz * (num_division[0] + 1) *
                                  (num_division[1] + 1));
        nodes_index.push_back(ix + 1 + (iy + 1) * (num_division[0] + 1) +
                              (iz + 1) * (num_division[0] + 1) *
                                  (num_division[1] + 1));
        nodes_index.push_back(ix + 1 + iy * (num_division[0] + 1) +
                              (iz + 1) * (num_division[0] + 1) *
                                  (num_division[1] + 1));
        nodes_index.push_back(ix + 2 + iy * (num_division[0] + 1) +
                              iz * (num_division[0] + 1) *
                                  (num_division[1] + 1));
        nodes_index.push_back(ix + 2 + (iy + 1) * (num_division[0] + 1) +
                              iz * (num_division[0] + 1) *
                                  (num_division[1] + 1));
        nodes_index.push_back(ix + 2 + (iy + 1) * (num_division[0] + 1) +
                              (iz + 1) * (num_division[0] + 1) *
                                  (num_division[1] + 1));
        nodes_index.push_back(ix + 2 + iy * (num_division[0] + 1) +
                              (iz + 1) * (num_division[0] + 1) *
                                  (num_division[1] + 1));
        outfile << ++index << " 5 2 1 0";
        for (auto&& node : nodes_index) outfile << " " << node;
        outfile << std::endl;
      }
    }
  }

  outfile << "$EndElements" << std::endl;
  MASTER_MESSAGE("\tElement data are written.\n");
  MASTER_MESSAGE("Gmsh3DHexa file is written.(Time: " +
                 std::to_string(STOP_TIMER()) + "s)\n");

  outfile.close();
}
}  // namespace deneb