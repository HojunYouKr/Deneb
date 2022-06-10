#include "deneb_contour.h"

#include <fstream>
#include <unordered_set>

#include "avocado.h"
#include "deneb_config_macro.h"
#include "deneb_data.h"
#include "deneb_equation.h"
#include "deneb_header.h"
#include "deneb_jacobian.h"

namespace deneb {
std::shared_ptr<Contour> DENEB_CONTOUR_NAME;
Contour::Contour(const int cell_post_order, const int face_post_order)
    : strandid_(-1),
      cell_post_order_(cell_post_order),
      face_post_order_(face_post_order) {
  MASTER_MESSAGE(avocado::GetTitle("Contour"));

  MASTER_MESSAGE("Cell post order = " + std::to_string(cell_post_order_) +
                 "\n");
  MASTER_MESSAGE("Face post order = " + std::to_string(face_post_order_) +
                 "\n");
  if (cell_post_order_ == 0 && face_post_order_ == 0)
    ERROR_MESSAGE("Both post orders are zero.\n");
}
void Contour::BuildData(void) {
  SYNCRO();
  MASTER_MESSAGE(avocado::GetTitle("Contour::BuildData()"));

  dimension_ = DENEB_EQUATION->GetDimension();
  cell_variables_ = DENEB_EQUATION->GetCellVariableNames();
  face_variables_ = DENEB_EQUATION->GetFaceVariableNames();

  // cell data
  START_TIMER();
  MASTER_MESSAGE("Building cell data... ");
  if (cell_post_order_ != 0)
    BuildCellData();
  else
    MASTER_MESSAGE("Nothing done... ");
  SYNCRO();
  MASTER_MESSAGE("Complete. (Time: " + std::to_string(STOP_TIMER()) + "s)\n");

  // face data
  START_TIMER();
  MASTER_MESSAGE("Building face data... ");
  if (face_post_order_ != 0)
    BuildFaceData();
  else
    MASTER_MESSAGE("Nothing done... ");
  SYNCRO();
  MASTER_MESSAGE("Complete. (Time: " + std::to_string(STOP_TIMER()) + "s)\n");

  {
    int max_num_points = 1;
    int num_variables =
        std::max(cell_variables_.size(), face_variables_.size());
    const int& num_states = DENEB_EQUATION->GetNumStates();

    {
      const std::vector<int>& num_points = cell_.num_nodes_;
      if (!num_points.empty()) {
        const int max_pt =
            *std::max_element(num_points.begin(), num_points.end());
        max_num_points = std::max(max_num_points, max_pt);
      }
    }
    {
      for (auto&& face : faces_) {
        const std::vector<int>& num_points = face.num_nodes_;
        if (!num_points.empty()) {
          const int max_pt =
              *std::max_element(num_points.begin(), num_points.end());
          max_num_points = std::max(max_num_points, max_pt);
        }
      }
    }
    sol_.resize(max_num_points * num_states);
    sol_grad_.resize(max_num_points * dimension_ * num_states);
    sol_post_.resize(max_num_points * num_variables);
  }
}
void Contour::CellSolution(const std::string& filename, const double* solution,
                           const double time) {
  if (cell_post_order_ == 0) return;
  if (strandid_ < 0) ERROR_MESSAGE("Strand-id is not initialized.\n");
  MASTER_MESSAGE("Exporting cell solution: " + filename + "\n");
  if (MYRANK == MASTER_NODE) avocado::MakeDirectory(filename);

  const FileType filetype = FileType::SOLUTION;
  const ZoneType zonetype =
      (dimension_ == 3) ? ZoneType::FETETS : ZoneType::FETRIS;
  const std::vector<std::string>& names = cell_variables_;
  const int num_variables = static_cast<int>(cell_variables_.size());
  const int& num_states = DENEB_EQUATION->GetNumStates();
  const int& num_bases = DENEB_DATA->GetNumBases();
  const int& num_cells = DENEB_DATA->GetNumCells();
  const int sb = num_states * num_bases;

  std::ofstream file;
  MPI_File mpifile = nullptr;
  int offset = 0;
  std::vector<int> offset_gather(NDOMAIN);

  std::vector<int> header;
  WriteTecplotVersion(header);
  WriteIntData(header, 1);  // default
  WriteIntData(header,
               static_cast<int>(filetype));  // file type
  const std::string title = "Tecplot file written by " DENEB_VERSION;
  WriteStringData(header, title);                           // title
  WriteIntData(header, num_variables);                      // num var
  for (auto&& name : names) WriteStringData(header, name);  // var names
  WriteReal32Data(header, 299.0);                           // zone marker
  const std::string cell_name = "cell";
  WriteStringData(header, cell_name);                // zone name
  WriteIntData(header, -1);                          // parent zone
  WriteIntData(header, strandid_++);                 // strand ID
  WriteReal64Data(header, time);                     // solution time
  WriteIntData(header, -1);                          // default = -1
  WriteIntData(header, static_cast<int>(zonetype));  // zone type
  WriteIntData(header, 0);                           // specify var location
  // Are raw local 1-to-1 face neighbors supplied?
  WriteIntData(header, 0);
  // number of miscellaneous user-defined face neighbor connections
  WriteIntData(header, 0);
  const int num_total_nodes =
      AVOCADO_MPI->Reduce(cell_.num_total_nodes_, avocado::MPI::Op::SUM);
  const int num_total_cells =
      AVOCADO_MPI->Reduce(cell_.num_total_cells_, avocado::MPI::Op::SUM);
  WriteIntData(header, num_total_nodes);  // number of nodes
  WriteIntData(header, num_total_cells);  // number of cells
  WriteIntData(header, 0);                // ICellDim
  WriteIntData(header, 0);                // JCellDim
  WriteIntData(header, 0);                // KCellDim
  WriteIntData(header, 0);                // auxiliary name
  WriteReal32Data(header, 357.0);         // EOH marker

  if (NDOMAIN == 1) {
    file.open(filename, std::ios::binary);
    file.write((char*)&header[0], sizeof(int) * header.size());
  } else {
    MPI_File_open(MPI_COMM_WORLD, filename.c_str(),
                  MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
    MPI_File_set_view(mpifile, offset, MPI_INT, MPI_INT, "native",
                      MPI_INFO_NULL);
    if (MYRANK == MASTER_NODE)
      MPI_File_write(mpifile, &header[0], header.size(), MPI_INT,
                     MPI_STATUS_IGNORE);
    offset += static_cast<int>(header.size()) * sizeof(int);
  }

  std::vector<int> zone;
  WriteReal32Data(zone, 299.0);  // zone marker
  for (int ivar = 0; ivar < num_variables; ivar++)
    WriteIntData(zone, 1);  // variable data format
  WriteIntData(zone, 0);    // has passive variables
  WriteIntData(zone, 0);    // has variable sharing
  // zeros based zone number to share connectivity list with (-1=no sharing)
  WriteIntData(zone, -1);

  std::vector<std::vector<float>> variables(num_variables);
  for (int ivar = 0; ivar < num_variables; ivar++)
    variables[ivar].reserve(cell_.num_total_nodes_);

  for (int icell = 0; icell < num_cells; icell++) {
    const int& num_points = cell_.num_nodes_[icell];

    gemmABT(1.0, &cell_.basis_[icell][0], &solution[icell * sb], 0.0, &sol_[0],
            num_points, num_bases, num_states);
    avocado::Kernel1::f54(&cell_.basis_grad_[icell][0], &solution[icell * sb],
                          &sol_grad_[0], dimension_, num_points, num_bases,
                          num_states, 1.0, 0.0);

    DENEB_EQUATION->GetCellPostSolution(num_points, sol_, sol_grad_, sol_post_);
    for (int ipoint = 0; ipoint < num_points; ipoint++)
      for (int ivar = 0; ivar < num_variables; ivar++)
        variables[ivar].push_back(sol_post_[ipoint * num_variables + ivar]);
  }

  for (int ivar = 0; ivar < num_variables; ivar++) {
    const double min_var =
        *std::min_element(variables[ivar].begin(), variables[ivar].end());
    const double max_var =
        *std::max_element(variables[ivar].begin(), variables[ivar].end());
    WriteReal64Data(zone, AVOCADO_MPI->Reduce(min_var, avocado::MPI::Op::MIN));
    WriteReal64Data(zone, AVOCADO_MPI->Reduce(max_var, avocado::MPI::Op::MAX));
  }

  if (NDOMAIN == 1) {
    file.write((char*)&zone[0], sizeof(int) * static_cast<int>(zone.size()));
    for (int ivar = 0; ivar < num_variables; ivar++)
      file.write((char*)&variables[ivar][0],
                 sizeof(float) * cell_.num_total_nodes_);
    file.close();
  } else {
    if (MYRANK == MASTER_NODE)
      MPI_File_write(mpifile, &zone[0], static_cast<int>(zone.size()), MPI_INT,
                     MPI_STATUS_IGNORE);
    offset += static_cast<int>(zone.size()) * sizeof(int);

    const int data_size = cell_.num_total_nodes_ * sizeof(float);
    const int total_data_size =
        AVOCADO_MPI->Reduce(data_size, avocado::MPI::Op::SUM);
    MPI_Allgather(&data_size, 1, MPI_INT, &offset_gather[0], 1, MPI_INT,
                  MPI_COMM_WORLD);
    for (int idomain = 0; idomain < MYRANK; idomain++)
      offset += offset_gather[idomain];

    for (int ivar = 0; ivar < num_variables; ivar++) {
      MPI_File_set_view(mpifile, offset, MPI_INT, MPI_INT, "native",
                        MPI_INFO_NULL);
      MPI_File_write(mpifile, &variables[ivar][0], cell_.num_total_nodes_,
                     MPI_FLOAT, MPI_STATUS_IGNORE);
      offset += total_data_size;
    }
    MPI_File_close(&mpifile);
  }
}
void Contour::CellGrid(const std::string& filename) {
  if (cell_post_order_ == 0) return;
  if (strandid_ < 0) ERROR_MESSAGE("Strand-id is not initialized.\n");
  MASTER_MESSAGE("Exporting cell grid: " + filename + "\n");
  if (MYRANK == MASTER_NODE) avocado::MakeDirectory(filename);

  const FileType filetype = FileType::GRIDONLY;
  const ZoneType zonetype =
      (dimension_ == 3) ? ZoneType::FETETS : ZoneType::FETRIS;
  const int csize = dimension_ + 1;
  std::vector<std::string> names;
  for (int idim = 0; idim < dimension_; idim++)
    names.push_back(DIRECTION(idim));

  std::ofstream file;
  MPI_File mpifile = nullptr;
  int offset = 0;
  std::vector<int> offset_gather(NDOMAIN);

  std::vector<int> header;
  WriteTecplotVersion(header);
  WriteIntData(header, 1);  // default
  WriteIntData(header,
               static_cast<int>(filetype));  // file type
  const std::string title = "Tecplot file written by " DENEB_VERSION;
  WriteStringData(header, title);                           // title
  WriteIntData(header, dimension_);                         // num var
  for (auto&& name : names) WriteStringData(header, name);  // var names
  WriteReal32Data(header, 299.0);                           // zone marker
  const std::string cell_name = "cell";
  WriteStringData(header, cell_name);                // zone name
  WriteIntData(header, -1);                          // parent zone
  WriteIntData(header, -1);                          // strand ID
  WriteReal64Data(header, 0.0);                      // solution time
  WriteIntData(header, -1);                          // default = -1
  WriteIntData(header, static_cast<int>(zonetype));  // zone type
  WriteIntData(header, 0);                           // specify var location
  // Are raw local 1-to-1 face neighbors supplied?
  WriteIntData(header, 0);
  // number of miscellaneous user-defined face neighbor connections
  WriteIntData(header, 0);
  const int num_total_nodes =
      AVOCADO_MPI->Reduce(cell_.num_total_nodes_, avocado::MPI::Op::SUM);
  const int num_total_cells =
      AVOCADO_MPI->Reduce(cell_.num_total_cells_, avocado::MPI::Op::SUM);
  WriteIntData(header, num_total_nodes);  // number of nodes
  WriteIntData(header, num_total_cells);  // number of cells
  WriteIntData(header, 0);                // ICellDim
  WriteIntData(header, 0);                // JCellDim
  WriteIntData(header, 0);                // KCellDim
  WriteIntData(header, 0);                // auxiliary name
  WriteReal32Data(header, 357.0);         // EOH marker

  if (NDOMAIN == 1) {
    file.open(filename, std::ios::binary);
    file.write((char*)&header[0], sizeof(int) * header.size());
  } else {
    MPI_File_open(MPI_COMM_WORLD, filename.c_str(),
                  MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
    MPI_File_set_view(mpifile, offset, MPI_INT, MPI_INT, "native",
                      MPI_INFO_NULL);
    if (MYRANK == MASTER_NODE)
      MPI_File_write(mpifile, &header[0], header.size(), MPI_INT,
                     MPI_STATUS_IGNORE);
    offset += static_cast<int>(header.size()) * sizeof(int);
  }

  std::vector<int> zone;
  WriteReal32Data(zone, 299.0);  // zone marker
  for (int idim = 0; idim < dimension_; idim++)
    WriteIntData(zone, 1);  // variable data format
  WriteIntData(zone, 0);    // has passive variables
  WriteIntData(zone, 0);    // has variable sharing
  // zeros based zone number to share connectivity list with (-1=no sharing)
  WriteIntData(zone, -1);

  for (int idim = 0; idim < dimension_; idim++) {
    const double min_var = *std::min_element(cell_.coords_[idim].begin(),
                                             cell_.coords_[idim].end());
    const double max_var = *std::max_element(cell_.coords_[idim].begin(),
                                             cell_.coords_[idim].end());
    WriteReal64Data(zone, AVOCADO_MPI->Reduce(min_var, avocado::MPI::Op::MIN));
    WriteReal64Data(zone, AVOCADO_MPI->Reduce(max_var, avocado::MPI::Op::MAX));
  }

  if (NDOMAIN == 1) {
    file.write((char*)&zone[0], sizeof(int) * static_cast<int>(zone.size()));
    for (int idim = 0; idim < dimension_; idim++)
      file.write((char*)&cell_.coords_[idim][0],
                 sizeof(float) * cell_.num_total_nodes_);
    file.write((char*)&cell_.connectivity_[0],
               sizeof(int) * csize * cell_.num_total_cells_);
    file.close();
  } else {
    if (MYRANK == MASTER_NODE)
      MPI_File_write(mpifile, &zone[0], static_cast<int>(zone.size()), MPI_INT,
                     MPI_STATUS_IGNORE);
    offset += static_cast<int>(zone.size()) * sizeof(int);

    const int data_size = cell_.num_total_nodes_ * sizeof(float);
    const int total_data_size =
        AVOCADO_MPI->Reduce(data_size, avocado::MPI::Op::SUM);
    const int offset_connect = offset + total_data_size * dimension_;
    MPI_Allgather(&data_size, 1, MPI_INT, &offset_gather[0], 1, MPI_INT,
                  MPI_COMM_WORLD);
    for (int idomain = 0; idomain < MYRANK; idomain++)
      offset += offset_gather[idomain];

    for (int idim = 0; idim < dimension_; idim++) {
      MPI_File_set_view(mpifile, offset, MPI_INT, MPI_INT, "native",
                        MPI_INFO_NULL);
      MPI_File_write(mpifile, &cell_.coords_[idim][0], cell_.num_total_nodes_,
                     MPI_FLOAT, MPI_STATUS_IGNORE);
      offset += total_data_size;
    }
    offset = offset_connect;

    const int connect_size = csize * cell_.num_total_cells_ * sizeof(int);
    MPI_Allgather(&connect_size, 1, MPI_INT, &offset_gather[0], 1, MPI_INT,
                  MPI_COMM_WORLD);
    for (int idomain = 0; idomain < MYRANK; idomain++)
      offset += offset_gather[idomain];

    MPI_File_set_view(mpifile, offset, MPI_INT, MPI_INT, "native",
                      MPI_INFO_NULL);
    MPI_File_write(mpifile, &cell_.connectivity_[0],
                   csize * cell_.num_total_cells_, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_close(&mpifile);
  }
}
void Contour::FaceSolution(const std::string& filename, const double* solution,
                           const double time) {
  if (face_post_order_ == 0) return;
  if (strandid_ < 0) ERROR_MESSAGE("Strand-id is not initialized.\n");
  MASTER_MESSAGE("Exporting face solution: " + filename + "\n");
  if (MYRANK == MASTER_NODE) avocado::MakeDirectory(filename);

  const FileType filetype = FileType::SOLUTION;
  const ZoneType zonetype =
      (dimension_ == 3) ? ZoneType::FETRIS : ZoneType::FELINE;
  const std::vector<std::string>& names = face_variables_;
  const int num_variables = static_cast<int>(face_variables_.size());
  const int& num_states = DENEB_EQUATION->GetNumStates();
  const int& num_bases = DENEB_DATA->GetNumBases();
  const int sb = num_states * num_bases;

  std::ofstream file;
  MPI_File mpifile = nullptr;
  int offset = 0;
  std::vector<int> offset_gather(NDOMAIN);

  std::vector<int> header;
  WriteTecplotVersion(header);
  WriteIntData(header, 1);  // default
  WriteIntData(header,
               static_cast<int>(filetype));  // file type
  const std::string title = "Tecplot file written by " DENEB_VERSION;
  WriteStringData(header, title);                           // title
  WriteIntData(header, num_variables);                      // num var
  for (auto&& name : names) WriteStringData(header, name);  // var names
  for (int i = 0; i < num_faces_; i++) {
    WriteReal32Data(header, 299.0);                    // zone marker
    WriteStringData(header, face_names_[i]);           // zone name
    WriteIntData(header, -1);                          // parent zone
    WriteIntData(header, strandid_++);                 // strand ID
    WriteReal64Data(header, time);                     // solution time
    WriteIntData(header, -1);                          // default = -1
    WriteIntData(header, static_cast<int>(zonetype));  // zone type
    WriteIntData(header, 0);                           // specify var location
    // Are raw local 1-to-1 face neighbors supplied?
    WriteIntData(header, 0);
    // number of miscellaneous user-defined face neighbor connections
    WriteIntData(header, 0);
    const int num_total_nodes =
        AVOCADO_MPI->Reduce(faces_[i].num_total_nodes_, avocado::MPI::Op::SUM);
    const int num_total_cells =
        AVOCADO_MPI->Reduce(faces_[i].num_total_cells_, avocado::MPI::Op::SUM);
    WriteIntData(header, num_total_nodes);  // number of nodes
    WriteIntData(header, num_total_cells);  // number of cells
    WriteIntData(header, 0);                // ICellDim
    WriteIntData(header, 0);                // JCellDim
    WriteIntData(header, 0);                // KCellDim
    WriteIntData(header, 0);                // auxiliary name
  }
  WriteReal32Data(header, 357.0);  // EOH marker

  if (NDOMAIN == 1) {
    file.open(filename, std::ios::binary);
    file.write((char*)&header[0], sizeof(int) * header.size());
  } else {
    MPI_File_open(MPI_COMM_WORLD, filename.c_str(),
                  MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
    MPI_File_set_view(mpifile, offset, MPI_INT, MPI_INT, "native",
                      MPI_INFO_NULL);
    if (MYRANK == MASTER_NODE)
      MPI_File_write(mpifile, &header[0], header.size(), MPI_INT,
                     MPI_STATUS_IGNORE);
    offset += static_cast<int>(header.size()) * sizeof(int);
  }

  for (int i = 0; i < num_faces_; i++) {
    std::vector<int> zone;
    WriteReal32Data(zone, 299.0);  // zone marker
    for (int ivar = 0; ivar < num_variables; ivar++)
      WriteIntData(zone, 1);  // variable data format
    WriteIntData(zone, 0);    // has passive variables
    WriteIntData(zone, 0);    // has variable sharing
    // zeros based zone number to share connectivity list with (-1=no sharing)
    WriteIntData(zone, -1);

    std::vector<std::vector<float>> variables(num_variables);
    for (int ivar = 0; ivar < num_variables; ivar++)
      variables[ivar].reserve(faces_[i].num_total_nodes_);

    for (int icell = 0; icell < faces_[i].num_target_cells_; icell++) {
      const int& cell = faces_[i].target_cells_[icell];
      const int& num_points = faces_[i].num_nodes_[icell];

      gemmABT(1.0, &faces_[i].basis_[icell][0], &solution[cell * sb], 0.0,
              &sol_[0], num_points, num_bases, num_states);
      avocado::Kernel1::f54(&faces_[i].basis_grad_[icell][0],
                            &solution[cell * sb], &sol_grad_[0], dimension_,
                            num_points, num_bases, num_states, 1.0, 0.0);

      DENEB_EQUATION->GetFacePostSolution(num_points, sol_, sol_grad_,
                                          faces_[i].normal_[icell], sol_post_);
      for (int ipoint = 0; ipoint < num_points; ipoint++)
        for (int ivar = 0; ivar < num_variables; ivar++)
          variables[ivar].push_back(sol_post_[ipoint * num_variables + ivar]);
    }

    for (int ivar = 0; ivar < num_variables; ivar++) {
      double min_var = 1.0e+100;
      double max_var = -1.0e+100;
      if (faces_[i].num_target_cells_ != 0) {
        min_var =
            *std::min_element(variables[ivar].begin(), variables[ivar].end());
        max_var =
            *std::max_element(variables[ivar].begin(), variables[ivar].end());
      }
      WriteReal64Data(zone,
                      AVOCADO_MPI->Reduce(min_var, avocado::MPI::Op::MIN));
      WriteReal64Data(zone,
                      AVOCADO_MPI->Reduce(max_var, avocado::MPI::Op::MAX));
    }

    if (NDOMAIN == 1) {
      file.write((char*)&zone[0], sizeof(int) * static_cast<int>(zone.size()));
      for (int ivar = 0; ivar < num_variables; ivar++)
        file.write((char*)&variables[ivar][0],
                   sizeof(float) * faces_[i].num_total_nodes_);
    } else {
      if (MYRANK == MASTER_NODE)
        MPI_File_write(mpifile, &zone[0], static_cast<int>(zone.size()),
                       MPI_INT, MPI_STATUS_IGNORE);
      offset += static_cast<int>(zone.size()) * sizeof(int);

      const int data_size = faces_[i].num_total_nodes_ * sizeof(float);
      const int total_data_size =
          AVOCADO_MPI->Reduce(data_size, avocado::MPI::Op::SUM);
      const int offset_connect = offset + total_data_size * num_variables;
      MPI_Allgather(&data_size, 1, MPI_INT, &offset_gather[0], 1, MPI_INT,
                    MPI_COMM_WORLD);
      for (int idomain = 0; idomain < MYRANK; idomain++)
        offset += offset_gather[idomain];

      for (int ivar = 0; ivar < num_variables; ivar++) {
        MPI_File_set_view(mpifile, offset, MPI_INT, MPI_INT, "native",
                          MPI_INFO_NULL);
        if (faces_[i].num_target_cells_ != 0)
          MPI_File_write(mpifile, &variables[ivar][0],
                         faces_[i].num_total_nodes_, MPI_FLOAT,
                         MPI_STATUS_IGNORE);
        offset += total_data_size;
      }
      offset = offset_connect;
      MPI_File_set_view(mpifile, offset, MPI_INT, MPI_INT, "native",
                        MPI_INFO_NULL);
    }
  }
  if (NDOMAIN == 1) {
    file.close();
  } else {
    MPI_File_close(&mpifile);
  }
}
void Contour::FaceGrid(const std::string& filename) {
  if (face_post_order_ == 0) return;
  if (strandid_ < 0) ERROR_MESSAGE("Strand-id is not initialized.\n");
  MASTER_MESSAGE("Exporting face grid: " + filename + "\n");
  if (MYRANK == MASTER_NODE) avocado::MakeDirectory(filename);

  const FileType filetype = FileType::GRIDONLY;
  const ZoneType zonetype =
      (dimension_ == 3) ? ZoneType::FETRIS : ZoneType::FELINE;
  const int csize = dimension_;
  std::vector<std::string> names;
  for (int idim = 0; idim < dimension_; idim++)
    names.push_back(DIRECTION(idim));

  std::ofstream file;
  MPI_File mpifile = nullptr;
  int offset = 0;
  std::vector<int> offset_gather(NDOMAIN);

  std::vector<int> header;
  WriteTecplotVersion(header);
  WriteIntData(header, 1);  // default
  WriteIntData(header,
               static_cast<int>(filetype));  // file type
  const std::string title = "Tecplot file written by " DENEB_VERSION;
  WriteStringData(header, title);                           // title
  WriteIntData(header, dimension_);                         // num var
  for (auto&& name : names) WriteStringData(header, name);  // var names
  for (int i = 0; i < num_faces_; i++) {
    WriteReal32Data(header, 299.0);                    // zone marker
    WriteStringData(header, face_names_[i]);           // zone name
    WriteIntData(header, -1);                          // parent zone
    WriteIntData(header, -1);                          // strand ID
    WriteReal64Data(header, 0.0);                      // solution time
    WriteIntData(header, -1);                          // default = -1
    WriteIntData(header, static_cast<int>(zonetype));  // zone type
    WriteIntData(header, 0);                           // specify var location
    // Are raw local 1-to-1 face neighbors supplied?
    WriteIntData(header, 0);
    // number of miscellaneous user-defined face neighbor connections
    WriteIntData(header, 0);
    const int num_total_nodes =
        AVOCADO_MPI->Reduce(faces_[i].num_total_nodes_, avocado::MPI::Op::SUM);
    const int num_total_cells =
        AVOCADO_MPI->Reduce(faces_[i].num_total_cells_, avocado::MPI::Op::SUM);
    WriteIntData(header, num_total_nodes);  // number of nodes
    WriteIntData(header, num_total_cells);  // number of cells
    WriteIntData(header, 0);                // ICellDim
    WriteIntData(header, 0);                // JCellDim
    WriteIntData(header, 0);                // KCellDim
    WriteIntData(header, 0);                // auxiliary name
  }
  WriteReal32Data(header, 357.0);  // EOH marker

  if (NDOMAIN == 1) {
    file.open(filename, std::ios::binary);
    file.write((char*)&header[0], sizeof(int) * header.size());
  } else {
    MPI_File_open(MPI_COMM_WORLD, filename.c_str(),
                  MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
    MPI_File_set_view(mpifile, offset, MPI_INT, MPI_INT, "native",
                      MPI_INFO_NULL);
    if (MYRANK == MASTER_NODE)
      MPI_File_write(mpifile, &header[0], header.size(), MPI_INT,
                     MPI_STATUS_IGNORE);
    offset += static_cast<int>(header.size()) * sizeof(int);
  }

  for (int i = 0; i < num_faces_; i++) {
    std::vector<int> zone;
    WriteReal32Data(zone, 299.0);  // zone marker
    for (int idim = 0; idim < dimension_; idim++)
      WriteIntData(zone, 1);  // variable data format
    WriteIntData(zone, 0);    // has passive variables
    WriteIntData(zone, 0);    // has variable sharing
    // zeros based zone number to share connectivity list with (-1=no sharing)
    WriteIntData(zone, -1);

    for (int idim = 0; idim < dimension_; idim++) {
      double min_var = 1.0e+100;
      double max_var = -1.0e+100;
      if (faces_[i].num_target_cells_ != 0) {
        min_var = *std::min_element(faces_[i].coords_[idim].begin(),
                                    faces_[i].coords_[idim].end());
        max_var = *std::max_element(faces_[i].coords_[idim].begin(),
                                    faces_[i].coords_[idim].end());
      }
      WriteReal64Data(zone,
                      AVOCADO_MPI->Reduce(min_var, avocado::MPI::Op::MIN));
      WriteReal64Data(zone,
                      AVOCADO_MPI->Reduce(max_var, avocado::MPI::Op::MAX));
    }

    if (NDOMAIN == 1) {
      file.write((char*)&zone[0], sizeof(int) * static_cast<int>(zone.size()));
      for (int idim = 0; idim < dimension_; idim++)
        file.write((char*)&faces_[i].coords_[idim][0],
                   sizeof(float) * faces_[i].num_total_nodes_);
      file.write((char*)&faces_[i].connectivity_[0],
                 sizeof(int) * csize * faces_[i].num_total_cells_);
    } else {
      if (MYRANK == MASTER_NODE)
        MPI_File_write(mpifile, &zone[0], static_cast<int>(zone.size()),
                       MPI_INT, MPI_STATUS_IGNORE);
      offset += static_cast<int>(zone.size()) * sizeof(int);

      const int data_size = faces_[i].num_total_nodes_ * sizeof(float);
      const int total_data_size =
          AVOCADO_MPI->Reduce(data_size, avocado::MPI::Op::SUM);
      const int offset_connect = offset + total_data_size * dimension_;
      MPI_Allgather(&data_size, 1, MPI_INT, &offset_gather[0], 1, MPI_INT,
                    MPI_COMM_WORLD);
      for (int idomain = 0; idomain < MYRANK; idomain++)
        offset += offset_gather[idomain];

      for (int idim = 0; idim < dimension_; idim++) {
        MPI_File_set_view(mpifile, offset, MPI_INT, MPI_INT, "native",
                          MPI_INFO_NULL);
        if (faces_[i].num_target_cells_ != 0)
          MPI_File_write(mpifile, &faces_[i].coords_[idim][0],
                         faces_[i].num_total_nodes_, MPI_FLOAT,
                         MPI_STATUS_IGNORE);
        offset += total_data_size;
      }
      offset = offset_connect;

      const int connect_size = csize * faces_[i].num_total_cells_ * sizeof(int);
      const int total_connect_size =
          AVOCADO_MPI->Reduce(connect_size, avocado::MPI::Op::SUM);
      const int offset_next_zone = offset + total_connect_size;

      MPI_Allgather(&connect_size, 1, MPI_INT, &offset_gather[0], 1, MPI_INT,
                    MPI_COMM_WORLD);
      for (int idomain = 0; idomain < MYRANK; idomain++)
        offset += offset_gather[idomain];

      MPI_File_set_view(mpifile, offset, MPI_INT, MPI_INT, "native",
                        MPI_INFO_NULL);
      if (faces_[i].num_target_cells_ != 0)
        MPI_File_write(mpifile, &faces_[i].connectivity_[0],
                       csize * faces_[i].num_total_cells_, MPI_INT,
                       MPI_STATUS_IGNORE);
      offset = offset_next_zone;
      MPI_File_set_view(mpifile, offset, MPI_INT, MPI_INT, "native",
                        MPI_INFO_NULL);
    }
  }
  if (NDOMAIN == 1) {
    file.close();
  } else {
    MPI_File_close(&mpifile);
  }
}
void Contour::BuildCellData(void) {
  const int& post_order = cell_post_order_;
  const int& num_cells = DENEB_DATA->GetNumCells();
  const std::vector<ElemType>& cell_elemtype = DENEB_DATA->GetCellElemtype();

  int num_total_nodes = 0;
  int num_total_cells = 0;
  std::vector<int> num_nodes(num_cells);
  std::vector<std::vector<double>> coords_temp(num_cells);
  std::vector<std::vector<int>> connectivity_temp(num_cells);
  for (int icell = 0; icell < num_cells; icell++) {
    std::vector<double> coords_local;
    DENEB_ELEMENT->GetSubcell(cell_elemtype[icell], post_order, coords_local,
                              connectivity_temp[icell]);
    for (auto&& connect : connectivity_temp[icell]) connect += num_total_nodes;
    DENEB_DATA->GetCellRefToPhyCoords(icell, coords_local, coords_temp[icell]);

    num_nodes[icell] = static_cast<int>(coords_temp[icell].size()) / dimension_;
    num_total_nodes += num_nodes[icell];
    num_total_cells +=
        static_cast<int>(connectivity_temp[icell].size()) / (dimension_ + 1);
  }

  std::vector<std::vector<double>> basis(num_cells);
  std::vector<std::vector<double>> basis_grad(num_cells);
  for (int icell = 0; icell < num_cells; icell++) {
    DENEB_DATA->GetCellBasisValues(icell, coords_temp[icell], basis[icell]);
    DENEB_DATA->GetCellBasisGradValues(icell, coords_temp[icell],
                                       basis_grad[icell]);
  }

  std::vector<std::vector<float>> coords(dimension_);
  for (int idim = 0; idim < dimension_; idim++) {
    coords[idim].reserve(num_total_nodes);
    for (int icell = 0; icell < num_cells; icell++)
      for (int inode = 0; inode < num_nodes[icell]; inode++)
        coords[idim].push_back(coords_temp[icell][inode * dimension_ + idim]);
  }
  coords_temp.clear();

  std::vector<int> connectivity;
  connectivity.reserve(num_total_cells * (dimension_ + 1));
  for (int icell = 0; icell < num_cells; icell++)
    for (auto&& i : connectivity_temp[icell]) connectivity.push_back(i);
  connectivity_temp.clear();

  if (NDOMAIN > 1) {
    std::vector<int> num_nodes_gather(NDOMAIN, 0);
    MPI_Allgather(&num_total_nodes, 1, MPI_INT, &num_nodes_gather[0], 1,
                  MPI_INT, MPI_COMM_WORLD);
    int offset = 0;
    for (int idomain = 0; idomain < MYRANK; idomain++)
      offset += num_nodes_gather[idomain];
    for (auto&& con : connectivity) con += offset;
  }

  cell_.num_total_nodes_ = std::move(num_total_nodes);
  cell_.num_total_cells_ = std::move(num_total_cells);
  cell_.num_nodes_ = std::move(num_nodes);
  cell_.connectivity_ = std::move(connectivity);
  cell_.coords_ = std::move(coords);
  cell_.basis_ = std::move(basis);
  cell_.basis_grad_ = std::move(basis_grad);
}
void Contour::BuildFaceData(void) {
  const int& post_order = face_post_order_;
  {
    std::vector<int> all_bdry_tag;
    {
      std::unordered_set<int> temp(DENEB_DATA->GetBdryTag().begin(),
                                   DENEB_DATA->GetBdryTag().end());
      std::vector<int> bdry_tag(temp.begin(), temp.end());
      temp.clear();

      std::vector<std::vector<int>> send_data(NDOMAIN, bdry_tag);
      bdry_tag.clear();
      std::vector<std::vector<int>> recv_data;
      AVOCADO_MPI->CommunicateData(send_data, recv_data);
      send_data.clear();

      for (auto&& data : recv_data) temp.insert(data.begin(), data.end());
      recv_data.clear();
      all_bdry_tag = std::vector<int>(temp.begin(), temp.end());
    }

    std::vector<int> all_peribdry_tag;
    {
      std::unordered_set<int> temp(DENEB_DATA->GetPeribdryTag().begin(),
                                   DENEB_DATA->GetPeribdryTag().end());
      std::vector<int> peribdry_tag(temp.begin(), temp.end());
      temp.clear();

      std::vector<std::vector<int>> send_data(NDOMAIN, peribdry_tag);
      peribdry_tag.clear();
      std::vector<std::vector<int>> recv_data;
      AVOCADO_MPI->CommunicateData(send_data, recv_data);
      send_data.clear();

      for (auto&& data : recv_data) temp.insert(data.begin(), data.end());
      recv_data.clear();
      all_peribdry_tag = std::vector<int>(temp.begin(), temp.end());
    }

    face_tags_.clear();
    for (auto&& tag : all_bdry_tag) face_tags_.push_back(tag);
    for (auto&& tag : all_peribdry_tag) face_tags_.push_back(tag);
  }

  auto& config = AVOCADO_CONFIG;
  num_faces_ = static_cast<int>(face_tags_.size());
  for (auto&& tag : face_tags_) {
    const std::string& type = config->GetConfigValue(BDRY_TYPE(tag));
    const std::string& name = config->GetConfigValue(BDRY_NAME(tag));
    face_names_.push_back(type + ", " + name + " (" + std::to_string(tag) +
                          ")");
  }

  faces_.clear();
  faces_.resize(num_faces_);

  const auto& cell_element = DENEB_DATA->GetCellElement();
  const auto& cell_jacobians = DENEB_DATA->GetCellJacobians();
  const auto& cell_order = DENEB_DATA->GetCellOrder();
  const auto& cell_subnode_ptr = DENEB_DATA->GetCellSubnodePtr();
  const auto& cell_subnode_ind = DENEB_DATA->GetCellSubnodeInd();

  for (int i = 0; i < num_faces_; i++) {
    const int& tag = face_tags_[i];
    const std::string& type = config->GetConfigValue(BDRY_TYPE(tag));
    const bool is_periodic = (!type.compare(PERIODIC)) ? true : false;

    const int& num_bdries = (is_periodic) ? DENEB_DATA->GetNumTotalPeribdries()
                                          : DENEB_DATA->GetNumBdries();
    const auto& bdry_owner_cell =
        (is_periodic) ? DENEB_DATA->GetPeribdryOwnerCell()
                      : DENEB_DATA->GetBdryOwnerCell();
    const auto& bdry_owner_type =
        (is_periodic) ? DENEB_DATA->GetPeribdryOwnerType()
                      : DENEB_DATA->GetBdryOwnerType();
    const auto& bdry_tag = (is_periodic) ? DENEB_DATA->GetPeribdryTag()
                                         : DENEB_DATA->GetBdryTag();

    int num_target_cells;
    std::vector<int> target_cells;
    std::vector<int> target_types;
    for (int ibdry = 0; ibdry < num_bdries; ibdry++) {
      if (bdry_tag[ibdry] != tag) continue;
      target_cells.push_back(bdry_owner_cell[ibdry]);
      target_types.push_back(bdry_owner_type[ibdry]);
    }
    num_target_cells = static_cast<int>(target_cells.size());

    int num_total_nodes = 0;
    int num_total_cells = 0;
    std::vector<int> num_nodes(num_target_cells);
    std::vector<std::vector<double>> coords_temp(num_target_cells);
    std::vector<std::vector<int>> connectivity_temp(num_target_cells);
    std::vector<std::vector<double>> normal(num_target_cells);
    for (int icell = 0; icell < num_target_cells; icell++) {
      const int& cell = target_cells[icell];
      const int& type = target_types[icell];
      const Element* element = cell_element[cell];
      const ElemType face_elemtype = element->GetFaceElemtype(type);

      std::vector<double> ref_face_coords;
      DENEB_ELEMENT->GetSubcell(face_elemtype, post_order, ref_face_coords,
                                connectivity_temp[icell]);
      for (auto&& connect : connectivity_temp[icell])
        connect += num_total_nodes;

      const int num_coords =
          static_cast<int>(ref_face_coords.size()) / (dimension_ - 1);
      std::vector<double> ref_coords(num_coords * dimension_);
      element->TransformToPhyCoords(type, num_coords, &ref_face_coords[0],
                                    &ref_coords[0]);
      DENEB_DATA->GetCellRefToPhyCoords(cell, ref_coords, coords_temp[icell]);

      num_nodes[icell] =
          static_cast<int>(coords_temp[icell].size()) / dimension_;
      num_total_nodes += num_nodes[icell];
      num_total_cells +=
          static_cast<int>(connectivity_temp[icell].size()) / dimension_;

      const double* ref_normal = element->GetFacetypeNormal(type);
      std::vector<double> cofactor(num_coords * dimension_ * dimension_);
      cell_jacobians[cell]->CalJacobianCofMat(num_coords, &ref_coords[0],
                                              &cofactor[0]);
      normal[icell].resize(num_coords * dimension_);
      for (int ipoint = 0; ipoint < num_coords; ipoint++) {
        gemvAx(1.0, &cofactor[ipoint * dimension_ * dimension_], dimension_,
               ref_normal, 1, 0.0, &normal[icell][ipoint * dimension_], 1,
               dimension_, dimension_);
        avocado::VecNormal(dimension_, &normal[icell][ipoint * dimension_]);
      }
    }

    std::vector<std::vector<double>> basis(num_target_cells);
    std::vector<std::vector<double>> basis_grad(num_target_cells);
    for (int icell = 0; icell < num_target_cells; icell++) {
      const int& cell = target_cells[icell];
      DENEB_DATA->GetCellBasisValues(cell, coords_temp[icell], basis[icell]);
      DENEB_DATA->GetCellBasisGradValues(cell, coords_temp[icell],
                                         basis_grad[icell]);
    }

    std::vector<std::vector<float>> coords(dimension_);
    for (int idim = 0; idim < dimension_; idim++) {
      coords[idim].reserve(num_total_nodes);
      for (int icell = 0; icell < num_target_cells; icell++)
        for (int inode = 0; inode < num_nodes[icell]; inode++)
          coords[idim].push_back(coords_temp[icell][inode * dimension_ + idim]);
    }
    coords_temp.clear();

    std::vector<int> connectivity;
    connectivity.reserve(num_total_cells * dimension_);
    for (int icell = 0; icell < num_target_cells; icell++)
      for (auto&& i : connectivity_temp[icell]) connectivity.push_back(i);
    connectivity_temp.clear();

    if (NDOMAIN > 1) {
      std::vector<int> num_nodes_gather(NDOMAIN, 0);
      MPI_Allgather(&num_total_nodes, 1, MPI_INT, &num_nodes_gather[0], 1,
                    MPI_INT, MPI_COMM_WORLD);
      int offset = 0;
      for (int idomain = 0; idomain < MYRANK; idomain++)
        offset += num_nodes_gather[idomain];
      for (auto&& con : connectivity) con += offset;
    }

    faces_[i].num_total_nodes_ = std::move(num_total_nodes);
    faces_[i].num_total_cells_ = std::move(num_total_cells);
    faces_[i].num_nodes_ = std::move(num_nodes);
    faces_[i].connectivity_ = std::move(connectivity);
    faces_[i].coords_ = std::move(coords);
    faces_[i].basis_ = std::move(basis);
    faces_[i].basis_grad_ = std::move(basis_grad);

    faces_[i].num_target_cells_ = std::move(num_target_cells);
    faces_[i].target_cells_ = std::move(target_cells);
    faces_[i].normal_ = std::move(normal);
  }
}
void Contour::WriteTecplotVersion(std::vector<int>& buffer) const {
  int ver1 = 0;
  int ver2 = 0;
  const char* version_ = "#!TDV112";
  ver1 = ver1 | (*(char*)&version_[3] << 24);
  ver1 = ver1 | (*(char*)&version_[2] << 16);
  ver1 = ver1 | (*(char*)&version_[1] << 8);
  ver1 = ver1 | *(char*)&version_[0];
  ver2 = ver2 | (*(char*)&version_[7] << 24);
  ver2 = ver2 | (*(char*)&version_[6] << 16);
  ver2 = ver2 | (*(char*)&version_[5] << 8);
  ver2 = ver2 | *(char*)&version_[4];
  buffer.push_back(ver1);
  buffer.push_back(ver2);
}
void Contour::WriteIntData(std::vector<int>& buffer, const int data) const {
  buffer.push_back(data);
}
void Contour::WriteReal32Data(std::vector<int>& buffer,
                              const float data) const {
  buffer.push_back(*(int*)&data);
}
void Contour::WriteReal64Data(std::vector<int>& buffer,
                              const double data) const {
  __int64 left = (*(__int64*)&data) >> 32;
  const __int64 mask = 0xffffffff00000000;
  __int64 right = ((*(__int64*)&data) & mask) >> 32;
  buffer.push_back((int)right);
  buffer.push_back((int)left);
}
void Contour::WriteStringData(std::vector<int>& buffer,
                              const std::string& data) const {
  const char* str = data.c_str();
  int ch;
  for (int i = 0; i < data.length(); i++) {
    ch = (int)str[i];
    buffer.push_back(ch);
  }
  buffer.push_back(0);
}
}  // namespace deneb