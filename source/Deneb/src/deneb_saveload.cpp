#include "deneb_saveload.h"

#include <ctime>
#include <unordered_map>
#include <sstream>

#include "avocado.h"
#include "deneb_data.h"
#include "deneb_equation.h"
#include "deneb_header.h"
#include "deneb_timescheme.h"

#define USE_PARALLEL_SAVE

namespace deneb {
std::shared_ptr<SaveLoad> DENEB_SAVELOAD_NAME;
void SaveLoad::Save(const std::string& filename, const double* solution,
                    const SaveData& data) const {
  std::stringstream ss;
  ss << "Saving solution: " << filename << "\n\tIteration: " << data.iteration_
     << "\n\tTime: " << DENEB_TIMESCHEME->ConvertTime(data.time_)
     << "\n\tStrandid: " << data.strandid_ << "\n";
  MASTER_MESSAGE(ss.str());
  if (MYRANK == MASTER_NODE) avocado::MakeDirectory(filename);

  std::vector<char> buffer;
  WriteData(buffer, DENEB_VERSION);
  WriteData(buffer, std::string(__DATE__) + " " + __TIME__);
  WriteData(buffer, avocado::GetLocalTime());

  const std::string comment = "Nothing";
  WriteData(buffer, comment);

  const int ndomain = NDOMAIN;
  WriteData(buffer, ndomain);

  const int order = DENEB_DATA->GetOrder();
  const int num_bases = DENEB_DATA->GetNumBases();
  const int num_global_cells = DENEB_DATA->GetNumGlobalCells();
  WriteData(buffer, order);
  WriteData(buffer, num_bases);
  WriteData(buffer, num_global_cells);

  const int dimension = DENEB_EQUATION->GetDimension();
  const int num_states = DENEB_EQUATION->GetNumStates();
  WriteData(buffer, dimension);
  WriteData(buffer, num_states);

  const int time_resolution = DENEB_TIMESCHEME->GetTimeResolution();
  WriteData(buffer, time_resolution);

  WriteData(buffer, data.iteration_);
  WriteData(buffer, data.strandid_);
  WriteData(buffer, data.time_);

  const int sb = num_states * num_bases;
  const int num_cells = DENEB_DATA->GetNumCells();
  const std::vector<int>& cell_global_index = DENEB_DATA->GetCellGlobalIndex();

  if (NDOMAIN == 1) {
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open())
      ERROR_MESSAGE("Save file is not opened: " + filename + "\n");

    file.write(&buffer[0], buffer.size());
    for (int icell = 0; icell < num_cells; icell++) {
      file.write((char*)&cell_global_index[icell], sizeof(int));
      file.write((char*)&solution[icell * sb], sizeof(double) * sb);
    }
    file.close();
  } else {
#ifdef USE_PARALLEL_SAVE
    MPI_File mpifile = nullptr;
    int offset = 0;
    std::vector<int> offset_gather(NDOMAIN);

    MPI_File_open(MPI_COMM_WORLD, filename.c_str(),
                  MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
    MPI_File_set_view(mpifile, offset, MPI_INT, MPI_INT, "native",
                      MPI_INFO_NULL);
    if (MYRANK == MASTER_NODE)
      MPI_File_write(mpifile, &buffer[0], buffer.size(), MPI_INT,
                     MPI_STATUS_IGNORE);
    offset += static_cast<int>(buffer.size()) * sizeof(char);

    const int data_size = num_cells * (sizeof(int) + sb * sizeof(double));
    MPI_Allgather(&data_size, 1, MPI_INT, &offset_gather[0], 1, MPI_INT,
                  MPI_COMM_WORLD);
    for (int idomain = 0; idomain < MYRANK; idomain++)
      offset += offset_gather[idomain];

    MPI_File_set_view(mpifile, offset, MPI_INT, MPI_INT, "native",
                      MPI_INFO_NULL);
    for (int icell = 0; icell < num_cells; icell++) {
      MPI_File_write(mpifile, &cell_global_index[icell], 1, MPI_INT,
                     MPI_STATUS_IGNORE);
      MPI_File_write(mpifile, &solution[icell * sb], sb, MPI_DOUBLE,
                     MPI_STATUS_IGNORE);
    }
    MPI_File_close(&mpifile);
#else
    std::vector<std::vector<int> > send_index(NDOMAIN);
    std::vector<std::vector<int> > recv_index;
    send_index[MASTER_NODE] = std::move(std::vector<int>(
        cell_global_index.begin(), cell_global_index.begin() + num_cells));
    AVOCADO_MPI->CommunicateData(send_index, recv_index);
    send_index.clear();

    std::vector<std::vector<double> > send_sol(NDOMAIN);
    std::vector<std::vector<double> > recv_sol(NDOMAIN);
    send_sol[MASTER_NODE] =
        std::move(std::vector<double>(solution, solution + num_cells * sb));
    AVOCADO_MPI->CommunicateData(send_sol, recv_sol);
    send_sol.clear();

    if (MYRANK == MASTER_NODE) {
      std::ofstream file(filename, std::ios::binary);
      if (!file.is_open())
        ERROR_MESSAGE("Save file is not opened: " + filename + "\n");

      file.write(&buffer[0], buffer.size());
      for (int i = 0; i < NDOMAIN; i++) {
        for (size_t icell = 0, len = recv_index[i].size(); icell < len;
             icell++) {
          file.write((char*)&recv_index[i][icell], sizeof(int));
          file.write((char*)&recv_sol[i][icell * sb], sizeof(double) * sb);
        }
      }
      file.close();
    }
#endif
  }
}
void SaveLoad::Load(const std::string& filename) {
  std::ifstream file(filename, std::ios::binary);
  if (!file.is_open())
    ERROR_MESSAGE("Load file is not opened: " + filename + "\n");

  std::string version, build_date, save_date;
  ReadData(file, version);
  ReadData(file, build_date);
  ReadData(file, save_date);

  std::string comment;
  ReadData(file, comment);

  int ndomain;
  ReadData(file, &ndomain);

  int order, num_bases, num_global_cells;
  ReadData(file, &order);
  ReadData(file, &num_bases);
  ReadData(file, &num_global_cells);

  Compare(file, "The number of global cells", num_global_cells,
          DENEB_DATA->GetNumGlobalCells());

  int dimension, num_states;
  ReadData(file, &dimension);
  ReadData(file, &num_states);

  Compare(file, "Dimension", dimension, DENEB_EQUATION->GetDimension());
  Compare(file, "The number of state variables", num_states,
          DENEB_EQUATION->GetNumStates());

  int time_resolution;
  ReadData(file, &time_resolution);

  unsigned __int64 time;
  ReadData(file, &data_.iteration_);
  ReadData(file, &data_.strandid_);
  ReadData(file, &time);
  {
    const int& current_time_resolution = DENEB_TIMESCHEME->GetTimeResolution();
    data_.time_ =
        time * std::pow(10, current_time_resolution - time_resolution);
  }

  std::stringstream ss;
  ss << "Loading solution: " << filename
     << "\n\tIteration: " << data_.iteration_
     << "\n\tTime: " << DENEB_TIMESCHEME->ConvertTime(data_.time_)
     << "\n\tStrandid: " << data_.strandid_ << "\n";
  MASTER_MESSAGE(ss.str());

  const int current_num_bases = DENEB_DATA->GetNumBases();
  const int num_cells = DENEB_DATA->GetNumCells();
  const std::vector<int>& cell_global_index = DENEB_DATA->GetCellGlobalIndex();
  std::unordered_map<int, int> global_to_local;
  for (int icell = 0; icell < num_cells; icell++)
    global_to_local[cell_global_index[icell]] = icell;

  solution_.resize(num_cells * num_states * current_num_bases);
  std::memset(&solution_[0], 0,
              num_cells * num_states * current_num_bases * sizeof(double));

  int global_index;
  const int min_num_bases = std::min(num_bases, current_num_bases);
  const int sb = num_states * num_bases;
  std::vector<double> buffer(sb);
  for (int icell = 0; icell < num_global_cells; icell++) {
    file.read((char*)&global_index, sizeof(int));
    auto&& iterator = global_to_local.find(global_index);
    if (iterator == global_to_local.end()) {
      file.seekg(sizeof(double) * sb, std::ios::cur);
      continue;
    }
    ReadData(file, &buffer[0], sb);
    for (int istate = 0; istate < num_states; istate++)
      for (int ibasis = 0; ibasis < min_num_bases; ibasis++)
        solution_[iterator->second * num_states * current_num_bases +
                  istate * current_num_bases + ibasis] =
            buffer[istate * num_bases + ibasis];
  }
  file.close();
}
void SaveLoad::Compare(std::ifstream& file, const std::string& name,
                       const int save, const int current) const {
  if (save != current) {
    file.close();
    ERROR_MESSAGE(name +
                  " does not match.\n\tSavefile: " + std::to_string(save) +
                  "\n\tCurrent: " + std::to_string(current) + "\n");
  }
}
}  // namespace deneb