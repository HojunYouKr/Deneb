#include "avocado_mpi.h"

#include <cstdio>

namespace avocado {
std::shared_ptr<MPI> AVOCADO_MPI_NAME = nullptr;

MPI::MPI(const std::string& codename, const std::string& filename)
    : codename_(codename), filename_(filename), myrank_(0), num_domains_(1) {
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank_);
  MPI_Comm_size(MPI_COMM_WORLD, &num_domains_);
}

void MPI::StopProgram(void) const {
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
}
void MPI::Syncro(void) const {
  MPI_Barrier(MPI_COMM_WORLD);
}

void MPI::PrintMessage(const bool only_master,
                       const std::string& message) const {
  if (!only_master)
    printf("(%4d)>> ", myrank_);
  else if (MASTER_NODE == myrank_)
    printf(">> ");
  if ((only_master && (MASTER_NODE == myrank_)) || !only_master)
    printf("%s", message.c_str());
  fflush(stdout);
}

void MPI::ErrorMessage(const std::string& message, const std::string& file,
                       const int line) const {
  printf(
      "(%4d)>> Program interruption"
      "\n\tFile: %s"
      "\n\tLine: %d"
      "\n\tMessage: %s",
      myrank_, file.c_str(), line, message.c_str());
  fflush(stdout);
  StopProgram();
}

Communicate::Communicate(const int size,
                         const std::vector<std::vector<int>>& send_data_list,
                         const std::vector<std::vector<int>>& recv_data_list)
    : size_(size),
      send_data_list_(send_data_list),
      recv_data_list_(recv_data_list) {
  if (NDOMAIN > 1) {
    recv_req_.resize(NDOMAIN);
    recv_sta_.resize(NDOMAIN);
    send_req_.resize(NDOMAIN);
    send_sta_.resize(NDOMAIN);

    num_send_data_.resize(NDOMAIN);
    num_recv_data_.resize(NDOMAIN);
    for (int idomain = 0; idomain < NDOMAIN; idomain++) {
      num_send_data_[idomain] =
          size_ * static_cast<int>(send_data_list_[idomain].size());
      num_recv_data_[idomain] =
          size_ * static_cast<int>(recv_data_list_[idomain].size());
    }

    send_data_.resize(NDOMAIN);
    recv_data_.resize(NDOMAIN);
    for (int idomain = 0; idomain < NDOMAIN; idomain++) {
      if (num_send_data_[idomain] == 0)
        send_data_[idomain].resize(1);
      else
        send_data_[idomain].resize(num_send_data_[idomain]);
      if (num_recv_data_[idomain] == 0)
        recv_data_[idomain].resize(1);
      else
        recv_data_[idomain].resize(num_recv_data_[idomain]);
    }
  }
}
void Communicate::CommunicateBegin(const double* send_data) {
  if (NDOMAIN > 1) {
    for (int idomain = 0; idomain < NDOMAIN; idomain++) {
      int loc = 0;
      for (auto&& cell : send_data_list_[idomain])
        for (int i = 0; i < size_; i++)
          send_data_[idomain][loc++] = send_data[cell * size_ + i];
    }
    nrecv_domain_ = 0;
    nsend_domain_ = 0;
    for (int idomain = 0; idomain < NDOMAIN; idomain++) {
      if (num_recv_data_[idomain] > 0)
        MPI_Irecv(&recv_data_[idomain][0], num_recv_data_[idomain], MPI_DOUBLE,
                  idomain, 22, MPI_COMM_WORLD, &recv_req_[nrecv_domain_++]);
      if (num_send_data_[idomain] > 0)
        MPI_Isend(&send_data_[idomain][0], num_send_data_[idomain], MPI_DOUBLE,
                  idomain, 22, MPI_COMM_WORLD, &send_req_[nsend_domain_++]);
    }
  }
}
void Communicate::CommunicateEnd(double* recv_data) {
  if (NDOMAIN > 1) {
    MPI_Waitall(nrecv_domain_, &recv_req_[0], &recv_sta_[0]);
    MPI_Waitall(nsend_domain_, &send_req_[0], &send_sta_[0]);

    for (int idomain = 0; idomain < NDOMAIN; idomain++) {
      int loc = 0;
      for (auto&& cell : recv_data_list_[idomain])
        for (int i = 0; i < size_; i++)
          recv_data[cell * size_ + i] = recv_data_[idomain][loc++];
    }
  }
}
}  // namespace avocado