#pragma once

#include <memory>
#include <string>
#include <vector>

#include <mpi.h>

#include "avocado_header.h"

#define AVOCADO_MPI_NAME mpi_global_ptr
#define AVOCADO_MPI avocado::AVOCADO_MPI_NAME
#define AVOCADO_MPI_INITIALIZE(codename, filename) \
  AVOCADO_MPI = std::make_shared<avocado::MPI>(codename, filename)
#define AVOCADO_MPI_FINALIZE() AVOCADO_MPI.reset()

#define MASTER_NODE 0
#define MYRANK AVOCADO_MPI->GetMyRank()
#define NDOMAIN AVOCADO_MPI->GetNumDomains()
#define SYNCRO() AVOCADO_MPI->Syncro()
#define MASTER_MESSAGE(message) AVOCADO_MPI->PrintMessage(true, (message))
#define MEMBER_MESSAGE(message) AVOCADO_MPI->PrintMessage(false, (message))
#define ERROR_MESSAGE(message) \
  AVOCADO_MPI->ErrorMessage(message, __FILE__, __LINE__)
#define ASSERT(expr)                                                        \
  (void)((!!(expr)) ||                                                      \
         (AVOCADO_MPI->ErrorMessage(                                        \
              "Invalid expression: " + std::string(#expr) + "\n", __FILE__, \
              __LINE__),                                                    \
          0))

namespace avocado {
class MPI {
 public:
  enum class Op : int { SUM = 0, MIN = 1, MAX = 2 };

 protected:
  std::string codename_;
  std::string filename_;
  int myrank_;
  int num_domains_;

 public:
  MPI(const std::string& codename, const std::string& filename);
  ~MPI(){};

  inline const std::string GetCodeName() const { return codename_; };
  inline int GetMyRank() const { return myrank_; };
  inline int GetNumDomains() const { return num_domains_; };

  void StopProgram(void) const;
  void Syncro(void) const;
  void PrintMessage(const bool only_master, const std::string& message) const;
  void ErrorMessage(const std::string& message, const std::string& file,
                    const int line) const;

  template <typename T>
  void CommunicateData(const std::vector<std::vector<T>>& send_data,
                       std::vector<std::vector<T>>& recv_data) const;
  template <typename T>
  T Reduce(const T data, const Op op) const;
};
extern std::shared_ptr<MPI> AVOCADO_MPI_NAME;

class Communicate {
 private:
  int size_;
  const std::vector<std::vector<int>>& send_data_list_;
  const std::vector<std::vector<int>>& recv_data_list_;

  std::vector<MPI_Request> recv_req_;
  std::vector<MPI_Status> recv_sta_;
  std::vector<MPI_Request> send_req_;
  std::vector<MPI_Status> send_sta_;
  std::vector<int> num_send_data_;
  std::vector<int> num_recv_data_;
  std::vector<std::vector<double>> send_data_;
  std::vector<std::vector<double>> recv_data_;

  int nrecv_domain_;
  int nsend_domain_;

 public:
  Communicate(const int size,
              const std::vector<std::vector<int>>& send_data_list,
              const std::vector<std::vector<int>>& recv_data_list);
  ~Communicate(){};

  void CommunicateBegin(const double* send_data);
  void CommunicateEnd(double* recv_data);
};

template <typename T>
void MPI::CommunicateData(const std::vector<std::vector<T>>& send_data,
                          std::vector<std::vector<T>>& recv_data) const {
  if (NDOMAIN > 1) {
    const int send_data_size = static_cast<int>(send_data.size());
    ASSERT(send_data_size == num_domains_);

    std::vector<int> num_send_data(num_domains_);
    std::vector<int> num_recv_data(num_domains_);

    for (int i = 0; i < num_domains_; i++)
      num_send_data[i] = static_cast<int>(send_data[i].size());
    MPI_Alltoall(&num_send_data[0], 1, MPI_INT, &num_recv_data[0], 1, MPI_INT,
                 MPI_COMM_WORLD);

    recv_data.clear();
    recv_data.resize(num_domains_);
    for (int i = 0; i < num_domains_; i++)
      recv_data[i].resize(num_recv_data[i]);

    std::vector<MPI_Request> recv_req(num_domains_);
    std::vector<MPI_Status> recv_sta(num_domains_);
    std::vector<MPI_Request> send_req(num_domains_);
    std::vector<MPI_Status> send_sta(num_domains_);
    int nrecv_domain = 0;
    int nsend_domain = 0;

    for (int i = 0; i < num_domains_; i++)
      if (num_recv_data[i] > 0 && i != myrank_)
        MPI_Irecv(&recv_data[i][0], num_recv_data[i] * sizeof(T), MPI_CHAR, i,
                  i, MPI_COMM_WORLD, &recv_req[nrecv_domain++]);
    for (int i = 0; i < num_domains_; i++)
      if (num_send_data[i] > 0 && i != myrank_)
        MPI_Isend(&send_data[i][0], num_send_data[i] * sizeof(T), MPI_CHAR, i,
                  myrank_, MPI_COMM_WORLD, &send_req[nsend_domain++]);
    recv_data[myrank_] = send_data[myrank_];
    MPI_Waitall(nrecv_domain, &recv_req[0], &recv_sta[0]);
    MPI_Waitall(nsend_domain, &send_req[0], &send_sta[0]);
  } else {
    recv_data = send_data;
  }
}
template <typename T>
T MPI::Reduce(const T data, const Op op) const {
  if (NDOMAIN > 1) {
    T reduce_data;
    MPI_Datatype datatype = MPI_INT;
    MPI_Op mpiop = MPI_SUM;
    if (typeid(T) == typeid(int))
      datatype = MPI_INT;
    else if (typeid(T) == typeid(double))
      datatype = MPI_DOUBLE;
    else
      ERROR_MESSAGE(
          "Not supported data type: " + std::string(typeid(T).name()) + "\n");
    if (op == Op::SUM)
      mpiop = MPI_SUM;
    else if (op == Op::MIN)
      mpiop = MPI_MIN;
    else if (op == Op::MAX)
      mpiop = MPI_MAX;
    else
      ERROR_MESSAGE("Not supported operation: " +
                    std::to_string(static_cast<int>(op)) + "\n");
    MPI_Allreduce(&data, &reduce_data, 1, datatype, mpiop, MPI_COMM_WORLD);
    return reduce_data;
  }
  return data;
}
}  // namespace avocado