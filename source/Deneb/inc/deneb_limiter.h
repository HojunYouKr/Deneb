#pragma once

#include <memory>
#include <string>
#include <vector>

#include "deneb_basis.h"
#include "deneb_element.h"
#include "deneb_jacobian.h"

#define DENEB_LIMITER_NAME limiter_global_ptr
#define DENEB_LIMITER deneb::DENEB_LIMITER_NAME
#define DENEB_LIMITER_INITIALIZE(name) \
  DENEB_LIMITER = deneb::Limiter::GetLimiter(name)
#define DENEB_LIMITER_FINALIZE() DENEB_LIMITER.reset()

namespace avocado {
class Communicate;
}

namespace deneb {
class Limiter {
 public:
  static std::shared_ptr<Limiter> GetLimiter(const std::string& name);

  Limiter(){};
  virtual ~Limiter(){};

  virtual void BuildData(void) = 0;
  virtual void Limiting(double* solution) = 0;
};
extern std::shared_ptr<Limiter> DENEB_LIMITER_NAME;

class NoLimiter : public Limiter {
 public:
  NoLimiter();
  virtual ~NoLimiter(){};

  virtual void BuildData(void);
  virtual void Limiting(double* solution){};
};

class hMLP : public Limiter {
 protected:
  enum class NODETYPE : char { NORMAL = 0, BOUNDARY = 1 };

  const double eps_;
  const double deact_eps_;

  std::vector<std::vector<int>> node_cells_;
  std::vector<std::vector<int>> node_vertices_;

  int target_state_;
  std::vector<int> num_bases_list_;
  std::vector<NODETYPE> nodetypes_;
  std::vector<double> foreign_solution_;
  std::shared_ptr<avocado::Communicate> communicate_;

  std::vector<std::vector<double>> cell_vertex_basis_value_;
  std::vector<double> foreign_cell_basis_value_;
  std::vector<double> cell_average_;
  std::vector<double> vertex_min_;
  std::vector<double> vertex_max_;

 public:
  hMLP();
  virtual ~hMLP(){};

  virtual void BuildData(void);
  virtual void Limiting(double* solution);

 protected:
  void ConstructNodeCells(void);
  void ComputeVertexMinMax(double* solution);
  bool Indicator(const double* solution_ptr, const int icell, const int order);
  bool MLPcondition(const double* solution_ptr, const int icell,
                    const int ivertex, const int istate, const double P1var,
                    const double vertex_value, const double average_value);
  bool Deactivation(const double* solution_ptr, const double volume,
                    const double vertex_value, const double average_value);
  bool SmoothDetect(const double* solution_ptr, const int icell,
                    const int ivertex, const int istate, const double P1var,
                    const double vertex_value, const double average_value);
  void Projection(double* solution_ptr, const int icell, const int order);
  void MLPu1(double* solution_ptr, const int icell);
};

class hMLP_BD : public Limiter {
 protected:
  enum class NODETYPE : char { NORMAL = 0, BOUNDARY = 1 };

  const double eps_;
  const double deact_eps_;

  std::vector<std::vector<int>> cell_cells_;
  std::vector<std::vector<int>> cell_faces_;
  std::vector<std::vector<int>> node_cells_;
  std::vector<std::vector<int>> node_vertices_;

  int target_state_;
  std::vector<int> num_bases_list_;
  std::vector<NODETYPE> nodetypes_;
  std::vector<double> foreign_solution_;
  std::shared_ptr<avocado::Communicate> communicate_;

  std::vector<int> cell_num_troubled_boundary_;
  std::vector<int> cell_num_for_type_one_;
  std::vector<int> cell_num_for_type_two_;
  std::vector<double> vertex_min_;
  std::vector<double> vertex_max_;
  std::vector<double> vertex_average_;
  std::vector<double> face_areas_;
  std::vector<double> face_characteristic_length_;
  std::vector<std::vector<double>> cell_vertex_basis_value_;
  std::vector<std::vector<double>> face_coefficients_;
  std::vector<std::vector<double>> simplex_average_coefficients_;
  std::vector<std::vector<double>> simplex_p1_coefficients_;
  std::vector<std::vector<double>> simplex_volume_;
  std::vector<std::vector<double>> face_difference_;

  typedef std::vector<std::vector<std::vector<int>>> vector3_t;
  const vector3_t SED_simplex_node_ind_ = {
      {{0, 1, 3}, {1, 2, 0}, {2, 3, 1}, {3, 0, 2}},  // Quad
      {{0, 1, 3, 4},
       {1, 2, 0, 5},
       {2, 3, 1, 6},
       {3, 0, 2, 7},
       {4, 7, 5, 0},
       {5, 4, 6, 1},
       {6, 5, 7, 2},
       {7, 6, 4, 3}},  // Hexa
      {{0, 1, 2, 3},
       {1, 2, 0, 4},
       {2, 0, 1, 5},
       {3, 5, 4, 0},
       {4, 3, 5, 1},
       {5, 4, 3, 2}},                                             // Pris
      {{0, 1, 3, 4}, {1, 2, 0, 4}, {2, 3, 1, 4}, {3, 0, 2, 4}}};  // Pyra

 public:
  hMLP_BD();
  virtual ~hMLP_BD(){};

  virtual void BuildData(void);
  virtual void Limiting(double* solution);

 protected:
  void ConstructCellCells(void);
  void ConstructNodeCells(void);
  void ConstructSimplex(void);
  void VertexMinMax(const double* solution);
  bool Indicator(const double* solution_ptr, const int icell, const int order);
  bool MLPcondition(const double* solution_ptr, const int icell,
                    const int ivertex, const int istate, const double P1var,
                    const double vertex_value, const double average_value);
  bool Deactivation(const double* solution_ptr, const double volume,
                    const double vertex_value, const double average_value);
  bool SmoothDetect(const double* solution_ptr, const int icell,
                    const int ivertex, const int istate, const double P1var,
                    const double vertex_value, const double average_value);
  void Projection(double* solution_ptr, const int icell, const int order);
  void MLPu1(double* solution_ptr, const int icell);
  bool BoundaryDetect(const int icell, const int num_troubled_boundary);
};
}  // namespace deneb