#pragma once

#include <vector>

#include "deneb_system_matrix.h"
#include "deneb_timescheme.h"

namespace deneb {
// -------------------- Rosenbrock R-K ------------------- //
class TimeschemeRosenbrock : public Timescheme {
 private:
  Mat sysmat_;
  GMRES solver_;

  Vec rhs_;
  Vec delta_;

  std::vector<Vec> stage_solution_;

 protected:
  int time_order_;
  int num_stages_;

  double coeff_gamma_;
  std::vector<double> coeff_m_;
  std::vector<double> coeff_r_;
  std::vector<double> coeff_alpha_;
  std::vector<std::vector<double>> coeff_a_;
  std::vector<std::vector<double>> coeff_c_;

 public:
  TimeschemeRosenbrock();
  ~TimeschemeRosenbrock();

  virtual void BuildData(void);
  virtual void Marching(void);
};

// ------------------------ ROS3P ------------------------ //
// Lang-Verwer ROS3P (RO3-3)
class TimeschemeROS3P : public TimeschemeRosenbrock {
 public:
  TimeschemeROS3P();
  ~TimeschemeROS3P();

  virtual void BuildData(void);
  virtual void Marching(void);
};
// ----------------------- RODAS3 ------------------------ //
// Hairer-Wanner RODAS3 (RO3-4)
class TimeschemeRODAS3 : public TimeschemeRosenbrock {
 public:
  TimeschemeRODAS3();
  ~TimeschemeRODAS3();

  virtual void BuildData(void);
  virtual void Marching(void);
};
// ------------------------ ROS4 ------------------------- //
// Shampine ROS4 (RO4-4)
class TimeschemeROS4 : public TimeschemeRosenbrock {
 public:
  TimeschemeROS4();
  ~TimeschemeROS4();

  virtual void BuildData(void);
  virtual void Marching(void);
};
// ----------------------- RODASP ------------------------ //
// Steinebach RODASP (RO4-6)
class TimeschemeRODASP : public TimeschemeRosenbrock {
 public:
  TimeschemeRODASP();
  ~TimeschemeRODASP();

  virtual void BuildData(void);
  virtual void Marching(void);
};
// ----------------------- RODAS5 ------------------------ //
// Di Marzo RODAS5-Rod5_1 (RO5-8)
class TimeschemeRODAS5 : public TimeschemeRosenbrock {
 public:
  TimeschemeRODAS5();
  ~TimeschemeRODAS5();

  virtual void BuildData(void);
  virtual void Marching(void);
};
// ----------------------- ROW6A ------------------------- //
// Kaps and Wanner ROW6A (RO6-6)
class TimeschemeROW6A : public TimeschemeRosenbrock {
 public:
  TimeschemeROW6A();
  ~TimeschemeROW6A();

  virtual void BuildData(void);
  virtual void Marching(void);
};
}  // namespace deneb