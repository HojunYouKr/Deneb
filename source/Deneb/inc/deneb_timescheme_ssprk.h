#pragma once

#include <vector>

#include "deneb_timescheme.h"

namespace deneb {
// ------------------------ SSPRK ------------------------ //
class TimeschemeSSPRK : public Timescheme {
 protected:
  double* ptr_;
  std::vector<double> c_;
  std::vector<double> rhs_;
  std::vector<double> u_;

 public:
  TimeschemeSSPRK();
  virtual ~TimeschemeSSPRK();

  virtual void BuildData(void);
  virtual void Marching(void);

protected:
  virtual void UpdateSolution(double* solution, const double t,
                              const double dt) = 0;
};

// ----------------------- SSPRK11 ----------------------- //
class TimeschemeSSPRK11 : public TimeschemeSSPRK {
 public:
  TimeschemeSSPRK11();
  virtual ~TimeschemeSSPRK11();

  virtual void BuildData(void);
  virtual void Marching(void);

 protected:
  virtual void UpdateSolution(double* solution, const double t,
                              const double dt);
};

// ----------------------- SSPRK33 ----------------------- //
class TimeschemeSSPRK33 : public TimeschemeSSPRK {
 private:
  const double c13_ = 1.0 / 3.0;

 public:
  TimeschemeSSPRK33();
  virtual ~TimeschemeSSPRK33();

  virtual void BuildData(void);
  virtual void Marching(void);

 protected:
  virtual void UpdateSolution(double* solution, const double t,
                              const double dt);
};

// ----------------------- SSPRK54 ----------------------- //
class TimeschemeSSPRK54 : public TimeschemeSSPRK {
 protected:
  double* ptr2_;
  std::vector<double> rhs2_;
  std::vector<double> u2_;

 public:
  TimeschemeSSPRK54();
  virtual ~TimeschemeSSPRK54();

  virtual void BuildData(void);
  virtual void Marching(void);

 protected:
  virtual void UpdateSolution(double* solution, const double t,
                              const double dt);
};
}  // namespace deneb