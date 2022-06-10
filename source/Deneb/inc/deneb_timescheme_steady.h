#pragma once

#include "deneb_system_matrix.h"
#include "deneb_timescheme.h"

namespace deneb {
// ------------------------ Steady ----------------------- //
class TimeschemeSteady : public Timescheme {
 private:
  Mat sysmat_;
  GMRES solver_;

  Vec rhs_;
  Vec delta_;

  // Options to increase CFL 
  double increase_amount_;
  double steady_stop_residual_;

 public:
  TimeschemeSteady();
  virtual ~TimeschemeSteady();

  virtual void BuildData(void);
  virtual void Marching(void);
};
}  // namespace deneb