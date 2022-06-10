#pragma once

#include "deneb_system_matrix.h"
#include "deneb_timescheme.h"

namespace deneb {
// -------------------- ImplicitEuler -------------------- //
class TimeschemeImpEuler : public Timescheme {
 private:
  Mat sysmat_;
  GMRES solver_;

  Vec rhs_;
  Vec delta_;

 public:
  TimeschemeImpEuler();
  virtual ~TimeschemeImpEuler();

  virtual void BuildData(void);
  virtual void Marching(void);
};
}  // namespace deneb