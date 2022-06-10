#pragma once

#include <string>
#include <vector>

#include "deneb_element.h"

#define DENEB_PRESSUREFIX_NAME pressurefix_global_ptr
#define DENEB_PRESSUREFIX deneb::DENEB_PRESSUREFIX_NAME
#define DENEB_PRESSUREFIX_INITIALIZE(name) \
  DENEB_PRESSUREFIX = deneb::PressureFix::GetPressureFix(name)
#define DENEB_PRESSUREFIX_FINALIZE() DENEB_PRESSUREFIX.reset()

namespace deneb {
class PressureFix {
 public:
  static std::shared_ptr<PressureFix> GetPressureFix(const std::string& name);

  PressureFix(){};
  virtual ~PressureFix(){};

  virtual void BuildData(void) = 0;
  virtual void Execute(double* solution) = 0;
};
extern std::shared_ptr<PressureFix> DENEB_PRESSUREFIX_NAME;

class PressureFixOff : public PressureFix {
 public:
  PressureFixOff();
  virtual ~PressureFixOff(){};

  virtual void BuildData(void);
  virtual void Execute(double* solution){};
};

class PressureFixOn : public PressureFix {
 protected:
  const int iteration_;

  std::vector<double> solution_;

 public:
  PressureFixOn();
  virtual ~PressureFixOn(){};

  virtual void BuildData(void);
  virtual void Execute(double* solution);
};
}  // namespace deneb