#pragma once

#include <memory>
#include <string>
#include <vector>

#define DENEB_TIMESCHEME_NAME timescheme_global_ptr
#define DENEB_TIMESCHEME deneb::DENEB_TIMESCHEME_NAME
#define DENEB_TIMESCHEME_INITIALIZE(name) \
  DENEB_TIMESCHEME = deneb::Timescheme::GetTimescheme(name)
#define DENEB_TIMESCHEME_FINALIZE() DENEB_TIMESCHEME.reset()

namespace deneb {
class Interrupt {
 public:
  enum class Control : int { ITER = 0, TIME = 1, NONE = 2 };

 private:
  int step_;
  Control control_;
  unsigned __int64 begin_;
  unsigned __int64 finish_;
  unsigned __int64 interval_;

 public:
  Interrupt()
      : step_(-1),
        control_(Control::NONE),
        begin_(-1),
        finish_(-1),
        interval_(-1){};
  ~Interrupt(){};

  inline void Initialize(const Control& control, const unsigned __int64& begin,
                         const unsigned __int64& finish,
                         const unsigned __int64& interval) {
    control_ = control;
    begin_ = begin;
    finish_ = finish;
    interval_ = interval;
  };
  bool CheckIterationFinish(const int iteration);
  bool CheckTimeFinish(const unsigned __int64& time,
                       unsigned __int64& timestep);
  bool CheckIterationEvent(const int iteration);
  bool CheckTimeEvent(const unsigned __int64& time, unsigned __int64& timestep);
};

class Timescheme {
 public:
  static std::shared_ptr<Timescheme> GetTimescheme(const std::string& name);

 protected:
  enum class TimestepControl : int { CFL = 0, DT = 1 };

 protected:
  bool is_implicit_;

  int order_;
  int max_iteration_;
  int iteration_;

  int time_resolution_;
  unsigned __int64 time_scale_;
  unsigned __int64 current_time_;  // (~ 1.8E+19)

  TimestepControl timestep_control_;
  double timestep_control_value_;

  Interrupt stop_;
  Interrupt post_;
  Interrupt save_;

  int length_;
  std::vector<double> solution_;
  std::vector<double> local_timestep_;
  double computing_cost_;

 public:
  Timescheme(const bool is_implicit);
  virtual ~Timescheme(){};

  inline bool IsImplicit(void) const { return is_implicit_; };
  inline int GetTimeResolution(void) const { return time_resolution_; };
  inline int GetCurrentIteration(void) const { return iteration_; };
  inline double GetCurrentTime(void) const {
    return ConvertTime(current_time_);
  };
  inline double GetComputingCost(void) const { return computing_cost_; };
  inline unsigned __int64 ConvertTime(const double time) const {
    return static_cast<unsigned __int64>(time * time_scale_);
  }
  inline double ConvertTime(const unsigned __int64 time) const {
    return time / static_cast<double>(time_scale_);
  };
  inline double* GetSolution(void) { return &solution_[0]; };

  virtual void BuildData(void) = 0;
  virtual void Marching(void) = 0;

 protected:
  void ComputeLocalTimestep(const double* solution,
                            std::vector<double>& local_timestep) const;
  double ComputeGlobalTimestep(const double* solution,
                               std::vector<double>& local_timestep) const;
  void InitSolution(void);
};
extern std::shared_ptr<Timescheme> DENEB_TIMESCHEME_NAME;
}  // namespace deneb