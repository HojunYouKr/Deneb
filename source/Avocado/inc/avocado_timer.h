#pragma once

#include <chrono>
#include <memory>
#include <string>
#include <unordered_map>

#define AVOCADO_TIMER_NAME timer_global_ptr
#define AVOCADO_TIMER avocado::AVOCADO_TIMER_NAME
#define AVOCADO_TIMER_INITIALIZE() \
  AVOCADO_TIMER = std::make_shared<avocado::Timer>()
#define AVOCADO_TIMER_FINALIZE() AVOCADO_TIMER.reset()

#define START_TIMER() AVOCADO_TIMER->StartTimer("")
#define STOP_TIMER() AVOCADO_TIMER->StopTimer("")
#define GET_RECORD() AVOCADO_TIMER->GetRecord("")
#define START_TIMER_TAG(tag) AVOCADO_TIMER->StartTimer(std::string(#tag))
#define STOP_TIMER_TAG(tag) AVOCADO_TIMER->StopTimer(std::string(#tag))
#define GET_RECORD_TAG(tag) AVOCADO_TIMER->GetRecord(std::string(#tag))

namespace avocado {
class Timer {
 private:
   struct watch {
     std::chrono::system_clock::time_point start_time;
     std::chrono::system_clock::time_point stop_time;
     double record;
   };

 private:
   std::unordered_map<std::string, watch> timer_;

 public:
  Timer(): timer_() {};
  ~Timer(){};

  double GetRecord(const std::string& tag = "") const;
  void StartTimer(const std::string& tag = "");
  double StopTimer(const std::string& tag = "");
};
extern std::shared_ptr<Timer> AVOCADO_TIMER_NAME;
}  // namespace avocado