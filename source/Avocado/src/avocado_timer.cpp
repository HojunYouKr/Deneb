#include "avocado_timer.h"

#include "avocado_mpi.h"

namespace avocado {
std::shared_ptr<Timer> AVOCADO_TIMER_NAME = nullptr;

double Timer::GetRecord(const std::string& tag) const {
  const auto& iterator = timer_.find(tag);
  if (iterator == timer_.end())
    ERROR_MESSAGE("No record for tag: " + tag + "\n");
  return iterator->second.record;
}

void Timer::StartTimer(const std::string& tag) {
  timer_[tag].start_time = std::chrono::system_clock::now();
}

double Timer::StopTimer(const std::string& tag) {
  const std::chrono::system_clock::time_point stop_time =
      std::chrono::system_clock::now();
  const auto& iterator = timer_.find(tag);
  if (iterator == timer_.end())
    ERROR_MESSAGE("No timer for tag: " + tag + "\n");
  iterator->second.stop_time = stop_time;
  const std::chrono::duration<double> gap = stop_time - iterator->second.start_time;
  iterator->second.record = gap.count();
  return iterator->second.record;
}
}  // namespace avocado