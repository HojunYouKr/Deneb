#pragma once

#include <memory>
#include <array>
#include <string>
#include <vector>

#include "avocado_header.h"

#ifdef AVOCADO_USE_WINDOWS
#include <psapi.h>
#include <tchar.h>
#include <windows.h>
#else
#include <limits.h>
#include <sys/param.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <unistd.h>
#endif

#define AVOCADO_MEMORY_PROFILER_NAME mp_global_ptr
#define AVOCADO_MEMORY_PROFILER avocado::AVOCADO_MEMORY_PROFILER_NAME
#define AVOCADO_MEMORY_PROFILER_INITIALIZE() \
  AVOCADO_MEMORY_PROFILER = std::make_shared<avocado::MemoryProfiler>()
#define AVOCADO_MEMORY_PROFILER_FINALIZE() AVOCADO_MEMORY_PROFILER.reset()

#define AVOCADO_MEMORY_CHECK() AVOCADO_MEMORY_PROFILER->MemoryCheck()
#define AVOCADO_MEMORY_CONSUMED() AVOCADO_MEMORY_PROFILER->GetMemoryConsumed()

namespace avocado {
class MemoryProfiler {
 private:
  std::string hostname_;
  int coretag_;
  std::array<double, 2> memory_used_;

  // only for MASTER_NODE
  std::vector<std::string> hostname_table_;
  std::vector<int> coretag_table_;
  std::vector<int> sort_mapping_;
  std::vector<int> host_ptr_;

 public:
  MemoryProfiler();
  ~MemoryProfiler(){};

  inline int GetMemoryConsumed(void) const {
    return static_cast<int>(memory_used_[1] - memory_used_[0]);
  };
  void PrintMemoryUsage(void);
  void MemoryCheck(void);

 private:
  const std::string GetHostName(void);
  size_t GetCurrentRSS(void);
  size_t GetMemorySize(void);
};
extern std::shared_ptr<MemoryProfiler> AVOCADO_MEMORY_PROFILER_NAME;
}  // namespace avocado