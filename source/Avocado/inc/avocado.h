#pragma once

#include "avocado_argument.h"
#include "avocado_blas.h"
#include "avocado_config.h"
#include "avocado_config_macro.h"
#include "avocado_dual_number.h"
#include "avocado_file.h"
#include "avocado_header.h"
#include "avocado_memory_profiler.h"
#include "avocado_mpi.h"
#include "avocado_string.h"
#include "avocado_timer.h"

namespace avocado {
void AVOCADO_Initialize(int argc, char* argv[], const char* codename,
                        const std::vector<std::string>& add_argv);
void AVOCADO_Finalize();
}  // namespace avocado