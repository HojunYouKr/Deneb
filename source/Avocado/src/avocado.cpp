#include "avocado.h"

#include <sstream>
#include <string>

namespace avocado {
void AVOCADO_Initialize(int argc, char* argv[], const char* codename,
                        const std::vector<std::string>& add_argv) {
  const char* filepath = argv[0];
  AVOCADO_MPI_INITIALIZE(codename, filepath);
  AVOCADO_ARGUMENT_INITIALIZE();
  AVOCADO_CONFIG_INITALIZE();
  AVOCADO_TIMER_INITIALIZE();
  AVOCADO_MEMORY_PROFILER_INITIALIZE();

  MASTER_MESSAGE(avocado::GetTitle("Avocado Initialization"));
  std::stringstream msg;
  msg << AVOCADO_VERSION;
  msg << " (Last update: " << AVOCADO_LAST_UPDATE << ")";
  msg << "\n\tCodename: " << codename;
  msg << "\n\tFilepath: " << filepath;
  msg << "\n\tBuild date: " << __DATE__ << " " << __TIME__;
  msg << "\n\tRun date: " << avocado::GetLocalTime();
  msg << "\n\tNumber of processors: " << NDOMAIN;
#ifdef AVOCADO_USE_WINDOWS
  msg << "\n\tRun on WINDOWS";
#elif defined(AVOCADO_USE_LINUX)
  msg << "\n\tRun on LINUX";
#endif
  MASTER_MESSAGE(msg.str() + "\n");

  // argument
  MASTER_MESSAGE(avocado::GetTitle("Avocado Argument Parsing"));
  auto& argument = AVOCADO_ARGUMENT;
  std::vector<std::string> total_argv;
  const int total_argc = argc + static_cast<int>(add_argv.size());
  for (int i = 0; i < argc; i++) total_argv.push_back(std::string(argv[i]));
  for (auto&& add : add_argv) total_argv.push_back(add);
  argument->ParseArguments(total_argc, total_argv);
  MASTER_MESSAGE(argument->PrintSummary() + "\n");

  // configuration
  MASTER_MESSAGE(avocado::GetTitle("Avocado Configuration"));
  auto& config = AVOCADO_CONFIG;
  for (auto&& option : argument->GetArguments(AVOCADO_PRE_OPTION))
    config->SetConfigValue(option);
  for (auto&& config_file_path : argument->GetArguments(AVOCADO_CONFIG_PATH))
    config->ReadConfigFile(config_file_path);
  for (auto&& option : argument->GetArguments(AVOCADO_POST_OPTION))
    config->SetConfigValue(option);
}
void AVOCADO_Finalize() {
  AVOCADO_MEMORY_PROFILER_FINALIZE();
  AVOCADO_TIMER_FINALIZE();
  AVOCADO_CONFIG_FINALIZE();
  AVOCADO_ARGUMENT_FINALIZE();
  AVOCADO_MPI_FINALIZE();
}
}  // namespace avocado