#pragma once

#include <memory>
#include <string>
#include <vector>

#define AVOCADO_ARGUMENT_NAME argument_global_ptr
#define AVOCADO_ARGUMENT avocado::AVOCADO_ARGUMENT_NAME
#define AVOCADO_ARGUMENT_INITIALIZE() \
  AVOCADO_ARGUMENT = std::make_shared<avocado::Argument>()
#define AVOCADO_ARGUMENT_FINALIZE() AVOCADO_ARGUMENT.reset()

#define AVOCADO_HELP "--help"
#define AVOCADO_CONFIG_PATH "-c"
#define AVOCADO_PRE_OPTION "-i"
#define AVOCADO_POST_OPTION "-f"
#define AVOCADO_OTHERS "-s"

namespace avocado {
class Argument {
 private:
  std::string help_option_;
  std::vector<std::string> options_;
  std::vector<std::string> comments_;
  std::vector<std::vector<std::string>> arguments_;

 public:
  Argument();
  ~Argument(){};

  const std::vector<std::string>& GetArguments(const std::string& option) const;
  void SetHelpOption(const std::string& help_option);
  void SetOption(const std::string& option, const std::string& comment);
  void SetArgument(const std::string& option, const std::string& argument);

  void ParseArguments(int argc, const std::vector<std::string>& argv);
  std::string PrintHelpMessage(void);
  std::string PrintSummary(void);
};
extern std::shared_ptr<Argument> AVOCADO_ARGUMENT_NAME;
}  // namespace avocado