#include "avocado_argument.h"

#include "avocado_mpi.h"

namespace avocado {
std::shared_ptr<Argument> AVOCADO_ARGUMENT_NAME = nullptr;

Argument::Argument() : help_option_(), options_(), comments_(), arguments_() {
  SetHelpOption(AVOCADO_HELP);
  SetOption(AVOCADO_CONFIG_PATH, "Configuration file path");
  SetOption(AVOCADO_PRE_OPTION, "Pre-processed configuration options");
  SetOption(AVOCADO_POST_OPTION, "Post-processed configuration options");
  SetOption(AVOCADO_OTHERS, "Other commands");
};

const std::vector<std::string>& Argument::GetArguments(
    const std::string& option) const {
  for (size_t loc = 0, len = options_.size(); loc < len; loc++)
    if (!options_[loc].compare(option)) return arguments_[loc];
  ERROR_MESSAGE("Illegal option (non-exist): " + option + "\n");
  return options_;
}

void Argument::SetHelpOption(const std::string& help_option) {
  for (size_t loc = 0, len = options_.size(); loc < len; loc++)
    if (!options_[loc].compare(help_option))
      ERROR_MESSAGE("Illegal help option (conflict with normal option): " + help_option + "\n");
  help_option_ = help_option;
}

void Argument::SetOption(const std::string& option,
                         const std::string& comment) {
  for (size_t loc = 0, len = options_.size(); loc < len; loc++)
    if (!options_[loc].compare(option))
      ERROR_MESSAGE("Illegal option (conflict): " + option + "\n");
  if (!help_option_.compare(option))
    ERROR_MESSAGE("Illegal option (Conflict with help option): " + option + "\n");
  options_.push_back(option);
  comments_.push_back(comment);
  arguments_.push_back(std::vector<std::string>());
}

void Argument::SetArgument(const std::string& option,
                           const std::string& argument) {
  if (!help_option_.compare(argument))
    ERROR_MESSAGE("Print help message\n" + PrintHelpMessage() + "\n");
  for (size_t loc = 0, len = options_.size(); loc < len; loc++) {
    if (!options_[loc].compare(option)) {
      arguments_[loc].push_back(argument);
      return;
    }
  }
  ERROR_MESSAGE("Illegal option (non-exist): " + option + "\n");
}

void Argument::ParseArguments(int argc, const std::vector<std::string>& argv) {
  const size_t len = options_.size();

  size_t argloc = -1;
  for (int i = 1; i < argc; i++) {
    const std::string temp = argv[i];
    if (!help_option_.compare(temp)) {
      MASTER_MESSAGE("Print help message\n" + PrintHelpMessage() + "\n");
      AVOCADO_MPI->StopProgram();
    }
    bool isoption = false;
    for (size_t loc = 0; loc < len; loc++) {
      if (!options_[loc].compare(temp)) {
        argloc = loc;
        isoption = true;
        break;
      }
    }
    if (isoption) continue;
    if (argloc != -1)
      arguments_[argloc].push_back(temp);
    else
      ERROR_MESSAGE("Wrong argument: " + temp + "\n\tPress " + help_option_ +
                       " for help\n");
  }
}

std::string Argument::PrintHelpMessage(void) {
  std::string message = "Help message";
  message += ("\n\t" + help_option_ + ",\t<Help>");
  for (size_t i = 0, len = options_.size(); i < len; i++)
    message += ("\n\t" + options_[i] + ",\t<" + comments_[i] + ">");
  return message;
}

std::string Argument::PrintSummary(void) {
  std::string summary = "Argument summary";
  for (size_t i = 0, len = options_.size(); i < len; i++) {
    summary += ("\n\t" + options_[i] + ",\t<" + comments_[i] + ">");
    for (auto&& arg : arguments_[i]) summary += ("\n\t\t" + arg);
  }
  return summary;
}
}  // namespace avocado