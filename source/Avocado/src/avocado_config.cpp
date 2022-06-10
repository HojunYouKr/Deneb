#include "avocado_config.h"

#include <fstream>

#include "avocado_config_macro.h"
#include "avocado_file.h"
#include "avocado_mpi.h"
#include "avocado_string.h"

namespace avocado {
std::shared_ptr<Config> AVOCADO_CONFIG_NAME = nullptr;

Config::Config() : configs_() {
  std::string pwd = "";
  GetPresentDirectory(pwd);
  configs_[PWD] = {pwd};
  configs_["MYRANK"] = {std::to_string(MYRANK)};
  configs_["NDOMAIN"] = {std::to_string(NDOMAIN)};
};

const std::string& Config::GetConfigValue(const std::string& config_type,
                                          const int n) {
  const auto& iterator = configs_.find(config_type);
  if (iterator == configs_.end())
    ERROR_MESSAGE("No configuration type: " + config_type + "\n");
  if (iterator->second.size() <= n)
    ERROR_MESSAGE("No configuration value: " + config_type + "." +
                  std::to_string(n) + "\n");
  return iterator->second[n];
}

const std::string& Config::GetConfigValue(const std::string& config_type) {
  std::string config_str = config_type;
  if (config_str.find('.') == std::string::npos) config_str += ".0";
  const std::vector<std::string> split = SplitString(config_str, '.');
  return GetConfigValue(split[0], std::stoi(split[1]));
}

void Config::SetConfigValue(const std::string& config_str) {
  std::string config_type;
  std::vector<std::string> config_values;
  {  // pre-process
    std::string temp = config_str;
    avocado::RemoveString(temp, std::vector<char>{' ', '\0', '\n', '\r', '\t'});
    const size_t comment_pos = temp.find("//");
    if (comment_pos != std::string::npos) temp.erase(comment_pos);
    if (temp.length() == 0) return;
    if (temp.find('=') == std::string::npos)
      ERROR_MESSAGE("Wrong configuration (no equal symbol): " + config_str +
                    "\n");
    const std::vector<std::string> split = avocado::SplitString(temp, '=');
    if (split.size() != 2)
      ERROR_MESSAGE("Wrong configuration (wrong usage of equal symbol): " +
                    config_str + "\n");
    config_type = split[0];
    config_values = avocado::SplitString(split[1], ',');
  }

  // variable
  for (auto&& config_value : config_values) {
    while (true) {
      const size_t start_pos = config_value.find("&(");
      if (start_pos == std::string::npos) break;
      const size_t end_pos = config_value.find(")", start_pos);
      if (end_pos == std::string::npos)
        ERROR_MESSAGE("Wrong configuration (wrong usage of &(variable)): " +
                      config_value + "\n");
      std::string variable =
          config_value.substr(start_pos + 2, end_pos - start_pos - 2);
      if (variable.find('.') == std::string::npos) variable += ".0";
      const std::vector<std::string> split = SplitString(variable, '.');
      variable = split[0];
      const int nvar = std::stoi(split[1]);

      const std::string value = GetConfigValue(variable, nvar);
      config_value.replace(start_pos, end_pos - start_pos + 1, value);
    }
  }

  // expression
  for (auto&& config_value : config_values) {
    while (true) {
      const size_t start_pos = config_value.find("[");
      if (start_pos == std::string::npos) break;
      const size_t end_pos = config_value.find("]", start_pos);
      if (end_pos == std::string::npos)
        ERROR_MESSAGE("Wrong configuration (wrong usage of [expr]): " +
                      config_value + "\n");
      std::string expression =
          config_value.substr(start_pos + 1, end_pos - start_pos - 1);
      avocado::EvaluationString(expression);
      config_value.replace(start_pos, end_pos - start_pos + 1, expression);
    }
  }

  // Command
  if (!config_type.compare(COMMAND)) {
    if (!config_values[0].compare(QUIT_PROGRAM))
      ERROR_MESSAGE("Quit program by configuration\n");
    else if (!config_values[0].compare(LOAD_CONFIG)) {
      if (config_values.size() == 1)
        ERROR_MESSAGE("No configuation file path for LoadConfig command\n");
      MASTER_MESSAGE("Load configuration file by config\n");
      ReadConfigFile(config_values[1]);
      return;
    } else
      ERROR_MESSAGE("Wrong command\n");
  }

  DeleteConfigValue(config_type);
  configs_[config_type] = config_values;
  std::string regist_message = config_type + " = ";
  for (auto&& config_value : config_values)
    regist_message += (config_value + ", ");
  regist_message.pop_back();
  regist_message.pop_back();
  MASTER_MESSAGE("Set configuration: " + regist_message + "\n");
}

void Config::DeleteConfigValue(const std::string& config_type) {
  std::unordered_map<std::string, std::vector<std::string>>::iterator iter;
  for (iter = configs_.begin(); iter != configs_.end();) {
    if (!iter->first.compare(config_type)) {
      MASTER_MESSAGE("Delete configuration: " + config_type + "\n");
      configs_.erase(iter++);
    } else
      ++iter;
  }
}

void Config::ReadConfigFile(const std::string& config_file_path) {
  std::string path = config_file_path;
  ToAbsolutePath(path);
  std::ifstream infile(path.c_str());
  if (infile.is_open())
    MASTER_MESSAGE("Open configuration file: " + path + "\n");
  else
    ERROR_MESSAGE("Illegal configuration file path (no-exist): " + path + "\n");
  std::string text;
  while (getline(infile, text)) SetConfigValue(text);
  infile.close();
  MASTER_MESSAGE("Close configuration file: " + path + "\n");
}
}  // namespace avocado