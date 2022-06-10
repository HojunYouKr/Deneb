#pragma once

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#define AVOCADO_CONFIG_NAME config_global_ptr
#define AVOCADO_CONFIG avocado::AVOCADO_CONFIG_NAME
#define AVOCADO_CONFIG_INITALIZE() \
  AVOCADO_CONFIG = std::make_shared<avocado::Config>()
#define AVOCADO_CONFIG_FINALIZE() AVOCADO_CONFIG.reset()

namespace avocado {
class Config {
 private:
  std::unordered_map<std::string, std::vector<std::string>> configs_;

 public:
  Config();
  ~Config(){};

  inline void CleanConfig(void) {
    configs_.clear();
  }

  const std::string& GetConfigValue(const std::string& config_type,
                                    const int n);
  const std::string& GetConfigValue(const std::string& config_type);
  void SetConfigValue(const std::string& config_str);
  void DeleteConfigValue(const std::string& config_type);
  void ReadConfigFile(const std::string& config_file_path);
};
extern std::shared_ptr<Config> AVOCADO_CONFIG_NAME;
}  // namespace deneb
