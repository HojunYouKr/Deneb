#pragma once

#include <memory>
#include <string>
#include <vector>
#include <fstream>

#define DENEB_SAVELOAD_NAME saveload_global_ptr
#define DENEB_SAVELOAD deneb::DENEB_SAVELOAD_NAME
#define DENEB_SAVELOAD_INITIALIZE() \
  DENEB_SAVELOAD = std::make_shared<deneb::SaveLoad>()
#define DENEB_SAVELOAD_FINALIZE() DENEB_SAVELOAD.reset()

namespace deneb {
class SaveLoad {
 public:
  struct SaveData {
    int iteration_;
    int strandid_;
    unsigned __int64 time_;
  };

 private:
  SaveData data_;
  std::vector<double> solution_;

 public:
  SaveLoad(){};
  ~SaveLoad(){};

  inline void ClearSolution(void) {
    solution_.clear();
  };
  inline void ClearSaveData(void) {
    data_ = SaveData({0, 0, 0});
  };
  inline const std::vector<double>& GetSolution(void) const {
    return solution_;
  };
  inline const SaveData& GetSaveData(void) const { return data_; };

  void Save(const std::string& filename, const double* solution,
            const SaveData& data) const;
  void Load(const std::string& filename);

 private:
  void Compare(std::ifstream& file, const std::string& name, const int save,
               const int current) const;

  template <typename T>
  void WriteData(std::vector<char>& buffer, const T data) const {
    for (size_t i = 0, s = sizeof(T); i < s; i++)
      buffer.push_back(*((char*)&data + i));
  };
  void WriteData(std::vector<char>& buffer, const std::string& data) const {
    for (size_t i = 0, s = data.length(); i < s; i++)
      buffer.push_back(data[i]);
    buffer.push_back('\0');
  };
  template <typename T>
  void ReadData(std::ifstream& file, T* data, const int n = 1) const {
    file.read((char*)data, sizeof(T) * n);
  };
  void ReadData(std::ifstream& file, std::string& data) const {
    char ch = '\0';
    do {
      file.read((char*)&ch, sizeof(char));
      if (ch != '\0') data.append(&ch, 1);
    } while (ch != '\0');
  };
};
extern std::shared_ptr<SaveLoad> DENEB_SAVELOAD_NAME;
}  // namespace deneb
