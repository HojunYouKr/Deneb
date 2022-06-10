#include "avocado_file.h"

#include "avocado_header.h"
#include "avocado_string.h"

namespace avocado {
void GetPresentDirectory(std::string& dir) {
  char* buffer;
#ifdef AVOCADO_USE_WINDOWS
  buffer = _getcwd(NULL, 0);
  dir = buffer;
  ReplaceString(dir, "\\", "/");
#elif defined(AVOCADO_USE_LINUX)
  buffer = getcwd(NULL, 0);
  dir = buffer;
#endif
  free(buffer);
  dir += "/";
}

bool ToAbsolutePath(std::string& path) {
  bool absolute_path = false;
#ifdef AVOCADO_USE_WINDOWS
  if (path.length() >= 1 && path.at(0) == '/') {
    path = "Wrong path for Windows: " + path + "\n";
    return false;
  }
  if (path.length() >= 2 && path.at(1) == ':' &&
      ('A' <= path.at(0) && path.at(0) <= 'Z'))
    absolute_path = true;
#elif defined(AVOCADO_USE_LINUX)
  if (path.length() >= 2 && path.at(1) == ':' &&
      ('A' <= path.at(0) && path.at(0) <= 'Z')) {
    path = "Wrong path for Linux: " + path + "\n";
    return false;
  }
  if (path.length() >= 1 && path.at(0) == '/') absolute_path = true;
#endif
  if (absolute_path == false) {
    std::string pwd;
    GetPresentDirectory(pwd);
    path = pwd + "/" + path;
  }
  ReplaceString(path, "//", "/");
  ReplaceString(path, "/./", "/");
  std::string::size_type offset = 0;
  while (true) {
    offset = path.find("/..");
    if (offset == std::string::npos) {
      break;
    } else {
      const std::string::size_type loc = path.rfind("/", offset - 1);
      if (loc == std::string::npos) {
        path = "Wrong path: " + path + "\n";
        return false;
      }
      path.replace(loc, offset - loc + 3, "");
    }
  }
  if (path.back() == '.') path.pop_back();
  return true;
}
pathattr GetPathAttribute(const std::string& path) {
  std::string abs_path = path;
  ToAbsolutePath(abs_path);
#ifdef AVOCADO_USE_WINDOWS
  DWORD attr = GetFileAttributes(abs_path.c_str());
  if (attr == INVALID_FILE_ATTRIBUTES)
    return pathattr::NOEXIST;
  else if (attr & FILE_ATTRIBUTE_DIRECTORY)
    return pathattr::DIR;
  else if (attr & FILE_ATTRIBUTE_ARCHIVE)
    return pathattr::FILE;
#elif defined(AVOCADO_USE_LINUX)
  struct stat st_buf;
  const int status = stat(abs_path.c_str(), &st_buf);
  if (status != 0)
    return pathattr::NOEXIST;
  else if (S_ISDIR(st_buf.st_mode))
    return pathattr::DIR;
  else if (S_ISREG(st_buf.st_mode))
    return pathattr::FILE;
#endif
  return pathattr::UNDEF;
}

bool MakeDirectory(const std::string& path) {
  std::string abs_path = path;
  if (ToAbsolutePath(abs_path) == false) return false;
  const std::vector<std::string> strtoks = SplitString(abs_path, '/');
  abs_path.clear();
  for (size_t i = 0, length = strtoks.size(); i < length - 1; i++) {
    abs_path = abs_path + strtoks[i] + "/";
    const pathattr dirattr = GetPathAttribute(abs_path);
    if (dirattr == pathattr::FILE) return false;
    if (dirattr == pathattr::NOEXIST) {
#ifdef AVOCADO_USE_WINDOWS
      if (_mkdir(abs_path.c_str()) == -1) return false;
#elif defined(AVOCADO_USE_LINUX)
      if (mkdir(abs_path.c_str(), 0755) == -1) return false;
#endif
    }
  }
  return true;
}

void CopyFileTo(const std::string& path, const std::string& dest_path) {
#ifdef AVOCADO_USE_WINDOWS
  std::string command = "copy " + path + " " + dest_path;
  ReplaceString(command, "/", "\\");
  system(command.c_str());
#elif defined(AVOCADO_USE_LINUX)
  std::string command = "cp " + path + " " + dest_path;
  system(command.c_str());
#endif
}
}  // namespace avocado