#pragma once

#include <string>

namespace avocado {
enum class pathattr {
  NOEXIST, DIR, FILE, UNDEF
};

// The return is an absolte path in a directory form
// that ends by '/'.
void GetPresentDirectory(std::string& dir);

// Both relative and absolute paths are allowable.
bool ToAbsolutePath(std::string& path);
pathattr GetPathAttribute(const std::string& path);
bool MakeDirectory(const std::string& path);
void CopyFileTo(const std::string& path, const std::string& dest_path);
}  // namespace avocado