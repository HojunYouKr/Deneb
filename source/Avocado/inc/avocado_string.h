#pragma once

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#define DEBUG_VECTOR_SHOW(vec) \
  std::string(#vec) + ": " + avocado::VecToString(vec) + "\n"
#define DEBUG_VAR_SHOW(var) \
  std::string(#var) + ": " + std::to_string(var) + "\n"

namespace avocado {
std::string GetTitle(const std::string& str);
std::string GetLocalTime(void);
void ReplaceString(std::string& str, const std::string& find_tok,
                   const std::string& replace_tok);
std::vector<std::string> SplitString(const std::string& str,
                                     const char delimiter);
void RemoveString(std::string& str, const std::vector<char>& drops);
void RemoveString(std::string& str, const std::string& remove_tok);
template <typename T>
const std::string VecToString(const std::vector<T>& vec) {
  const int size = static_cast<int>(vec.size());
  if (size == 0) return std::string("");
  std::string str = std::to_string(vec[0]);
  for (int i = 1; i < size; i++) str += (", " + std::to_string(vec[i]));
  return str;
}

template <typename T>
const std::string FloatToString(const T value, const int n = 15) {
  std::ostringstream out;
  out.precision(n);
  out << std::scientific << value;
  return out.str();
}

void EvaluationString(std::string& str);
bool GetFrontNumber(std::string& str, std::string& number);
void EvaluationStringFloat(std::string& str);
void EvaluationStringInt(std::string& str);

template <typename T>
T Bioperation(T a, T b, char oper) {
  if (oper == '^')
    return static_cast<T>(std::pow(a, b));
  else if (oper == '*')
    return a * b;
  else if (oper == '/')
    return a / b;
  else if (oper == '+')
    return a + b;
  else if (oper == '-')
    return a - b;
  return -1;
}
}  // namespace avocado