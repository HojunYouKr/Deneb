#include "avocado_string.h"

#include <sstream>
#include <algorithm>
#include <ctime>
#include <map>
#include <stack>


namespace avocado {
std::string GetTitle(const std::string& str) {
  const int size = 80;
  std::string title;
  title.push_back('\n');
  title.append(std::string(size, '='));
  title.push_back('\n');
  const int size_str = static_cast<int>(str.length());
  const int offset = (size - size_str) / 2;
  title.append(std::string(offset, ' '));
  title.append(str);
  title.push_back('\n');
  title.append(std::string(size, '='));
  title.push_back('\n');
  return title;
}
std::string GetLocalTime(void) {
  const std::time_t timer = std::time(nullptr);
  const struct tm* t = std::localtime(&timer);
  const std::vector<std::string> month = {"Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};

  const int curr_year = t->tm_year + 1900;
  const std::string& curr_month = month[t->tm_mon];
  const int curr_day = t->tm_mday;
  const int curr_hour = t->tm_hour;
  const int curr_minute = t->tm_min;
  const int curr_second = t->tm_sec;

  std::stringstream now;
  now << curr_month << "  " << curr_day << " " << curr_year << " " << curr_hour
     << ":" << curr_minute << ":" << curr_second;
  return std::string(now.str());
}
void ReplaceString(std::string& str, const std::string& find_tok,
                   const std::string& replace_tok) {
  std::string::size_type offset = 0;
  while (true) {
    offset = str.find(find_tok);
    if (offset == std::string::npos)
      break;
    else
      str.replace(offset, find_tok.length(), replace_tok);
  }
}

std::vector<std::string> SplitString(const std::string& str,
                                     const char delimiter) {
  std::vector<std::string> toks;
  std::stringstream ss(str);
  std::string tok;
  while (std::getline(ss, tok, delimiter)) toks.push_back(tok);
  if (str.back() == delimiter) toks.push_back("");
  return toks;
}

void RemoveString(std::string& str, const std::vector<char>& drops) {
  for (auto&& it = drops.begin(); it < drops.end(); it++)
    str.erase(std::remove(str.begin(), str.end(), *it), str.end());
}

void RemoveString(std::string& str, const std::string& remove_tok) {
  ReplaceString(str, remove_tok, "");
}

void EvaluationString(std::string& str) {
  if (str.find('.') == std::string::npos)
    EvaluationStringInt(str);
  else
    EvaluationStringFloat(str);
}

bool GetFrontNumber(std::string& str, std::string& number) {
  number = "";
  bool point = false;
  bool exp = false;
  for (size_t i = 0, e = str.length(); i < e; i++) {
    const char& ch = str.at(i);
    if (i == 0 && ch == '-') {
      number.push_back(ch);
    } else if ('0' <= ch && ch <= '9') {
      number.push_back(ch);
    } else if (point == false && exp == false && ch == '.') {
      point = true;
      number.push_back(ch);
    } else if (exp == false && (ch == 'e' || ch == 'E')) {
      exp = true;
      number.push_back(ch);
      const char& exp_sign = str.at(i + 1);
      if (exp_sign == '-' || exp_sign == '+') {
        number.push_back(exp_sign);
        i++;
      }
    } else {
      break;
    }
  }
  const size_t size = number.length();
  if (size > 0) {
    str.erase(0, size);
    return true;
  }
  return false;
}

void EvaluationStringFloat(std::string& str) {
  const std::map<char, int> priority = {
      {'^', 3}, {'*', 2}, {'/', 2}, {'+', 1}, {'-', 1}};
  struct oper {
    int p;
    char o;
  };
  std::stack<double> numbers;
  std::stack<oper> operations;

  str.push_back('\n');
  std::string tok;
  while (true) {
    const bool isnumber = GetFrontNumber(str, tok);
    const char oper = str.at(0);
    str.erase(0, 1);

    if (isnumber) numbers.push(std::stod(tok));
    if (oper == '\n') break;

    auto&& iterator = priority.find(oper);
    if (oper == '(') {
      operations.push({0, oper});
    } else if (oper == ')') {
      while (operations.top().o != '(') {
        const double b = numbers.top();
        numbers.pop();
        const double a = numbers.top();
        numbers.pop();
        const char oper = operations.top().o;
        operations.pop();
        numbers.push(Bioperation(a, b, oper));
      }
      operations.pop();
    } else if (iterator != priority.end()) {
      const int prior = iterator->second;
      while (!operations.empty() && prior <= operations.top().p) {
        const double b = numbers.top();
        numbers.pop();
        const double a = numbers.top();
        numbers.pop();
        const char oper = operations.top().o;
        operations.pop();
        numbers.push(Bioperation(a, b, oper));
      }
      operations.push({prior, oper});
    }
  }
  while (!operations.empty()) {
    const double b = numbers.top();
    numbers.pop();
    const double a = numbers.top();
    numbers.pop();
    const char oper = operations.top().o;
    operations.pop();
    numbers.push(Bioperation(a, b, oper));
  }
  str = FloatToString(numbers.top());
}

void EvaluationStringInt(std::string& str) {
  const std::map<char, int> priority = {
      {'^', 3}, {'*', 2}, {'/', 2}, {'+', 1}, {'-', 1}};
  struct oper {
    int p;
    char o;
  };
  std::stack<int> numbers;
  std::stack<oper> operations;

  str.push_back('\n');
  std::string tok;
  while (true) {
    const bool isnumber = GetFrontNumber(str, tok);
    const char oper = str.at(0);
    str.erase(0, 1);

    if (isnumber) numbers.push(std::stoi(tok));
    if (oper == '\n') break;

    auto&& iterator = priority.find(oper);
    if (oper == '(') {
      operations.push({0, oper});
    } else if (oper == ')') {
      while (operations.top().o != '(') {
        const int b = numbers.top();
        numbers.pop();
        const int a = numbers.top();
        numbers.pop();
        const char oper = operations.top().o;
        operations.pop();
        numbers.push(Bioperation(a, b, oper));
      }
      operations.pop();
    } else if (iterator != priority.end()) {
      const int prior = iterator->second;
      while (!operations.empty() && prior <= operations.top().p) {
        const int b = numbers.top();
        numbers.pop();
        const int a = numbers.top();
        numbers.pop();
        const char oper = operations.top().o;
        operations.pop();
        numbers.push(Bioperation(a, b, oper));
      }
      operations.push({prior, oper});
    }
  }
  while (!operations.empty()) {
    const int b = numbers.top();
    numbers.pop();
    const int a = numbers.top();
    numbers.pop();
    const char oper = operations.top().o;
    operations.pop();
    numbers.push(Bioperation(a, b, oper));
  }
  str = std::to_string(numbers.top());
}
}  // namespace avocado