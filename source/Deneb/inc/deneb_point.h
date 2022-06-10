#pragma once

#include <functional>
#include <string>
#include <unordered_map>
#include <vector>
#include <memory>

#define DENEB_POINT_NAME point_global_ptr
#define DENEB_POINT deneb::DENEB_POINT_NAME
#define DENEB_POINT_INITIALIZE(directory) \
  DENEB_POINT = std::make_shared<deneb::Point>(directory)
#define DENEB_POINT_FINALIZE() DENEB_POINT.reset()

namespace deneb {
class Point {
 public:
  enum class Option : bool { DEGREE, NUM_POINTS };
  enum class Type : int {
    LINE_GL,    // Gauss-Legendre (Quadrature)
    LINE_GLL,   // Gauss-Legendre-Lobatto
    TRIS_WV,    // Witherden-Vincent (Quadrature)
    TRIS_ALPHA, // Alpha-optimized
    QUAD_WV,    // Witherden-Vincent (Quadrature)
    TETS_WV,    // Witherden-Vincent (Quadrature)
    TETS_ALPHA, // Alpha-optimized
    PYRA_SFP    // Shape function points
  };

 private:
  struct Tag {
    Type type;
    int degree;
  };
  friend bool operator==(const Tag& left, const Tag& right);
  friend size_t Hash_Point_Tag(const Tag& tag);
  struct Rule {
    std::vector<double> points;
    std::vector<double> weights;
  };

 private:
  std::string directory_;
  std::unordered_map<Tag, Rule, std::function<size_t(const Tag&)>> registry_;
  bool set_regist_;

 public:
  Point(const std::string& directory);
  ~Point();

  bool GetPoint(const Type type, const Option opt, const int opt_value,
                std::vector<double>& pts, std::vector<double>& wts);

 private:
  inline void RegisterPoint(const Tag& tag, const Rule& rule) {
    registry_[tag] = rule;
  };
  bool GetPoint(const Tag& tag, Rule& rule) const;
};
extern std::shared_ptr<Point> DENEB_POINT_NAME;
bool operator==(const Point::Tag& left, const Point::Tag& right);
size_t Hash_Point_Tag(const Point::Tag& tag);
}  // namespace deneb