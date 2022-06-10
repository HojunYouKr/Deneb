#include "deneb_point.h"

#include <fstream>

#include "avocado.h"

namespace deneb {
std::shared_ptr<Point> DENEB_POINT_NAME = nullptr;
Point::Point(const std::string& directory)
    : directory_(directory),
      registry_(100, Hash_Point_Tag),
      set_regist_(true) {
};
Point::~Point() {
  registry_.clear();
};
bool Point::GetPoint(const Type type, const Option opt, const int opt_value,
              std::vector<double>& pts, std::vector<double>& wts) {
  pts.clear();
  wts.clear();
  bool is_wts = true; // default
  int dimension;
  std::vector<int> support_deg;
  std::vector<int> support_npt;
  std::function<const std::string(const int, const int)> filepath;
  switch (type) {
    case Type::LINE_GL:
      dimension = 1;
      for (int ipt = 1; ipt <= 200; ipt++) {
        support_deg.push_back(ipt * 2 - 1);
        support_npt.push_back(ipt);
      }
      filepath = [](const int deg, const int npt) -> const std::string {
        return "Points/line/GaussLegendre_n" + std::to_string(npt) +".dat";
      };
      break;
    case Type::LINE_GLL:
      dimension = 1;
      for (int ipt = 2; ipt <= 200; ipt++) {
        support_deg.push_back(ipt * 2 - 3);
        support_npt.push_back(ipt);
      }
      filepath = [](const int deg, const int npt) -> const std::string {
        return "Points/line/GaussLegendreLobatto_n" + std::to_string(npt) + ".dat";
      };
      break;
    case Type::TRIS_WV:
      dimension = 2;
      support_deg = {1,  2,  4,  5,  6,  7,  8,  9,  10, 11,
                     12, 13, 14, 15, 16, 17, 18, 19, 20};
      support_npt = {1,  3,  6,  7,  12, 15, 16, 19, 25, 28,
                     33, 37, 42, 49, 55, 60, 67, 73, 79};
      filepath = [](const int deg, const int npt) -> const std::string {
        return "Points/tris/WitherdenVincent_n" + std::to_string(npt) +
               "_d" + std::to_string(deg) + ".dat";
      };
      break;
    case Type::TRIS_ALPHA:
      is_wts = false;
      dimension = 2;
      support_deg = {0, 1,  2,  3,  4,  5,  6,  7,  8,
                     9, 10, 11, 12, 13, 14, 15, 16, 17};
      support_npt = {1,  2,  6,  10, 15,  21,  28,  36,  45,
                     55, 66, 78, 91, 105, 120, 136, 153, 171};
      filepath = [](const int deg, const int npt) -> const std::string {
        return "Points/tris/AlphaOptimized_n" + std::to_string(npt) + ".dat";
      };
      break;
    case Type::QUAD_WV:
      dimension = 2;
      support_deg = {1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21};
      support_npt = {1, 4, 8, 12, 20, 28, 37, 48, 60, 72, 85};
      filepath = [](const int deg, const int npt) -> const std::string {
        return "Points/quad/WitherdenVincent_n" + std::to_string(npt) +
               "_d" + std::to_string(deg) + ".dat";
      };
      break;
    case Type::TETS_WV:
      dimension = 3;
      support_deg = {1, 2, 3, 5, 6, 7, 8, 9, 10};
      support_npt = {1, 4, 8, 14, 24, 35, 46, 59, 81};
      filepath = [](const int deg, const int npt) -> const std::string {
        return "Points/tets/WitherdenVincent_n" + std::to_string(npt) +
               "_d" + std::to_string(deg) + ".dat";
      };
      break;
    case Type::TETS_ALPHA:
      is_wts = false;
      dimension = 3;
      support_deg = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
      support_npt = {1, 4, 10, 20, 35, 56, 84, 120, 165, 220};
      filepath = [](const int deg, const int npt) -> const std::string {
        return "Points/tets/AlphaOptimized_n" + std::to_string(npt) + ".dat";
      };
      break;
    case Type::PYRA_SFP:
      is_wts = false;
      dimension = 3;
      support_deg = {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                     10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
      support_npt = {1,   5,   14,  30,   55,   91,   140,  204,  285,  385,
                     506, 650, 819, 1015, 1240, 1496, 1785, 2109, 2470, 2870};
      filepath = [](const int deg, const int npt) -> const std::string {
        return "Points/pyra/ShapeFunctionPoints_n" + std::to_string(npt) + "_d" +
               std::to_string(deg) + ".dat";
      };
      break;
    default:
      ERROR_MESSAGE("Invalid point type.\n");
  }

  const std::vector<int>* support_value = &support_deg;
  if (opt == Option::NUM_POINTS) support_value = &support_npt;
  if (opt_value > (*support_value).back()) return false;
  for (size_t i = 0, len = (*support_value).size(); i < len; i++) {
    if ((*support_value)[i] >= opt_value) {
      const int& deg = support_deg[i];
      const int& npt = support_npt[i];
      const Tag tag{type, deg};
      Rule rule;
      if (!GetPoint(tag, rule)) {
        const std::string filename = directory_ + filepath(deg, npt);
        std::ifstream infile(filename.c_str());
        if (!infile.is_open())
          ERROR_MESSAGE("Invalid filename (no-exist):" + filename + "\n");

        rule.points.resize(npt * dimension);
        if (is_wts) rule.weights.resize(npt);
        for (int ipt = 0; ipt < npt; ipt++) {
          for (int idim = 0; idim < dimension; idim++)
            infile >> rule.points[ipt * dimension + idim];
          if (is_wts) infile >> rule.weights[ipt];
        }

        infile.close();
      }
      if (set_regist_) RegisterPoint(tag, rule);
      pts = std::move(rule.points);
      if (is_wts) wts = std::move(rule.weights);
      break;
    }
  }
  return true;
}
bool Point::GetPoint(const Tag& tag, Rule& rule) const {
  auto&& iterator = registry_.find(tag);
  if (iterator == registry_.end()) return false;
  rule = iterator->second;
  return true;
}

bool operator==(const Point::Tag& left, const Point::Tag& right) {
  return left.type == right.type && left.degree == right.degree;
};
size_t Hash_Point_Tag(const Point::Tag& tag) {
  return std::hash<int>()(static_cast<int>(tag.type)) ^
         std::hash<int>()(tag.degree) << 1;
};
}  // namespace deneb