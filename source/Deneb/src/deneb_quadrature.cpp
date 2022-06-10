#include "deneb_quadrature.h"

#include "avocado.h"
#include "deneb_point.h"

namespace deneb {
namespace quadrature {
void Line_Poly1D(const int degree, std::vector<double>& points,
                 std::vector<double>& weights) {
  points.clear();
  weights.clear();
  if (!DENEB_POINT->GetPoint(Point::Type::LINE_GL, Point::Option::DEGREE, degree,
                       points, weights))
    ERROR_MESSAGE("Fail to achieve Line_Poly1D quadrature (degree=" +
                  std::to_string(degree) + ")\n");
}
void Tris_Poly2D(const int degree, std::vector<double>& points,
                 std::vector<double>& weights) {
  points.clear();
  weights.clear();
  if (!DENEB_POINT->GetPoint(Point::Type::TRIS_WV, Point::Option::DEGREE, degree,
                       points, weights)) {
    std::vector<double> line_points, line_weights;
    if (!DENEB_POINT->GetPoint(Point::Type::LINE_GL, Point::Option::DEGREE,
                         degree + 1, line_points, line_weights))
      ERROR_MESSAGE("Fail to achieve Tris_Poly2D quadrature (degree=" +
                    std::to_string(degree) + ")\n");
    const size_t npt_line = line_weights.size();
    points.reserve(npt_line * npt_line * 2);
    weights.reserve(npt_line * npt_line);
    for (size_t i = 0; i < npt_line; i++) {
      for (size_t j = 0; j < npt_line; j++) {
        const double& r = line_points[i];
        const double& s = line_points[j];
        points.push_back(0.5 * (r * (1.0 - s) - (1.0 + s)));
        points.push_back(s);
        weights.push_back(line_weights[i] * line_weights[j] * (1.0 - s) * 0.5);
      }
    }
  }
}
void Quad_QuadShape(const int degree, std::vector<double>& points,
                    std::vector<double>& weights) {
  points.clear();
  weights.clear();
  std::vector<double> line_points, line_weights;
  if (!DENEB_POINT->GetPoint(Point::Type::LINE_GL, Point::Option::DEGREE, degree,
                       line_points, line_weights))
    ERROR_MESSAGE("Fail to achieve Quad_QuadShape quadrature (degree=" +
                  std::to_string(degree) + ")\n");
  const size_t npt_line = line_weights.size();
  points.reserve(npt_line * npt_line * 2);
  weights.reserve(npt_line * npt_line);
  for (size_t i = 0; i < npt_line; i++) {
    for (size_t j = 0; j < npt_line; j++) {
      points.push_back(line_points[i]);
      points.push_back(line_points[j]);
      weights.push_back(line_weights[i] * line_weights[j]);
    }
  }
}
void Tets_Poly3D(const int degree, std::vector<double>& points,
                 std::vector<double>& weights) {
  points.clear();
  weights.clear();
  if (!DENEB_POINT->GetPoint(Point::Type::TETS_WV, Point::Option::DEGREE, degree,
                       points, weights)) {
    std::vector<double> line_points, line_weights;
    if (!DENEB_POINT->GetPoint(Point::Type::LINE_GL, Point::Option::DEGREE,
                         degree + 2, line_points, line_weights))
      ERROR_MESSAGE("Fail to achieve Tets_Poly3D quadrature (degree=" +
                    std::to_string(degree) + ")\n");
    std::vector<double> tris_points, tris_weights;
    Tris_Poly2D(degree, tris_points, tris_weights);

    const size_t npt_line = line_weights.size();
    const size_t npt_tris = tris_weights.size();
    points.reserve(npt_line * npt_tris * 3);
    weights.reserve(npt_line * npt_tris);
    for (size_t i = 0; i < npt_line; i++) {
      for (size_t j = 0; j < npt_tris; j++) {
        const double& r = tris_points[2 * j];
        const double& s = tris_points[2 * j + 1];
        const double& t = line_points[i];
        points.push_back(0.5 * (r * (1.0 - t) - (1.0 + t)));
        points.push_back(0.5 * (s * (1.0 - t) - (1.0 + t)));
        points.push_back(t);
        weights.push_back(line_weights[i] * tris_weights[j] * (1.0 - t) *
                          (1.0 - t) * 0.25);
      }
    }
  }
}
void Hexa_HexaShape(const int degree, std::vector<double>& points,
                    std::vector<double>& weights) {
  points.clear();
  weights.clear();
  std::vector<double> line_points, line_weights;
  if (!DENEB_POINT->GetPoint(Point::Type::LINE_GL, Point::Option::DEGREE, degree,
                       line_points, line_weights))
    ERROR_MESSAGE("Fail to achieve Hexa_HexaShape quadrature (degree=" +
                  std::to_string(degree) + ")\n");
  const size_t npt_line = line_weights.size();
  points.reserve(npt_line * npt_line * npt_line * 3);
  weights.reserve(npt_line * npt_line * npt_line);
  for (size_t i = 0; i < npt_line; i++) {
    for (size_t j = 0; j < npt_line; j++) {
      for (size_t k = 0; k < npt_line; k++) {
        points.push_back(line_points[i]);
        points.push_back(line_points[j]);
        points.push_back(line_points[k]);
        weights.push_back(line_weights[i] * line_weights[j] * line_weights[k]);
      }
    }
  }
}
void Pris_PrisShape(const int degree_tris, const int degree_line,
                    std::vector<double>& points, std::vector<double>& weights) {
  points.clear();
  weights.clear();
  std::vector<double> line_points, line_weights;
  if (!DENEB_POINT->GetPoint(Point::Type::LINE_GL, Point::Option::DEGREE, degree_line,
                       line_points, line_weights))
    ERROR_MESSAGE("Fail to achieve Pris_PrisShape quadrature (degree=" +
                  std::to_string(degree_tris) + "," +
                  std::to_string(degree_line) + ")\n");
  std::vector<double> tris_points, tris_weights;
  Tris_Poly2D(degree_tris, tris_points, tris_weights);

  const size_t npt_line = line_weights.size();
  const size_t npt_tris = tris_weights.size();
  points.reserve(npt_line * npt_tris * 3);
  weights.reserve(npt_line * npt_tris);
  for (size_t i = 0; i < npt_line; i++) {
    for (size_t j = 0; j < npt_tris; j++) {
      points.push_back(tris_points[2 * j]);
      points.push_back(tris_points[2 * j + 1]);
      points.push_back(line_points[i]);
      weights.push_back(tris_weights[j] * line_weights[i]);
    }
  }
}
void Pyra_PyraShape(const int degree, std::vector<double>& points,
                    std::vector<double>& weights) {
  points.clear();
  weights.clear();
  std::vector<double> line_points, line_weights;
  if (!DENEB_POINT->GetPoint(Point::Type::LINE_GL, Point::Option::DEGREE, degree,
                       line_points, line_weights))
    ERROR_MESSAGE("Fail to achieve Pyra_PyraShape quadrature (degree=" +
                  std::to_string(degree) + ")\n");
  const size_t npt_line = line_weights.size();
  points.reserve(npt_line * npt_line * npt_line * 3);
  weights.reserve(npt_line * npt_line * npt_line);
  for (size_t i = 0; i < npt_line; i++) {
    for (size_t j = 0; j < npt_line; j++) {
      for (size_t k = 0; k < npt_line; k++) {
        const double& r = line_points[i];
        const double& s = line_points[j];
        const double& t = line_points[k];
        points.push_back(0.5 * r * (1.0 - t));
        points.push_back(0.5 * s * (1.0 - t));
        points.push_back(t);
        weights.push_back(line_weights[i] * line_weights[j] * line_weights[k] *
                          0.25 * (1.0 - t) * (1.0 - t));
      }
    }
  }
}
}  // namespace quadrature
}  // namespace deneb