#include "deneb_DRM.h"

#include "avocado.h"
#include "deneb_basis.h"
#include "deneb_point.h"

namespace deneb {
std::shared_ptr<DRMVolume> DRMVolume::GetDRM(const ElemType elemtype) {
  switch (elemtype) {
    case ElemType::TRIS:
      return std::make_shared<DRMVolumeTris>();
    case ElemType::QUAD:
      return std::make_shared<DRMVolumeQuad>();
    case ElemType::TETS:
      return std::make_shared<DRMVolumeTets>();
    case ElemType::HEXA:
      return std::make_shared<DRMVolumeHexa>();
    case ElemType::PRIS:
      return std::make_shared<DRMVolumePris>();
    case ElemType::PYRA:
      return std::make_shared<DRMVolumePyra>();
  }
  ERROR_MESSAGE("Invalid DRMVolume element yype: " +
                std::to_string(static_cast<int>(elemtype)) + "\n");
  return nullptr;
}

// DRMVolumeTris: Poly2D, alpha-optimized points
DRMVolumeTris::DRMVolumeTris() {
  approx_basis_ =
      Basis::GetAllocatedBasis(ElemType::TRIS, std::make_shared<Poly2D>());
}
void DRMVolumeTris::GetDRMPoints(const int num_points,
                                 std::vector<double>& points) const {
  std::vector<double> weights;
  if (!DENEB_POINT->GetPoint(Point::Type::TRIS_ALPHA, Point::Option::NUM_POINTS,
                             num_points, points, weights))
    ERROR_MESSAGE("Fail to achieve DRM points for DRMVolumeTris."
                  "\n\tNum points: " + std::to_string(num_points) + "\n");
}

// DRMVolumeQuad: PolyQuadShape, Gauss-Legendre points
DRMVolumeQuad::DRMVolumeQuad() {
  approx_basis_ =
      Basis::GetAllocatedBasis(ElemType::QUAD, std::make_shared<PolyQuadShape>());
}
void DRMVolumeQuad::GetDRMPoints(const int num_points,
                                 std::vector<double>& points) const {
  std::vector<double> line_weights, line_points;
  const int np = static_cast<int>(
      std::ceil(std::sqrt(static_cast<double>(num_points) - 0.5)) + 0.5);
  if (!DENEB_POINT->GetPoint(Point::Type::LINE_GL, Point::Option::NUM_POINTS,
                             np, line_points, line_weights))
    ERROR_MESSAGE(
        "Fail to achieve DRM points for DRMVolumeQuad."
        "\n\tNum points: " +
        std::to_string(num_points) + "\n");

  const size_t npt_line = line_weights.size();
  points.clear();
  points.reserve(npt_line * npt_line * 2);
  for (size_t i = 0; i < npt_line; i++) {
    for (size_t j = 0; j < npt_line; j++) {
      points.push_back(line_points[i]);
      points.push_back(line_points[j]);
    }
  }
}

// DRMVolumeTets: Poly3D, alpha-optimized points
DRMVolumeTets::DRMVolumeTets() {
  approx_basis_ = Basis::GetAllocatedBasis(ElemType::TETS,
                                           std::make_shared<Poly3D>());
}
void DRMVolumeTets::GetDRMPoints(const int num_points,
                                 std::vector<double>& points) const {
  std::vector<double> weights;
  if (!DENEB_POINT->GetPoint(Point::Type::TETS_ALPHA, Point::Option::NUM_POINTS,
                             num_points, points, weights))
    ERROR_MESSAGE(
        "Fail to achieve DRM points for DRMVolumeTets."
        "\n\tNum points: " +
        std::to_string(num_points) + "\n");
}

// DRMVolumeHexa: PolyHexaShape, Gauss-Legendre points
DRMVolumeHexa::DRMVolumeHexa() {
  approx_basis_ = Basis::GetAllocatedBasis(ElemType::HEXA,
                                           std::make_shared<PolyHexaShape>());
}
void DRMVolumeHexa::GetDRMPoints(const int num_points,
                                 std::vector<double>& points) const {
  std::vector<double> line_weights, line_points;
  const int np = static_cast<int>(
      std::ceil(std::pow(static_cast<double>(num_points) - 0.5, 1.0 / 3.0)) +
      0.5);
  if (!DENEB_POINT->GetPoint(Point::Type::LINE_GL, Point::Option::NUM_POINTS,
                             np, line_points, line_weights))
    ERROR_MESSAGE(
        "Fail to achieve DRM points for DRMVolumeHexa."
        "\n\tNum points: " +
        std::to_string(num_points) + "\n");

  const size_t npt_line = line_weights.size();
  points.clear();
  points.reserve(npt_line * npt_line * npt_line * 3);
  for (size_t i = 0; i < npt_line; i++) {
    for (size_t j = 0; j < npt_line; j++) {
      for (size_t k = 0; k < npt_line; k++) {
        points.push_back(line_points[i]);
        points.push_back(line_points[j]);
        points.push_back(line_points[k]);
      }
    }
  }
}

// DRMVolumePris: PolyPrisShape, Gauss-Legendre points * alpha-optimized points
DRMVolumePris::DRMVolumePris() {
  approx_basis_ = Basis::GetAllocatedBasis(ElemType::PRIS,
                                           std::make_shared<PolyPrisShape>());
}
void DRMVolumePris::GetDRMPoints(const int num_points,
                                 std::vector<double>& points) const {
  std::vector<double> line_weights, line_points;
  std::vector<double> tris_weights, tris_points;
  int np = 0;
  for (np = 1; np < 100; np++)
    if (np * np * (np + 1) / 2 >= num_points) break;
  if (!DENEB_POINT->GetPoint(Point::Type::LINE_GL, Point::Option::NUM_POINTS,
                             np, line_points, line_weights) ||
      !DENEB_POINT->GetPoint(Point::Type::TRIS_ALPHA, Point::Option::NUM_POINTS,
                             np * (np + 1) / 2, tris_points, tris_weights))
    ERROR_MESSAGE(
        "Fail to achieve DRM points for DRMVolumePris."
        "\n\tNum points: " +
        std::to_string(num_points) + "\n");

  const size_t npt_line = line_weights.size();
  const size_t npt_tris = tris_weights.size();
  points.clear();
  points.reserve(npt_line * npt_tris * 3);
  for (size_t i = 0; i < npt_line; i++) {
    for (size_t j = 0; j < npt_tris; j++) {
      points.push_back(tris_points[2 * j]);
      points.push_back(tris_points[2 * j + 1]);
      points.push_back(line_points[i]);
    }
  }
}

// DRMVolumePyra: PolyPyra, pyramid DRM-SFP points
DRMVolumePyra::DRMVolumePyra() {
  approx_basis_ =
      Basis::GetAllocatedBasis(ElemType::PYRA, std::make_shared<PolyPyra>());
}
void DRMVolumePyra::GetDRMPoints(const int num_points,
                                 std::vector<double>& points) const {
  std::vector<double> weights;
  if (!DENEB_POINT->GetPoint(Point::Type::PYRA_SFP, Point::Option::NUM_POINTS,
                             num_points, points, weights))
    ERROR_MESSAGE(
        "Fail to achieve DRM points for DRMVolumePyra."
        "\n\tNum points: " +
        std::to_string(num_points) + "\n");
}

std::shared_ptr<DRMSurface> DRMSurface::GetDRM(const ElemType elemtype) {
  switch (elemtype) {
    case ElemType::LINE:
      return std::make_shared<DRMSurfaceLine>();
    case ElemType::TRIS:
      return std::make_shared<DRMSurfaceTris>();
    case ElemType::QUAD:
      return std::make_shared<DRMSurfaceQuad>();
  }
  ERROR_MESSAGE("Invalid DRMSurface element yype: " +
                std::to_string(static_cast<int>(elemtype)) + "\n");
  return nullptr;
}

// DRMSurfaceLine: Poly2D, Gauss-Legendre points
DRMSurfaceLine::DRMSurfaceLine() {
  approx_basis_ =
      BasisSurface::GetAllocatedBasis(ElemType::LINE, std::make_shared<Poly2DSorted>());
}
void DRMSurfaceLine::GetDRMPoints(const int num_points,
                                 std::vector<double>& points) const {
  std::vector<double> weights;
  if (!DENEB_POINT->GetPoint(Point::Type::LINE_GL, Point::Option::NUM_POINTS,
                             num_points, points, weights))
    ERROR_MESSAGE(
        "Fail to achieve DRM points for DRMSurfaceLine."
        "\n\tNum points: " +
        std::to_string(num_points) + "\n");
}
// DRMSurfaceTris: Poly3D, Witherden-Vincent points
DRMSurfaceTris::DRMSurfaceTris() {
  approx_basis_ = BasisSurface::GetAllocatedBasis(ElemType::TRIS, std::make_shared<Poly3DSorted>());
}
void DRMSurfaceTris::GetDRMPoints(const int num_points,
                                  std::vector<double>& points) const {
  std::vector<double> weights;
  if (!DENEB_POINT->GetPoint(Point::Type::TRIS_WV, Point::Option::NUM_POINTS,
                             num_points, points, weights))
    ERROR_MESSAGE(
        "Fail to achieve DRM points for DRMSurfaceTris."
        "\n\tNum points: " +
        std::to_string(num_points) + "\n");
}
// DRMSurfaceQuad: Poly3D, Witherden-Vincent points
DRMSurfaceQuad::DRMSurfaceQuad() {
  approx_basis_ = BasisSurface::GetAllocatedBasis(ElemType::QUAD, std::make_shared<Poly3DSorted>());
}
void DRMSurfaceQuad::GetDRMPoints(const int num_points,
                                  std::vector<double>& points) const {
  std::vector<double> weights;
  if (!DENEB_POINT->GetPoint(Point::Type::QUAD_WV, Point::Option::NUM_POINTS,
                             num_points, points, weights))
    ERROR_MESSAGE(
        "Fail to achieve DRM points for DRMSurfaceQuad."
        "\n\tNum points: " +
        std::to_string(num_points) + "\n");
}
}  // namespace deneb