#include "deneb_polynomial.h"

#include <vector>

#include "avocado.h"
#include "deneb_quadrature.h"

namespace deneb {
// Polynomial
int Polynomial::NumTerms(const int order) const {
  ERROR_MESSAGE(poly_ + " has no NumTerms method.\n");
  return 0;
}
void Polynomial::F(const int order, const double* coord, double* result) const {
  ERROR_MESSAGE(poly_ + " has no F method.\n");
}
void Polynomial::Fdr(const int order, const double* coord,
                     double* result) const {
  ERROR_MESSAGE(poly_ + " has no Fdr method.\n");
}
void Polynomial::Fds(const int order, const double* coord,
                     double* result) const {
  ERROR_MESSAGE(poly_ + " has no Fds method.\n");
}
void Polynomial::Fdt(const int order, const double* coord,
                     double* result) const {
  ERROR_MESSAGE(poly_ + " has no Fdt method.\n");
}
void Polynomial::Fdrr(const int order, const double* coord,
                      double* result) const {
  ERROR_MESSAGE(poly_ + " has no Fdrr method.\n");
}
void Polynomial::Fdrs(const int order, const double* coord,
                      double* result) const {
  ERROR_MESSAGE(poly_ + " has no Fdrs method.\n");
}
void Polynomial::Fdss(const int order, const double* coord,
                      double* result) const {
  ERROR_MESSAGE(poly_ + " has no Fdss method.\n");
}
void Polynomial::GetVolumeQuadrature(const int order, const ElemType elemtype,
                                     const int elemorder,
                                     std::vector<double>& points,
                                     std::vector<double>& weights) const {
  ERROR_MESSAGE(poly_ + " has no GetVolumeQuadrature method.\n");
}
void Polynomial::GetSurfaceQuadrature(const int order, const ElemType elemtype,
                                  const int elemorder,
                                  std::vector<double>& points,
                                      std::vector<double>& weights) const {
  ERROR_MESSAGE(poly_ + " has no GetSurfaceQuadrature method.\n");
}

// Poly1D: 1-D standard polynomial of r
int Poly1D::NumTerms(const int order) const { return order + 1; }
void Poly1D::F(const int order, const double* coord, double* result) const {
  const double& r = coord[0];
  result[0] = 1.0;
  for (int i = 1; i < order + 1; i++) result[i] = result[i - 1] * r;
}
void Poly1D::Fdr(const int order, const double* coord, double* result) const {
  const double& r = coord[0];
  result[0] = 0.0;
  if (order == 0) return;
  result[1] = 1.0;
  for (int i = 2; i < order + 1; i++)
    result[i] = double(i) / double(i - 1) * result[i - 1] * r;
}
void Poly1D::Fdrr(const int order, const double* coord, double* result) const {
  const double& r = coord[0];
  result[0] = 0.0;
  if (order == 0) return;
  result[1] = 0.0;
  if (order == 1) return;
  result[2] = 2.0;
  for (int i = 3; i < order + 1; i++)
    result[i] = double(i) / double(i - 2) * result[i - 1] * r;
}

// Poly2D: 2-D standard polynomial of r, s
int Poly2D::NumTerms(const int order) const {
  return (order + 1) * (order + 2) / 2;
}
void Poly2D::F(const int order, const double* coord, double* result) const {
  const double& r = coord[0];
  const double& s = coord[1];
  result[0] = 1.0;
  int index = 1;
  for (int i = 1; i <= order; i++) {
    result[index] = result[index - i] * r;
    index++;
    for (int j = 0; j < i; j++) {
      result[index] = result[index - i - 1] * s;
      index++;
    }
  }
}
void Poly2D::Fdr(const int order, const double* coord, double* result) const {
  const double& s = coord[1];
  const Poly1D poly1d;
  std::vector<double> poly1d_dr(poly1d.NumTerms(order));
  poly1d.Fdr(order, coord, &poly1d_dr[0]);
  int index = 0;
  for (int i = 0; i <= order; i++) {
    result[index] = poly1d_dr[i];
    index++;
    for (int j = 0; j < i; j++) {
      result[index] = result[index - i - 1] * s;
      index++;
    }
  }
}
void Poly2D::Fds(const int order, const double* coord, double* result) const {
  const double& r = coord[0];
  const Poly1D poly1d;
  std::vector<double> poly1d_ds(poly1d.NumTerms(order));
  poly1d.Fdr(order, coord + 1, &poly1d_ds[0]);
  int index = 0;
  for (int i = 0; i <= order; i++) {
    for (int j = 0; j < i; j++) {
      result[index] = result[index - i] * r;
      index++;
    }
    result[index] = poly1d_ds[i];
    index++;
  }
}
void Poly2D::Fdrr(const int order, const double* coord, double* result) const {
  const double& s = coord[1];
  const Poly1D poly1d;
  std::vector<double> poly1d_drr(poly1d.NumTerms(order));
  poly1d.Fdrr(order, coord, &poly1d_drr[0]);
  int index = 0;
  for (int i = 0; i <= order; i++) {
    result[index] = poly1d_drr[i];
    index++;
    for (int j = 0; j < i; j++) {
      result[index] = result[index - i - 1] * s;
      index++;
    }
  }
}
void Poly2D::Fdrs(const int order, const double* coord, double* result) const {
  const Poly1D poly1d;
  std::vector<double> poly1d_dr(poly1d.NumTerms(order));
  std::vector<double> poly1d_ds(poly1d.NumTerms(order));
  poly1d.Fdr(order, coord, &poly1d_dr[0]);
  poly1d.Fdr(order, coord + 1, &poly1d_ds[0]);
  int index = 0;
  for (int i = 0; i <= order; i++)
    for (int j = 0; j <= i; j++)
      result[index++] = poly1d_dr[i - j] * poly1d_ds[j];
}
void Poly2D::Fdss(const int order, const double* coord, double* result) const {
  const double& r = coord[0];
  const Poly1D poly1d;
  std::vector<double> poly1d_dss(poly1d.NumTerms(order));
  poly1d.Fdrr(order, coord + 1, &poly1d_dss[0]);
  int index = 0;
  for (int i = 0; i <= order; i++) {
    for (int j = 0; j < i; j++) {
      result[index] = result[index - i] * r;
      index++;
    }
    result[index] = poly1d_dss[i];
    index++;
  }
}
void Poly2D::GetVolumeQuadrature(const int order, const ElemType elemtype,
                                 const int elemorder,
                                 std::vector<double>& points,
                                 std::vector<double>& weights) const {
  int deg = 0;
  switch (elemtype) {
    case ElemType::TRIS:
      deg = order * elemorder + 2 * (elemorder - 1);
      quadrature::Tris_Poly2D(deg, points, weights);
      return;
    case ElemType::QUAD:
      deg = order * elemorder + 2 * elemorder - 1;
      quadrature::Quad_QuadShape(deg, points, weights);
      return;
    default:
      ERROR_MESSAGE(poly_ + " has no GetVolumeQuadrature method for " +
                    Element::GetElementName(elemtype) + ".\n");
  }
}

// Poly2DSorted: 2-D standard polynomial of r, s
// sorted by the powers of r
int Poly2DSorted::NumTerms(const int order) const {
  return (order + 1) * (order + 2) / 2;
}
void Poly2DSorted::F(const int order, const double* coord,
                     double* result) const {
  const Poly1D poly1d;
  std::vector<double> poly1d_r(poly1d.NumTerms(order));
  std::vector<double> poly1d_s(poly1d.NumTerms(order));
  poly1d.F(order, coord, &poly1d_r[0]);
  poly1d.F(order, coord + 1, &poly1d_s[0]);
  int index = 0;
  for (int i = 0; i <= order; i++)
    for (int j = 0; j <= order - i; j++)
      result[index++] = poly1d_s[i] * poly1d_r[j];
}
void Poly2DSorted::GetSurfaceQuadrature(const int order,
                                        const ElemType elemtype,
                                        const int elemorder,
                                        std::vector<double>& points,
                                        std::vector<double>& weights) const {
  int deg = 0;
  switch (elemtype) {
    case ElemType::LINE:
      deg = order * elemorder + elemorder - 1;
      quadrature::Line_Poly1D(deg, points, weights);
      return;
    default:
      ERROR_MESSAGE(poly_ + " has no GetSurfaceQuadrature method for " +
                    Element::GetElementName(elemtype) + ".\n");
  }
}

// PolyQuadShape: Shape polynomial of quadrilateral of r, s
int PolyQuadShape::NumTerms(const int order) const {
  return (order + 1) * (order + 1);
}
void PolyQuadShape::F(const int order, const double* coord,
                      double* result) const {
  const double& r = coord[0];
  const double& s = coord[1];
  result[0] = 1.0;
  int index = 1;
  for (int i = 1; i <= order; i++) {
    result[index] = result[index - 2 * i + 1] * r;
    index++;
    for (int j = 0; j < i; j++) {
      result[index] = result[index - 1] * s;
      index++;
    }
    for (int j = 0; j < i; j++) {
      result[index] = result[index - 2 * i - 1] * s;
      index++;
    }
  }
}
void PolyQuadShape::Fdr(const int order, const double* coord,
                        double* result) const {
  const double& s = coord[1];
  const Poly1D poly1d;
  std::vector<double> poly1d_dr(poly1d.NumTerms(order));
  poly1d.Fdr(order, coord, &poly1d_dr[0]);
  int index = 0;
  for (int i = 0; i <= order; i++) {
    result[index] = poly1d_dr[i];
    index++;
    for (int j = 0; j < i; j++) {
      result[index] = result[index - 1] * s;
      index++;
    }
    for (int j = 0; j < i; j++) {
      result[index] = result[index - 2 * i - 1] * s;
      index++;
    }
  }
}
void PolyQuadShape::Fds(const int order, const double* coord,
                        double* result) const {
  const double& r = coord[0];
  const Poly1D poly1d;
  std::vector<double> poly1d_ds(poly1d.NumTerms(order));
  poly1d.Fdr(order, coord + 1, &poly1d_ds[0]);
  result[0] = 0.0;
  int index = 1;
  for (int i = 1; i <= order; i++) {
    for (int j = 0; j < i; j++) {
      result[index] = result[index - 2 * i + 1] * r;
      index++;
    }
    index = (i + 1) * (i + 1) - 1;
    result[index] = poly1d_ds[i];
    index--;
    for (int j = 0; j < i; j++) {
      result[index] = result[index + 1] * r;
      index--;
    }
    index = (i + 1) * (i + 1);
  }
}
void PolyQuadShape::Fdrr(const int order, const double* coord,
                         double* result) const {
  const double& s = coord[1];
  const Poly1D poly1d;
  std::vector<double> poly1d_drr(poly1d.NumTerms(order));
  poly1d.Fdrr(order, coord, &poly1d_drr[0]);
  int index = 0;
  for (int i = 0; i <= order; i++) {
    result[index] = poly1d_drr[i];
    index++;
    for (int j = 0; j < i; j++) {
      result[index] = result[index - 1] * s;
      index++;
    }
    for (int j = 0; j < i; j++) {
      result[index] = result[index - 2 * i - 1] * s;
      index++;
    }
  }
}
void PolyQuadShape::Fdrs(const int order, const double* coord,
                         double* result) const {
  const Poly1D poly1d;
  std::vector<double> poly1d_dr(poly1d.NumTerms(order));
  std::vector<double> poly1d_ds(poly1d.NumTerms(order));
  poly1d.Fdr(order, coord, &poly1d_dr[0]);
  poly1d.Fdr(order, coord + 1, &poly1d_ds[0]);
  int index = 0;
  for (int i = 0; i <= order; i++) {
    for (int j = 0; j < i; j++) result[index++] = poly1d_dr[i] * poly1d_ds[j];
    for (int j = i - 1; j >= 0; j--)
      result[index++] = poly1d_dr[j] * poly1d_ds[i];
  }
}
void PolyQuadShape::Fdss(const int order, const double* coord,
                         double* result) const {
  const double& r = coord[0];
  const Poly1D poly1d;
  std::vector<double> poly1d_dss(poly1d.NumTerms(order));
  poly1d.Fdrr(order, coord + 1, &poly1d_dss[0]);
  result[0] = 0.0;
  int index = 1;
  for (int i = 1; i <= order; i++) {
    for (int j = 0; j < i; j++) {
      result[index] = result[index - 2 * i + 1] * r;
      index++;
    }
    index = (i + 1) * (i + 1) - 1;
    result[index] = poly1d_dss[i];
    index--;
    for (int j = 0; j < i; j++) {
      result[index] = result[index + 1] * r;
      index--;
    }
    index = (i + 1) * (i + 1);
  }
}
void PolyQuadShape::GetVolumeQuadrature(const int order,
                                        const ElemType elemtype,
                                        const int elemorder,
                                        std::vector<double>& points,
                                        std::vector<double>& weights) const {
  int deg = 0;
  switch (elemtype) {
    case ElemType::QUAD:
      deg = 2 * order * elemorder + 2 * elemorder - 1;
      quadrature::Quad_QuadShape(deg, points, weights);
      return;
    default:
      ERROR_MESSAGE(poly_ + " has no GetVolumeQuadrature method for " +
                    Element::GetElementName(elemtype) + ".\n");
  }
}

// PolyQuadShapeSorted: Shape polynomial of quadrilateral of r, s
// sorted by the powers of r
int PolyQuadShapeSorted::NumTerms(const int order) const {
  return (order + 1) * (order + 1);
}
void PolyQuadShapeSorted::F(const int order, const double* coord,
                            double* result) const {
  const Poly1D poly1d;
  std::vector<double> poly1d_r(poly1d.NumTerms(order));
  std::vector<double> poly1d_s(poly1d.NumTerms(order));
  poly1d.F(order, coord, &poly1d_r[0]);
  poly1d.F(order, coord + 1, &poly1d_s[0]);
  int index = 0;
  for (int i = 0; i <= order; i++)
    for (int j = 0; j <= order; j++)
      result[index++] = poly1d_s[i] * poly1d_r[j];
}

// Poly3D: 3-D standard polynomial of r, s, t
int Poly3D::NumTerms(const int order) const {
  return (order + 1) * (order + 2) * (order + 3) / 6;
}
void Poly3D::F(const int order, const double* coord, double* result) const {
  const double& r = coord[0];
  const double& s = coord[1];
  const double& t = coord[2];
  result[0] = 1.0;
  int index = 1;
  int bef_index = 0;
  int def = 0;
  for (int i = 1; i <= order; i++) {
    result[index] = result[bef_index] * r;
    index++;
    def = index - bef_index;
    bef_index = index - 1;
    for (int j = 1; j <= i; j++) {
      result[index] = result[index - def] * s;
      index++;
      def++;
      for (int k = 1; k <= j; k++) {
        result[index] = result[index - def] * t;
        index++;
      }
    }
  }
}
void Poly3D::Fdr(const int order, const double* coord, double* result) const {
  const double& s = coord[1];
  const double& t = coord[2];
  const Poly1D poly1d;
  std::vector<double> poly1d_dr(poly1d.NumTerms(order));
  poly1d.Fdr(order, coord, &poly1d_dr[0]);
  result[0] = 0.0;
  int index = 1;
  int bef_index = 0;
  int def = 0;
  for (int i = 1; i <= order; i++) {
    result[index] = poly1d_dr[i];
    index++;
    def = index - bef_index;
    bef_index = index - 1;
    for (int j = 1; j <= i; j++) {
      result[index] = result[index - def] * s;
      index++;
      def++;
      for (int k = 1; k <= j; k++) {
        result[index] = result[index - def] * t;
        index++;
      }
    }
  }
}
void Poly3D::Fds(const int order, const double* coord, double* result) const {
  const double& r = coord[0];
  const double& t = coord[2];
  const Poly1D poly1d;
  std::vector<double> poly1d_ds(poly1d.NumTerms(order));
  poly1d.Fdr(order, coord + 1, &poly1d_ds[0]);
  result[0] = 0.0;
  int index = 1;
  int bef_index = 0;
  int def = 0;
  for (int i = 1; i <= order; i++) {
    def = index - bef_index;
    bef_index = index;
    for (int j = 0; j < i; j++) {
      for (int k = 0; k <= j; k++) {
        result[index] = result[index - def] * r;
        index++;
      }
    }
    result[index] = poly1d_ds[i];
    index++;
    def += (i + 1);
    for (int k = 1; k <= i; k++) {
      result[index] = result[index - def] * t;
      index++;
    }
  }
}
void Poly3D::Fdt(const int order, const double* coord, double* result) const {
  const double& r = coord[0];
  const double& s = coord[1];
  const Poly1D poly1d;
  std::vector<double> poly1d_dt(poly1d.NumTerms(order));
  poly1d.Fdr(order, coord + 2, &poly1d_dt[0]);
  result[0] = 0.0;
  int index = 1;
  int bef_index = 0;
  int def = 0;
  for (int i = 1; i <= order; i++) {
    def = index - bef_index;
    bef_index = index;
    for (int j = 0; j < i; j++) {
      for (int k = 0; k <= j; k++) {
        result[index] = result[index - def] * r;
        index++;
      }
    }
    def += i;
    for (int k = 0; k < i; k++) {
      result[index] = result[index - def] * s;
      index++;
    }
    result[index] = poly1d_dt[i];
    index++;
  }
}
void Poly3D::GetVolumeQuadrature(const int order, 
                                 const ElemType elemtype, const int elemorder,
                                 std::vector<double>& points,
                                 std::vector<double>& weights) const {
  int deg = 0;
  int deg_tris = 0;
  int deg_line = 0;
  switch (elemtype) {
    case ElemType::TETS:
      deg = order * elemorder + 3 * (elemorder - 1);
      quadrature::Tets_Poly3D(deg, points, weights);
      return;
    case ElemType::HEXA:
      deg = order * elemorder + 3 * elemorder - 1;
      quadrature::Hexa_HexaShape(deg, points, weights);
      return;
    case ElemType::PRIS:
      deg_tris = order * elemorder + 3 * elemorder - 2;
      deg_line = order * elemorder + 3 * elemorder - 1;
      quadrature::Pris_PrisShape(deg_tris, deg_line, points, weights);
      return;
    case ElemType::PYRA:
      deg = order * elemorder + 3 * elemorder - 1;
      quadrature::Pyra_PyraShape(deg, points, weights);
      return;
    default:
      ERROR_MESSAGE(poly_ + " has no GetVolumeQuadrature method for " +
                    Element::GetElementName(elemtype) + ".\n");
  }
}

// Poly3DSorted: 3-D standard polynomial of r, s, t
// sorted by the powers of r, s
int Poly3DSorted::NumTerms(const int order) const {
  return (order + 1) * (order + 2) * (order + 3) / 6;
}
void Poly3DSorted::F(const int order, const double* coord,
                     double* result) const {
  const Poly1D poly1d;
  const Poly2D poly2d;
  std::vector<double> poly2d_rs(poly2d.NumTerms(order));
  std::vector<double> poly1d_t(poly1d.NumTerms(order));
  poly2d.F(order, coord, &poly2d_rs[0]);
  poly1d.F(order, coord + 2, &poly1d_t[0]);
  int index = 0;
  for (int i = 0; i <= order; i++)
    for (int j = 0, len = poly2d.NumTerms(order - i); j < len; j++)
      result[index++] = poly2d_rs[j] * poly1d_t[i];
}
void Poly3DSorted::GetSurfaceQuadrature(const int order,
                                        const ElemType elemtype,
                                        const int elemorder,
                                        std::vector<double>& points,
                                        std::vector<double>& weights) const {
  int deg = 0;
  switch (elemtype) {
    case ElemType::TRIS:
      deg = order * elemorder + 2*(elemorder - 1);
      quadrature::Tris_Poly2D(deg, points, weights);
      return;
    case ElemType::QUAD:
      deg = order * elemorder + 2 * elemorder - 1;
      quadrature::Quad_QuadShape(deg, points, weights);
      return;
    default:
      ERROR_MESSAGE(poly_ + " has no GetSurfaceQuadrature method for " +
                    Element::GetElementName(elemtype) + ".\n");
  }
}

// Poly3DQuad: 3-D standard polynomial of r, s, t
// with PolyQuadShape of r, s
int Poly3DQuad::NumTerms(const int order) const {
  return (order + 1) * ((order + 2) * (order + 3) + 3 * order) / 6;
}
void Poly3DQuad::F(const int order, const double* coord, double* result) const {
  const Poly1D poly1d;
  const Poly2D poly2d;
  const PolyQuadShape quad;
  std::vector<double> poly1d_t(poly1d.NumTerms(order));
  std::vector<double> poly2d_rs(poly2d.NumTerms(order));
  std::vector<double> quad_rs(quad.NumTerms(order));
  poly1d.F(order, coord + 2, &poly1d_t[0]);
  poly2d.F(order, coord, &poly2d_rs[0]);
  quad.F(order, coord, &quad_rs[0]);

  int index = 0;
  for (int i = 0, num_quad = quad.NumTerms(order); i < num_quad; i++)
    result[index++] = quad_rs[i] * poly1d_t[0];
  for (int i = 1; i <= order; i++)
    for (int j = 0, num_tris = poly2d.NumTerms(order - i); j < num_tris; j++)
      result[index++] = poly2d_rs[j] * poly1d_t[i];
}

// PolyHexaShape: Shape polynomial of Hexahedron of r, s, t
int PolyHexaShape::NumTerms(const int order) const {
  return (order + 1) * (order + 1) * (order + 1);
}
void PolyHexaShape::F(const int order, const double* coord,
                      double* result) const {
  const double& r = coord[0];
  const double& s = coord[1];
  const double& t = coord[2];
  result[0] = 1.0;
  int index = 1;
  int bef_index = 0;
  int def = 0;
  for (int i = 1; i <= order; i++) {
    def = index - bef_index;
    result[index] = result[bef_index] * r;
    bef_index = index;
    index++;
    for (int j = 1; j <= i; j++) {  // bottom
      result[index] = result[index - 1] * s;
      index++;
    }
    def += 2;
    for (int j = 1; j <= i; j++) {
      result[index] = result[index - def] * s;
      index++;
    }
    def = 2 * i + 1;
    for (int j = 1; j < i; j++) {  // wall
      for (int k = 0; k < def; k++) {
        result[index] = result[index - def] * t;
        index++;
      }
    }
    // top
    def = i * (3 * i + 1);
    for (int j = 0; j < i; j++) {
      for (int k = 0; k < 2 * j + 1; k++) {
        result[index] = result[index - def] * t;
        index++;
      }
    }
    def = (i + 1) * (i + 1);
    for (int j = 0; j < 2 * i + 1; j++) {
      result[index] = result[index - def] * t;
      index++;
    }
  }
}
void PolyHexaShape::Fdr(const int order, const double* coord,
                        double* result) const {
  const double& s = coord[1];
  const double& t = coord[2];
  const Poly1D poly1d;
  std::vector<double> poly1d_dr(poly1d.NumTerms(order));
  poly1d.Fdr(order, coord, &poly1d_dr[0]);
  result[0] = 0.0;
  int index = 1;
  int bef_index = 0;
  int def = 0;
  for (int i = 1; i <= order; i++) {
    result[index] = poly1d_dr[i];
    def = index - bef_index;
    bef_index = index;
    index++;
    for (int j = 1; j <= i; j++) {  // bottom
      result[index] = result[index - 1] * s;
      index++;
    }
    def += 2;
    for (int j = 1; j < i; j++) {
      result[index] = result[index - def] * s;
      index++;
    }
    result[index++] = 0.0;
    def = 2 * i + 1;
    for (int j = 1; j < i; j++) {  // wall
      for (int k = 1; k < def; k++) {
        result[index] = result[index - def] * t;
        index++;
      }
      result[index++] = 0.0;
    }
    // top
    def = i * (3 * i + 1);
    for (int j = 0; j < i; j++) {
      for (int k = 1; k < 2 * j + 1; k++) {
        result[index] = result[index - def] * t;
        index++;
      }
      result[index++] = 0.0;
    }
    def = (i + 1) * (i + 1);
    for (int j = 1; j < 2 * i + 1; j++) {
      result[index] = result[index - def] * t;
      index++;
    }
    result[index++] = 0.0;
  }
}
void PolyHexaShape::Fds(const int order, const double* coord,
                        double* result) const {
  const double& r = coord[0];
  const double& t = coord[2];
  const Poly1D poly1d;
  std::vector<double> poly1d_ds(poly1d.NumTerms(order));
  poly1d.Fdr(order, coord + 1, &poly1d_ds[0]);
  result[0] = 0.0;
  int index = 1;
  int bef_index = 0;
  int def = 0;
  for (int i = 1; i <= order; i++) {
    def = index - bef_index;
    bef_index = index;
    result[index++] = 0.0;
    for (int j = 1; j < i; j++) {  // bottom
      result[index] = result[index - def] * r;
      index++;
    }
    index += i;
    result[index] = poly1d_ds[i];
    for (int j = 1; j <= i; j++) {
      index--;
      result[index] = result[index + 1] * r;
    }
    index += (i + 1);
    def = 2 * i + 1;
    for (int j = 1; j < i; j++) {  // wall
      result[index++] = 0.0;
      for (int k = 1; k < def; k++) {
        result[index] = result[index - def] * t;
        index++;
      }
    }
    // top
    def = i * (3 * i + 1);
    for (int j = 0; j < i; j++) {
      result[index++] = 0.0;
      for (int k = 1; k < 2 * j + 1; k++) {
        result[index] = result[index - def] * t;
        index++;
      }
    }
    def = (i + 1) * (i + 1);
    result[index++] = 0.0;
    for (int j = 1; j < 2 * i + 1; j++) {
      result[index] = result[index - def] * t;
      index++;
    }
  }
}
void PolyHexaShape::Fdt(const int order, const double* coord,
                        double* result) const {
  const double& r = coord[0];
  const double& s = coord[1];
  const Poly1D poly1d;
  std::vector<double> poly1d_dt(poly1d.NumTerms(order));
  poly1d.Fdr(order, coord + 2, &poly1d_dt[0]);
  result[0] = 0.0;
  int index = 1;
  int bef_index = 0;
  int def = 0;
  for (int i = 1; i <= order; i++) {
    def = index - bef_index + 2;
    bef_index = index;
    for (int j = 0; j < 2 * i + 1; j++) {  // bottom
      result[index++] = 0.0;
    }
    for (int j = 1; j < i; j++) {  // wall
      if (j == i - 1) {
        def -= (i - 1) * (i - 1);
      }
      result[index] = result[index - def] * r;
      index++;
      for (int k = 1; k <= i; k++) {
        result[index] = result[index - 1] * s;
        index++;
      }
      def += 2;
      for (int k = 1; k <= i; k++) {
        result[index] = result[index - def] * s;
        index++;
      }
    }
    // top
    result[index++] = poly1d_dt[i];
    def = 1;
    for (int j = 1; j <= i; j++) {
      result[index] = result[index - def] * r;
      index++;
      for (int k = 1; k <= j; k++) {
        result[index] = result[index - 1] * s;
        index++;
      }
      def += 2;
      for (int k = 1; k <= j; k++) {
        result[index] = result[index - def] * s;
        index++;
      }
    }
  }
}
void PolyHexaShape::GetVolumeQuadrature(const int order,
                                 const ElemType elemtype, const int elemorder,
                                 std::vector<double>& points,
                                 std::vector<double>& weights) const {
  int deg = 0;
  switch (elemtype) {
    case ElemType::HEXA:
      deg = 3 * order * elemorder + 3 * elemorder - 1;
      quadrature::Hexa_HexaShape(deg, points, weights);
      return;
    default:
      ERROR_MESSAGE(poly_ + " has no GetVolumeQuadrature method for " +
                    Element::GetElementName(elemtype) + ".\n");
  }
}

// PolyHexaShapeSorted: Shape polynomial of Hexahedron of r, s, t
// sorted by the powers of r, s
int PolyHexaShapeSorted::NumTerms(const int order) const {
  return (order + 1) * (order + 1) * (order + 1);
}
void PolyHexaShapeSorted::F(const int order, const double* coord,
                            double* result) const {
  const Poly1D poly1d;
  PolyQuadShape quad;
  std::vector<double> quad_rs(quad.NumTerms(order));
  std::vector<double> poly1d_t(poly1d.NumTerms(order));
  quad.F(order, coord, &quad_rs[0]);
  poly1d.F(order, coord + 2, &poly1d_t[0]);
  int index = 0;
  for (int i = 0; i <= order; i++)
    for (int j = 0, len = quad.NumTerms(order); j < len; j++)
      result[index++] = quad_rs[j] * poly1d_t[i];
}

// PolyPrisShape: Shape polynomial of prism of r, s, t
int PolyPrisShape::NumTerms(const int order) const {
  return (order + 1) * (order + 1) * (order + 2) / 2;
}
void PolyPrisShape::F(const int order, const double* coord,
                      double* result) const {
  const double& r = coord[0];
  const double& s = coord[1];
  const double& t = coord[2];
  result[0] = 1.0;
  int index = 1;
  int bef_index = 0;
  int def = 0;
  for (int i = 1; i <= order; i++) {
    result[index] = result[bef_index] * r;
    index++;
    def = index - bef_index;
    bef_index = index - 1;
    for (int j = 1; j <= i; j++) {  // bottom
      result[index] = result[index - def] * s;
      index++;
    }
    def = i + 1;
    for (int j = 1; j < i; j++) {  // wall
      for (int k = 0; k < def; k++) {
        result[index] = result[index - def] * t;
        index++;
      }
    }
    // top
    def = i * (i + 1) * 3 / 2;
    for (int j = 1; j <= i; j++) {
      for (int k = 1; k <= j; k++) {
        result[index] = result[index - def] * t;
        index++;
      }
    }
    def = (i + 1) * (i + 2) / 2;
    for (int j = 0; j < i + 1; j++) {
      result[index] = result[index - def] * t;
      index++;
    }
  }
}
void PolyPrisShape::Fdr(const int order, const double* coord,
                        double* result) const {
  const double& s = coord[1];
  const double& t = coord[2];
  const Poly1D poly1d;
  std::vector<double> poly1d_dr(poly1d.NumTerms(order));
  poly1d.Fdr(order, coord, &poly1d_dr[0]);
  result[0] = 0.0;
  int index = 1;
  int bef_index = 0;
  int def = 0;
  for (int i = 1; i <= order; i++) {
    result[index] = poly1d_dr[i];
    index++;
    def = index - bef_index;
    bef_index = index - 1;
    for (int j = 1; j <= i; j++) {  // bottom
      result[index] = result[index - def] * s;
      index++;
    }
    def = i + 1;
    for (int j = 1; j < i; j++) {  // wall
      for (int k = 0; k < def; k++) {
        result[index] = result[index - def] * t;
        index++;
      }
    }
    // top
    def = i * (i + 1) * 3 / 2;
    for (int j = 1; j <= i; j++) {
      for (int k = 1; k <= j; k++) {
        result[index] = result[index - def] * t;
        index++;
      }
    }
    def = (i + 1) * (i + 2) / 2;
    for (int j = 0; j < i + 1; j++) {
      result[index] = result[index - def] * t;
      index++;
    }
  }
}
void PolyPrisShape::Fds(const int order, const double* coord,
                        double* result) const {
  const double& r = coord[0];
  const double& t = coord[2];
  const Poly1D poly1d;
  std::vector<double> poly1d_ds(poly1d.NumTerms(order));
  poly1d.Fdr(order, coord + 1, &poly1d_ds[0]);
  result[0] = 0.0;
  int index = 1;
  int bef_index = 0;
  int def = 0;
  for (int i = 1; i <= order; i++) {
    def = index - bef_index;
    bef_index = index;
    result[index++] = 0.0;
    for (int j = 1; j < i; j++) {  // bottom
      result[index] = result[index - def] * r;
      index++;
    }
    result[index++] = poly1d_ds[i];
    def = i + 1;
    for (int j = 1; j < i; j++) {  // wall
      for (int k = 0; k < def; k++) {
        result[index] = result[index - def] * t;
        index++;
      }
    }
    // top
    def = i * (i + 1) * 3 / 2;
    for (int j = 1; j <= i; j++) {
      for (int k = 1; k <= j; k++) {
        result[index] = result[index - def] * t;
        index++;
      }
    }
    def = (i + 1) * (i + 2) / 2;
    for (int j = 0; j < i + 1; j++) {
      result[index] = result[index - def] * t;
      index++;
    }
  }
}
void PolyPrisShape::Fdt(const int order, const double* coord,
                        double* result) const {
  const double& r = coord[0];
  const double& s = coord[1];
  const Poly1D poly1d;
  std::vector<double> poly1d_dt(poly1d.NumTerms(order));
  poly1d.Fdr(order, coord + 2, &poly1d_dt[0]);
  result[0] = 0.0;
  int index = 1;
  int bef_index = 0;
  int def = 0;
  for (int i = 1; i <= order; i++) {
    def = index - bef_index;
    bef_index = index;

    for (int j = 1; j <= i + 1; j++)  // bottom
      result[index++] = 0.0;

    def++;
    for (int j = 1; j < i; j++) {  // wall
      if (j == i - 1) def = i * (i + 1) - 1;
      result[index] = result[index - def] * r;
      index++;
      def++;
      for (int k = 1; k <= i; k++) {
        result[index] = result[index - def] * s;
        index++;
      }
    }

    // top
    result[index++] = poly1d_dt[i];
    def = 1;
    for (int j = 1; j <= i; j++) {
      result[index] = result[index - def] * r;
      index++;
      def++;
      for (int k = 1; k <= j; k++) {
        result[index] = result[index - def] * s;
        index++;
      }
    }
  }
}
void PolyPrisShape::GetVolumeQuadrature(const int order,
                                        const ElemType elemtype,
                                        const int elemorder,
                                        std::vector<double>& points,
                                        std::vector<double>& weights) const {
  int deg_tris = 0;
  int deg_line = 0;
  switch (elemtype) {
    case ElemType::PRIS:
      deg_tris = 2 * order * elemorder + 3 * elemorder - 2;
      deg_line = 2 * order * elemorder + 3 * elemorder - 1;
      quadrature::Pris_PrisShape(deg_tris, deg_line, points, weights);
      return;
    default:
      ERROR_MESSAGE(poly_ + " has no GetVolumeQuadrature method for " +
                    Element::GetElementName(elemtype) + ".\n");
  }
}

// PolyPyraShape: Shape (rational) polynomial of pyramid of r, s, t
int PolyPyraShape::NumTerms(const int order) const {
  return (order + 1) * (order + 2) * (2 * order + 3) / 6;
}
void PolyPyraShape::F(const int order, const double* coord,
                      double* result) const {
  const double& r = coord[0];
  const double& s = coord[1];
  const double& t = coord[2];
  result[0] = 1.0;
  if (order == 0) return;

  double mix = 0.0;
  if (std::abs(1.0 - t) > 1.0E-10) mix = r * s / (1.0 - t);
  const Poly1D poly1d;
  const Poly2D poly2d;
  const Poly3D poly3d;
  std::vector<double> poly1d_mix(poly1d.NumTerms(order));
  std::vector<double> poly2d_rs(poly2d.NumTerms(order - 1));
  std::vector<double> poly3d_rst(poly3d.NumTerms(order));
  poly1d.F(order, &mix, &poly1d_mix[0]);
  poly2d.F(order - 1, coord, &poly2d_rs[0]);
  poly3d.F(order, coord, &poly3d_rst[0]);

  std::vector<int> poly2d_index(order + 2);
  std::vector<int> poly3d_index(order + 2);
  for (int i = 0; i <= order + 1; i++) {
    poly2d_index[i] = poly2d.NumTerms(i - 1);
    poly3d_index[i] = poly3d.NumTerms(i - 1);
  }

  int index = 1;
  for (int i = 1; i <= order; i++) {
    for (int j = poly3d_index[i]; j < poly3d_index[i + 1]; j++)
      result[index++] = poly3d_rst[j];
    for (int j = 0; j < i; j++)
      for (int k = poly2d_index[i - j - 1]; k < poly2d_index[i - j]; k++)
        result[index++] = poly1d_mix[j + 1] * poly2d_rs[k];
  }
}
void PolyPyraShape::Fdr(const int order, const double* coord,
                        double* result) const {
  const double& r = coord[0];
  const double& s = coord[1];
  const double& t = coord[2];
  result[0] = 0.0;
  if (order == 0) return;

  double mix = 0.0;
  double mix_dr = 0.0;
  if (std::abs(1.0 - t) > 1.0E-10) {
    mix = r * s / (1.0 - t);
    mix_dr = s / (1.0 - t);
  }
  const Poly1D poly1d;
  const Poly2D poly2d;
  const Poly3D poly3d;
  std::vector<double> poly1d_mix(poly1d.NumTerms(order));
  std::vector<double> poly2d_rs(poly2d.NumTerms(order - 1));
  std::vector<double> poly3d_rst(poly3d.NumTerms(order));
  poly1d.F(order, &mix, &poly1d_mix[0]);
  poly2d.F(order - 1, coord, &poly2d_rs[0]);
  poly3d.F(order, coord, &poly3d_rst[0]);

  std::vector<double> poly1d_mix_dr(poly1d.NumTerms(order));
  std::vector<double> poly2d_rs_dr(poly2d.NumTerms(order - 1));
  std::vector<double> poly3d_rst_dr(poly3d.NumTerms(order));
  poly1d.Fdr(order, &mix, &poly1d_mix_dr[0]);
  poly2d.Fdr(order - 1, coord, &poly2d_rs_dr[0]);
  poly3d.Fdr(order, coord, &poly3d_rst_dr[0]);

  std::vector<int> poly2d_index(order + 2);
  std::vector<int> poly3d_index(order + 2);
  for (int i = 0; i <= order + 1; i++) {
    poly2d_index[i] = poly2d.NumTerms(i - 1);
    poly3d_index[i] = poly3d.NumTerms(i - 1);
  }

  int index = 1;
  for (int i = 1; i <= order; i++) {
    for (int j = poly3d_index[i]; j < poly3d_index[i + 1]; j++)
      result[index++] = poly3d_rst_dr[j];
    for (int j = 0; j < i; j++)
      for (int k = poly2d_index[i - j - 1]; k < poly2d_index[i - j]; k++)
        result[index++] = poly1d_mix_dr[j + 1] * mix_dr * poly2d_rs[k] +
                          poly1d_mix[j + 1] * poly2d_rs_dr[k];
  }
}
void PolyPyraShape::Fds(const int order, const double* coord,
                        double* result) const {
  const double& r = coord[0];
  const double& s = coord[1];
  const double& t = coord[2];
  result[0] = 0.0;
  if (order == 0) return;

  double mix = 0.0;
  double mix_ds = 0.0;
  if (std::abs(1.0 - t) > 1.0E-10) {
    mix = r * s / (1.0 - t);
    mix_ds = r / (1.0 - t);
  }
  const Poly1D poly1d;
  const Poly2D poly2d;
  const Poly3D poly3d;
  std::vector<double> poly1d_mix(poly1d.NumTerms(order));
  std::vector<double> poly2d_rs(poly2d.NumTerms(order - 1));
  std::vector<double> poly3d_rst(poly3d.NumTerms(order));
  poly1d.F(order, &mix, &poly1d_mix[0]);
  poly2d.F(order - 1, coord, &poly2d_rs[0]);
  poly3d.F(order, coord, &poly3d_rst[0]);

  std::vector<double> poly1d_mix_ds(poly1d.NumTerms(order));
  std::vector<double> poly2d_rs_ds(poly2d.NumTerms(order - 1));
  std::vector<double> poly3d_rst_ds(poly3d.NumTerms(order));
  poly1d.Fdr(order, &mix, &poly1d_mix_ds[0]);
  poly2d.Fds(order - 1, coord, &poly2d_rs_ds[0]);
  poly3d.Fds(order, coord, &poly3d_rst_ds[0]);

  std::vector<int> poly2d_index(order + 2);
  std::vector<int> poly3d_index(order + 2);
  for (int i = 0; i <= order + 1; i++) {
    poly2d_index[i] = poly2d.NumTerms(i - 1);
    poly3d_index[i] = poly3d.NumTerms(i - 1);
  }

  int index = 1;
  for (int i = 1; i <= order; i++) {
    for (int j = poly3d_index[i]; j < poly3d_index[i + 1]; j++)
      result[index++] = poly3d_rst_ds[j];
    for (int j = 0; j < i; j++)
      for (int k = poly2d_index[i - j - 1]; k < poly2d_index[i - j]; k++)
        result[index++] = poly1d_mix_ds[j + 1] * mix_ds * poly2d_rs[k] +
                          poly1d_mix[j + 1] * poly2d_rs_ds[k];
  }
}
void PolyPyraShape::Fdt(const int order, const double* coord,
                        double* result) const {
  const double& r = coord[0];
  const double& s = coord[1];
  const double& t = coord[2];
  result[0] = 0.0;
  if (order == 0) return;

  double mix = 0.0;
  double mix_dt = 0.0;
  if (std::abs(1.0 - t) > 1.0E-10) {
    mix = r * s / (1.0 - t);
    mix_dt = mix / (1.0 - t);
  }
  const Poly1D poly1d;
  const Poly2D poly2d;
  const Poly3D poly3d;
  std::vector<double> poly1d_mix(poly1d.NumTerms(order));
  std::vector<double> poly2d_rs(poly2d.NumTerms(order - 1));
  std::vector<double> poly3d_rst(poly3d.NumTerms(order));
  poly1d.F(order, &mix, &poly1d_mix[0]);
  poly2d.F(order - 1, coord, &poly2d_rs[0]);
  poly3d.F(order, coord, &poly3d_rst[0]);

  std::vector<double> poly1d_mix_dt(poly1d.NumTerms(order));
  std::vector<double> poly3d_rst_dt(poly3d.NumTerms(order));
  poly1d.Fdr(order, &mix, &poly1d_mix_dt[0]);
  poly3d.Fdt(order, coord, &poly3d_rst_dt[0]);

  std::vector<int> poly2d_index(order + 2);
  std::vector<int> poly3d_index(order + 2);
  for (int i = 0; i <= order + 1; i++) {
    poly2d_index[i] = poly2d.NumTerms(i - 1);
    poly3d_index[i] = poly3d.NumTerms(i - 1);
  }

  int index = 1;
  for (int i = 1; i <= order; i++) {
    for (int j = poly3d_index[i]; j < poly3d_index[i + 1]; j++)
      result[index++] = poly3d_rst_dt[j];
    for (int j = 0; j < i; j++)
      for (int k = poly2d_index[i - j - 1]; k < poly2d_index[i - j]; k++)
        result[index++] = poly1d_mix_dt[j + 1] * mix_dt * poly2d_rs[k];
  }
}

// PolyPyra: Polynomial of r, s, t
// defined by sum_{i=0}^{order}{PolyQuadShape(i;r,s)*t**i}
int PolyPyra::NumTerms(const int order) const {
  return (order + 1) * (order + 2) * (2 * order + 3) / 6;
}
void PolyPyra::F(const int order, const double* coord, double* result) const {
  result[0] = 1.0;
  if (order == 0) return;

  const Poly1D poly1d;
  const PolyQuadShape quad;
  std::vector<double> poly1d_t(poly1d.NumTerms(order));
  std::vector<double> quad_rs(quad.NumTerms(order));
  poly1d.F(order, coord + 2, &poly1d_t[0]);
  quad.F(order, coord, &quad_rs[0]);

  int index = 1;
  for (int i = 1; i <= order; i++)
    for (int j = 0; j <= i; j++)
      for (int k = j * j; k < (j + 1) * (j + 1); k++)
        result[index++] = quad_rs[k] * poly1d_t[i - j];
}
void PolyPyra::Fdr(const int order, const double* coord, double* result) const {
  result[0] = 0.0;
  if (order == 0) return;

  const Poly1D poly1d;
  const PolyQuadShape quad;
  std::vector<double> poly1d_t(poly1d.NumTerms(order));
  std::vector<double> quad_rs_dr(quad.NumTerms(order));
  poly1d.F(order, coord + 2, &poly1d_t[0]);
  quad.Fdr(order, coord, &quad_rs_dr[0]);

  int index = 1;
  for (int i = 1; i <= order; i++)
    for (int j = 0; j <= i; j++)
      for (int k = j * j; k < (j + 1) * (j + 1); k++)
        result[index++] = quad_rs_dr[k] * poly1d_t[i - j];
}
void PolyPyra::Fds(const int order, const double* coord, double* result) const {
  result[0] = 0.0;
  if (order == 0) return;

  const Poly1D poly1d;
  const PolyQuadShape quad;
  std::vector<double> poly1d_t(poly1d.NumTerms(order));
  std::vector<double> quad_rs_ds(quad.NumTerms(order));
  poly1d.F(order, coord + 2, &poly1d_t[0]);
  quad.Fds(order, coord, &quad_rs_ds[0]);

  int index = 1;
  for (int i = 1; i <= order; i++)
    for (int j = 0; j <= i; j++)
      for (int k = j * j; k < (j + 1) * (j + 1); k++)
        result[index++] = quad_rs_ds[k] * poly1d_t[i - j];
}
void PolyPyra::Fdt(const int order, const double* coord, double* result) const {
  result[0] = 0.0;
  if (order == 0) return;

  const Poly1D poly1d;
  const PolyQuadShape quad;
  std::vector<double> poly1d_dt(poly1d.NumTerms(order));
  std::vector<double> quad_rs(quad.NumTerms(order));
  poly1d.Fdr(order, coord + 2, &poly1d_dt[0]);
  quad.F(order, coord, &quad_rs[0]);

  int index = 1;
  for (int i = 1; i <= order; i++)
    for (int j = 0; j <= i; j++)
      for (int k = j * j; k < (j + 1) * (j + 1); k++)
        result[index++] = quad_rs[k] * poly1d_dt[i - j];
}
void PolyPyra::GetVolumeQuadrature(const int order,
                                        const ElemType elemtype,
                                        const int elemorder,
                                        std::vector<double>& points,
                                        std::vector<double>& weights) const {
  int deg = 0;
  switch (elemtype) {
    case ElemType::PYRA:
      deg = 2 * order * elemorder + 3 * elemorder - 1;
      quadrature::Pyra_PyraShape(deg, points, weights);
      return;
    default:
      ERROR_MESSAGE(poly_ + " has no GetVolumeQuadrature method for " +
                    Element::GetElementName(elemtype) + ".\n");
  }
}
}  // namespace deneb