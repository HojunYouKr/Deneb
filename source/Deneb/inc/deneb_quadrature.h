#pragma once

#include <vector>

namespace deneb {
namespace quadrature {
// Gauss-Legendre
void Line_Poly1D(const int degree, std::vector<double>& points,
                 std::vector<double>& weights);
// Witherden-Vincent
void Tris_Poly2D(const int degree, std::vector<double>& points,
                 std::vector<double>& weights);
// Gauss-Legendre
void Quad_QuadShape(const int degree, std::vector<double>& points,
                    std::vector<double>& weights);
// Witherden-Vincent
void Tets_Poly3D(const int degree, std::vector<double>& points,
                 std::vector<double>& weights);
// Gauss-Legendre
void Hexa_HexaShape(const int degree, std::vector<double>& points,
                    std::vector<double>& weights);
// Gauss-Legendre (height) * Witherden-Vincent (base)
void Pris_PrisShape(const int degree_tris, const int degree_line,
                    std::vector<double>& points, std::vector<double>& weights);
// Gauss-Legendre
void Pyra_PyraShape(const int degree, std::vector<double>& points,
                    std::vector<double>& weights);
}  // namespace quadrature
}  // namespace deneb
