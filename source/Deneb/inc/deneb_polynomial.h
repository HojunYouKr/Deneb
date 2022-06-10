#pragma once

#include <string>
#include <vector>

#include "deneb_element.h"

namespace deneb {
class Polynomial {
 protected:
  const std::string poly_;
  void (Polynomial::*grad_[3])(const int order, const double* coord,
                               double* result) const;

 public:
  Polynomial(const std::string poly) : poly_(poly) {
    grad_[0] = &Polynomial::Fdr;
    grad_[1] = &Polynomial::Fds;
    grad_[2] = &Polynomial::Fdt;
  };
  virtual ~Polynomial(){};

  virtual int NumTerms(const int order) const;
  virtual void F(const int order, const double* coord, double* result) const;

  virtual void Fdr(const int order, const double* coord, double* result) const;
  virtual void Fds(const int order, const double* coord, double* result) const;
  virtual void Fdt(const int order, const double* coord, double* result) const;

  virtual void Fdrr(const int order, const double* coord, double* result) const;
  virtual void Fdrs(const int order, const double* coord, double* result) const;
  virtual void Fdss(const int order, const double* coord, double* result) const;

  virtual void GetVolumeQuadrature(const int order, const ElemType elemtype,
                                   const int elemorder,
                                   std::vector<double>& points,
                                   std::vector<double>& weights) const;
  virtual void GetSurfaceQuadrature(const int order, const ElemType elemtype,
                                    const int elemorder,
                                    std::vector<double>& points,
                                    std::vector<double>& weights) const;

  inline void dF(const int dim, const int order, const double* coord,
                 double* result) const {
    (this->*grad_[dim])(order, coord, result);
  }
};

// Poly1D: 1-D standard polynomial of r
class Poly1D : public Polynomial {
 public:
  Poly1D() : Polynomial("Poly1D"){};
  virtual ~Poly1D(){};

  virtual int NumTerms(const int order) const;
  virtual void F(const int order, const double* coord, double* result) const;
  virtual void Fdr(const int order, const double* coord, double* result) const;
  virtual void Fdrr(const int order, const double* coord, double* result) const;
};

// Poly2D: 2-D standard polynomial of r, s
class Poly2D : public Polynomial {
 public:
  Poly2D() : Polynomial("Poly2D"){};
  virtual ~Poly2D(){};

  virtual int NumTerms(const int order) const;
  virtual void F(const int order, const double* coord, double* result) const;
  virtual void Fdr(const int order, const double* coord, double* result) const;
  virtual void Fds(const int order, const double* coord, double* result) const;
  virtual void Fdrr(const int order, const double* coord, double* result) const;
  virtual void Fdrs(const int order, const double* coord, double* result) const;
  virtual void Fdss(const int order, const double* coord, double* result) const;

  virtual void GetVolumeQuadrature(const int order, const ElemType elemtype,
                                   const int elemorder,
                                   std::vector<double>& points,
                                   std::vector<double>& weights) const;
};

// Poly2DSorted: 2-D standard polynomial of r, s
// sorted by the powers of r
class Poly2DSorted : public Polynomial {
 public:
  Poly2DSorted() : Polynomial("Poly2DSorted"){};
  virtual ~Poly2DSorted(){};

  virtual int NumTerms(const int order) const;
  virtual void F(const int order, const double* coord, double* result) const;

  virtual void GetSurfaceQuadrature(const int order, const ElemType elemtype,
                                    const int elemorder,
                                    std::vector<double>& points,
                                    std::vector<double>& weights) const;
};

// PolyQuadShape: Shape polynomial of quadrilateral of r, s
class PolyQuadShape : public Polynomial {
 public:
  PolyQuadShape() : Polynomial("PolyQuadShape"){};
  virtual ~PolyQuadShape(){};

  virtual int NumTerms(const int order) const;
  virtual void F(const int order, const double* coord, double* result) const;
  virtual void Fdr(const int order, const double* coord, double* result) const;
  virtual void Fds(const int order, const double* coord, double* result) const;
  virtual void Fdrr(const int order, const double* coord, double* result) const;
  virtual void Fdrs(const int order, const double* coord, double* result) const;
  virtual void Fdss(const int order, const double* coord, double* result) const;

  virtual void GetVolumeQuadrature(const int order, const ElemType elemtype,
                                   const int elemorder,
                                   std::vector<double>& points,
                                   std::vector<double>& weights) const;
};

// PolyQuadShapeSorted: Shape polynomial of quadrilateral of r, s
// sorted by the powers of r
class PolyQuadShapeSorted : public Polynomial {
 public:
  PolyQuadShapeSorted() : Polynomial("PolyQuadShapeSorted"){};
  virtual ~PolyQuadShapeSorted(){};

  virtual int NumTerms(const int order) const;
  virtual void F(const int order, const double* coord, double* result) const;
};

// Poly3D: 3-D standard polynomial of r, s, t
class Poly3D : public Polynomial {
 public:
  Poly3D() : Polynomial("Poly3D"){};
  virtual ~Poly3D(){};

  virtual int NumTerms(const int order) const;
  virtual void F(const int order, const double* coord, double* result) const;
  virtual void Fdr(const int order, const double* coord, double* result) const;
  virtual void Fds(const int order, const double* coord, double* result) const;
  virtual void Fdt(const int order, const double* coord, double* result) const;

  virtual void GetVolumeQuadrature(const int order, const ElemType elemtype,
                                   const int elemorder,
                                   std::vector<double>& points,
                                   std::vector<double>& weights) const;
};

// Poly3DSorted: 3-D standard polynomial of r, s, t
// sorted by the powers of r, s
class Poly3DSorted : public Polynomial {
 public:
  Poly3DSorted() : Polynomial("Poly3DSorted"){};
  virtual ~Poly3DSorted(){};

  virtual int NumTerms(const int order) const;
  virtual void F(const int order, const double* coord, double* result) const;

  virtual void GetSurfaceQuadrature(const int order, const ElemType elemtype,
                                    const int elemorder,
                                    std::vector<double>& points,
                                    std::vector<double>& weights) const;
};

// Poly3DQuad: 3-D standard polynomial of r, s, t
// with PolyQuadShape of r, s
class Poly3DQuad : public Polynomial {
 public:
  Poly3DQuad() : Polynomial("Poly3DQuad"){};
  virtual ~Poly3DQuad(){};

  virtual int NumTerms(const int order) const;
  virtual void F(const int order, const double* coord, double* result) const;
};

// PolyHexaShape: Shape polynomial of Hexahedron of r, s, t
class PolyHexaShape : public Polynomial {
 public:
  PolyHexaShape() : Polynomial("PolyHexaShape"){};
  virtual ~PolyHexaShape(){};

  virtual int NumTerms(const int order) const;
  virtual void F(const int order, const double* coord, double* result) const;
  virtual void Fdr(const int order, const double* coord, double* result) const;
  virtual void Fds(const int order, const double* coord, double* result) const;
  virtual void Fdt(const int order, const double* coord, double* result) const;

  virtual void GetVolumeQuadrature(const int order, const ElemType elemtype,
                                   const int elemorder,
                                   std::vector<double>& points,
                                   std::vector<double>& weights) const;
};

// PolyHexaShapeSorted: Shape polynomial of Hexahedron of r, s, t
// sorted by the powers of r, s
class PolyHexaShapeSorted : public Polynomial {
 public:
  PolyHexaShapeSorted() : Polynomial("PolyHexaShapeSorted"){};
  virtual ~PolyHexaShapeSorted(){};

  virtual int NumTerms(const int order) const;
  virtual void F(const int order, const double* coord, double* result) const;
};

// PolyPrisShape: Shape polynomial of prism of r, s, t
class PolyPrisShape : public Polynomial {
 public:
  PolyPrisShape() : Polynomial("PolyPrisShape"){};
  virtual ~PolyPrisShape(){};

  virtual int NumTerms(const int order) const;
  virtual void F(const int order, const double* coord, double* result) const;
  virtual void Fdr(const int order, const double* coord, double* result) const;
  virtual void Fds(const int order, const double* coord, double* result) const;
  virtual void Fdt(const int order, const double* coord, double* result) const;

  virtual void GetVolumeQuadrature(const int order, const ElemType elemtype,
                                   const int elemorder,
                                   std::vector<double>& points,
                                   std::vector<double>& weights) const;
};

// PolyPyraShape: Shape (rational) polynomial of pyramid of r, s, t
class PolyPyraShape : public Polynomial {
 public:
  PolyPyraShape() : Polynomial("PolyPyraShape"){};
  virtual ~PolyPyraShape(){};

  virtual int NumTerms(const int order) const;
  virtual void F(const int order, const double* coord, double* result) const;
  virtual void Fdr(const int order, const double* coord, double* result) const;
  virtual void Fds(const int order, const double* coord, double* result) const;
  virtual void Fdt(const int order, const double* coord, double* result) const;
};

// PolyPyra: Polynomial of r, s, t
// defined by sum_{i=0}^{order}{PolyQuadShape(i;r,s)*t**i}
class PolyPyra : public Polynomial {
 public:
  PolyPyra() : Polynomial("PolyPyra"){};
  virtual ~PolyPyra(){};

  virtual int NumTerms(const int order) const;
  virtual void F(const int order, const double* coord, double* result) const;
  virtual void Fdr(const int order, const double* coord, double* result) const;
  virtual void Fds(const int order, const double* coord, double* result) const;
  virtual void Fdt(const int order, const double* coord, double* result) const;

  virtual void GetVolumeQuadrature(const int order, const ElemType elemtype,
                                   const int elemorder,
                                   std::vector<double>& points,
                                   std::vector<double>& weights) const;
};
}  // namespace deneb