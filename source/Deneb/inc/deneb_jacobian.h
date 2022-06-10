#pragma once

#include <memory>
#include <vector>

#include "avocado.h"
#include "deneb_polynomial.h"

namespace deneb {
class Jacobian {
 public:
  static std::shared_ptr<Jacobian> GetJacobian(const ElemType elemtype);

 protected:
  ElemType elemtype_;
  int dimension_;

  int elemorder_;
  int num_bases_;
  std::vector<double> transform_;

  std::shared_ptr<Polynomial> shape_;

 public:
  Jacobian(const ElemType elemtype, const int dimension)
      : elemtype_(elemtype),
        dimension_(dimension) {};
  virtual ~Jacobian(){};

  inline ElemType GetElemType() const { return elemtype_; };
  inline int GetElemOrder() const { return elemorder_; };
  inline const std::vector<double>& GetTransform() const { return transform_; };
  inline void SetTransform(const int elemorder,
                           const std::vector<double>& transform) {
    elemorder_ = elemorder;
    num_bases_ = shape_->NumTerms(elemorder);
    transform_ = transform;
  }

  void SetTopology(const int elemorder, const double* phy_coords);

  void TransformToPhyCoords(const int num_coords, const double* ref_coords,
                            double* phy_coords) const;

  void CalJacobianMatT(const int num_coords, const double* ref_coords,
                       double* results) const;
  void CalJacobianDet(const int num_coords, const double* ref_coords,
                      double* results) const;
  void CalJacobianCofMat(const int num_coords, const double* ref_coords,
                         double* results) const;
};

class JacobianTris : public Jacobian {
 public:
  JacobianTris() : Jacobian(ElemType::TRIS, 2) {
    shape_ = std::make_shared<Poly2D>();
  };
  virtual ~JacobianTris(){};
};
class JacobianQuad : public Jacobian {
 public:
  JacobianQuad() : Jacobian(ElemType::QUAD, 2) {
    shape_ = std::make_shared<PolyQuadShape>();
  };
  virtual ~JacobianQuad(){};
};
class JacobianTets : public Jacobian {
 public:
  JacobianTets() : Jacobian(ElemType::TETS, 3) {
    shape_ = std::make_shared<Poly3D>();
  };
  virtual ~JacobianTets(){};
};
class JacobianHexa : public Jacobian {
 public:
  JacobianHexa() : Jacobian(ElemType::HEXA, 3) {
    shape_ = std::make_shared<PolyHexaShape>();
  };
  virtual ~JacobianHexa(){};
};
class JacobianPris : public Jacobian {
 public:
  JacobianPris() : Jacobian(ElemType::PRIS, 3) {
    shape_ = std::make_shared<PolyPrisShape>();
  };
  virtual ~JacobianPris(){};
};
class JacobianPyra : public Jacobian {
 public:
  JacobianPyra() : Jacobian(ElemType::PYRA, 3) {
    shape_ = std::make_shared<PolyPyraShape>();
  };
  virtual ~JacobianPyra(){};
};
}  // namespace deneb