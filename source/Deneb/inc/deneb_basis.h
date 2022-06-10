#pragma once

#include <memory>
#include <vector>

#include "deneb_polynomial.h"

namespace deneb {
class Polynomial;

class Basis {
 public:
  static std::shared_ptr<Basis> GetStandardBasis(const ElemType elemtype);
  static std::shared_ptr<Basis> GetAllocatedBasis(
      const ElemType elemtype, std::shared_ptr<Polynomial> polynomial);

 protected:
  static const int iteration_ = 3;
  ElemType elemtype_;
  int dimension_;

  int order_;
  int num_bases_;
  // transform_ = [A, b]**T = [A**T;b**T]
  // Transformation does A*x+b,
  // where A, b are least square solutions of r=A*x+b
  // x: physical coordinate
  // r: reference coordinate
  std::vector<double> transform_;
  std::vector<double> connecting_matrix_;

  std::shared_ptr<Polynomial> basis_;

 public:
  Basis(const ElemType elemtype, const int dimension,
        std::shared_ptr<Polynomial> basis)
      : elemtype_(elemtype),
        dimension_(dimension),
        order_(-1),
        num_bases_(-1),
        transform_((dimension + 1) * dimension, 0.0),
        connecting_matrix_(),
        basis_(basis) {
    for (int idim = 0; idim < dimension_; idim++)
      transform_[idim * (dimension_ + 1)] = 1.0;
  };
  virtual ~Basis(){};

  inline int GetNumBases(void) const { return num_bases_; };
  inline const std::vector<double>& GetTransform(void) const {
    return transform_;
  };
  inline const std::vector<double>& GetConnectingMatrix(void) const {
    return connecting_matrix_;
  };
  inline std::shared_ptr<Polynomial> GetBasisPolynomial(void) const {
    return basis_;
  };
  inline void SetBasis(const int order, const std::vector<double>& transform,
                       const std::vector<double>& connecting_matrix) {
    order_ = order;
    num_bases_ = basis_->NumTerms(order);
    transform_ = transform;
    connecting_matrix_ = connecting_matrix;
  }

  virtual void SetTransform(const std::vector<double>& coords);
  // Data order of results: p*b
  virtual void GetBasis(const int num_coords, const double* coords,
                        double* results) const;
  // Data order of results: p*b*d
  void GetBasisGrad(const int num_coords, const double* coords,
                    double* results) const;
  void ComputeConnectingMatrix(const int order,
                               const std::vector<double>& points,
                               const std::vector<double>& weights);

 protected:
  void TransformCoord(const double* coord, double* result) const;
};

class BasisSurface : public Basis {
 public:
  static std::shared_ptr<BasisSurface> GetAllocatedBasis(
      const ElemType elemtype, std::shared_ptr<Polynomial> polynomial);
  
 protected:
  std::vector<int> index_order_;

 public:
  BasisSurface(const ElemType elemtype, std::shared_ptr<Polynomial> basis)
      : Basis(elemtype, Element::GetDimension(elemtype) + 1, basis) {};
  virtual ~BasisSurface(){};

  virtual void SetTransform(const std::vector<double>& coords);
  // Data order of results: p*b
  virtual void GetBasis(const int num_coords, const double* coords,
                        double* results) const;
  void ComputeConnectingMatrix(const int order,
                               const std::vector<double>& points,
                               const std::vector<double>& weights,
                               const double eps);
};
}  // namespace deneb