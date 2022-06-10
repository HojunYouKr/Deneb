#include "deneb_jacobian.h"

namespace deneb {
std::shared_ptr<Jacobian> Jacobian::GetJacobian(const ElemType elemtype) {
  switch (elemtype) {
    case ElemType::TRIS:
      return std::make_shared<JacobianTris>();
    case ElemType::QUAD:
      return std::make_shared<JacobianQuad>();
    case ElemType::TETS:
      return std::make_shared<JacobianTets>();
    case ElemType::HEXA:
      return std::make_shared<JacobianHexa>();
    case ElemType::PRIS:
      return std::make_shared<JacobianPris>();
    case ElemType::PYRA:
      return std::make_shared<JacobianPyra>();
  }
  ERROR_MESSAGE("Illegal element type: " +
               std::to_string(static_cast<int>(elemtype)) + "\n");
  return std::shared_ptr<Jacobian>();
}

void Jacobian::SetTopology(const int elemorder, const double* phy_coords) {
  elemorder_ = elemorder;
  num_bases_ = shape_->NumTerms(elemorder_);

  std::vector<double> rs(num_bases_ * num_bases_);
  const Element* element = DENEB_ELEMENT->GetElement(elemtype_);
  const std::vector<double>& ref_coords = element->GetNodesCoord(elemorder_);
  for (int i = 0; i < num_bases_; i++)
    shape_->F(elemorder_, &ref_coords[i * dimension_], &rs[i * num_bases_]);

  std::vector<double> transform(dimension_ * num_bases_);
  avocado::MatInv(num_bases_, &rs[0]);
  gemmATBT(1.0, phy_coords, &rs[0], 0.0, &transform[0], num_bases_, dimension_,
           num_bases_);
  transform_ = std::move(transform);
}

void Jacobian::TransformToPhyCoords(const int num_coords,
                                    const double* ref_coords,
                                    double* phy_coords) const {
  std::vector<double> rs(num_coords * num_bases_);
  for (int i = 0; i < num_coords; i++)
    shape_->F(elemorder_, &ref_coords[i * dimension_], &rs[i * num_bases_]);
  gemmABT(1.0, &rs[0], &transform_[0], 0.0, phy_coords, num_coords, num_bases_,
          dimension_);
}

void Jacobian::CalJacobianMatT(const int num_coords, const double* ref_coords,
                               double* results) const {
  std::vector<double> rs(num_coords * num_bases_);
  for (int idim = 0; idim < dimension_; idim++) {
    for (int i = 0; i < num_coords; i++)
      shape_->dF(idim, elemorder_, &ref_coords[i * dimension_],
                          &rs[i * num_bases_]);
    gemmABTn(1.0, &rs[0], num_bases_, &transform_[0], num_bases_, 0.0,
             &results[idim * dimension_], dimension_ * dimension_,
             num_coords, num_bases_, dimension_);
  }
}
void Jacobian::CalJacobianDet(const int num_coords, const double* ref_coords,
                              double* results) const {
  const int dd = dimension_ * dimension_;
  std::vector<double> jacmat(num_coords * dd);
  CalJacobianMatT(num_coords, ref_coords, &jacmat[0]);
  for (int i = 0; i < num_coords; i++)
    results[i] = avocado::MatDet(dimension_, &jacmat[i * dd]);
}
void Jacobian::CalJacobianCofMat(const int num_coords, const double* ref_coords,
                                 double* results) const {
  const int dd = dimension_ * dimension_;
  std::vector<double> jacmat(num_coords * dd);
  CalJacobianMatT(num_coords, ref_coords, &jacmat[0]);
  int ind = 0;
  if (dimension_ == 2) {
    for (int i = 0; i < num_coords; i++) {
      ind = i * dd;
      const double& xr = jacmat[ind++];
      const double& yr = jacmat[ind++];
      const double& xs = jacmat[ind++];
      const double& ys = jacmat[ind];
      ind = i * dd;
      results[ind++] = ys;
      results[ind++] = -yr;
      results[ind++] = -xs;
      results[ind] = xr;
    }
  } else if (dimension_ == 3) {
    for (int i = 0; i < num_coords; i++) {
      ind = i * dd;
      const double& xr = jacmat[ind++];
      const double& yr = jacmat[ind++];
      const double& zr = jacmat[ind++];
      const double& xs = jacmat[ind++];
      const double& ys = jacmat[ind++];
      const double& zs = jacmat[ind++];
      const double& xt = jacmat[ind++];
      const double& yt = jacmat[ind++];
      const double& zt = jacmat[ind];
      ind = i * dd;
      results[ind++] = ys * zt - yt * zs;
      results[ind++] = yt * zr - yr * zt;
      results[ind++] = yr * zs - ys * zr;
      results[ind++] = xt * zs - xs * zt;
      results[ind++] = xr * zt - xt * zr;
      results[ind++] = xs * zr - xr * zs;
      results[ind++] = xs * yt - xt * ys;
      results[ind++] = xt * yr - xr * yt;
      results[ind] = xr * ys - xs * yr;
    }
  } else
    ERROR_MESSAGE("Invalid dimension for " +
                 Element::GetElementName(elemtype_) + ": " +
                 std::to_string(dimension_) + "\n");
}
};  // namespace deneb