#include "deneb_basis.h"

#include <algorithm>
#include <cmath>

#include "avocado.h"

namespace deneb {
std::shared_ptr<Basis> Basis::GetStandardBasis(const ElemType elemtype) {
  switch (elemtype) {
    case ElemType::TRIS:
    case ElemType::QUAD:
      return GetAllocatedBasis(elemtype, std::make_shared<Poly2D>());
    case ElemType::TETS:
    case ElemType::HEXA:
    case ElemType::PRIS:
    case ElemType::PYRA:
      return GetAllocatedBasis(elemtype, std::make_shared<Poly3D>());
  }
  ERROR_MESSAGE("Invalid element type: " +
                std::to_string(static_cast<int>(elemtype)) + "\n");
  return nullptr;
}

std::shared_ptr<Basis> Basis::GetAllocatedBasis(
    const ElemType elemtype, std::shared_ptr<Polynomial> polynomial) {
  const int dimension = Element::GetDimension(elemtype);
  return std::make_shared<Basis>(elemtype, dimension, polynomial);
}

void Basis::SetTransform(const std::vector<double>& coords) {
  const Element* element = DENEB_ELEMENT->GetElement(elemtype_);
  const int num_nodes = element->GetNumNodes(1);
  const std::vector<double>& Rt = element->GetNodesCoord(1);

  std::vector<double> Xt;
  Xt.reserve(num_nodes * (dimension_ + 1));
  for (int inode = 0; inode < num_nodes; inode++) {
    for (int idim = 0; idim < dimension_; idim++)
      Xt.push_back(coords[inode * dimension_ + idim]);
    Xt.push_back(1.0);
  }

  std::vector<double> XRt((dimension_ + 1) * dimension_);
  std::vector<double> XXt((dimension_ + 1) * (dimension_ + 1));
  gemmATB(1.0, &Xt[0], &Rt[0], 0.0, &XRt[0], num_nodes, dimension_ + 1,
          dimension_);
  gemmATB(1.0, &Xt[0], &Xt[0], 0.0, &XXt[0], num_nodes, dimension_ + 1,
          dimension_ + 1);
  avocado::MatInv(dimension_ + 1, &XXt[0]);
  gemmAB(1.0, &XXt[0], &XRt[0], 0.0, &transform_[0], dimension_ + 1,
         dimension_ + 1, dimension_);
}
void Basis::GetBasis(const int num_coords, const double* coords,
                     double* results) const {
  std::vector<double> trans_coord(dimension_);
  std::vector<double> func_vars(num_bases_);
  for (int ipt = 0; ipt < num_coords; ipt++) {
    TransformCoord(&coords[ipt * dimension_], &trans_coord[0]);
    basis_->F(order_, &trans_coord[0], &func_vars[0]);
    for (int ibasis = 0, ind = 0; ibasis < num_bases_; ind += (++ibasis))
      results[ipt * num_bases_ + ibasis] = avocado::VecInnerProd(
          ibasis + 1, &connecting_matrix_[ind], &func_vars[0]);
  }
}
void Basis::GetBasisGrad(const int num_coords, const double* coords,
                         double* results) const {
  std::vector<double> trans_coord(dimension_);
  std::vector<double> mono_vars(dimension_ * num_bases_);
  std::vector<double> func_vars(dimension_ * num_bases_);
  for (int ipt = 0; ipt < num_coords; ipt++) {
    TransformCoord(&coords[ipt * dimension_], &trans_coord[0]);
    for (int idim = 0; idim < dimension_; idim++)
      basis_->dF(idim, order_, &trans_coord[0], &mono_vars[idim * num_bases_]);
    gemmAB(1.0, &transform_[0], &mono_vars[0], 0.0, &func_vars[0], dimension_,
           dimension_, num_bases_);
    for (int ibasis = 0, ind = 0; ibasis < num_bases_; ind += (++ibasis))
      gemvAx(1.0, &func_vars[0], num_bases_, &connecting_matrix_[ind], 1, 0.0,
             &results[(ipt * num_bases_ + ibasis) * dimension_], 1, dimension_,
             ibasis + 1);
  }
}
void Basis::ComputeConnectingMatrix(const int order,
                                    const std::vector<double>& points,
                                    const std::vector<double>& weights) {
  order_ = order;
  num_bases_ = basis_->NumTerms(order);

  const int num_points = static_cast<int>(weights.size());
  std::vector<double> q_var(num_bases_ * num_points);
  {
    std::vector<double> trans_coord(dimension_);
    std::vector<double> func_vars(num_bases_);
    for (int ipt = 0; ipt < num_points; ipt++) {
      TransformCoord(&points[ipt * dimension_], &trans_coord[0]);
      basis_->F(order_, &trans_coord[0], &func_vars[0]);
      for (int ibasis = 0; ibasis < num_bases_; ibasis++)
        q_var[ibasis * num_points + ipt] = func_vars[ibasis];
    }
  }

  std::vector<int> ind(num_bases_ + 1, 0);
  std::vector<double> conmat(num_bases_ * (num_bases_ + 1) / 2, 0.0);
  for (int i = 0; i < num_bases_; i++) {
    ind[i + 1] = ind[i] + (i + 1);
    conmat[ind[i + 1] - 1] = 1.0;
  }

  std::vector<double> r(num_bases_);
  for (int i = 0; i < num_bases_; i++) {
    for (int iter = 0; iter < iteration_; iter++) {
      for (int j = 0; j < i; j++) {
        r[j] = avocado::VecInnerProd(num_points, &q_var[i * num_points],
                                     &q_var[j * num_points], &weights[0]);
        cblas_daxpy(num_points, -r[j], &q_var[j * num_points], 1,
                    &q_var[i * num_points], 1);
      }
      const double norm2 =
          avocado::VecInnerProd(num_points, &q_var[i * num_points],
                                &q_var[i * num_points], &weights[0]);
      r[i] = 1.0 / std::sqrt(norm2);
      avocado::VecScale(num_points, r[i], &q_var[i * num_points]);
      avocado::VecScale(i, -r[i], &r[0]);
      for (int j = 0; j <= i; j++) {
        double sum = 0.0;
        for (int k = j; k <= i; k++) sum += r[k] * conmat[ind[k] + j];
        conmat[ind[i] + j] = sum;
      }
    }
  }
  connecting_matrix_ = std::move(conmat);
}
void Basis::TransformCoord(const double* coord, double* result) const {
  std::vector<double> temp(coord, coord + dimension_);
  temp.push_back(1);
  gemvATx(1.0, &transform_[0], dimension_, &temp[0], 1, 0.0, result, 1,
          dimension_ + 1, dimension_);
}

// BasisSurface
std::shared_ptr<BasisSurface> BasisSurface::GetAllocatedBasis(
    const ElemType elemtype, std::shared_ptr<Polynomial> polynomial) {
  return std::make_shared<BasisSurface>(elemtype, polynomial);
}
void BasisSurface::SetTransform(const std::vector<double>& coords) {
  const Element* element = DENEB_ELEMENT->GetElement(elemtype_);
  const int num_nodes = element->GetNumNodes(1);

  std::vector<double> rotate(dimension_ * dimension_, 0.0);
  {
    for (int idim = 0; idim < dimension_; idim++)
      rotate[idim * (dimension_ + 1)] = 1.0;
    std::vector<double> stdev(dimension_);
    for (int idim = 0; idim < dimension_; idim++)
      stdev[idim] = avocado::VecStdDev(num_nodes, &coords[idim], dimension_);
    const int dim = static_cast<int>(
        std::min_element(stdev.begin(), stdev.end()) - stdev.begin());
    std::vector<double> A, b;
    A.reserve(num_nodes * dimension_);
    b.reserve(num_nodes);
    for (int inode = 0; inode < num_nodes; inode++) {
      for (int idim = 0; idim < dimension_; idim++)
        if (idim != dim) A.push_back(coords[inode * dimension_ + idim]);
      A.push_back(1.0);
      b.push_back(coords[inode * dimension_ + dim]);
    }
    int num_ranks;
    avocado::MatPseudoInv(&A[0], num_nodes, dimension_, num_ranks);
    std::vector<double> plane(dimension_);
    gemvAx(1.0, &A[0], num_nodes, &b[0], 1, 0.0, &plane[0], 1, dimension_,
           num_nodes);
    std::vector<double> normal;
    int ind = 0;
    for (int idim = 0; idim < dimension_; idim++) {
      if (idim != dim)
        normal.push_back(plane[ind++]);
      else
        normal.push_back(-1.0);
    }
    avocado::VecNormal(dimension_, &normal[0]);

    if (dimension_ == 2) {
      const double length = std::abs(normal[0]);
      if (length > 1.0E-12) {
        const double uz = normal[0] / length;
        const double theta = std::acos(normal[1]);
        const double c = std::cos(theta);
        const double s = std::sin(theta);
        rotate = {c, -uz * s, uz * s, c};
      }
    } else if (dimension_ == 3) {
      const double length =
          std::sqrt(normal[0] * normal[0] + normal[1] * normal[1]);
      if (length > 1.0E-12) {
        const double ux = normal[1] / length;
        const double uy = -normal[0] / length;
        const double theta = std::acos(normal[2]);
        const double c = std::cos(theta);
        const double s = std::sin(theta);
        rotate = {c + ux * ux * (1.0 - c),
                  ux * uy * (1.0 - c),
                  uy * s,
                  ux * uy * (1.0 - c),
                  c + uy * uy * (1.0 - c),
                  -ux * s,
                  -uy * s,
                  ux * s,
                  c};
      }
    }
  }

  std::vector<double> rotated_coords(num_nodes * dimension_);
  gemmABT(1.0, &coords[0], &rotate[0], 0.0, &rotated_coords[0], num_nodes,
          dimension_, dimension_);

  std::vector<double> transform;
  const std::vector<double>& R = element->GetNodesCoord(1);
  {
    std::vector<double> avg(dimension_);
    for (int idim = 0; idim < dimension_; idim++)
      avg[idim] =
          avocado::VecAverage(num_nodes, &rotated_coords[idim], dimension_);
    if (dimension_ == 2) {
      // [r] = [a 0 b][x]
      // [0]   [0 1 c][y]
      //              [1]
      const double c = -avg[dimension_ - 1];
      std::vector<double> A, B;
      A.reserve(num_nodes * 2);
      B.reserve(num_nodes);
      for (int inode = 0; inode < num_nodes; inode++) {
        A.push_back(rotated_coords[inode * dimension_]);
        A.push_back(1.0);
        B.push_back(R[inode]);
      }
      int num_ranks;
      avocado::MatPseudoInv(&A[0], num_nodes, 2, num_ranks);
      std::vector<double> ab(2);
      gemvAx(1.0, &A[0], num_nodes, &B[0], 1, 0.0, &ab[0], 1, 2, num_nodes);
      const double& a = ab[0];
      const double& b = ab[1];
      transform = {a, 0.0, 0.0, 1.0, b, c};
    } else if (dimension_ == 3) {
      // [r] = [a  b 0 c][x]       [r] = [a b c][x]
      // [s]   [d  e 0 f][y]  =>   [s]   [d e f][y]
      // [t]   [0  0 1 g][z]                    [1]
      //                 [1]
      const double g = -avg[dimension_ - 1];
      std::vector<double> r(2 * num_nodes);
      std::vector<double> X(3 * num_nodes);
      for (int inode = 0; inode < num_nodes; inode++) {
        const double& x = rotated_coords[inode * dimension_];
        const double& y = rotated_coords[inode * dimension_ + 1];
        r[inode] = R[inode * 2];
        r[num_nodes + inode] = R[inode * 2 + 1];
        X[inode] = x;
        X[num_nodes + inode] = y;
        X[2 * num_nodes + inode] = 1.0;
      }
      std::vector<double> coef(2 * 3, 0.0);
      {
        std::vector<double> rXt(2 * 3);
        std::vector<double> XXt(3 * 3);
        gemmABT(1.0, &r[0], &X[0], 0.0, &rXt[0], 2, num_nodes, 3);
        gemmABT(1.0, &X[0], &X[0], 0.0, &XXt[0], 3, num_nodes, 3);
        avocado::MatInv(3, &XXt[0]);
        gemmAB(1.0, &rXt[0], &XXt[0], 0.0, &coef[0], 2, 3, 3);
      }
      transform = {coef[0], coef[3], 0.0, coef[1], coef[4], 0.0,
                   0.0,     0.0,     1.0, coef[2], coef[5], g};
    }
  }

  transform_ = transform;
  gemmATB(1.0, &rotate[0], &transform[0], 0.0, &transform_[0], dimension_,
          dimension_, dimension_);
}
void BasisSurface::GetBasis(const int num_coords, const double* coords,
                            double* results) const {
  const int num_original_bases = basis_->NumTerms(order_);
  std::vector<double> trans_coord(dimension_);
  std::vector<double> func_vars(num_original_bases);
  for (int ipt = 0; ipt < num_coords; ipt++) {
    TransformCoord(&coords[ipt * dimension_], &trans_coord[0]);
    basis_->F(order_, &trans_coord[0], &func_vars[0]);
    for (int ibasis = 0, ind = 0; ibasis < num_bases_; ind += (++ibasis)) {
      double& sum = results[ipt * num_bases_ + ibasis];
      sum = 0.0;
      for (int i = 0; i <= ibasis; i++)
        sum += connecting_matrix_[ind + i] * func_vars[index_order_[i]];
    }
  }
}
void BasisSurface::ComputeConnectingMatrix(const int order,
                                           const std::vector<double>& points,
                                           const std::vector<double>& weights,
                                           const double eps) {
  order_ = order;
  const int num_original_bases = basis_->NumTerms(order);

  double length_ratio = 1.0;
  {
    std::vector<double> rotate(transform_.begin(),
                               transform_.begin() + dimension_ * dimension_);
    length_ratio = std::abs(avocado::MatDet(dimension_, &rotate[0]));
  }

  double maxh = 0.0;
  std::vector<double> normalize;
  const int num_points = static_cast<int>(weights.size());
  std::vector<double> q_var(num_original_bases * num_points);
  {
    std::vector<double> trans_coord(dimension_);
    std::vector<double> func_vars(num_original_bases);
    for (int ipt = 0; ipt < num_points; ipt++) {
      TransformCoord(&points[ipt * dimension_], &trans_coord[0]);
      maxh = std::max(std::abs(trans_coord[dimension_ - 1]), maxh);
      basis_->F(order_, &trans_coord[0], &func_vars[0]);
      for (int ibasis = 0; ibasis < num_original_bases; ibasis++)
        q_var[ibasis * num_points + ipt] = func_vars[ibasis];
    }

    std::vector<double> maxh_vec(dimension_, 1.0);
    maxh_vec[dimension_ - 1] = maxh;
    normalize.resize(num_original_bases, 1.0);
    if (maxh >= 100.0 * eps * length_ratio)
      basis_->F(order_, &maxh_vec[0], &normalize[0]);
  }

  std::vector<int> ind(num_original_bases + 1, 0);
  for (int i = 0; i < num_original_bases; i++) ind[i + 1] = ind[i] + (i + 1);

  int index = 0;
  int jndex = 0;
  double bef_var = 1.0;
  std::vector<char> stiff(num_original_bases, false);
  std::vector<double> conmat;

  for (int i = 0; i < num_original_bases; i++) {
    const std::vector<double> q_var_store(q_var.begin() + i * num_points,
                                          q_var.begin() + (i + 1) * num_points);
    for (int iter = 0; iter < iteration_; iter++) {
      std::vector<double> r(index + 1, 0.0);
      jndex = 0;
      for (int j = 0; j < i; j++) {
        if (static_cast<bool>(stiff[j]) == true) continue;
        r[jndex] = avocado::VecInnerProd(num_points, &q_var[i * num_points],
                                         &q_var[j * num_points], &weights[0]);
        cblas_daxpy(num_points, -r[jndex], &q_var[j * num_points], 1,
                    &q_var[i * num_points], 1);
        jndex++;
      }
      const double norm =
          std::sqrt(avocado::VecInnerProd(num_points, &q_var[i * num_points],
                                          &q_var[i * num_points], &weights[0]));

      if (iter == 0) {
        const double var = norm / normalize[i];
        if (var < 100.0 * eps * length_ratio || var < 0.1 * bef_var ||
            std::isnan(var)) {
          stiff[i] = true;
          cblas_dcopy(num_points, &q_var_store[0], 1, &q_var[i * num_points],
                      1);
          break;
        }
        bef_var = var;
        conmat.resize(conmat.size() + index + 1, 0.0);
        conmat.back() = 1.0;
      }

      r[index] = 1.0 / norm;
      avocado::VecScale(num_points, r[index], &q_var[i * num_points]);
      avocado::VecScale(index, -r[index], &r[0]);
      for (int j = 0; j <= index; j++) {
        double sum = 0.0;
        for (int k = j; k <= index; k++) sum += r[k] * conmat[ind[k] + j];
        conmat[ind[index] + j] = sum;
      }
    }
    if (!stiff[i]) index++;
  }

  std::vector<char> second_stiff(num_original_bases, true);
  for (int i = 0; i < num_original_bases; i++) {
    if (!stiff[i]) continue;
    for (int iter = 0; iter < iteration_; iter++) {
      std::vector<double> r(index + 1, 0.0);
      jndex = 0;
      for (int j = 0; j < num_original_bases; j++) {
        if (static_cast<bool>(stiff[j]) == true) continue;
        r[jndex] = avocado::VecInnerProd(num_points, &q_var[i * num_points],
                                         &q_var[j * num_points], &weights[0]);
        cblas_daxpy(num_points, -r[jndex], &q_var[j * num_points], 1,
                    &q_var[i * num_points], 1);
        jndex++;
      }
      for (int j = 0; j < i; j++) {
        if (static_cast<bool>(second_stiff[j]) == true) continue;
        r[jndex] = avocado::VecInnerProd(num_points, &q_var[i * num_points],
                                         &q_var[j * num_points], &weights[0]);
        cblas_daxpy(num_points, -r[jndex], &q_var[j * num_points], 1,
                    &q_var[i * num_points], 1);
        jndex++;
      }
      const double norm =
          std::sqrt(avocado::VecInnerProd(num_points, &q_var[i * num_points],
                                          &q_var[i * num_points], &weights[0]));

      if (iter == 0) {
        const double var = norm / normalize[i];
        if (var < 100.0 * eps * length_ratio || std::isnan(var)) break;
        second_stiff[i] = false;
        conmat.resize(conmat.size() + index + 1, 0.0);
        conmat.back() = 1.0;
      }

      r[index] = 1.0 / norm;
      avocado::VecScale(num_points, r[index], &q_var[i * num_points]);
      avocado::VecScale(index, -r[index], &r[0]);
      for (int j = 0; j <= index; j++) {
        double sum = 0.0;
        for (int k = j; k <= index; k++) sum += r[k] * conmat[ind[k] + j];
        conmat[ind[index] + j] = sum;
      }
    }
    if (!second_stiff[i]) index++;
  }

  std::vector<int> match(num_original_bases, -1);
  jndex = 0;
  for (int i = 0; i < num_original_bases; i++) {
    if (stiff[i]) continue;
    match[i] = jndex;
    jndex++;
  }
  for (int i = 0; i < num_original_bases; i++) {
    if (second_stiff[i]) continue;
    match[i] = jndex;
    jndex++;
  }
  num_bases_ = jndex;

  std::vector<int> index_order(num_bases_);
  for (int ibasis = 0; ibasis < num_original_bases; ibasis++) {
    if (match[ibasis] == -1) continue;
    index_order[match[ibasis]] = ibasis;
  }
  index_order_ = move(index_order);
  connecting_matrix_ = std::move(conmat);
}
}  // namespace deneb