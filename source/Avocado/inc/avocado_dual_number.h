#pragma once

#include <algorithm>
#include <cstring>

template <int dim>
struct Dual {
  double f;
  double df[dim];

  Dual() {
    f = 0.0;
    memset(df, 0, dim * sizeof(double));
  }
  Dual(const double fvar) {
    f = fvar;
    memset(df, 0, dim * sizeof(double));
  }
  Dual(const double fvar, const int dimvar) {
    f = fvar;
    memset(df, 0, dim * sizeof(double));
    df[dimvar] = 1.0;
  }
  Dual<dim> operator-() const {
    Dual<dim> var;
    var.f = -f;
    for (int idim = 0; idim < dim; idim++) var.df[idim] = -df[idim];
    return var;
  }
};

struct DualSystem {
  double f;
  double fx;
  double fy;
  double fxx;
  double fxy;
  double fyy;

  double ft;
  double ftt;
  double fxt;
  double fyt;
  double fxxt;
  double fxyt;
  double fyyt;

  void Init(void) {
    f = fx = fy = fxx = fxy = fyy = 0.0;
    ft = ftt = fxt = fyt = fxxt = fxyt = fyyt = 0.0;
  }
  DualSystem() { Init(); }
  DualSystem(const double fvar) {
    Init();
    f = fvar;
  }
  DualSystem(const double fvar, const double fxvar, const double fyvar,
             const double ftvar) {
    Init();
    f = fvar;
    fx = fxvar;
    fy = fyvar;
    ft = ftvar;
  }
};

DualSystem operator+(const DualSystem& left, const DualSystem& right);
DualSystem operator-(const DualSystem& left, const DualSystem& right);
DualSystem operator*(const DualSystem& left, const DualSystem& right);
DualSystem operator/(const DualSystem& left, const DualSystem& right);

template <int dim>
Dual<dim> operator+(const Dual<dim>& left, const Dual<dim>& right) {
  Dual<dim> var(left.f + right.f);
  for (int idim = 0; idim < dim; idim++)
    var.df[idim] = left.df[idim] + right.df[idim];
  return var;
}
template <int dim>
Dual<dim> operator-(const Dual<dim>& left, const Dual<dim>& right) {
  Dual<dim> var(left.f - right.f);
  for (int idim = 0; idim < dim; idim++)
    var.df[idim] = left.df[idim] - right.df[idim];
  return var;
}
template <int dim>
Dual<dim> operator*(const Dual<dim>& left, const Dual<dim>& right) {
  Dual<dim> var(left.f * right.f);
  for (int idim = 0; idim < dim; idim++)
    var.df[idim] = left.f * right.df[idim] + left.df[idim] * right.f;
  return var;
}
template <int dim>
Dual<dim> operator/(const Dual<dim>& left, const Dual<dim>& right) {
  const double t1 = 1.0 / right.f;
  Dual<dim> var(left.f * t1);
  const double t3 = var.f * t1;
  for (int idim = 0; idim < dim; idim++)
    var.df[idim] = t1 * left.df[idim] - t3 * right.df[idim];
  return var;
}
template <int dim>
bool operator<(const Dual<dim>& left, const Dual<dim>& right) {
  return (left.f < right.f);
}
template <int dim>
bool operator>(const Dual<dim>& left, const Dual<dim>& right) {
  return (left.f > right.f);
}
template <int dim>
bool operator<=(const Dual<dim>& left, const Dual<dim>& right) {
  return (left.f <= right.f);
}
template <int dim>
bool operator>=(const Dual<dim>& left, const Dual<dim>& right) {
  return (left.f >= right.f);
}

template <int dim>
Dual<dim> operator+(const double& left, const Dual<dim>& right) {
  Dual<dim> var = right;
  var.f += left;
  return var;
}
template <int dim>
Dual<dim> operator-(const double& left, const Dual<dim>& right) {
  Dual<dim> var = -right;
  var.f += left;
  return var;
}
template <int dim>
Dual<dim> operator*(const double& left, const Dual<dim>& right) {
  Dual<dim> var(left * right.f);
  for (int idim = 0; idim < dim; idim++) var.df[idim] = left * right.df[idim];
  return var;
}
template <int dim>
Dual<dim> operator/(const double& left, const Dual<dim>& right) {
  const double t1 = 1.0 / right.f;
  Dual<dim> var(left * t1);
  const double t3 = var.f * t1;
  for (int idim = 0; idim < dim; idim++) var.df[idim] = -t3 * right.df[idim];
  return var;
}
template <int dim>
bool operator<(const double& left, const Dual<dim>& right) {
  return (left < right.f);
}
template <int dim>
bool operator>(const double& left, const Dual<dim>& right) {
  return (left > right.f);
}
template <int dim>
bool operator<=(const double& left, const Dual<dim>& right) {
  return (left <= right.f);
}
template <int dim>
bool operator>=(const double& left, const Dual<dim>& right) {
  return (left >= right.f);
}

template <int dim>
Dual<dim> operator+(const Dual<dim>& left, const double& right) {
  Dual<dim> var = left;
  var.f += right;
  return var;
}
template <int dim>
Dual<dim> operator-(const Dual<dim>& left, const double& right) {
  Dual<dim> var = left;
  var.f -= right;
  return var;
}
template <int dim>
Dual<dim> operator*(const Dual<dim>& left, const double& right) {
  Dual<dim> var(left.f * right);
  for (int idim = 0; idim < dim; idim++) var.df[idim] = left.df[idim] * right;
  return var;
}
template <int dim>
Dual<dim> operator/(const Dual<dim>& left, const double& right) {
  const double t1 = 1.0 / right;
  Dual<dim> var(left.f * t1);
  for (int idim = 0; idim < dim; idim++) var.df[idim] = t1 * left.df[idim];
  return var;
}
template <int dim>
bool operator<(const Dual<dim>& left, const double& right) {
  return (left.f < right);
}
template <int dim>
bool operator>(const Dual<dim>& left, const double& right) {
  return (left.f > right);
}
template <int dim>
bool operator<=(const Dual<dim>& left, const double& right) {
  return (left.f <= right);
}
template <int dim>
bool operator>=(const Dual<dim>& left, const double& right) {
  return (left.f >= right);
}

template <typename T>
int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

namespace std {
DualSystem cos(const DualSystem& var);
DualSystem sin(const DualSystem& var);
DualSystem exp(const DualSystem& var);

template <int dim>
Dual<dim> abs(const Dual<dim>& var) {
  const int sign = sgn(var.f);
  Dual<dim> result(std::abs(var.f));
  for (int idim = 0; idim < dim; idim++) result.df[idim] = sign * var.df[idim];
  return result;
}
template <int dim>
Dual<dim> sqrt(const Dual<dim>& var) {
  Dual<dim> result(std::sqrt(var.f));
  const double t = 0.5 / result.f;
  for (int idim = 0; idim < dim; idim++) result.df[idim] = t * var.df[idim];
  return result;
}
template <int dim>
Dual<dim> cos(const Dual<dim>& var) {
  Dual<dim> result(std::cos(var.f));
  const double t = -std::sin(var.f);
  for (int idim = 0; idim < dim; idim++) result.df[idim] = t * var.df[idim];
  return result;
}
template <int dim>
Dual<dim> sin(const Dual<dim>& var) {
  Dual<dim> result(std::sin(var.f));
  const double t = std::cos(var.f);
  for (int idim = 0; idim < dim; idim++) result.df[idim] = t * var.df[idim];
  return result;
}
template <int dim>
Dual<dim> max(const Dual<dim>& left, const Dual<dim>& right) {
  if (left.f > right.f)
    return left;
  else
    return right;
}
template <int dim>
Dual<dim> max(const double& left, const Dual<dim>& right) {
  if (left > right.f)
    return Dual<dim>(left);
  else
    return right;
}
template <int dim>
Dual<dim> max(const Dual<dim>& left, const double& right) {
  if (left.f > right)
    return left;
  else
    return Dual<dim>(right);
}
template <int dim>
Dual<dim> min(const Dual<dim>& left, const Dual<dim>& right) {
  if (left.f < right.f)
    return left;
  else
    return right;
}
template <int dim>
Dual<dim> min(const double& left, const Dual<dim>& right) {
  if (left < right.f)
    return Dual<dim>(left);
  else
    return right;
}
template <int dim>
Dual<dim> min(const Dual<dim>& left, const double& right) {
  if (left.f < right)
    return left;
  else
    return Dual<dim>(right);
}
template <int dim>
Dual<dim> pow(const Dual<dim>& var, const double& exp) {
  Dual<dim> result(std::pow(var.f, exp));
  const double t = exp * std::pow(var.f, exp - 1.0);
  for (int idim = 0; idim < dim; idim++) result.df[idim] = t * var.df[idim];
  return result;
}
template <int dim>
Dual<dim> pow(const Dual<dim>& var, const int& exp) {
  Dual<dim> result(std::pow(var.f, exp));
  const double t = exp * std::pow(var.f, exp - 1.0);
  for (int idim = 0; idim < dim; idim++) result.df[idim] = t * var.df[idim];
  return result;
}
template <int dim>
Dual<dim> exp(const Dual<dim>& var) {
  const double t = std::exp(var.f);
  Dual<dim> result(t);
  for (int idim = 0; idim < dim; idim++) result.df[idim] = t * var.df[idim];
  return result;
}
template <int dim>
Dual<dim> atan(const Dual<dim>& var) {
  const double t = 1.0 / (1.0 + var.f * var.f);
  Dual<dim> result(std::atan(var.f));
  for (int idim = 0; idim < dim; idim++) result.df[idim] = t * var.df[idim];
  return result;
}
}  // namespace std