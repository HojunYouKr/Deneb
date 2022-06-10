#include "avocado_dual_number.h"

DualSystem operator+(const DualSystem& left, const DualSystem& right) {
  DualSystem var;
  var.f = left.f + right.f;
  var.fx = left.fx + right.fx;
  var.fy = left.fy + right.fy;
  var.fxx = left.fxx + right.fxx;
  var.fxy = left.fxy + right.fxy;
  var.fyy = left.fyy + right.fyy;

  var.ft = left.ft + right.ft;
  var.ftt = left.ftt + right.ftt;
  var.fxt = left.fxt + right.fxt;
  var.fyt = left.fyt + right.fyt;
  var.fxxt = left.fxxt + right.fxxt;
  var.fxyt = left.fxyt + right.fxyt;
  var.fyyt = left.fyyt + right.fyyt;
  return var;
}
DualSystem operator-(const DualSystem& left, const DualSystem& right) {
  DualSystem var;
  var.f = left.f - right.f;
  var.fx = left.fx - right.fx;
  var.fy = left.fy - right.fy;
  var.fxx = left.fxx - right.fxx;
  var.fxy = left.fxy - right.fxy;
  var.fyy = left.fyy - right.fyy;

  var.ft = left.ft - right.ft;
  var.ftt = left.ftt - right.ftt;
  var.fxt = left.fxt - right.fxt;
  var.fyt = left.fyt - right.fyt;
  var.fxxt = left.fxxt - right.fxxt;
  var.fxyt = left.fxyt - right.fxyt;
  var.fyyt = left.fyyt - right.fyyt;
  return var;
}
DualSystem operator*(const DualSystem& left, const DualSystem& right) {
  DualSystem var;
  var.f = left.f * right.f;
  var.fx = left.fx * right.f + left.f * right.fx;
  var.fy = left.fy * right.f + left.f * right.fy;
  var.fxx = left.fxx * right.f + 2.0 * left.fx * right.fx + left.f * right.fxx;
  var.fxy = left.fxy * right.f + left.fx * right.fy + left.fy * right.fx +
            left.f * right.fxy;
  var.fyy = left.fyy * right.f + 2.0 * left.fy * right.fy + left.f * right.fyy;

  var.ft = left.ft * right.f + left.f * right.ft;
  var.ftt = left.ftt * right.f + 2.0 * left.ft * right.ft + left.f * right.ftt;
  var.fxt = left.fxt * right.f + left.fx * right.ft + left.ft * right.fx +
            left.f * right.fxt;
  var.fyt = left.fyt * right.f + left.fy * right.ft + left.ft * right.fy +
            left.f * right.fyt;
  var.fxxt = left.fxxt * right.f + left.fxx * right.ft +
             4.0 * left.fx * right.fxt + left.ft * right.fxx +
             left.f * right.fxxt;
  var.fxyt = left.fxyt * right.f + left.fxy * right.ft + left.fxt * right.fy +
             left.fx * right.fyt + left.fyt * right.fx + left.fy * right.fxt +
             left.ft * right.fxy + left.f * right.fxyt;
  var.fyyt = left.fyyt * right.f + left.fyy * right.ft +
             4.0 * left.fy * right.fyt + left.ft * right.fyy +
             left.f * right.fyyt;
  return var;
}
DualSystem operator/(const DualSystem& left, const DualSystem& right) {
  DualSystem var;
  var.f = left.f / right.f;
  var.fx = left.fx / right.f - (right.fx * left.f) / std::pow(right.f, 2);
  var.fy = left.fy / right.f - (right.fy * left.f) / std::pow(right.f, 2);
  var.fxx =
      left.fxx / right.f -
      (left.f * right.fxx + 2.0 * left.fx * right.fx) / std::pow(right.f, 2) +
      (2.0 * std::pow(right.fx, 2) * left.f) / std::pow(right.f, 3);
  var.fxy = left.fxy / right.f -
            (left.f * right.fxy + left.fx * right.fy + left.fy * right.fx) /
                std::pow(right.f, 2) +
            (2.0 * right.fx * right.fy * left.f) / std::pow(right.f, 3);
  var.fyy =
      left.fyy / right.f -
      (left.f * right.fyy + 2.0 * left.fy * right.fy) / std::pow(right.f, 2) +
      (2.0 * std::pow(right.fy, 2) * left.f) / std::pow(right.f, 3);

  var.ft = left.ft / right.f - (right.ft * left.f) / std::pow(right.f, 2);
  var.ftt =
      left.ftt / right.f -
      (left.f * right.ftt + 2.0 * left.ft * right.ft) / std::pow(right.f, 2) +
      (2.0 * std::pow(right.ft, 2) * left.f) / std::pow(right.f, 3);
  var.fxt = left.fxt / right.f -
            (right.fxt * left.f + left.ft * right.fx + left.fx * right.ft) /
                std::pow(right.f, 2) +
            (2.0 * right.ft * right.fx * left.f) / std::pow(right.f, 3);
  var.fyt = left.fyt / right.f -
            (right.fyt * left.f + left.ft * right.fy + left.fy * right.ft) /
                std::pow(right.f, 2) +
            (2.0 * right.ft * right.fy * left.f) / std::pow(right.f, 3);
  var.fxxt =
      left.fxxt / right.f -
      (left.ft * right.fxx + 2.0 * left.fx * right.fxt + right.ft * left.fxx +
       2.0 * right.fx * left.fxt + left.f * right.fxxt) /
          std::pow(right.f, 2) +
      (2.0 * left.ft * std::pow(right.fx, 2) +
       4.0 * left.fx * right.ft * right.fx +
       2.0 * right.ft * right.fxx * left.f +
       4.0 * right.fx * right.fxt * left.f) /
          std::pow(right.f, 3) -
      (6.0 * right.ft * std::pow(right.fx, 2) * left.f) / std::pow(right.f, 4);
  var.fxyt =
      left.fxyt / right.f -
      (left.ft * right.fxy + left.fx * right.fyt + left.fy * right.fxt +
       right.ft * left.fxy + right.fx * left.fyt + right.fy * left.fxt +
       left.f * right.fxyt) /
          std::pow(right.f, 2) +
      (2.0 * left.ft * right.fx * right.fy +
       2.0 * left.fx * right.ft * right.fy +
       2.0 * left.fy * right.ft * right.fx + 2 * right.ft * right.fxy * left.f +
       2.0 * right.fx * right.fyt * left.f +
       2.0 * right.fy * right.fxt * left.f) /
          std::pow(right.f, 3) -
      (6.0 * right.ft * right.fx * right.fy * left.f) / std::pow(right.f, 4);
  var.fyyt =
      left.fyyt / right.f -
      (left.ft * right.fyy + 2.0 * left.fy * right.fyt + right.ft * left.fyy +
       2.0 * right.fy * left.fyt + left.f * right.fyyt) /
          std::pow(right.f, 2) +
      (2.0 * left.ft * std::pow(right.fy, 2) +
       4.0 * left.fy * right.ft * right.fy +
       2.0 * right.ft * right.fyy * left.f +
       4.0 * right.fy * right.fyt * left.f) /
          std::pow(right.f, 3) -
      (6.0 * right.ft * std::pow(right.fy, 2) * left.f) / std::pow(right.f, 4);
  return var;
}

DualSystem std::cos(const DualSystem& var) {
  const double cosvar = std::cos(var.f);
  const double sinvar = std::sin(var.f);
  DualSystem result;
  result.f = cosvar;
  result.fx = -sinvar * var.fx;
  result.fy = -sinvar * var.fy;
  result.fxx = -cosvar * std::pow(var.fx, 2) - sinvar * var.fxx;
  result.fxy = -cosvar * var.fx * var.fy - sinvar * var.fxy;
  result.fyy = -cosvar * std::pow(var.fy, 2) - sinvar * var.fyy;

  result.ft = -sinvar * var.ft;
  result.ftt = -cosvar * std::pow(var.ft, 2) - sinvar * var.ftt;
  result.fxt = -cosvar * var.fx * var.ft - sinvar * var.fxt;
  result.fyt = -cosvar * var.fy * var.ft - sinvar * var.fyt;
  result.fxxt = sinvar * var.ft * std::pow(var.fx, 2) -
                2.0 * cosvar * var.fx * var.fxt - cosvar * var.ft * var.fxx -
                sinvar * var.fxxt;
  result.fxyt = sinvar * var.fx * var.fy * var.ft -
                cosvar * (var.fxt * var.fy + var.fyt * var.fx) -
                cosvar * var.fxy * var.ft - sinvar * var.fxyt;
  result.fyyt = sinvar * var.ft * std::pow(var.fy, 2) -
                2.0 * cosvar * var.fy * var.fyt - cosvar * var.ft * var.fyy -
                sinvar * var.fyyt;
  return result;
}
DualSystem std::sin(const DualSystem& var) {
  const double cosvar = std::cos(var.f);
  const double sinvar = std::sin(var.f);
  DualSystem result;
  result.f = sinvar;
  result.fx = cosvar * var.fx;
  result.fy = cosvar * var.fy;
  result.fxx = -sinvar * std::pow(var.fx, 2) + cosvar * var.fxx;
  result.fxy = -sinvar * var.fx * var.fy + cosvar * var.fxy;
  result.fyy = -sinvar * std::pow(var.fy, 2) + cosvar * var.fyy;

  result.ft = cosvar * var.ft;
  result.ftt = -sinvar * std::pow(var.ft, 2) + cosvar * var.ftt;
  result.fxt = -sinvar * var.fx * var.ft + cosvar * var.fxt;
  result.fyt = -sinvar * var.fy * var.ft + cosvar * var.fyt;
  result.fxxt = -cosvar * var.ft * std::pow(var.fx, 2) -
                2.0 * sinvar * var.fx * var.fxt - sinvar * var.ft * var.fxx +
                cosvar * var.fxxt;
  result.fxyt = -cosvar * var.fx * var.fy * var.ft -
                sinvar * (var.fxt * var.fy + var.fyt * var.fx) -
                sinvar * var.fxy * var.ft + cosvar * var.fxyt;
  result.fyyt = -cosvar * var.ft * std::pow(var.fy, 2) -
                2.0 * sinvar * var.fy * var.fyt - sinvar * var.ft * var.fyy +
                cosvar * var.fyyt;
  return result;
}
DualSystem std::exp(const DualSystem& var) {
  const double expvar = std::cos(var.f);
  DualSystem result;
  result.f = expvar;
  result.fx = expvar * var.fx;
  result.fy = expvar * var.fy;
  result.fxx = expvar * (var.fx * var.fx + var.fxx);
  result.fxy = expvar * (var.fx * var.fy + var.fxy);
  result.fyy = expvar * (var.fy * var.fy + var.fyy);

  result.ft = expvar * var.ft;
  result.ftt = expvar * (var.ft * var.ft + var.ftt);
  result.fxt = expvar * (var.fx * var.ft + var.fxt);
  result.fyt = expvar * (var.fy * var.ft + var.fyt);
  result.fxxt = expvar * var.ft * (var.fx * var.fx + var.fxx) +
                expvar * (2.0 * var.fx * var.fxt + var.fxxt);
  result.fxyt = expvar * var.ft * (var.fx * var.fy + var.fxy) +
                expvar * (var.fxt * var.fy + var.fx * var.fyt + var.fxyt);
  result.fyyt = expvar * var.ft * (var.fy * var.fy + var.fyy) +
                expvar * (2.0 * var.fy * var.fyt + var.fyyt);
  return result;
}