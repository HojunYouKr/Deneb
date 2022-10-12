#include "avocado_blas.h"

#include <cmath>
#include <vector>

#include "avocado_mpi.h"

#ifndef USE_MKL
void vdSub(const int n, const double* a, const double* b, double* y) {
  for (int i = 0; i < n; i++) y[i] = a[i] - b[i];
}
#endif

namespace avocado {
void VecCopy(const int num_index, const int* index, const double* vec1,
             const int stride1, double* vec2, const int stride2) {
  for (int i = 0; i < num_index; i++)
    for (int s = 0; s < stride2; s++)
      vec2[i * stride2 + s] = vec1[index[i] * stride1 + s];
}
void VecNormal(const int n, double* vec) {
  const double length_inv = 1.0 / VecLength(n, vec);
  VecScale(n, length_inv, vec);
}
void VecCrossProd(double* result, const double* vec1, const double* vec2) {
  result[0] = (vec1[1] * vec2[2] - vec1[2] * vec2[1]);
  result[1] = (vec1[2] * vec2[0] - vec1[0] * vec2[2]);
  result[2] = (vec1[0] * vec2[1] - vec1[1] * vec2[0]);
}
bool VecIsSame(const int n, const double* vec1, const double* vec2,
               const double err) {
  for (int i = 0; i < n; i++)
    if (std::abs(vec1[i] - vec2[i]) > err) return false;
  return true;
}
void VecScale(const int n, const double alpha, double* vec) {
  cblas_dscal(n, alpha, vec, 1);
}
double VecInnerProd(const int n, const double* vec1, const double* vec2) {
  return cblas_ddot(n, vec1, 1, vec2, 1);
}
double VecInnerProd(const int n, const double* vec1, const double* vec2,
                    const double* weights) {
  double result = 0;
  for (int i = 0; i < n; i++) result += (vec1[i] * vec2[i] * weights[i]);
  return result;
}
double VecLength(const int n, const double* vec) {
  return std::sqrt(VecInnerProd(n, vec, vec));
}
double VecDistance(const int n, const double* vec1, const double* vec2) {
  double result = 0;
  for (int i = 0; i < n; i++)
    result += (vec2[i] - vec1[i]) * (vec2[i] - vec1[i]);
  return std::sqrt(result);
}
double VecAverage(const int n, const double* vec, const int stride) {
  double avg = 0.0;
  for (int i = 0; i < n; i++) avg += vec[i * stride];
  return avg / static_cast<double>(n);
}
double VecStdDev(const int n, const double* vec, const int stride) {
  const double avg = VecAverage(n, vec, stride);
  double avg2 = 0.0;
  for (int i = 0; i < n; i++) avg2 += (vec[i * stride] * vec[i * stride]);
  return std::sqrt(avg2 / static_cast<double>(n) - avg * avg);
}
double ParVecInnerProd(const int n, const double* vec1, const double* vec2) {
  const double value = VecInnerProd(n, vec1, vec2);
  return AVOCADO_MPI->Reduce(value, MPI::Op::SUM);
}
double ParVecLength(const int n, const double* vec) {
  return std::sqrt(ParVecInnerProd(n, vec, vec));
}

void dgetrf(const int n, double* mat, int* ipiv) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < i; j++) {
      for (int k = 0; k < j; k++)
        mat[i * n + k] -= mat[i * n + j] * mat[j * n + k];
      for (int k = j + 1; k < n; k++)
        mat[i * n + k] -= mat[i * n + j] * mat[j * n + k];
      mat[i * n + j] *= (-mat[j * n + j]);
    }

    double maxvar = std::abs(mat[i * n + i]);
    ipiv[i] = i;
    for (int j = i + 1; j < n; j++) {
      if (maxvar < std::abs(mat[i * n + j])) {
        maxvar = std::abs(mat[i * n + j]);
        ipiv[i] = j;
      }
    }
    if (ipiv[i] != i) {
      for (int j = 0; j < n; j++) {
        const double temp = mat[j * n + i];
        mat[j * n + i] = mat[j * n + ipiv[i]];
        mat[j * n + ipiv[i]] = temp;
      }
    }

    mat[i * n + i] = 1.0 / mat[i * n + i];
    for (int j = 0; j < i; j++) mat[i * n + j] *= mat[i * n + i];
    for (int j = i + 1; j < n; j++) mat[i * n + j] *= mat[i * n + i];
  }
}
void dgetri(const int n, double* mat, const int* ipiv) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < i; j++) {
      for (int k = 0; k < i; k++)
        mat[j * n + k] -= mat[j * n + i] * mat[i * n + k];
      for (int k = i + 1; k < n; k++)
        mat[j * n + k] -= mat[j * n + i] * mat[i * n + k];
      mat[j * n + i] *= (-mat[i * n + i]);
    }
  }
  for (int i = n - 1; i >= 0; i--) {
    if (ipiv[i] != i) {
      for (int j = 0; j < n; j++) {
        const double temp = mat[i * n + j];
        mat[i * n + j] = mat[ipiv[i] * n + j];
        mat[ipiv[i] * n + j] = temp;
      }
    }
  }
}
void MatInv(const int n, double* mat) {
  std::vector<int> ipiv(n, 0);
  dgetrf(n, mat, &ipiv[0]);
  dgetri(n, mat, &ipiv[0]);
}
double MatDet(const int n, double* mat) {
  std::vector<int> ipiv(n, 0);
  dgetrf(n, mat, &ipiv[0]);

  double det = 1.0;
  for (int i = 0; i < n; i++) {
    det *= mat[i * n + i];
    if (ipiv[i] != i) det = -det;
  };
  return 1.0 / det;
}
double MatInvDet(const int n, double* mat) {
  std::vector<int> ipiv(n, 0);
  dgetrf(n, mat, &ipiv[0]);

  double det = 1.0;
  for (int i = 0; i < n; i++) {
    det *= mat[i * n + i];
    if (ipiv[i] != i) det = -det;
  };
  dgetri(n, mat, &ipiv[0]);
  return 1.0 / det;
}

void MatPseudoInv(double* A, const int nar, const int nac, int& num_ranks,
                  const double eps) {
  num_ranks = -1;
  std::vector<double> vt(nar * nac);
  std::vector<double> u(nar * nar);
  std::vector<double> s(nar);
  std::vector<double> superb(nar, 0.0);
  int info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'S', nar, nac, &A[0], nac,
                            &s[0], &u[0], nar, &vt[0], nac, &superb[0]);

  num_ranks = 0;
  for (auto&& si : s) {
    if (std::abs(si) > eps) {
      si = 1.0 / si;
      num_ranks++;
    } else {
      si = 0.0;
    }
  }
  for (int i = 0; i < num_ranks; i++) {
    for (int j = 0; j < nac; j++) vt[i * nac + j] *= s[i];
  }
  for (int i = num_ranks * nac; i < nar * nac; i++) vt[i] = 0.0;
  gemmATBT(1.0, &vt[0], &u[0], 0.0, &A[0], nar, nac, nar);
}

void TensorTranspose(double* tensor, const int rank, const char* index_from,
                     const int* index_size, const char* index_to) {
  std::vector<int> mapping(rank, -1);  // to -> from
  for (int i = 0; i < rank; i++) {
    for (int j = 0; j < rank; j++) {
      if (index_to[i] == index_from[j]) {
        mapping[i] = j;
        break;
      }
    }
    if (mapping[i] == -1) {
      ERROR_MESSAGE("Index for tensor transpose is not matched.\n");
      return;
    }
  }

  std::vector<int> strands_from(rank, 1);
  for (int i = 1; i < rank; i++)
    strands_from[i] = strands_from[i - 1] * index_size[rank - i];
  std::vector<int> strands_to(rank, 1);
  for (int i = 1; i < rank; i++)
    strands_to[i] = strands_to[i - 1] * index_size[mapping[rank - i]];

  int total_size = 1;
  for (int i = 0; i < rank; i++) total_size *= index_size[i];
  std::vector<double> temp(total_size);
  for (int i = 0; i < total_size; i++) temp[i] = tensor[i];

  std::vector<int> ijk(rank);
  for (int i = 0; i < total_size; i++) {
    int var_temp = i;
    for (int i = rank - 1; i >= 0; i--) {
      ijk[i] = var_temp / strands_from[i];
      var_temp -= ijk[i] * strands_from[i];
    }
    int ito = 0;
    for (int i = 0; i < rank; i++)
      ito += strands_to[i] * ijk[rank - 1 - mapping[rank - 1 - i]];
    tensor[ito] = temp[i];
  }
}

// ----------------------------------------- Kernel 0
// ----------------------------------------- //
void Kernel0::f1(const double* A, const double* B, double* C, const int a,
                 const int b, const int c, const double alpha,
                 const double beta) {
  gemmAB(alpha, &A[0], &B[0], beta, &C[0], a, b, c);
}
void Kernel0::f2(const double* A, const double* B, double* C, const int a,
                 const int b, const int c, const double alpha,
                 const double beta) {
  gemmATBT(alpha, &B[0], &A[0], beta, &C[0], b, c, a);
}
void Kernel0::f3(const double* A, const double* B, double* C, const int a,
                 const int b, const int c, const double alpha,
                 const double beta) {
  gemmABT(alpha, &A[0], &B[0], beta, &C[0], a, b, c);
}
void Kernel0::f4(const double* A, const double* B, double* C, const int a,
                 const int b, const int c, const double alpha,
                 const double beta) {
  gemmABT(alpha, &B[0], &A[0], beta, &C[0], c, b, a);
}
void Kernel0::f5(const double* A, const double* B, double* C, const int a,
                 const int b, const int c, const double alpha,
                 const double beta) {
  gemmATB(alpha, &A[0], &B[0], beta, &C[0], b, a, c);
}
void Kernel0::f6(const double* A, const double* B, double* C, const int a,
                 const int b, const int c, const double alpha,
                 const double beta) {
  gemmATB(alpha, &B[0], &A[0], beta, &C[0], b, c, a);
}
void Kernel0::f7(const double* A, const double* B, double* C, const int a,
                 const int b, const int c, const double alpha,
                 const double beta) {
  gemmATBT(alpha, &A[0], &B[0], beta, &C[0], b, a, c);
}
void Kernel0::f8(const double* A, const double* B, double* C, const int a,
                 const int b, const int c, const double alpha,
                 const double beta) {
  gemmAB(alpha, &B[0], &A[0], beta, &C[0], c, b, a);
}

// ----------------------------------------- Kernel 1
// ----------------------------------------- //
void Kernel1::f1(const double* A, const double* B, double* C, const int a,
                 const int b, const int c, const int d, const double alpha,
                 const double beta) {
  const int ab = a * b;
  gemmAB(alpha, &A[0], &B[0], beta, &C[0], ab, c, d);
}
void Kernel1::f2(const double* A, const double* B, double* C, const int a,
                 const int b, const int c, const int d, const double alpha,
                 const double beta) {
  for (int i = 0; i < a; i++)
    gemmAB(alpha, &A[i * b * c], &B[0], beta, &C[i * b * d], b, c, d);
}
void Kernel1::f3(const double* A, const double* B, double* C, const int a,
                 const int b, const int c, const int d, const double alpha,
                 const double beta) {
  const int bc = b * c;
  const int bd = b * d;
  for (int i = 0; i < b; i++)
    gemmABn(alpha, &A[i * c], bc, &B[0], d, beta, &C[i * d], bd, a, c, d);
}
void Kernel1::f4(const double* A, const double* B, double* C, const int a,
                 const int b, const int c, const int d, const double alpha,
                 const double beta) {
  for (int i = 0; i < a; i++)
    gemmATBT(alpha, &B[0], &A[i * b * c], beta, &C[i * d * b], c, d, b);
}
void Kernel1::f5(const double* A, const double* B, double* C, const int a,
                 const int b, const int c, const int d, const double alpha,
                 const double beta) {
  const int ad = a * d;
  for (int i = 0; i < a; i++)
    gemmABn(alpha, &A[i * b * c], c, &B[0], d, beta, &C[i * d], ad, b, c, d);
}
void Kernel1::f6(const double* A, const double* B, double* C, const int a,
                 const int b, const int c, const int d, const double alpha,
                 const double beta) {
  const int bc = b * c;
  for (int i = 0; i < b; i++)
    gemmABn(alpha, &A[i * c], bc, &B[0], d, beta, &C[i * a * d], d, a, c, d);
}
void Kernel1::f7(const double* A, const double* B, double* C, const int a,
                 const int b, const int c, const int d, const double alpha,
                 const double beta) {
  const int bc = b * c;
  for (int i = 0; i < b; i++)
    gemmATBTn(alpha, &B[0], d, &A[i * c], bc, beta, &C[i * d * a], a, c, d, a);
}
void Kernel1::f8(const double* A, const double* B, double* C, const int a,
                 const int b, const int c, const int d, const double alpha,
                 const double beta) {
  const int ab = a * b;
  gemmATBT(alpha, &B[0], &A[0], beta, &C[0], c, d, ab);
}
void Kernel1::f9(const double* A, const double* B, double* C, const int a,
                 const int b, const int c, const int d, const double alpha,
                 const double beta) {
  const int ab = a * b;
  for (int i = 0; i < a; i++)
    gemmATBTn(alpha, &B[0], d, &A[i * b * c], c, beta, &C[i * b], ab, c, d, b);
}
void Kernel1::f10(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bc = b * c;
  const int ba = b * a;
  for (int i = 0; i < b; i++)
    gemmATBTn(alpha, &B[0], d, &A[i * c], bc, beta, &C[i * a], ba, c, d, a);
}
void Kernel1::f11(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ab = a * b;
  gemmABT(alpha, &A[0], &B[0], beta, &C[0], ab, c, d);
}
void Kernel1::f12(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  for (int i = 0; i < a; i++)
    gemmABT(alpha, &A[i * b * c], &B[0], beta, &C[i * b * d], b, c, d);
}
void Kernel1::f13(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bc = b * c;
  const int bd = b * d;
  for (int i = 0; i < b; i++)
    gemmABTn(alpha, &A[i * c], bc, &B[0], c, beta, &C[i * d], bd, a, c, d);
}
void Kernel1::f14(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  for (int i = 0; i < a; i++)
    gemmABT(alpha, &B[0], &A[i * b * c], beta, &C[i * d * b], d, c, b);
}
void Kernel1::f15(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ad = a * d;
  for (int i = 0; i < a; i++)
    gemmABTn(alpha, &A[i * b * c], c, &B[0], c, beta, &C[i * d], ad, b, c, d);
}
void Kernel1::f16(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bc = b * c;
  for (int i = 0; i < b; i++)
    gemmABTn(alpha, &A[i * c], bc, &B[0], c, beta, &C[i * a * d], d, a, c, d);
}
void Kernel1::f17(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bc = b * c;
  for (int i = 0; i < b; i++)
    gemmABTn(alpha, &B[0], c, &A[i * c], bc, beta, &C[i * d * a], a, d, c, a);
}
void Kernel1::f18(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ab = a * b;
  gemmABT(alpha, &B[0], &A[0], beta, &C[0], d, c, ab);
}
void Kernel1::f19(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ab = a * b;
  for (int i = 0; i < a; i++)
    gemmABTn(alpha, &B[0], c, &A[i * b * c], c, beta, &C[i * b], ab, d, c, b);
}
void Kernel1::f20(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bc = b * c;
  const int ba = b * a;
  for (int i = 0; i < b; i++)
    gemmABTn(alpha, &B[0], c, &A[i * c], bc, beta, &C[i * a], ba, d, c, a);
}
void Kernel1::f21(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  for (int i = 0; i < a; i++)
    gemmATB(alpha, &A[i * c * b], &B[0], beta, &C[i * b * d], c, b, d);
}
void Kernel1::f22(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  for (int i = 0; i < a; i++)
    gemmATB(alpha, &B[0], &A[i * c * b], beta, &C[i * d * b], c, d, b);
}
void Kernel1::f23(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ad = a * d;
  for (int i = 0; i < a; i++)
    gemmATBn(alpha, &A[i * c * b], b, &B[0], d, beta, &C[i * d], ad, c, b, d);
}
void Kernel1::f24(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ab = a * b;
  for (int i = 0; i < a; i++)
    gemmATBn(alpha, &B[0], d, &A[i * c * b], b, beta, &C[i * b], ab, c, d, b);
}
void Kernel1::f25(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  for (int i = 0; i < a; i++)
    gemmATBT(alpha, &A[i * c * b], &B[0], beta, &C[i * b * d], c, b, d);
}
void Kernel1::f26(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  for (int i = 0; i < a; i++)
    gemmAB(alpha, &B[0], &A[i * c * b], beta, &C[i * d * b], d, c, b);
}
void Kernel1::f27(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ad = a * d;
  for (int i = 0; i < a; i++)
    gemmATBTn(alpha, &A[i * c * b], b, &B[0], c, beta, &C[i * d], ad, c, b, d);
}
void Kernel1::f28(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ab = a * b;
  for (int i = 0; i < a; i++)
    gemmABn(alpha, &B[0], c, &A[i * c * b], b, beta, &C[i * b], ab, d, c, b);
}
void Kernel1::f29(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ac = a * c;
  for (int i = 0; i < a; i++)
    gemmABn(alpha, &A[i * c], ac, &B[0], d, beta, &C[i * b * d], d, b, c, d);
}
void Kernel1::f30(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bd = b * d;
  for (int i = 0; i < b; i++)
    gemmABn(alpha, &A[i * a * c], c, &B[0], d, beta, &C[i * d], bd, a, c, d);
}
void Kernel1::f31(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ac = a * c;
  for (int i = 0; i < a; i++)
    gemmATBTn(alpha, &B[0], d, &A[i * c], ac, beta, &C[i * d * b], b, c, d, b);
}
void Kernel1::f32(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ab = a * b;
  gemmAB(alpha, &A[0], &B[0], beta, &C[0], ab, c, d);
}
void Kernel1::f33(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ac = a * c;
  const int ad = a * d;
  for (int i = 0; i < a; i++)
    gemmABn(alpha, &A[i * c], ac, &B[0], d, beta, &C[i * d], ad, b, c, d);
}
void Kernel1::f34(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  for (int i = 0; i < b; i++)
    gemmAB(alpha, &A[i * a * c], &B[0], beta, &C[i * a * d], a, c, d);
}
void Kernel1::f35(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  for (int i = 0; i < b; i++)
    gemmATBT(alpha, &B[0], &A[i * a * c], beta, &C[i * d * a], c, d, a);
}
void Kernel1::f36(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ac = a * c;
  const int ab = a * b;
  for (int i = 0; i < a; i++)
    gemmATBTn(alpha, &B[0], d, &A[i * c], ac, beta, &C[i * b], ab, c, d, b);
}
void Kernel1::f37(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ab = a * b;
  gemmATBT(alpha, &B[0], &A[0], beta, &C[0], c, d, ab);
}
void Kernel1::f38(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ba = b * a;
  for (int i = 0; i < b; i++)
    gemmATBTn(alpha, &B[0], d, &A[i * a * c], c, beta, &C[i * a], ba, c, d, a);
}
void Kernel1::f39(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ac = a * c;
  for (int i = 0; i < a; i++)
    gemmABTn(alpha, &A[i * c], ac, &B[0], c, beta, &C[i * b * d], d, b, c, d);
}
void Kernel1::f40(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bd = b * d;
  for (int i = 0; i < b; i++)
    gemmABTn(alpha, &A[i * a * c], c, &B[0], c, beta, &C[i * d], bd, a, c, d);
}
void Kernel1::f41(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ac = a * c;
  for (int i = 0; i < a; i++)
    gemmABTn(alpha, &B[0], c, &A[i * c], ac, beta, &C[i * d * b], b, d, c, b);
}
void Kernel1::f42(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ab = a * b;
  gemmABT(alpha, &A[0], &B[0], beta, &C[0], ab, c, d);
}
void Kernel1::f43(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ac = a * c;
  const int ad = a * d;
  for (int i = 0; i < a; i++)
    gemmABTn(alpha, &A[i * c], ac, &B[0], c, beta, &C[i * d], ad, b, c, d);
}
void Kernel1::f44(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  for (int i = 0; i < b; i++)
    gemmABT(alpha, &A[i * a * c], &B[0], beta, &C[i * a * d], a, c, d);
}
void Kernel1::f45(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  for (int i = 0; i < b; i++)
    gemmABT(alpha, &B[0], &A[i * a * c], beta, &C[i * d * a], d, c, a);
}
void Kernel1::f46(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ac = a * c;
  const int ab = a * b;
  for (int i = 0; i < a; i++)
    gemmABTn(alpha, &B[0], c, &A[i * c], ac, beta, &C[i * b], ab, d, c, b);
}
void Kernel1::f47(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ab = a * b;
  gemmABT(alpha, &B[0], &A[0], beta, &C[0], d, c, ab);
}
void Kernel1::f48(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ba = b * a;
  for (int i = 0; i < b; i++)
    gemmABTn(alpha, &B[0], c, &A[i * a * c], c, beta, &C[i * a], ba, d, c, a);
}
void Kernel1::f49(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bd = b * d;
  for (int i = 0; i < b; i++)
    gemmATBn(alpha, &A[i * c * a], a, &B[0], d, beta, &C[i * d], bd, c, a, d);
}
void Kernel1::f50(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  for (int i = 0; i < b; i++)
    gemmATB(alpha, &A[i * c * a], &B[0], beta, &C[i * a * d], c, a, d);
}
void Kernel1::f51(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  for (int i = 0; i < b; i++)
    gemmATB(alpha, &B[0], &A[i * c * a], beta, &C[i * d * a], c, d, a);
}
void Kernel1::f52(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ba = b * a;
  for (int i = 0; i < b; i++)
    gemmATBn(alpha, &B[0], d, &A[i * c * a], a, beta, &C[i * a], ba, c, d, a);
}
void Kernel1::f53(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bd = b * d;
  for (int i = 0; i < b; i++)
    gemmATBTn(alpha, &A[i * c * a], a, &B[0], c, beta, &C[i * d], bd, c, a, d);
}
void Kernel1::f54(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  for (int i = 0; i < b; i++)
    gemmATBT(alpha, &A[i * c * a], &B[0], beta, &C[i * a * d], c, a, d);
}
void Kernel1::f55(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  for (int i = 0; i < b; i++)
    gemmAB(alpha, &B[0], &A[i * c * a], beta, &C[i * d * a], d, c, a);
}
void Kernel1::f56(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ba = b * a;
  for (int i = 0; i < b; i++)
    gemmABn(alpha, &B[0], c, &A[i * c * a], a, beta, &C[i * a], ba, d, c, a);
}
void Kernel1::f57(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ab = a * b;
  gemmATB(alpha, &A[0], &B[0], beta, &C[0], c, ab, d);
}
void Kernel1::f58(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ab = a * b;
  for (int i = 0; i < a; i++)
    gemmATBn(alpha, &A[i * b], ab, &B[0], d, beta, &C[i * b * d], d, c, b, d);
}
void Kernel1::f59(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ab = a * b;
  for (int i = 0; i < a; i++)
    gemmATBn(alpha, &B[0], d, &A[i * b], ab, beta, &C[i * d * b], b, c, d, b);
}
void Kernel1::f60(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ab = a * b;
  const int ad = a * d;
  for (int i = 0; i < a; i++)
    gemmATBn(alpha, &A[i * b], ab, &B[0], d, beta, &C[i * d], ad, c, b, d);
}
void Kernel1::f61(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ab = a * b;
  gemmATB(alpha, &B[0], &A[0], beta, &C[0], c, d, ab);
}
void Kernel1::f62(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ab = a * b;
  for (int i = 0; i < a; i++)
    gemmATBn(alpha, &B[0], d, &A[i * b], ab, beta, &C[i * b], ab, c, d, b);
}
void Kernel1::f63(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ab = a * b;
  gemmATBT(alpha, &A[0], &B[0], beta, &C[0], c, ab, d);
}
void Kernel1::f64(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ab = a * b;
  for (int i = 0; i < a; i++)
    gemmATBTn(alpha, &A[i * b], ab, &B[0], c, beta, &C[i * b * d], d, c, b, d);
}
void Kernel1::f65(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ab = a * b;
  for (int i = 0; i < a; i++)
    gemmABn(alpha, &B[0], c, &A[i * b], ab, beta, &C[i * d * b], b, d, c, b);
}
void Kernel1::f66(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ab = a * b;
  const int ad = a * d;
  for (int i = 0; i < a; i++)
    gemmATBTn(alpha, &A[i * b], ab, &B[0], c, beta, &C[i * d], ad, c, b, d);
}
void Kernel1::f67(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ab = a * b;
  gemmAB(alpha, &B[0], &A[0], beta, &C[0], d, c, ab);
}
void Kernel1::f68(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ab = a * b;
  for (int i = 0; i < a; i++)
    gemmABn(alpha, &B[0], c, &A[i * b], ab, beta, &C[i * b], ab, d, c, b);
}
void Kernel1::f69(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ba = b * a;
  const int bd = b * d;
  for (int i = 0; i < b; i++)
    gemmATBn(alpha, &A[i * a], ba, &B[0], d, beta, &C[i * d], bd, c, a, d);
}
void Kernel1::f70(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ab = a * b;
  gemmATB(alpha, &A[0], &B[0], beta, &C[0], c, ab, d);
}
void Kernel1::f71(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ba = b * a;
  for (int i = 0; i < b; i++)
    gemmATBn(alpha, &A[i * a], ba, &B[0], d, beta, &C[i * a * d], d, c, a, d);
}
void Kernel1::f72(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ba = b * a;
  for (int i = 0; i < b; i++)
    gemmATBn(alpha, &B[0], d, &A[i * a], ba, beta, &C[i * d * a], a, c, d, a);
}
void Kernel1::f73(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ab = a * b;
  gemmATB(alpha, &B[0], &A[0], beta, &C[0], c, d, ab);
}
void Kernel1::f74(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ba = b * a;
  for (int i = 0; i < b; i++)
    gemmATBn(alpha, &B[0], d, &A[i * a], ba, beta, &C[i * a], ba, c, d, a);
}
void Kernel1::f75(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ba = b * a;
  const int bd = b * d;
  for (int i = 0; i < b; i++)
    gemmATBTn(alpha, &A[i * a], ba, &B[0], c, beta, &C[i * d], bd, c, a, d);
}
void Kernel1::f76(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ab = a * b;
  gemmATBT(alpha, &A[0], &B[0], beta, &C[0], c, ab, d);
}
void Kernel1::f77(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ba = b * a;
  for (int i = 0; i < b; i++)
    gemmATBTn(alpha, &A[i * a], ba, &B[0], c, beta, &C[i * a * d], d, c, a, d);
}
void Kernel1::f78(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ba = b * a;
  for (int i = 0; i < b; i++)
    gemmABn(alpha, &B[0], c, &A[i * a], ba, beta, &C[i * d * a], a, d, c, a);
}
void Kernel1::f79(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ab = a * b;
  gemmAB(alpha, &B[0], &A[0], beta, &C[0], d, c, ab);
}
void Kernel1::f80(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ba = b * a;
  for (int i = 0; i < b; i++)
    gemmABn(alpha, &B[0], c, &A[i * a], ba, beta, &C[i * a], ba, d, c, a);
}

// ----------------------------------------- Kernel 2
// ----------------------------------------- //
void Kernel2::f1(const double* A, const double* B, double* C, const int a,
                 const int b, const int c, const int d, const double alpha,
                 const double beta) {
  const int bc = b * c;
  gemmAB(alpha, &A[0], &B[0], beta, &C[0], a, bc, d);
}
void Kernel2::f2(const double* A, const double* B, double* C, const int a,
                 const int b, const int c, const int d, const double alpha,
                 const double beta) {
  const int bc = b * c;
  gemmABn(alpha, &A[0], bc, &B[0], d, beta, &C[0], d, a, c, d);
  for (int i = 1; i < b; i++)
    gemmABn(alpha, &A[i * c], bc, &B[i * c * d], d, 1.0, &C[0], d, a, c, d);
}
void Kernel2::f3(const double* A, const double* B, double* C, const int a,
                 const int b, const int c, const int d, const double alpha,
                 const double beta) {
  const int bc = b * c;
  gemmATBT(alpha, &B[0], &A[0], beta, &C[0], bc, d, a);
}
void Kernel2::f4(const double* A, const double* B, double* C, const int a,
                 const int b, const int c, const int d, const double alpha,
                 const double beta) {
  const int bc = b * c;
  gemmATBTn(alpha, &B[0], d, &A[0], bc, beta, &C[0], a, c, d, a);
  for (int i = 1; i < b; i++)
    gemmATBTn(alpha, &B[i * c * d], d, &A[i * c], bc, 1.0, &C[0], a, c, d, a);
}
void Kernel2::f5(const double* A, const double* B, double* C, const int a,
                 const int b, const int c, const int d, const double alpha,
                 const double beta) {
  const int bc = b * c;
  gemmABTn(alpha, &A[0], bc, &B[0], c, beta, &C[0], d, a, c, d);
  for (int i = 1; i < b; i++)
    gemmABTn(alpha, &A[i * c], bc, &B[i * d * c], c, 1.0, &C[0], d, a, c, d);
}
void Kernel2::f6(const double* A, const double* B, double* C, const int a,
                 const int b, const int c, const int d, const double alpha,
                 const double beta) {
  const int bc = b * c;
  gemmABTn(alpha, &B[0], c, &A[0], bc, beta, &C[0], a, d, c, a);
  for (int i = 1; i < b; i++)
    gemmABTn(alpha, &B[i * d * c], c, &A[i * c], bc, 1.0, &C[0], a, d, c, a);
}
void Kernel2::f7(const double* A, const double* B, double* C, const int a,
                 const int b, const int c, const int d, const double alpha,
                 const double beta) {
  const int bc = b * c;
  const int bd = b * d;
  gemmABn(alpha, &A[0], bc, &B[0], bd, beta, &C[0], d, a, c, d);
  for (int i = 1; i < b; i++)
    gemmABn(alpha, &A[i * c], bc, &B[i * d], bd, 1.0, &C[0], d, a, c, d);
}
void Kernel2::f8(const double* A, const double* B, double* C, const int a,
                 const int b, const int c, const int d, const double alpha,
                 const double beta) {
  const int bc = b * c;
  const int bd = b * d;
  gemmATBTn(alpha, &B[0], bd, &A[0], bc, beta, &C[0], a, c, d, a);
  for (int i = 1; i < b; i++)
    gemmATBTn(alpha, &B[i * d], bd, &A[i * c], bc, 1.0, &C[0], a, c, d, a);
}
void Kernel2::f9(const double* A, const double* B, double* C, const int a,
                 const int b, const int c, const int d, const double alpha,
                 const double beta) {
  const int bc = b * c;
  gemmABT(alpha, &A[0], &B[0], beta, &C[0], a, bc, d);
}
void Kernel2::f10(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bc = b * c;
  gemmABTn(alpha, &A[0], bc, &B[0], bc, beta, &C[0], d, a, c, d);
  for (int i = 1; i < b; i++)
    gemmABTn(alpha, &A[i * c], bc, &B[i * c], bc, 1.0, &C[0], d, a, c, d);
}
void Kernel2::f11(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bc = b * c;
  gemmABT(alpha, &B[0], &A[0], beta, &C[0], d, bc, a);
}
void Kernel2::f12(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bc = b * c;
  gemmABTn(alpha, &B[0], bc, &A[0], bc, beta, &C[0], a, d, c, a);
  for (int i = 1; i < b; i++)
    gemmABTn(alpha, &B[i * c], bc, &A[i * c], bc, 1.0, &C[0], a, d, c, a);
}
void Kernel2::f13(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int cb = c * b;
  const int cd = c * d;
  gemmABn(alpha, &A[0], cb, &B[0], cd, beta, &C[0], d, a, b, d);
  for (int i = 1; i < c; i++)
    gemmABn(alpha, &A[i * b], cb, &B[i * d], cd, 1.0, &C[0], d, a, b, d);
}
void Kernel2::f14(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int cb = c * b;
  const int cd = c * d;
  gemmATBTn(alpha, &B[0], cd, &A[0], cb, beta, &C[0], a, b, d, a);
  for (int i = 1; i < c; i++)
    gemmATBTn(alpha, &B[i * d], cd, &A[i * b], cb, 1.0, &C[0], a, b, d, a);
}
void Kernel2::f15(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bc = b * c;
  gemmAB(alpha, &A[0], &B[0], beta, &C[0], a, bc, d);
}
void Kernel2::f16(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int cb = c * b;
  gemmABn(alpha, &A[0], cb, &B[0], d, beta, &C[0], d, a, b, d);
  for (int i = 1; i < c; i++)
    gemmABn(alpha, &A[i * b], cb, &B[i * b * d], d, 1.0, &C[0], d, a, b, d);
}
void Kernel2::f17(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bc = b * c;
  gemmATBT(alpha, &B[0], &A[0], beta, &C[0], bc, d, a);
}
void Kernel2::f18(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int cb = c * b;
  gemmATBTn(alpha, &B[0], d, &A[0], cb, beta, &C[0], a, b, d, a);
  for (int i = 1; i < c; i++)
    gemmATBTn(alpha, &B[i * b * d], d, &A[i * b], cb, 1.0, &C[0], a, b, d, a);
}
void Kernel2::f19(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int cb = c * b;
  gemmABTn(alpha, &A[0], cb, &B[0], b, beta, &C[0], d, a, b, d);
  for (int i = 1; i < c; i++)
    gemmABTn(alpha, &A[i * b], cb, &B[i * d * b], b, 1.0, &C[0], d, a, b, d);
}
void Kernel2::f20(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int cb = c * b;
  gemmABTn(alpha, &B[0], b, &A[0], cb, beta, &C[0], a, d, b, a);
  for (int i = 1; i < c; i++)
    gemmABTn(alpha, &B[i * d * b], b, &A[i * b], cb, 1.0, &C[0], a, d, b, a);
}
void Kernel2::f21(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bc = b * c;
  gemmABT(alpha, &A[0], &B[0], beta, &C[0], a, bc, d);
}
void Kernel2::f22(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int cb = c * b;
  gemmABTn(alpha, &A[0], cb, &B[0], cb, beta, &C[0], d, a, b, d);
  for (int i = 1; i < c; i++)
    gemmABTn(alpha, &A[i * b], cb, &B[i * b], cb, 1.0, &C[0], d, a, b, d);
}
void Kernel2::f23(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bc = b * c;
  gemmABT(alpha, &B[0], &A[0], beta, &C[0], d, bc, a);
}
void Kernel2::f24(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int cb = c * b;
  gemmABTn(alpha, &B[0], cb, &A[0], cb, beta, &C[0], a, d, b, a);
  for (int i = 1; i < c; i++)
    gemmABTn(alpha, &B[i * b], cb, &A[i * b], cb, 1.0, &C[0], a, d, b, a);
}
void Kernel2::f25(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  gemmAB(alpha, &A[0], &B[0], beta, &C[0], a, c, d);
  for (int i = 1; i < b; i++)
    gemmAB(alpha, &A[i * a * c], &B[i * c * d], 1.0, &C[0], a, c, d);
}
void Kernel2::f26(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  gemmATBT(alpha, &B[0], &A[0], beta, &C[0], c, d, a);
  for (int i = 1; i < b; i++)
    gemmATBT(alpha, &B[i * c * d], &A[i * a * c], 1.0, &C[0], c, d, a);
}
void Kernel2::f27(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  gemmABT(alpha, &A[0], &B[0], beta, &C[0], a, c, d);
  for (int i = 1; i < b; i++)
    gemmABT(alpha, &A[i * a * c], &B[i * d * c], 1.0, &C[0], a, c, d);
}
void Kernel2::f28(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  gemmABT(alpha, &B[0], &A[0], beta, &C[0], d, c, a);
  for (int i = 1; i < b; i++)
    gemmABT(alpha, &B[i * d * c], &A[i * a * c], 1.0, &C[0], d, c, a);
}
void Kernel2::f29(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bd = b * d;
  gemmABn(alpha, &A[0], c, &B[0], bd, beta, &C[0], d, a, c, d);
  for (int i = 1; i < b; i++)
    gemmABn(alpha, &A[i * a * c], c, &B[i * d], bd, 1.0, &C[0], d, a, c, d);
}
void Kernel2::f30(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bd = b * d;
  gemmATBTn(alpha, &B[0], bd, &A[0], c, beta, &C[0], a, c, d, a);
  for (int i = 1; i < b; i++)
    gemmATBTn(alpha, &B[i * d], bd, &A[i * a * c], c, 1.0, &C[0], a, c, d, a);
}
void Kernel2::f31(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bc = b * c;
  gemmABTn(alpha, &A[0], c, &B[0], bc, beta, &C[0], d, a, c, d);
  for (int i = 1; i < b; i++)
    gemmABTn(alpha, &A[i * a * c], c, &B[i * c], bc, 1.0, &C[0], d, a, c, d);
}
void Kernel2::f32(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bc = b * c;
  gemmABTn(alpha, &B[0], bc, &A[0], c, beta, &C[0], a, d, c, a);
  for (int i = 1; i < b; i++)
    gemmABTn(alpha, &B[i * c], bc, &A[i * a * c], c, 1.0, &C[0], a, d, c, a);
}
void Kernel2::f33(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bc = b * c;
  gemmATB(alpha, &A[0], &B[0], beta, &C[0], bc, a, d);
}
void Kernel2::f34(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  gemmATB(alpha, &A[0], &B[0], beta, &C[0], c, a, d);
  for (int i = 1; i < b; i++)
    gemmATB(alpha, &A[i * c * a], &B[i * c * d], 1.0, &C[0], c, a, d);
}
void Kernel2::f35(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ca = c * a;
  const int cd = c * d;
  gemmATBn(alpha, &A[0], ca, &B[0], cd, beta, &C[0], d, b, a, d);
  for (int i = 1; i < c; i++)
    gemmATBn(alpha, &A[i * a], ca, &B[i * d], cd, 1.0, &C[0], d, b, a, d);
}
void Kernel2::f36(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bc = b * c;
  gemmATB(alpha, &B[0], &A[0], beta, &C[0], bc, d, a);
}
void Kernel2::f37(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  gemmATB(alpha, &B[0], &A[0], beta, &C[0], c, d, a);
  for (int i = 1; i < b; i++)
    gemmATB(alpha, &B[i * c * d], &A[i * c * a], 1.0, &C[0], c, d, a);
}
void Kernel2::f38(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ca = c * a;
  const int cd = c * d;
  gemmATBn(alpha, &B[0], cd, &A[0], ca, beta, &C[0], a, b, d, a);
  for (int i = 1; i < c; i++)
    gemmATBn(alpha, &B[i * d], cd, &A[i * a], ca, 1.0, &C[0], a, b, d, a);
}
void Kernel2::f39(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  gemmATBT(alpha, &A[0], &B[0], beta, &C[0], c, a, d);
  for (int i = 1; i < b; i++)
    gemmATBT(alpha, &A[i * c * a], &B[i * d * c], 1.0, &C[0], c, a, d);
}
void Kernel2::f40(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  gemmAB(alpha, &B[0], &A[0], beta, &C[0], d, c, a);
  for (int i = 1; i < b; i++)
    gemmAB(alpha, &B[i * d * c], &A[i * c * a], 1.0, &C[0], d, c, a);
}
void Kernel2::f41(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bd = b * d;
  gemmATBn(alpha, &A[0], a, &B[0], bd, beta, &C[0], d, c, a, d);
  for (int i = 1; i < b; i++)
    gemmATBn(alpha, &A[i * c * a], a, &B[i * d], bd, 1.0, &C[0], d, c, a, d);
}
void Kernel2::f42(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ca = c * a;
  gemmATBn(alpha, &A[0], ca, &B[0], d, beta, &C[0], d, b, a, d);
  for (int i = 1; i < c; i++)
    gemmATBn(alpha, &A[i * a], ca, &B[i * b * d], d, 1.0, &C[0], d, b, a, d);
}
void Kernel2::f43(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bd = b * d;
  gemmATBn(alpha, &B[0], bd, &A[0], a, beta, &C[0], a, c, d, a);
  for (int i = 1; i < b; i++)
    gemmATBn(alpha, &B[i * d], bd, &A[i * c * a], a, 1.0, &C[0], a, c, d, a);
}
void Kernel2::f44(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ca = c * a;
  gemmATBn(alpha, &B[0], d, &A[0], ca, beta, &C[0], a, b, d, a);
  for (int i = 1; i < c; i++)
    gemmATBn(alpha, &B[i * b * d], d, &A[i * a], ca, 1.0, &C[0], a, b, d, a);
}
void Kernel2::f45(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ca = c * a;
  gemmATBTn(alpha, &A[0], ca, &B[0], b, beta, &C[0], d, b, a, d);
  for (int i = 1; i < c; i++)
    gemmATBTn(alpha, &A[i * a], ca, &B[i * d * b], b, 1.0, &C[0], d, b, a, d);
}
void Kernel2::f46(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ca = c * a;
  gemmABn(alpha, &B[0], b, &A[0], ca, beta, &C[0], a, d, b, a);
  for (int i = 1; i < c; i++)
    gemmABn(alpha, &B[i * d * b], b, &A[i * a], ca, 1.0, &C[0], a, d, b, a);
}
void Kernel2::f47(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bc = b * c;
  gemmATBT(alpha, &A[0], &B[0], beta, &C[0], bc, a, d);
}
void Kernel2::f48(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bc = b * c;
  gemmATBTn(alpha, &A[0], a, &B[0], bc, beta, &C[0], d, c, a, d);
  for (int i = 1; i < b; i++)
    gemmATBTn(alpha, &A[i * c * a], a, &B[i * c], bc, 1.0, &C[0], d, c, a, d);
}
void Kernel2::f49(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bc = b * c;
  gemmAB(alpha, &B[0], &A[0], beta, &C[0], d, bc, a);
}
void Kernel2::f50(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bc = b * c;
  gemmABn(alpha, &B[0], bc, &A[0], a, beta, &C[0], a, d, c, a);
  for (int i = 1; i < b; i++)
    gemmABn(alpha, &B[i * c], bc, &A[i * c * a], a, 1.0, &C[0], a, d, c, a);
}
void Kernel2::f51(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ca = c * a;
  const int cb = c * b;
  gemmATBTn(alpha, &A[0], ca, &B[0], cb, beta, &C[0], d, b, a, d);
  for (int i = 1; i < c; i++)
    gemmATBTn(alpha, &A[i * a], ca, &B[i * b], cb, 1.0, &C[0], d, b, a, d);
}
void Kernel2::f52(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ca = c * a;
  const int cb = c * b;
  gemmABn(alpha, &B[0], cb, &A[0], ca, beta, &C[0], a, d, b, a);
  for (int i = 1; i < c; i++)
    gemmABn(alpha, &B[i * b], cb, &A[i * a], ca, 1.0, &C[0], a, d, b, a);
}
void Kernel2::f53(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int cd = c * d;
  gemmABn(alpha, &A[0], b, &B[0], cd, beta, &C[0], d, a, b, d);
  for (int i = 1; i < c; i++)
    gemmABn(alpha, &A[i * a * b], b, &B[i * d], cd, 1.0, &C[0], d, a, b, d);
}
void Kernel2::f54(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int cd = c * d;
  gemmATBTn(alpha, &B[0], cd, &A[0], b, beta, &C[0], a, b, d, a);
  for (int i = 1; i < c; i++)
    gemmATBTn(alpha, &B[i * d], cd, &A[i * a * b], b, 1.0, &C[0], a, b, d, a);
}
void Kernel2::f55(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  gemmAB(alpha, &A[0], &B[0], beta, &C[0], a, b, d);
  for (int i = 1; i < c; i++)
    gemmAB(alpha, &A[i * a * b], &B[i * b * d], 1.0, &C[0], a, b, d);
}
void Kernel2::f56(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  gemmATBT(alpha, &B[0], &A[0], beta, &C[0], b, d, a);
  for (int i = 1; i < c; i++)
    gemmATBT(alpha, &B[i * b * d], &A[i * a * b], 1.0, &C[0], b, d, a);
}
void Kernel2::f57(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  gemmABT(alpha, &A[0], &B[0], beta, &C[0], a, b, d);
  for (int i = 1; i < c; i++)
    gemmABT(alpha, &A[i * a * b], &B[i * d * b], 1.0, &C[0], a, b, d);
}
void Kernel2::f58(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  gemmABT(alpha, &B[0], &A[0], beta, &C[0], d, b, a);
  for (int i = 1; i < c; i++)
    gemmABT(alpha, &B[i * d * b], &A[i * a * b], 1.0, &C[0], d, b, a);
}
void Kernel2::f59(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int cb = c * b;
  gemmABTn(alpha, &A[0], b, &B[0], cb, beta, &C[0], d, a, b, d);
  for (int i = 1; i < c; i++)
    gemmABTn(alpha, &A[i * a * b], b, &B[i * b], cb, 1.0, &C[0], d, a, b, d);
}
void Kernel2::f60(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int cb = c * b;
  gemmABTn(alpha, &B[0], cb, &A[0], b, beta, &C[0], a, d, b, a);
  for (int i = 1; i < c; i++)
    gemmABTn(alpha, &B[i * b], cb, &A[i * a * b], b, 1.0, &C[0], a, d, b, a);
}
void Kernel2::f61(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ba = b * a;
  gemmATBn(alpha, &A[0], ba, &B[0], d, beta, &C[0], d, c, a, d);
  for (int i = 1; i < b; i++)
    gemmATBn(alpha, &A[i * a], ba, &B[i * c * d], d, 1.0, &C[0], d, c, a, d);
}
void Kernel2::f62(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int cd = c * d;
  gemmATBn(alpha, &A[0], a, &B[0], cd, beta, &C[0], d, b, a, d);
  for (int i = 1; i < c; i++)
    gemmATBn(alpha, &A[i * b * a], a, &B[i * d], cd, 1.0, &C[0], d, b, a, d);
}
void Kernel2::f63(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ba = b * a;
  gemmATBn(alpha, &B[0], d, &A[0], ba, beta, &C[0], a, c, d, a);
  for (int i = 1; i < b; i++)
    gemmATBn(alpha, &B[i * c * d], d, &A[i * a], ba, 1.0, &C[0], a, c, d, a);
}
void Kernel2::f64(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int cd = c * d;
  gemmATBn(alpha, &B[0], cd, &A[0], a, beta, &C[0], a, b, d, a);
  for (int i = 1; i < c; i++)
    gemmATBn(alpha, &B[i * d], cd, &A[i * b * a], a, 1.0, &C[0], a, b, d, a);
}
void Kernel2::f65(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ba = b * a;
  gemmATBTn(alpha, &A[0], ba, &B[0], c, beta, &C[0], d, c, a, d);
  for (int i = 1; i < b; i++)
    gemmATBTn(alpha, &A[i * a], ba, &B[i * d * c], c, 1.0, &C[0], d, c, a, d);
}
void Kernel2::f66(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ba = b * a;
  gemmABn(alpha, &B[0], c, &A[0], ba, beta, &C[0], a, d, c, a);
  for (int i = 1; i < b; i++)
    gemmABn(alpha, &B[i * d * c], c, &A[i * a], ba, 1.0, &C[0], a, d, c, a);
}
void Kernel2::f67(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bc = b * c;
  gemmATB(alpha, &A[0], &B[0], beta, &C[0], bc, a, d);
}
void Kernel2::f68(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ba = b * a;
  const int bd = b * d;
  gemmATBn(alpha, &A[0], ba, &B[0], bd, beta, &C[0], d, c, a, d);
  for (int i = 1; i < b; i++)
    gemmATBn(alpha, &A[i * a], ba, &B[i * d], bd, 1.0, &C[0], d, c, a, d);
}
void Kernel2::f69(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  gemmATB(alpha, &A[0], &B[0], beta, &C[0], b, a, d);
  for (int i = 1; i < c; i++)
    gemmATB(alpha, &A[i * b * a], &B[i * b * d], 1.0, &C[0], b, a, d);
}
void Kernel2::f70(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bc = b * c;
  gemmATB(alpha, &B[0], &A[0], beta, &C[0], bc, d, a);
}
void Kernel2::f71(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ba = b * a;
  const int bd = b * d;
  gemmATBn(alpha, &B[0], bd, &A[0], ba, beta, &C[0], a, c, d, a);
  for (int i = 1; i < b; i++)
    gemmATBn(alpha, &B[i * d], bd, &A[i * a], ba, 1.0, &C[0], a, c, d, a);
}
void Kernel2::f72(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  gemmATB(alpha, &B[0], &A[0], beta, &C[0], b, d, a);
  for (int i = 1; i < c; i++)
    gemmATB(alpha, &B[i * b * d], &A[i * b * a], 1.0, &C[0], b, d, a);
}
void Kernel2::f73(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  gemmATBT(alpha, &A[0], &B[0], beta, &C[0], b, a, d);
  for (int i = 1; i < c; i++)
    gemmATBT(alpha, &A[i * b * a], &B[i * d * b], 1.0, &C[0], b, a, d);
}
void Kernel2::f74(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  gemmAB(alpha, &B[0], &A[0], beta, &C[0], d, b, a);
  for (int i = 1; i < c; i++)
    gemmAB(alpha, &B[i * d * b], &A[i * b * a], 1.0, &C[0], d, b, a);
}
void Kernel2::f75(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ba = b * a;
  const int bc = b * c;
  gemmATBTn(alpha, &A[0], ba, &B[0], bc, beta, &C[0], d, c, a, d);
  for (int i = 1; i < b; i++)
    gemmATBTn(alpha, &A[i * a], ba, &B[i * c], bc, 1.0, &C[0], d, c, a, d);
}
void Kernel2::f76(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int ba = b * a;
  const int bc = b * c;
  gemmABn(alpha, &B[0], bc, &A[0], ba, beta, &C[0], a, d, c, a);
  for (int i = 1; i < b; i++)
    gemmABn(alpha, &B[i * c], bc, &A[i * a], ba, 1.0, &C[0], a, d, c, a);
}
void Kernel2::f77(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bc = b * c;
  gemmATBT(alpha, &A[0], &B[0], beta, &C[0], bc, a, d);
}
void Kernel2::f78(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int cb = c * b;
  gemmATBTn(alpha, &A[0], a, &B[0], cb, beta, &C[0], d, b, a, d);
  for (int i = 1; i < c; i++)
    gemmATBTn(alpha, &A[i * b * a], a, &B[i * b], cb, 1.0, &C[0], d, b, a, d);
}
void Kernel2::f79(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int bc = b * c;
  gemmAB(alpha, &B[0], &A[0], beta, &C[0], d, bc, a);
}
void Kernel2::f80(const double* A, const double* B, double* C, const int a,
                  const int b, const int c, const int d, const double alpha,
                  const double beta) {
  const int cb = c * b;
  gemmABn(alpha, &B[0], cb, &A[0], a, beta, &C[0], a, d, b, a);
  for (int i = 1; i < c; i++)
    gemmABn(alpha, &B[i * b], cb, &A[i * b * a], a, 1.0, &C[0], a, d, b, a);
}
}  // namespace avocado