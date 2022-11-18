#pragma once

#ifdef USE_OPENBLAS
#include <cblas.h>
#include <lapacke.h>
#else
#define USE_MKL
#include <mkl.h>
#endif

#define gemmAB(alpha, A, B, beta, C, nar, nac, nbc)                           \
  cblas_dgemm(CBLAS_LAYOUT::CblasRowMajor, CBLAS_TRANSPOSE::CblasNoTrans,     \
              CBLAS_TRANSPOSE::CblasNoTrans, nar, nbc, nac, alpha, A, nac, B, \
              nbc, beta, C, nbc)
#define gemmATB(alpha, A, B, beta, C, nar, nac, nbc)                          \
  cblas_dgemm(CBLAS_LAYOUT::CblasRowMajor, CBLAS_TRANSPOSE::CblasTrans,       \
              CBLAS_TRANSPOSE::CblasNoTrans, nac, nbc, nar, alpha, A, nac, B, \
              nbc, beta, C, nbc)
#define gemmABT(alpha, A, B, beta, C, nar, nac, nbr)                        \
  cblas_dgemm(CBLAS_LAYOUT::CblasRowMajor, CBLAS_TRANSPOSE::CblasNoTrans,   \
              CBLAS_TRANSPOSE::CblasTrans, nar, nbr, nac, alpha, A, nac, B, \
              nac, beta, C, nbr)
#define gemmATBT(alpha, A, B, beta, C, nar, nac, nbr)                       \
  cblas_dgemm(CBLAS_LAYOUT::CblasRowMajor, CBLAS_TRANSPOSE::CblasTrans,     \
              CBLAS_TRANSPOSE::CblasTrans, nac, nbr, nar, alpha, A, nac, B, \
              nar, beta, C, nbr)
#define gemmABn(alpha, A, na, B, nb, beta, C, nc, nar, nac, nbc)             \
  cblas_dgemm(CBLAS_LAYOUT::CblasRowMajor, CBLAS_TRANSPOSE::CblasNoTrans,    \
              CBLAS_TRANSPOSE::CblasNoTrans, nar, nbc, nac, alpha, A, na, B, \
              nb, beta, C, nc)
#define gemmATBn(alpha, A, na, B, nb, beta, C, nc, nar, nac, nbc)            \
  cblas_dgemm(CBLAS_LAYOUT::CblasRowMajor, CBLAS_TRANSPOSE::CblasTrans,      \
              CBLAS_TRANSPOSE::CblasNoTrans, nac, nbc, nar, alpha, A, na, B, \
              nb, beta, C, nc)
#define gemmABTn(alpha, A, na, B, nb, beta, C, nc, nar, nac, nbr)              \
  cblas_dgemm(CBLAS_LAYOUT::CblasRowMajor, CBLAS_TRANSPOSE::CblasNoTrans,      \
              CBLAS_TRANSPOSE::CblasTrans, nar, nbr, nac, alpha, A, na, B, nb, \
              beta, C, nc)
#define gemmATBTn(alpha, A, na, B, nb, beta, C, nc, nar, nac, nbr)             \
  cblas_dgemm(CBLAS_LAYOUT::CblasRowMajor, CBLAS_TRANSPOSE::CblasTrans,        \
              CBLAS_TRANSPOSE::CblasTrans, nac, nbr, nar, alpha, A, na, B, nb, \
              beta, C, nc)
#define gemvAx(alpha, A, lda, x, incx, beta, y, incy, nar, nac)                \
  cblas_dgemv(CBLAS_LAYOUT::CblasRowMajor, CBLAS_TRANSPOSE::CblasNoTrans, nar, \
              nac, alpha, A, lda, x, incx, beta, y, incy)
#define gemvATx(alpha, A, lda, x, incx, beta, y, incy, nar, nac)             \
  cblas_dgemv(CBLAS_LAYOUT::CblasRowMajor, CBLAS_TRANSPOSE::CblasTrans, nar, \
              nac, alpha, A, lda, x, incx, beta, y, incy)

#ifndef USE_MKL
void vdSub(const int n, const double* a, const double* b, double* y);
#endif

namespace avocado {
void VecCopy(const int num_index, const int* index, const double* vec1,
             const int stride1, double* vec2, const int stride2);
void VecNormal(const int n, double* vec);
void VecCrossProd(double* result, const double* vec1, const double* vec2);
bool VecIsSame(const int n, const double* vec1, const double* vec2,
               const double err = 1.0e-10);
void VecScale(const int n, const double alpha, double* vec);
double VecInnerProd(const int n, const double* vec1, const double* vec2);
double VecInnerProd(const int n, const double* vec1, const double* vec2,
                    const double* weights);
double VecLength(const int n, const double* vec);
double VecDistance(const int n, const double* vec1, const double* vec2);
double VecAverage(const int n, const double* vec, const int stride);
double VecStdDev(const int n, const double* vec, const int stride);

double ParVecInnerProd(const int n, const double* vec1, const double* vec2);
double ParVecLength(const int n, const double* vec);

void dgetrf(const int n, double* mat, int* ipiv);
void dgetri(const int n, double* mat, const int* ipiv);

void MatInv(const int n, double* mat);
void MatPseudoInv(double* A, int nar, int nac, int& num_ranks,
                  const double eps = 1.0E-9);
double MatDet(const int n, double* mat);
double MatInvDet(const int n, double* mat);

void TensorTranspose(double* tensor, const int rank, const char* index_from,
                     const int* index_size, const char* index_to);

namespace Kernel0 {
void f1(const double* A, const double* B, double* C, const int a, const int b,
        const int c, const double alpha, const double beta);
void f2(const double* A, const double* B, double* C, const int a, const int b,
        const int c, const double alpha, const double beta);
void f3(const double* A, const double* B, double* C, const int a, const int b,
        const int c, const double alpha, const double beta);
void f4(const double* A, const double* B, double* C, const int a, const int b,
        const int c, const double alpha, const double beta);
void f5(const double* A, const double* B, double* C, const int a, const int b,
        const int c, const double alpha, const double beta);
void f6(const double* A, const double* B, double* C, const int a, const int b,
        const int c, const double alpha, const double beta);
void f7(const double* A, const double* B, double* C, const int a, const int b,
        const int c, const double alpha, const double beta);
void f8(const double* A, const double* B, double* C, const int a, const int b,
        const int c, const double alpha, const double beta);
}  // namespace Kernel0

namespace Kernel1 {
void f1(const double* A, const double* B, double* C, const int a, const int b,
        const int c, const int d, const double alpha, const double beta);
void f2(const double* A, const double* B, double* C, const int a, const int b,
        const int c, const int d, const double alpha, const double beta);
void f3(const double* A, const double* B, double* C, const int a, const int b,
        const int c, const int d, const double alpha, const double beta);
void f4(const double* A, const double* B, double* C, const int a, const int b,
        const int c, const int d, const double alpha, const double beta);
void f5(const double* A, const double* B, double* C, const int a, const int b,
        const int c, const int d, const double alpha, const double beta);
void f6(const double* A, const double* B, double* C, const int a, const int b,
        const int c, const int d, const double alpha, const double beta);
void f7(const double* A, const double* B, double* C, const int a, const int b,
        const int c, const int d, const double alpha, const double beta);
void f8(const double* A, const double* B, double* C, const int a, const int b,
        const int c, const int d, const double alpha, const double beta);
void f9(const double* A, const double* B, double* C, const int a, const int b,
        const int c, const int d, const double alpha, const double beta);
void f10(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f11(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f12(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f13(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f14(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f15(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f16(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f17(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f18(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f19(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f20(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f21(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f22(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f23(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f24(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f25(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f26(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f27(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f28(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f29(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f30(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f31(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f32(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f33(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f34(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f35(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f36(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f37(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f38(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f39(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f40(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f41(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f42(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f43(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f44(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f45(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f46(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f47(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f48(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f49(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f50(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f51(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f52(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f53(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f54(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f55(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f56(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f57(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f58(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f59(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f60(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f61(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f62(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f63(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f64(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f65(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f66(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f67(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f68(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f69(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f70(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f71(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f72(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f73(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f74(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f75(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f76(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f77(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f78(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f79(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f80(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
}  // namespace Kernel1

namespace Kernel2 {
void f1(const double* A, const double* B, double* C, const int a, const int b,
        const int c, const int d, const double alpha, const double beta);
void f2(const double* A, const double* B, double* C, const int a, const int b,
        const int c, const int d, const double alpha, const double beta);
void f3(const double* A, const double* B, double* C, const int a, const int b,
        const int c, const int d, const double alpha, const double beta);
void f4(const double* A, const double* B, double* C, const int a, const int b,
        const int c, const int d, const double alpha, const double beta);
void f5(const double* A, const double* B, double* C, const int a, const int b,
        const int c, const int d, const double alpha, const double beta);
void f6(const double* A, const double* B, double* C, const int a, const int b,
        const int c, const int d, const double alpha, const double beta);
void f7(const double* A, const double* B, double* C, const int a, const int b,
        const int c, const int d, const double alpha, const double beta);
void f8(const double* A, const double* B, double* C, const int a, const int b,
        const int c, const int d, const double alpha, const double beta);
void f9(const double* A, const double* B, double* C, const int a, const int b,
        const int c, const int d, const double alpha, const double beta);
void f10(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f11(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f12(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f13(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f14(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f15(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f16(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f17(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f18(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f19(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f20(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f21(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f22(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f23(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f24(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f25(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f26(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f27(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f28(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f29(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f30(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f31(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f32(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f33(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f34(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f35(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f36(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f37(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f38(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f39(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f40(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f41(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f42(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f43(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f44(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f45(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f46(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f47(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f48(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f49(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f50(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f51(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f52(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f53(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f54(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f55(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f56(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f57(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f58(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f59(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f60(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f61(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f62(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f63(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f64(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f65(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f66(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f67(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f68(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f69(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f70(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f71(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f72(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f73(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f74(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f75(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f76(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f77(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f78(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f79(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
void f80(const double* A, const double* B, double* C, const int a, const int b,
         const int c, const int d, const double alpha, const double beta);
}  // namespace Kernel2
}  // namespace avocado