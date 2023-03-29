#include "mnblas.h"
#include <stdio.h>
#include "complexe.h"

void mncblas_sgemv(const MNCBLAS_LAYOUT layout,
                   const MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                   const float alpha, const float *A, const int lda,
                   const float *X, const int incX, const float beta,
                   float *Y, const int incY)
{

    // Y := alpha*A*x + beta*y

    int i, j;
    float d;

    #pragma omp parallel for private(d, i, j)
    // alpha * A * x
    for (i = 0; i < M; ++i)
    {
        d = 0.0;
        for (j = 0; j < N; ++j){
            d += A[j+i*N] * X[j] * alpha;
        }
        // beta * y
        Y[i] = Y[i]*beta + d;
    }
}

void mncblas_dgemv(MNCBLAS_LAYOUT layout,
                   MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                   const double alpha, const double *A, const int lda,
                   const double *X, const int incX, const double beta,
                   double *Y, const int incY)
{
    int i, j;
    double d;

    for (i = 0; i < N; ++i){
        Y[i] *= beta;
    }

    #pragma omp parallel for private(d, i, j)
    for (i = 0; i < M; ++i)
    {
        d = 0.0;
        // alpha * A * x
        for (j = 0; j < N; ++j){
            d += A[j+i*M] * X[j] * alpha;
        }
        // beta * y
        Y[i] += d;
    }
}

void mncblas_cgemv(MNCBLAS_LAYOUT layout,
                   MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                   const void *alpha, const void *A, const int lda,
                   const void *X, const int incX, const void *beta,
                   void *Y, const int incY)
{
    register const complexe_float_t *a = A;
    register const complexe_float_t *x = X;
    register complexe_float_t *y = Y;
    register complexe_float_t *al = (complexe_float_t *)alpha;
    register complexe_float_t *be = (complexe_float_t *)beta;

    int i, j;
    complexe_float_t d;


    for (i = 0; i < N; ++i){
        y[i] = mult_complexe_float(*be, y[i]);
    }


    #pragma omp parallel for private(d, i, j)
    for (i = 0; i < M; ++i)
    {
        d.real = 0.0;
        d.imaginary = 0.0;
        // alpha * A * x
        for (j = 0; j < N; ++j){
            d = add_complexe_float(d,mult_complexe_float(mult_complexe_float(a[j+i*M], x[j]), *al));
        }
        // beta * y
        y[i] = add_complexe_float(y[i],d);
    }

    Y = y;

}

void mncblas_zgemv(MNCBLAS_LAYOUT layout,
                   MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                   const void *alpha, const void *A, const int lda,
                   const void *X, const int incX, const void *beta,
                   void *Y, const int incY)
{

    register const complexe_double_t *a = A;
    register const complexe_double_t *x = X;
    register complexe_double_t *y = Y;
    register complexe_double_t *al = (complexe_double_t *)alpha;
    register complexe_double_t *be = (complexe_double_t *)beta;

    int i, j;
    complexe_double_t d;


    for (i = 0; i < N; ++i){
        y[i] = mult_complexe_double(*be, y[i]);
    }


    #pragma omp parallel for private(d, i, j)
    for (i = 0; i < M; ++i)
    {
        d.real = 0.0;
        d.imaginary = 0.0;
        // alpha * A * x
        for (j = 0; j < N; ++j){
            d = add_complexe_double(d,mult_complexe_double(mult_complexe_double(a[j+i*M], x[j]), *al));
        }
        // beta * y
        y[i] = add_complexe_double(y[i],d);
    }

    Y = y;
}