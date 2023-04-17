#include "mnblas.h"
#include <stdio.h>
#include "complexe.h"

void mncblas_sgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                   MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                   const int K, const float alpha, const float *A,
                   const int lda, const float *B, const int ldb,
                   const float beta, float *C, const int ldc)
{
    // C = A*B*alpha + C*Beta
    int i, j, k;
    float d;

    #pragma omp parallel for private(d, i, j, k)
    for (i = 0; i < M; i++)
    {
        for (j = 0; j < M; j++)
        {
            d = 0.0;
            for (k = 0; k < M; k++)
                d += A[i * M + k] * B[ k * M + j] * alpha;
            C[i * M + j] = d + beta * C[i * M + j];
        }
    }
}

void mncblas_dgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                   MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                   const int K, const double alpha, const double *A,
                   const int lda, const double *B, const int ldb,
                   const double beta, double *C, const int ldc)
{

    int i, j, k;
    double d;

    #pragma omp parallel for private(d, i, j, k)
    for (i = 0; i < M; ++i)
    {
        for (j = 0; j < M; j++)
        {
            d = 0.0;
            for (k = 0; k < M; k++)
                d += A[i * M + k] * B[k * M + j] * alpha;
            C[i * M + j] = d + beta * C[i * M + j];
        }
    }
}

void mncblas_cgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                   MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                   const int K, const void *alpha, const void *A,
                   const int lda, const void *B, const int ldb,
                   const void *beta, void *C, const int ldc)
{
    int i, j, k;
    complexe_float_t d;
    const complexe_float_t *a = (const complexe_float_t*) A;
    register const complexe_float_t *b = (const complexe_float_t*)B;
    register complexe_float_t *c = ( complexe_float_t*)C;

    register complexe_float_t *al = (complexe_float_t *)alpha;
    register complexe_float_t *be = (complexe_float_t *)beta;


    #pragma omp parallel for private(d, i, j, k)
    for (i = 0; i < M; i++)
    {
        for (j = 0; j < M; j++)
        {
            d .imaginary= 0.0;
            d .real= 0.0;
            
            for (k = 0; k < M; k++){
                d =add_complexe_float(d,mult_complexe_float(mult_complexe_float(a[i * M + k], b[k * M + j]), *al));
            }
            c[i * M + j] = add_complexe_float(d,mult_complexe_float(*be,c[i * M + j]));
        }
    }

    C = c;
}

void mncblas_zgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                   MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                   const int K, const void *alpha, const void *A,
                   const int lda, const void *B, const int ldb,
                   const void *beta, void *C, const int ldc)
{

    int i, j, k;
    complexe_double_t d;
    const complexe_double_t *a = (const complexe_double_t*) A;
    register const complexe_double_t *b = (const complexe_double_t*)B;
    register complexe_double_t *c = ( complexe_double_t*)C;

    register complexe_double_t *al = (complexe_double_t *)alpha;
    register complexe_double_t *be = (complexe_double_t *)beta;


    #pragma omp parallel for private(d, i, j, k)
    for (i = 0; i < M; i++)
    {
        for (j = 0; j < M; j++)
        {
            d .imaginary= 0.0;
            d .real= 0.0;
            
            for (k = 0; k < M; k++){
                d =add_complexe_double(d,mult_complexe_double(mult_complexe_double(a[i * M + k], b[k * M + j]), *al));
            }
            c[i * M + j] = add_complexe_double(d,mult_complexe_double(*be,c[i * M + j]));
        }
    }

    C = c;
}
