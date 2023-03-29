#include "mnblas.h"
#include "complexe.h"

void mnblas_saxpy(const int N, const float alpha, const float *X, const int incX, float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
    Y[i] = alpha * X[i] + Y[i];
    }
}

void mnblas_daxpy(const int N, const double alpha, const double *X,
                  const int incX, double *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
    Y[i] = alpha * X[i] + Y[i];
    }
}

void mnblas_caxpy(const int N, const void *alpha, const void *X,
                  const int incX, void *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  float Falpha = (*(float*)alpha);
  
  for (; ((i <= N) && (j <= N)) ; i += incX, j+=incY)
    {
      ((float*)Y )[j] = ((float*)X )[i] * Falpha + ((float*)Y )[j];
      ((float*)Y )[j+1] = ((float*)X )[i+1] * Falpha + ((float*)Y )[j+1];
    }
  return ;
}

void mnblas_zaxpy(const int N, const void *alpha, const void *X,
                  const int incX, void *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  complexe_double_t Dalpha = *((complexe_double_t *)alpha);
  
  for (; ((i <= N) && (j <= N)) ; i += incX, j+=incY)
    {
      ((complexe_double_t*)Y )[j] = mult_complexe_double((*(complexe_double_t *)X[i]) ,alpha) + ((double*)Y )[j];
    }
  return;
}