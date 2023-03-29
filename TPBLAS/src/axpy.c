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
  register complexe_float_t Dalpha = *((complexe_float_t *)alpha);
  register complexe_float_t* X2 = X;
  register complexe_float_t* Y2 = Y;
  
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      Y2[i] = add_complexe_float(mult_complexe_float(Dalpha, X2[i]), Y2[i]);
    }

    Y = Y2;
}

void mnblas_zaxpy(const int N, const void *alpha, const void *X,
                  const int incX, void *Y, const int incY)
{

  register complexe_double_t Dalpha = *((complexe_double_t *)alpha);
  register complexe_double_t* X2 = X;
  register complexe_double_t* Y2 = Y;
  
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      Y2[i] = add_complexe_double(mult_complexe_double(Dalpha, X2[i]), Y2[i]);
    }

    Y = Y2;
}