#include "mnblas.h"
#include <stdio.h>
#include "complexe.h"

/*
float mncblas_sdot(const int N, const float *X, const int incX, 
                 const float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register float dot = 0.0 ;
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      dot = dot + X [i] * Y [j] ;
    }

  return dot ;
}
*/

float mncblas_sdot(const int N, const float *X, const int incX, 
                 const float *Y, const int incY)
{
  float dot = 0;
  int i ;
  #pragma parallel omp for private (i) reduction (+:dot)  
  for (i = 0;i < N; i++)
    {
      dot += Y[i]*X[i];
    }
  return dot;
}

double mncblas_ddot(const int N, const double *X, const int incX, 
                 const double *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  double dot = 0;
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      dot += Y[j]*X[i];
    }
  return dot;
}

void   mncblas_cdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu) // complex float
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  const register complexe_float_t * Ytab = Y;
  const register complexe_float_t * Xtab = X;

  complexe_float_t *res = dotu;
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      *res = add_complexe_float(*res, mult_complexe_float(Ytab[j], Xtab[i]));
    }
  return;
}

void   mncblas_cdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  const register complexe_float_t * Ytab = Y;
  const register complexe_float_t * Xtab = X;

  complexe_float_t *res = dotc;
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      *res = add_complexe_float(*res, mult_complexe_float(Ytab[j], conjuguer_complexe_float(Xtab[i])));
    }
  return;
}

void   mncblas_zdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  const register complexe_double_t * Ytab = Y;
  const register complexe_double_t * Xtab = X;

  complexe_double_t *res = dotu;
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      *res = add_complexe_double(*res, mult_complexe_double(Ytab[j], Xtab[i]));
    }
  return;
}
  
void   mncblas_zdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  const register complexe_double_t * Ytab = Y;
  const register complexe_double_t * Xtab = X;

  complexe_double_t *res = dotc;
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      *res = add_complexe_double(*res, mult_complexe_double(Ytab[j], conjuguer_complexe_double(Xtab[i])));
    }
  return;
}




