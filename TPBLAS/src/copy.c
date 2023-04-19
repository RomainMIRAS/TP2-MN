#include "mnblas.h"
#include "complexe.h"
#include "omp.h"
#include <stdio.h>

void mncblas_scopy(const int N, const float *X, const int incX, 
                 float *Y, const int incY)
{
  register unsigned int i; 
  #pragma omp parallel for
  for (i = 0;i < N; i ++)
    {
      Y [i] = X [i] ;
      printf("Indice %d\n",i);
    }

  return ;
}

void mncblas_dcopy(const int N, const double *X, const int incX, 
                 double *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      Y [j] = X [i] ;
    }

  return ;
}

void mncblas_ccopy(const int N, const void *X, const int incX,
                   void *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  register complexe_double_t * Ytab = Y;
  register const complexe_double_t * Xtab = X;
  
  for (; i < N && j < N ; i += incX, j+=incY)
    {
      Ytab[j] = Xtab[i] ;
    }
  return ;
}

void mncblas_zcopy(const int N, const void *X, const int incX, 
		                    void *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  register complexe_double_t * Ytab = Y;
  register const complexe_double_t * Xtab = X;
  
  for (; i < N && j < N ; i += incX, j+=incY)
    {
      Ytab[j] = Xtab[i] ;
    }
  return ;
}

