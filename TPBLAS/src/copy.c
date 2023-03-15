#include "mnblas.h"

void mncblas_scopy(const int N, const float *X, const int incX, 
                 float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      Y [j] = X [i] ;
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
  
  for (; ((i < N*2) && (j < N*2)) ; i += incX, j+=incY)
    {
      ((float*)Y )[j] = ((float*)X )[i] ;
      ((float*)Y )[j+1] = ((float*)X )[i+1] ;
    }
  return ;
}

void mncblas_zcopy(const int N, const void *X, const int incX, 
		                    void *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  
  for (; ((i < N*2) && (j < N*2)) ; i += incX, j+=incY)
    {
      ((double*)Y )[j] = ((double*)X )[i] ;
      ((double*)Y )[j+1] = ((double*)X )[i+1] ;
    }
  return ;
}

