#include "mnblas.h"

void mncblas_sswap(const int N, float *X, const int incX, 
                 float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register float save ;
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = Y [j] ;
      Y [j] = X [i] ;
      X [i] = save ;
    }

  return ;
}

void mncblas_dswap(const int N, double *X, const int incX, 
                 double *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register double save ;
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = Y [j] ;
      Y [j] = X [i] ;
      X [i] = save ;
    }

  return ;
}

void mncblas_cswap(const int N, void *X, const int incX, 
		                    void *Y, const int incY) // float // Réel first iM second
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register float save;

      
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = ((float*)Y )[j] ;
      ((float*)Y )[j] = ((float*)X )[i] ;
      ((float*)X )[i]  = save ;

      save = ((float*)Y )[j+1] ;
      ((float*)Y )[j+1] = ((float*)X )[i+1] ;
      ((float*)X )[i+1] = save ;
    }

      printf("VALUE FIRST OF (3, 6) : %f\n", ((float*)Y)[0]);
      printf("VALUE FIRST OF (1, 2) : %f\n", ((float*)X)[0]);

      printf("VALUE SECOND OF (3, 6) : %f\n", ((float*)Y)[1]);
      printf("VALUE SECOND OF (1, 2) : %f\n", ((float*)X)[1]);
  return ;
}

void mncblas_zswap(const int N, void *X, const int incX, 
		                    void *Y, const int incY) // double
{

  return ;
}

