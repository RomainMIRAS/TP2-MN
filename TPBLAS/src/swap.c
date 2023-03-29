#include "mnblas.h"
#include "complexe.h"

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
  
  for (; ((i <= N) && (j <= N)) ; i += incX, j+=incY)
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

  register complexe_float_t * Ytab = Y;
  register complexe_float_t * Xtab = X;

  register complexe_float_t saveC;
  
  for (; i < N && j < N ; i += incX, j+=incY)
    {
      saveC = Ytab[j];
      Ytab[j] = Xtab[i] ;
      Xtab[i]  = saveC ;
    }
  return ;
}

void mncblas_zswap(const int N, void *X, const int incX, 
		                    void *Y, const int incY) // double
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  register complexe_double_t * Ytab = Y;
  register complexe_double_t * Xtab = X;

  register complexe_double_t saveC;
  
  for (; i < N && j < N ; i += incX, j+=incY)
    {
      saveC = Ytab[j];
      Ytab[j] = Xtab[i] ;
      Xtab[i]  = saveC ;
    }
  return ;
}

