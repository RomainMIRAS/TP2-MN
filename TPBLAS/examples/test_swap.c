#include <stdio.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define VECSIZE 65536

#define NB_FOIS 10

typedef float vfloat[VECSIZE];

typedef double vdouble[VECSIZE];

void vectorF_init(vfloat V, float x)
{
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++)
    V[i] = x;

  return;
}

int vectorF_equal(vfloat V1, float val)
{
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++)
  {
    if (V1[i] != val)
      return 0;
  }
  return 1;
}

void vectorD_init(vdouble V, double x)
{
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++)
    V[i] = x;

  return;
}

int vectorD_equal(vdouble V1, double val)
{
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++)
  {
    if (V1[i] != val)
      return 0;
  }
  return 1;
}

// void vector_print (vfloat V)
// {
//   register unsigned int i ;

//   for (i = 0; i < VECSIZE; i++)
//     printf ("%f ", V[i]) ;
//   printf ("\n") ;

//   return ;
// }

int main(int argc, char **argv)
{
  struct timeval start, end;
  unsigned long long int start_tsc, end_tsc;

  float res;
  int i;

  init_flop_tsc();
  printf("==========================================================\n");
  printf("TEST SWAPS\n");
  printf("==========================================================\n");

  vfloat vec1, vec2;

  init_flop_tsc();
  float val1 = 1.0;
  float val2 = 3.0;

  for (i = 0; i < NB_FOIS; i++)
  {
    vectorF_init(vec1, val1);
    
    vectorF_init(vec2, val2);

    start_tsc = _rdtsc();
    mncblas_sswap(VECSIZE, vec1, 1, vec2, 1);
    end_tsc = _rdtsc();

    calcul_flop_nano("sdot nano ", 2 * VECSIZE, end_tsc - start_tsc);
  }

  printf("WORKING ? %s\n", (vectorF_equal(vec1, val2) == 0) ? "NON" : "OUI");
  printf("WORKING ? %s\n", (vectorF_equal(vec2, val1) == 0) ? "NON" : "OUI");

  printf("==========================================================\n");
  printf("TEST SWAPD\n");
  printf("==========================================================\n");

  vdouble vec3, vec4;

  init_flop_tsc();
  double val3 = 3.0;
  double val4 = 3.0;

  for (i = 0; i < NB_FOIS; i++)
  {
    vectorF_init(vec3, val3);
    
    vectorF_init(vec4, val4);

    start_tsc = _rdtsc();
    mncblas_sswap(VECSIZE, vec3, 3, vec4, 3);
    end_tsc = _rdtsc();

    calcul_flop_nano("sdot nano ", 4 * VECSIZE, end_tsc - start_tsc);
  }

  printf("WORKING ? %s\n", (vectorF_equal(vec3, val4) == 0) ? "NON" : "OUI");
  printf("WORKING ? %s\n", (vectorF_equal(vec4, val3) == 0) ? "NON" : "OUI");
}