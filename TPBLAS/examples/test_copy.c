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

int vectorF_equal(vfloat V1, vfloat V2)
{
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++)
  {
    if (V1[i] != V2[i])
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

int vectorD_equal(vdouble V1, vdouble V2)
{
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++)
  {
    if (V1[i] != V2[i])
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
  printf("TEST COPYS\n");
  printf("==========================================================\n");

  vfloat vec1, vec2;

  init_flop_tsc();

  for (i = 0; i < NB_FOIS; i++)
  {
    vectorF_init(vec1, 1.0);
    res = 0.0;

    start_tsc = _rdtsc();
    mncblas_scopy(VECSIZE, vec1, 1, vec2, 1);
    end_tsc = _rdtsc();

    calcul_flop_nano("sdot nano ", 2 * VECSIZE, end_tsc - start_tsc);
  }

  printf("WORKING ? %s\n", (vectorF_equal(vec1, vec2) == 0) ? "NON" : "OUI");
  printf("==========================================================\n");

  printf("==========================================================\n");
  printf("TEST COPYD\n");
  printf("==========================================================\n");

  vdouble vec3, vec4;

  init_flop_tsc();

  for (i = 0; i < NB_FOIS; i++)
  {
    vectorD_init(vec3, 1.0);

    res = 0.0;

    start_tsc = _rdtsc();
    mncblas_dcopy(VECSIZE, vec3, 1, vec4, 1);
    end_tsc = _rdtsc();

    calcul_flop_nano("sdot nano ", 2 * VECSIZE, end_tsc - start_tsc);
  }

  printf("WORKING ? %s\n", (vectorD_equal(vec3, vec4) == 0) ? "NON" : "OUI");
}