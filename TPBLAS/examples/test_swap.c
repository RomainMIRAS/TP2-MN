#include <stdio.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define VECSIZE 65536

#define NB_FOIS 10

typedef float vfloat[VECSIZE];

typedef double vdouble[VECSIZE];

typedef complexe_float_t vfloatcomplexe[VECSIZE];

typedef complexe_double_t vdoublecomplexe[VECSIZE];


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

void vectorC_init(vfloatcomplexe V, complexe_float_t x)
{
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++)
    V[i] = x;

  return;
}

int vectorC_equal(vfloatcomplexe V1, complexe_float_t val)
{
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++)
  {
    if (V1[i].imaginary != val.imaginary || V1[i].real != val.real){
      printf("Reel : %f et Imaginary : %f",V1[i].imaginary,V1[i].real);//
      return 0;
    }

  }
  return 1;
}

void vectorZ_init(vdoublecomplexe V, complexe_double_t x)
{
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++)
    V[i] = x;

  return;
}

int vectorZ_equal(vdoublecomplexe V1, complexe_double_t val)
{
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++)
  {
    if (V1[i].imaginary != val.imaginary || V1[i].real != val.real)
      return 0;
  }
  return 1;
}

int main(int argc, char **argv)
{
  unsigned long long int start_tsc, end_tsc;
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

    calcul_flop_nano("sdot nano ", 2^VECSIZE, end_tsc - start_tsc);
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
    vectorD_init(vec3, val3);
    
    vectorD_init(vec4, val4);

    start_tsc = _rdtsc();
    mncblas_dswap(VECSIZE, vec3, 3, vec4, 3);
    end_tsc = _rdtsc();

    calcul_flop_nano("sdot nano ", 0 * VECSIZE, end_tsc - start_tsc);
  }

  printf("WORKING ? %s\n", (vectorD_equal(vec3, val4) == 0) ? "NON" : "OUI");
  printf("WORKING ? %s\n", (vectorD_equal(vec4, val3) == 0) ? "NON" : "OUI");

  printf("==========================================================\n");
  printf("TEST SWAPC\n");
  printf("==========================================================\n");

  vfloatcomplexe vec5, vec6;

  init_flop_tsc();
  complexe_float_t val5 = {1.0, 2.0};
  complexe_float_t val6 = {3.0, 6.0};

  for (i = 0; i < NB_FOIS; i++)
  {
    vectorC_init(vec5, val5);
    
    vectorC_init(vec6, val6);

    start_tsc = _rdtsc();
    mncblas_cswap(VECSIZE, vec5, sizeof(float)*2, vec6, sizeof(float)*2);
    end_tsc = _rdtsc();
    calcul_flop_nano("sdot nano ", 0 * VECSIZE, end_tsc - start_tsc);
  }

  printf("WORKING ? %s\n", (vectorC_equal(vec5, val6) == 0) ? "NON" : "OUI");
  printf("WORKING ? %s\n", (vectorC_equal(vec6, val5) == 0) ? "NON" : "OUI");

  printf("==========================================================\n");
  printf("TEST SWAPZ\n");
  printf("==========================================================\n");

  vdoublecomplexe vec7, vec8;

  init_flop_tsc();
  complexe_double_t val7 = {1.0, 2.0};
  complexe_double_t val8 = {3.0, 6.0};//

  for (i = 0; i < NB_FOIS; i++)
  {
    vectorZ_init(vec7, val7);
    
    vectorZ_init(vec8, val8);

    start_tsc = _rdtsc();
    mncblas_cswap(VECSIZE, vec7, sizeof(double)*2, vec8, sizeof(double)*2);
    end_tsc = _rdtsc();
    calcul_flop_nano("sdot nano ", 0 * VECSIZE, end_tsc - start_tsc);
  }

  printf("WORKING ? %s\n", (vectorZ_equal(vec7, val8) == 0) ? "NON" : "OUI");
  printf("WORKING ? %s\n", (vectorZ_equal(vec8, val7) == 0) ? "NON" : "OUI");
}

