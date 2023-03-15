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

void vectorC_init(vfloatcomplexe V, complexe_float_t x)
{
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++)
    V[i] = x;

  return;
}

int vectorC_equal(vfloatcomplexe V1, vfloatcomplexe V2)
{
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++)
  {
    if (V1[i].imaginary != V2[i].imaginary || V1[i].real != V2[i].real)
      return 0;
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

int vectorZ_equal(vdoublecomplexe V1, vdoublecomplexe V2)
{
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++)
  {
    if (V1[i].imaginary != V2[i].imaginary || V1[i].real != V2[i].real)
      return 0;
  }
  return 1;
}

int main(int argc, char **argv)
{
  struct timeval start, end;
  unsigned long long int start_tsc, end_tsc;

  float res;
  int i;

  init_flop_tsc();
  printf("==========================================================\n");
  printf("TEST copyS\n");
  printf("==========================================================\n");

  vfloat vec1, vec2;

  init_flop_tsc();
  float val1 = 1.0;
  float val2 = 3.0;

  for (i = 0; i < NB_FOIS; i++)
  {
    vectorF_init(vec1, val1);
    
    vectorF_init(vec2, val2);

     TOP_MICRO(start) ;
        mncblas_scopy(VECSIZE, vec1, 1, vec2, 1);
     TOP_MICRO(end) ;

     calcul_flop_micro ("sdot micro", 2 * VECSIZE, tdiff_micro (&start, &end)) ;
  }

  printf("WORKING ? %s\n", (vectorF_equal(vec1, vec2) == 0) ? "NON" : "OUI");

  printf("==========================================================\n");
  printf("TEST copyD\n");
  printf("==========================================================\n");

  vdouble vec3, vec4;

  init_flop_tsc();
  double val3 = 3.0;
  double val4 = 3.0;

  for (i = 0; i < NB_FOIS; i++)
  {
    vectorF_init(vec3, val3);
    
    vectorF_init(vec4, val4);

     TOP_MICRO(start) ;
        mncblas_dcopy(VECSIZE, vec3, 1, vec4, 1);
     TOP_MICRO(end) ;

     calcul_flop_micro ("sdot micro", 2 * VECSIZE, tdiff_micro (&start, &end)) ;
  }

  printf("WORKING ? %s\n", (vectorF_equal(vec3, vec4) == 0) ? "NON" : "OUI");

  printf("==========================================================\n");
  printf("TEST copyC\n");
  printf("==========================================================\n");

  vfloatcomplexe vec5, vec6;

  init_flop_tsc();
  complexe_float_t val5 = {1.0, 2.0};
  complexe_float_t val6 = {3.0, 6.0};

  for (i = 0; i < NB_FOIS; i++)
  {
    vectorC_init(vec5, val5);
    
    vectorC_init(vec6, val6);

     TOP_MICRO(start) ;
        mncblas_ccopy(VECSIZE, vec5, 1, vec6, 1);
     TOP_MICRO(end) ;

     calcul_flop_micro ("sdot micro", 4 * VECSIZE, tdiff_micro (&start, &end)) ;
  }

  printf("WORKING ? %s\n", (vectorF_equal(vec5, vec6) == 0) ? "NON" : "OUI");

  printf("==========================================================\n");
  printf("TEST copyZ\n");
  printf("==========================================================\n");

  vdoublecomplexe vec7, vec8;

  init_flop_tsc();
  complexe_double_t val7 = {1.0, 2.0};
  complexe_double_t val8 = {3.0, 6.0};

  for (i = 0; i < NB_FOIS; i++)
  {
    vectorZ_init(vec7, val7);
    vectorZ_init(vec8, val8);

     TOP_MICRO(start) ;
        mncblas_zcopy(VECSIZE, vec7, 1, vec8, 1);
     TOP_MICRO(end) ;

     calcul_flop_micro ("sdot micro", 4 * VECSIZE, tdiff_micro (&start, &end)) ;
  }

  printf("WORKING ? %s\n", (vectorF_equal(vec7, vec8) == 0) ? "NON" : "OUI");
}

