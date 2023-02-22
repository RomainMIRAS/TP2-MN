#include <stdio.h>
#include <stdlib.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define NB_FOIS 512

int main(int argc, char **argv)
{
  complexe_float_t c1 = {1.0, 2.0};
  complexe_float_t c2 = {3.0, 6.0};

  complexe_double_t cd1;
  complexe_double_t cd2;

  struct timeval start, end;

  int i;

  init_flop_micro();

  c1 = add_complexe_float(c1, c2);

  printf("c1.r : %f ET c1.i : %f\n", c1.real, c1.imaginary);
  printf("c2.r : %f ET c2.i : %f\n", c2.real, c2.imaginary);

  cd1 = (complexe_double_t){10.0, 7.0};
  cd2 = (complexe_double_t){25.0, 32.0};

  printf("-----------------------\n");
  printf("TEST ADD COMPLEXE DOUBLE\n");
  printf("-----------------------\n");

  TOP_MICRO(start);

  for (i = 0; i < NB_FOIS; i++)
  {
    cd1 = add_complexe_double(cd1, cd2);
  }

  TOP_MICRO(end);

  printf("boucle ADD : cd1.real %f cd1.imaginary %f duree %f \n", cd1.real, cd1.imaginary, tdiff_micro(&start, &end));

  calcul_flop_micro("Calcul complexe ADD : ", NB_FOIS * 2, tdiff_micro(&start, &end));

  printf("-----------------------\n");
  printf("TEST MULT COMPLEXE DOUBLE\n");
  printf("-----------------------\n");

  TOP_MICRO(start);

  for (i = 0; i < NB_FOIS; i++)
  {
    cd1 = mult_complexe_double(cd1, cd2);
  }

  TOP_MICRO(end);

  printf("boucle MULT : cd1.real %f cd1.imaginary %f duree %f \n", cd1.real, cd1.imaginary, tdiff_micro(&start, &end));

  calcul_flop_micro("Calcul complexe MULT : ", NB_FOIS * 2, tdiff_micro(&start, &end));

  printf("-----------------------\n");
  printf("TEST DIV COMPLEXE DOUBLE\n");
  printf("-----------------------\n");

  TOP_MICRO(start);

  for (i = 0; i < NB_FOIS; i++)
  {
    cd1 = div_complexe_double(cd1, cd2);
  }

  TOP_MICRO(end);

  printf("boucle DIV : cd1.real %f cd1.imaginary %f duree %f \n", cd1.real, cd1.imaginary, tdiff_micro(&start, &end));

  calcul_flop_micro("Calcul complexe DIV : ", NB_FOIS * 2, tdiff_micro(&start, &end));

  exit(0);
}