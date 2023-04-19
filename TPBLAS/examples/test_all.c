#include <stdio.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define VECSIZE 512
#define NB_FOIS 1024

typedef float vfloat[VECSIZE];

typedef double vdouble[VECSIZE];

typedef float mfloat[VECSIZE][VECSIZE];

typedef double mdouble[VECSIZE][VECSIZE];

typedef complexe_float_t vfloatcomplexe[VECSIZE];

typedef complexe_double_t vdoublecomplexe[VECSIZE];

typedef complexe_float_t mfloatcomplexe[VECSIZE][VECSIZE];

typedef complexe_double_t mdoublecomplexe[VECSIZE][VECSIZE];

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

void matrixF_init(mfloat M, float x)
{
    register unsigned int i, j;

    for (i = 0; i < VECSIZE; i++)
        for (j = 0; j < VECSIZE; j++)
            M[i][j] = x;

    return;
}

void matrixD_init(mdouble M, double x)
{
    register unsigned int i, j;

    for (i = 0; i < VECSIZE; i++)
        for (j = 0; j < VECSIZE; j++)
            M[i][j] = x;

    return;
}

void matrixC_init(mfloatcomplexe M, complexe_float_t x)
{
    register unsigned int i, j;

    for (i = 0; i < VECSIZE; i++)
        for (j = 0; j < VECSIZE; j++)
            M[i][j] = x;

    return;
}

void matrixZ_init(mdoublecomplexe M, complexe_double_t x)
{
    register unsigned int i, j;

    for (i = 0; i < VECSIZE; i++)
        for (j = 0; j < VECSIZE; j++)
            M[i][j] = x;

    return;
}

void show_matrixF(mfloat M)
{
    register unsigned int i, j;

    for (i = 0; i < VECSIZE; i++)
    {
        for (j = 0; j < VECSIZE; j++)
        {
            printf("%f ", M[i][j]);
        }
        printf("\n");
    }
}

int main(int argc, char **argv)
{
    struct timeval start, end;

    init_flop_tsc();
    printf("==========================================================\n");
    printf("TEST F BLAS1 \n");
    printf("==========================================================\n");

    vfloat vec1, vec2;

    init_flop_tsc();
    float val1 = 1.0;
    float val2 = 3.0;

    vectorF_init(vec1, val1);

    vectorF_init(vec2, val2);

    TOP_MICRO(start);
    for(size_t i = 0; i < NB_FOIS; i++)
        mncblas_scopy(VECSIZE, vec1, 1, vec2, 1);
    TOP_MICRO(end);

    // printf("SHOW F VEC1\n");
    // for (size_t i = 0; i < VECSIZE; i++)
    // {
    //     printf("VEC 1:%ld : %f\n", i, vec1[i]);
    // }

    // printf("SHOW F VEC2\n");
    // for (size_t i = 0; i < VECSIZE; i++)
    // {
    //     printf("VEC 2:%ld : %f\n", i, vec2[i]);
    // }

    calcul_flop_micro("scopy micro", 2 * VECSIZE*NB_FOIS, tdiff_micro(&start, &end));

    // printf("==========================================================\n");
    // printf("TEST D\n");
    // printf("==========================================================\n");

    // vdouble vec3, vec4;

    // init_flop_tsc();
    // double val3 = 3.0;
    // double val4 = 5.0;

    // vectorD_init(vec3, val3);

    // vectorD_init(vec4, val4);

    // TOP_MICRO(start);
    // mncblas_dcopy(VECSIZE, vec3, 1, vec4, 1);
    // TOP_MICRO(end);

    // printf("SHOW D VEC1\n");
    // for (size_t i = 0; i < VECSIZE; i++)
    // {
    //     printf("VEC 1:%ld : %f\n", i, vec3[i]);
    // }

    // printf("SHOW D VEC2\n");
    // for (size_t i = 0; i < VECSIZE; i++)
    // {
    //     printf("VEC 2:%ld : %f\n", i, vec4[i]);
    // }

    // calcul_flop_micro("sdot micro", 2 * VECSIZE, tdiff_micro(&start, &end));

    // printf("==========================================================\n");
    // printf("TEST C\n");
    // printf("==========================================================\n");

    // vfloatcomplexe vec5, vec6;

    // init_flop_tsc();
    // complexe_float_t val5 = {1.0, 2.0};
    // complexe_float_t val6 = {3.0, 6.0};

    // vectorC_init(vec5, val5);

    // vectorC_init(vec6, val6);

    // TOP_MICRO(start);
    // mncblas_ccopy(VECSIZE, vec5, 1, vec6, 1);
    // TOP_MICRO(end);

    // printf("SHOW C VEC1\n");
    // for (size_t i = 0; i < VECSIZE; i++)
    // {
    //     printf("VEC 1:%ld : %f ; %fi\n", i, vec5[i].real, vec5[i].imaginary);
    // }

    // printf("SHOW C VEC2\n");
    // for (size_t i = 0; i < VECSIZE; i++)
    // {
    //     printf("VEC 2:%ld : %f ; %fi\n", i, vec6[i].real, vec6[i].imaginary);
    // }

    // printf("==========================================================\n");
    // printf("TEST Z\n");
    // printf("==========================================================\n");

    // vdoublecomplexe vec7, vec8;

    // init_flop_tsc();
    // complexe_double_t val7 = {1.0, 2.0};
    // complexe_double_t val8 = {3.0, 6.0};

    // vectorZ_init(vec7, val7);
    // vectorZ_init(vec8, val8);

    // TOP_MICRO(start);
    // mncblas_zcopy(VECSIZE, vec7, 1, vec8, 1);
    // TOP_MICRO(end);

    // printf("SHOW Z VEC1\n");
    // for (size_t i = 0; i < VECSIZE; i++)
    // {
    //     printf("VEC 1:%ld : %f ; %fi\n", i, vec7[i].real, vec7[i].imaginary);
    // }

    // printf("SHOW Z VEC2\n");
    // for (size_t i = 0; i < VECSIZE; i++)
    // {
    //     printf("VEC 2:%ld : %f ; %fi\n", i, vec8[i].real, vec8[i].imaginary);
    // }


    // init_flop_tsc();
    // printf("==========================================================\n");
    // printf("TEST SWAP FLOAT\n");
    // printf("==========================================================\n");

    // init_flop_tsc();
    // val1 = 1.0;
    // val2 = 3.0;

    // vectorF_init(vec1, val1);

    // vectorF_init(vec2, val2);

    // TOP_MICRO(start);
    // mncblas_sswap(VECSIZE, vec1, 1, vec2, 1);
    // TOP_MICRO(end);

    // printf("SHOW F VEC1\n");
    // for (size_t i = 0; i < VECSIZE; i++)
    // {
    //     printf("VEC 1:%ld : %f\n", i, vec1[i]);
    // }

    // printf("SHOW F VEC2\n");
    // for (size_t i = 0; i < VECSIZE; i++)
    // {
    //     printf("VEC 2:%ld : %f\n", i, vec2[i]);
    // }

    // calcul_flop_micro("sdot micro", 2 * VECSIZE, tdiff_micro(&start, &end));

    // printf("==========================================================\n");
    // printf("TEST SWAP DOUBLE\n");
    // printf("==========================================================\n");

    // init_flop_tsc();
    // val3 = 1.0;
    // val4 = 3.0;

    // vectorD_init(vec3, val3);

    // vectorD_init(vec4, val4);

    // TOP_MICRO(start);
    // mncblas_dswap(VECSIZE, vec3, 1, vec4, 1);
    // TOP_MICRO(end);

    // printf("SHOW D VEC1\n");
    // for (size_t i = 0; i < VECSIZE; i++)
    // {
    //     printf("VEC 1:%ld : %f\n", i, vec3[i]);
    // }

    // printf("SHOW D VEC2\n");
    // for (size_t i = 0; i < VECSIZE; i++)
    // {
    //     printf("VEC 2:%ld : %f\n", i, vec4[i]);
    // }

    // calcul_flop_micro("sdot micro", 2 * VECSIZE, tdiff_micro(&start, &end));

    // printf("==========================================================\n");
    // printf("TEST SWAP COMPLEX FLOAT\n");
    // printf("==========================================================\n");

    // init_flop_tsc();

    // vectorC_init(vec5, val5);

    // vectorC_init(vec6, val6);

    // TOP_MICRO(start);
    // mncblas_cswap(VECSIZE, vec5, 1, vec6, 1);
    // TOP_MICRO(end);

    // printf("SHOW C VEC1\n");
    // for (size_t i = 0; i < VECSIZE; i++)
    // {
    //     printf("VEC 1:%ld : %f ; %fi\n", i, vec5[i].real, vec5[i].imaginary);
    // }

    // printf("SHOW C VEC2\n");
    // for (size_t i = 0; i < VECSIZE; i++)
    // {
    //     printf("VEC 2:%ld : %f ; %fi\n", i, vec6[i].real, vec6[i].imaginary);
    // }

    // calcul_flop_micro("sdot micro", 2 * VECSIZE, tdiff_micro(&start, &end));

    // printf("==========================================================\n");
    // printf("TEST SWAP COMPLEX DOUBLE\n");
    // printf("==========================================================\n");

    // init_flop_tsc();

    // val7.real = 1.0;
    // val7.imaginary = 2.0;
    // val8.real = 3.0;
    // val8.imaginary = 6.0;

    // vectorZ_init(vec7, val7);
    // vectorZ_init(vec8, val8);

    // TOP_MICRO(start);
    // mncblas_zswap(VECSIZE, vec7, 1, vec8, 1);
    // TOP_MICRO(end);

    // printf("SHOW Z VEC1\n");
    // for (size_t i = 0; i < VECSIZE; i++)
    // {
    //     printf("VEC 1:%ld : %f ; %fi\n", i, vec7[i].real, vec7[i].imaginary);
    // }

    // printf("SHOW Z VEC2\n");
    // for (size_t i = 0; i < VECSIZE; i++)
    // {
    //     printf("VEC 2:%ld : %f ; %fi\n", i, vec8[i].real, vec8[i].imaginary);
    // }

    // printf("==========================================================\n");
    // printf("TEST cdotu COMPLEX FLOAT\n");
    // printf("==========================================================\n");

    // init_flop_tsc();

    // val5.real = 1.0;
    // val5.imaginary = 2.0;
    // val6.real = 3.0;
    // val6.imaginary = 6.0;

    // void * res = malloc(sizeof(complexe_float_t));

    // vectorC_init(vec5, val5);
    // vectorC_init(vec6, val6);

    // TOP_MICRO(start);
    // mncblas_cdotu_sub(VECSIZE, vec5, 1, vec6, 1, res);
    // TOP_MICRO(end);

    // printf("SHOW Z VEC1\n");
    // for (size_t i = 0; i < VECSIZE; i++)
    // {
    //     printf("VEC 1:%ld : %f ; %fi\n", i, vec5[i].real, vec5[i].imaginary);
    // }

    // printf("SHOW Z VEC2\n");
    // for (size_t i = 0; i < VECSIZE; i++)
    // {
    //     printf("VEC 2:%ld : %f ; %fi\n", i, vec6[i].real, vec6[i].imaginary);
    // }
    // complexe_float_t *RES = res;

    // printf("Result : (%f, %f)\n", RES[0].real, RES[0].imaginary);

    // printf("==========================================================\n");
    // printf("TEST cdotc COMPLEX FLOAT\n");
    // printf("==========================================================\n");

    // init_flop_tsc();

    // val5.real = 1.0;
    // val5.imaginary = 2.0;
    // val6.real = 3.0;
    // val6.imaginary = 6.0;

    // void * res2 = malloc(sizeof(complexe_float_t));

    // vectorC_init(vec5, val5);
    // vectorC_init(vec6, val6);

    // TOP_MICRO(start);
    // mncblas_cdotc_sub(VECSIZE, vec5, 1, vec6, 1, res2);
    // TOP_MICRO(end);

    // printf("SHOW Z VEC1\n");
    // for (size_t i = 0; i < VECSIZE; i++)
    // {
    //     printf("VEC 1:%ld : %f ; %fi\n", i, vec5[i].real, vec5[i].imaginary);
    // }

    // printf("SHOW Z VEC2\n");
    // for (size_t i = 0; i < VECSIZE; i++)
    // {
    //     printf("VEC 2:%ld : %f ; %fi\n", i, vec6[i].real, vec6[i].imaginary);
    // }
    // complexe_float_t *RES2 = res2;

    // printf("Result : (%f, %f)\n", RES2[0].real, RES2[0].imaginary);














    // printf("==========================================================\n");
    // printf("TEST zdotu COMPLEX FLOAT\n");
    // printf("==========================================================\n");

    // init_flop_tsc();

    // val7.real = 1.0;
    // val7.imaginary = 2.0;
    // val8.real = 3.0;
    // val8.imaginary = 6.0;

    // void * res3 = malloc(sizeof(complexe_double_t));

    // vectorZ_init(vec7, val7);
    // vectorZ_init(vec8, val8);

    // TOP_MICRO(start);
    // mncblas_zdotu_sub(VECSIZE, vec7, 1, vec8, 1, res3);
    // TOP_MICRO(end);

    // printf("SHOW Z VEC1\n");
    // for (size_t i = 0; i < VECSIZE; i++)
    // {
    //     printf("VEC 1:%ld : %f ; %fi\n", i, vec7[i].real, vec7[i].imaginary);
    // }

    // printf("SHOW Z VEC2\n");
    // for (size_t i = 0; i < VECSIZE; i++)
    // {
    //     printf("VEC 2:%ld : %f ; %fi\n", i, vec8[i].real, vec8[i].imaginary);
    // }
    // complexe_double_t *RES3 = res3;

    // printf("Result : (%f, %f)\n", RES3[0].real, RES3[0].imaginary);

    // printf("==========================================================\n");
    // printf("TEST zdotc COMPLEX FLOAT\n");
    // printf("==========================================================\n");

    // init_flop_tsc();

    // val7.real = 1.0;
    // val7.imaginary = 2.0;
    // val8.real = 3.0;
    // val8.imaginary = 6.0;

    // void * res4 = malloc(sizeof(complexe_double_t));

    // vectorZ_init(vec7, val7);
    // vectorZ_init(vec8, val8);

    // TOP_MICRO(start);
    // mncblas_zdotc_sub(VECSIZE, vec7, 1, vec8, 1, res4);
    // TOP_MICRO(end);

    // printf("SHOW Z VEC1\n");
    // for (size_t i = 0; i < VECSIZE; i++)
    // {
    //     printf("VEC 1:%ld : %f ; %fi\n", i, vec7[i].real, vec7[i].imaginary);
    // }

    // printf("SHOW Z VEC2\n");
    // for (size_t i = 0; i < VECSIZE; i++)
    // {
    //     printf("VEC 2:%ld : %f ; %fi\n", i, vec8[i].real, vec8[i].imaginary);
    // }
    // complexe_double_t *RES4 = res4;

    // printf("Result : (%f, %f)\n", RES4[0].real, RES4[0].imaginary);

    printf("==========================================================\n");
    printf("TEST GEMV FLOAT\n");
    printf("==========================================================\n");

    vfloat vec9, vec10;
    mfloat mat4;

    val1 = 3.0;
    val2 = 4.0;


    matrixF_init(mat4, val1); // INIT A
    vectorF_init(vec9, val2); // INIT X
    vectorF_init(vec10, 1.0); // INIT Y

    TOP_MICRO(start);
    // AMPHA 2 et BETA 5
           for(size_t i = 0; i < NB_FOIS; i++)
                mncblas_sgemv(MNCblasRowMajor,MNCblasNoTrans,VECSIZE,VECSIZE,2,mat4, 1, vec9, 1,5,vec10,0); 
    TOP_MICRO(end);

    calcul_flop_micro("sgemv micro", NB_FOIS*(5*(VECSIZE^2) + 4*VECSIZE), tdiff_micro(&start, &end));


    // printf("==========================================================\n");
    // printf("TEST GEMM FLOAT\n");
    // printf("==========================================================\n");

    // mfloat mat1, mat2, mat3;

    // init_flop_tsc();
    // float valM3 = 1.0;
    // val1 = 2.0;
    // val2 = 3.0;
    // // vectorF_init(vec1, val1);
    // // vectorF_init(vec2, val2);
    // matrixF_init(mat1, val1);
    // matrixF_init(mat2, val2);
    // matrixF_init(mat3, valM3);

    // TOP_MICRO(start);
    // mncblas_sgemm(MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,VECSIZE,VECSIZE,VECSIZE,2,mat1, 1, mat2, 1,5,mat3,0); // AMPHA 2 et BETA 5
    // TOP_MICRO(end);

    // printf("SHOW F MAT3\n");
    // show_matrixF(mat3);
}
