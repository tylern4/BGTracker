#ifndef GUARD_BGCuda_h
#define GUARD_BGCuda_h


#include <cuda.h>
#include <cuda_runtime_api.h>

__global__ void mysum (int N, double *m, double *y);
__global__ void drif (double *m, double *x);

// Drif
__global__ void xeqmdotx (double *m, double *x);
__global__ void xupdate13 (double *x, double *k);


// Drif
__global__ void dsympass4 (double *X, double *TM, double L, int N);

// Quad
__global__ void qsympass4 (double *x, double *Ma, double *Mb, double K1Lg, double K1Ld, double dL, int NKick, int N);
__device__ void dxeqmdotx (double *m, double *x);

// Bend
__global__ void bsympass4 (double *X, double *Ma, double *Mb, double *M1, double *M2, double K1Ld, double K1Lg, double K2Ld, double K2Lg, double Lg, double Ld,  double dL, double R, int NKick, int N);

// Sext
__global__ void ssympass4 (double *X, double *Ma, double *Mb, double K2Ld, double K2Lg, double dL, int NKick, int N);

#endif
