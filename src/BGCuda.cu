#include "BGCuda.h"
#include "cuda_runtime.h"

__global__ void mysum (int N, double *m, double *y)
{
  for (int i = 0; i != N; ++i) {
    *y += m[i];
  }
}



__global__ void drif (double *m, double *x)
{
  int IParticle = blockDim.x * blockIdx.x / 6 + threadIdx.x / 6;
  int IRow = threadIdx.x % 6;

  double value = 0;

  for (int i = 0; i != 6; ++i) {
    value += m[IRow * 6 + i] * x[IParticle * 6 + i];
  }

  __syncthreads();
  x[IParticle * 6 + IRow] = value;

  __syncthreads();
  if (IRow == 4) {
    x[IParticle * 6 + IRow] = value + (sqrt(1 + x[IParticle*6+1]*x[IParticle*6+1] + x[IParticle*6+3]*x[IParticle*6+3]) - 1) * m[1];
  }
  return;
}




__global__ void xeqmdotx (double *m, double *x)
{
  int IParticle = blockDim.x * blockIdx.x / 6 + threadIdx.x / 6;
  int IRow = threadIdx.x % 6;
  
  double value = 0;
  for (int i = 0; i != 6; ++i) {
    value += m[IRow * 6 + i] * x[IParticle * 6 + i];
  }
  __syncthreads();
  x[IParticle * 6 + IRow] = value;

  return;
}

__device__ void dXeqMdotX (double *m, double *x)
{
  double value[6] = { 0 };
  for (int iy = 0; iy != 6; ++iy) {
    for (int i = 0; i != 6; ++i) {
      value[iy] += m[iy * 6 + i] * x[i];
    }
  }
  for (int i = 0; i != 6; ++i) {
    x[i] = value[i];
  }

  return;
}


__global__ void xupdate13 (double *x, double *k)
{
  int IParticle = blockDim.x * blockIdx.x / 2 + threadIdx.x / 2;

  int i = (threadIdx.x % 2) * 2 + 1;
  double sign = i - 2;
  x[IParticle * 6 + i] += sign * *k * x[IParticle * 6 + i - 1] / (1 + x[IParticle * 6 + 5]);
}

__global__ void addtos (double *s, double *x)
{
  int IParticle = blockDim.x * blockIdx.x + threadIdx.x;

  double x2p = x[IParticle*6 + 1];
  double y2p = x[IParticle*6 + 3];
  //double xp = (x1p + x2p);
  *s += 0;
}




__global__ void dsympass4 (double *X, double *TM, double L, int N)
{
  int IParticle = blockDim.x * blockIdx.x + threadIdx.x;

  if (IParticle >= N) {
    return;
  }

  double *x = &(X[IParticle * 6]);

  dXeqMdotX(TM, x);
  x[4] += (sqrt(1 + x[1]*x[1] + x[3]*x[3]) - 1) * L;

  return;
}







__global__ void qsympass4 (double *X, double *Ma, double *Mb, double K1Lg, double K1Ld, double dL, int NKick, int N)
{
  int IParticle = blockDim.x * blockIdx.x + threadIdx.x;

  if (IParticle >= N) {
    return;
  }

  double *x = &(X[IParticle * 6]);

  double x1p;
  double y1p;
  double xp;
  double yp;

  double S = 0;
  for (int i = 0; i != NKick; ++i) {
    x1p = x[1];
    y1p = x[3];

    dXeqMdotX(Ma, x);
    x[1] -= K1Lg * x[0] / (1 + x[5]);
    x[3] += K1Lg * x[2] / (1 + x[5]);
    dXeqMdotX(Mb, x);
    x[1] -= K1Ld * x[0] / (1 + x[5]);
    x[3] += K1Ld * x[2] / (1 + x[5]);
    dXeqMdotX(Mb, x);
    x[1] -= K1Lg * x[0] / (1 + x[5]);
    x[3] += K1Lg * x[2] / (1 + x[5]);
    dXeqMdotX(Ma, x);

    xp = (x1p + x[1])/2;
    yp = (y1p + x[3])/2;

    S += sqrt(1 + (xp*xp) + (yp*yp)) * (dL);

  }

  x[4] += S - NKick * (dL);

  return;
}




__global__ void bsympass4 (double *X, double *Ma, double *Mb, double *M1, double *M2, double K1Ld, double K1Lg, double K2Ld, double K2Lg, double Lg, double Ld,  double dL, double R, int NKick, int N)
{
  int IParticle = blockDim.x * blockIdx.x + threadIdx.x;

  if (IParticle >= N) {
    return;
  }

  double *x = &(X[IParticle * 6]);

  double x1p;
  double y1p;
  double x1;
  double xav;
  double xp;
  double yp;

  // skipped tilt, dx dy
  

  // M1.  include this as if??
  dXeqMdotX(M1, x);

  double S = 0;
  for (int i = 0; i != NKick; ++i) {
    x1p = x[1];
    y1p = x[3];
    x1 = x[0];

    dXeqMdotX(Ma, x);
    x[1] -= K2Lg / 2 * (x[0]*x[0] - x[2]*x[2])/(1 + x[5]) + K1Lg * x[0] / (1 + x[5]) - Lg * x[5] / R + Lg * x[0] / (R * R);
    x[3] += K2Lg * x[0] * x[2] / (1 + x[5]) + K1Lg * x[2] / (1 + x[5]);
    dXeqMdotX(Mb, x);
    x[1] -= K2Ld / 2 * (x[0]*x[0] - x[2]*x[2])/(1 + x[5]) + K1Ld * x[0] / (1 + x[5]) - Ld * x[5] / R + Ld * x[0] / (R * R);
    x[3] += K2Ld * x[0] * x[2] / (1 + x[5]) + K1Ld * x[2] / (1 + x[5]);
    dXeqMdotX(Mb, x);
    x[1] -= K2Lg / 2 * (x[0]*x[0] - x[2]*x[2])/(1 + x[5]) + K1Lg * x[0] / (1 + x[5]) - Lg * x[5] / R + Lg * x[0] / (R * R);
    x[3] += K2Lg * x[0] * x[2] / (1 + x[5]) + K1Lg * x[2] / (1 + x[5]);
    dXeqMdotX(Ma, x);

    xp = (x1p + x[1]) / 2;
    yp = (y1p + x[3]) / 2;
    xav = (x1 + x[0]) / 2;

    S += sqrt(1 + xp*xp + yp*yp) * (1 + xav/R) * dL;
  }

  dXeqMdotX(M2, x);

  // skipping for now tilt dx dy

  x[4] += S - dL * NKick;

  return;
}







__global__ void ssympass4 (double *X, double *Ma, double *Mb, double K2Ld, double K2Lg, double dL, int NKick, int N)
{
  int IParticle = blockDim.x * blockIdx.x + threadIdx.x;

  if (IParticle >= N) {
    return;
  }

  double *x = &(X[IParticle * 6]);

  double x1p;
  double y1p;
  double x1;
  double xav;
  double xp;
  double yp;

  // skipped tilt, dx dy

  double S = 0;
  for (int i = 0; i != NKick; ++i) {
    x1p = x[1];
    y1p = x[3];

    dXeqMdotX(Ma, x);
    x[1] -= K2Lg / 2 * (x[0]*x[0] - x[2]*x[2])/(1 + x[5]);
    x[3] += K2Lg * x[0] * x[2] / (1 + x[5]);
    dXeqMdotX(Mb, x);
    x[1] -= K2Ld / 2 * (x[0]*x[0] - x[2]*x[2])/(1 + x[5]);
    x[3] += K2Ld * x[0] * x[2] / (1 + x[5]);
    dXeqMdotX(Mb, x);
    x[1] -= K2Lg / 2 * (x[0]*x[0] - x[2]*x[2])/(1 + x[5]);
    x[3] += K2Lg * x[0] * x[2] / (1 + x[5]);
    dXeqMdotX(Ma, x);

    xp = (x1p + x[1]) / 2;
    yp = (y1p + x[3]) / 2;
    xav = (x1 + x[0]) / 2;

    S += sqrt(1 + xp*xp + yp*yp) * dL;
  }


  // skipping for now tilt dx dy

  x[4] += S - dL * NKick;

  return;
}
