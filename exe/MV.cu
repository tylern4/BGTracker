////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed Aug 17 14:29:41 EDT 2022
//
////////////////////////////////////////////////////////////////////


#include <iostream>


__global__ void mdotv (double *m, double *v)
{
  int IV = blockDim.x * blockIdx.x / 6 + threadIdx.x / 6;
  int IR = threadIdx.x % 6;

  double value = 0;

  for (int icol = 0; icol != 6; ++icol) {
    value += m[IR * 6 + icol] * v[IV * 6 + icol];
  }

  __syncthreads();
  v[IV * 6 + IR] = value;
}


void printp (double* p)
{
  for (int ix = 0; ix != 6; ++ix) {
    printf("%1.6e ", p[ix]);
  }
  std::cout << std::endl;
  return;
}

void printm (double* m)
{
  for (int iy = 0; iy != 6; ++iy) {
    for (int ix = 0; ix != 6; ++ix) {
      printf("%1.6e ", m[iy*6 + ix]);
    }
    std::cout << std::endl;
  }
  return;
}

int main (int argc, char* argv[])
{

  double* p;
  double* m;
  cudaMallocManaged(&p, 6*sizeof(double));
  cudaMallocManaged(&m, 6*6*sizeof(double));

  p[0] = 1;
  p[1] = 0;
  p[2] = 1;
  p[3] = 1;
  p[4] = 1;
  p[5] = 0.1;

  for (int i = 0; i != 6; ++i) {
    m[i*6 +i] = 1;
  }
  m[0*6 + 1] = 1.234;
  m[2*6 + 3] = 1.234;
  cudaDeviceSynchronize();

  printp(p);
  std::cout << std::endl;
  printm(m);

  mdotv<<<1, 6>>>(m, p);
  cudaDeviceSynchronize();
  std::cout << std::endl;
  std::cout << std::endl;

  printp(p);
  std::cout << std::endl;
  printm(m);


  return 0;
}
