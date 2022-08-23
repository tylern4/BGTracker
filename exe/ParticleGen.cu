////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Mon Aug 15 14:27:55 EDT 2022
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "Drif.h"
#include "Quad.h"
#include "Bend.h"
#include "Sext.h"
#include "BGCuda.h"

int ParticleGen (int const N)
{
  int const NParticles = N;
  int const NThreadsPerBlock = 192;
  int const NBlocks = N % NThreadsPerBlock == 0 ? ((int) N / NThreadsPerBlock) : ((int) N / NThreadsPerBlock) + 1;

  std::cout << "Blocks: " << NBlocks << std::endl;

  std::srand(1234987);
  std::cout << NParticles*6*sizeof(double) / 1024/1024 << " MB"  << std::endl;


  double* particles;
  cudaMallocManaged(&particles, NParticles*6*sizeof(double));
  double* S;
  cudaMallocManaged(&S, NParticles*sizeof(double));
  double* XPYP;
  cudaMallocManaged(&XPYP, NParticles*4*sizeof(double));

  for (int i = 0; i != NParticles; ++i) {
      particles[i*6 + 0] = 1;
      particles[i*6 + 1] = 0;
      particles[i*6 + 2] = 1;
      particles[i*6 + 3] = 1;
      particles[i*6 + 4] = 1;
      particles[i*6 + 5] = 0.1;
    //for (int j = 0; j != 6; ++j) {
    //  particles[i*6 + j] = (double)std::rand()/(float) RAND_MAX;
    //}
  }


  Drif d1("D1", 1.234);
  for (int iy = 0; iy != 6; ++iy) {
    for (int ix = 0; ix != 6; ++ix) {
      printf("%+1.3e ", d1.fTM[iy * 6 + ix]);
    }
    std::cout << std::endl;
  }
  cudaDeviceSynchronize();
  dsympass4<<<NBlocks, NThreadsPerBlock>>>(particles, d1.fTM, d1.fL, NParticles);
  cudaDeviceSynchronize();
  if (true)
  for (int i = 0; i != NParticles; ++i) {
    for (int j = 0; j != 6; ++j) {
      printf("%.6e ", particles[i*6 + j]);
    }
    std::cout << std::endl;
  }
  exit(0);

  Sext s1("sh1", 0.2, 20.0, 4, 0, 0, 0);
  for (int iy = 0; iy != 6; ++iy) {
    for (int ix = 0; ix != 6; ++ix) {
      printf("%+1.3e ", s1.fMb[iy * 6 + ix]);
    }
    std::cout << std::endl;
  }

  cudaDeviceSynchronize();
  ssympass4<<<NBlocks, NThreadsPerBlock>>> (particles, s1.fMa, s1.fMb, s1.fK2Ld, s1.fK2Lg, s1.fdL, s1.fNKick, NParticles);
  cudaDeviceSynchronize();
  if (true)
  for (int i = 0; i != NParticles; ++i) {
    for (int j = 0; j != 6; ++j) {
      printf("%.6e ", particles[i*6 + j]);
    }
    std::cout << std::endl;
  }
  exit(0);




  Bend b1("B01", 1, 0.02, 0.1, 0.1, 0.5, 0.0, 20, 0, 0.5, 0, 0, 0);
  for (int iy = 0; iy != 6; ++iy) {
    for (int ix = 0; ix != 6; ++ix) {
      printf("%+1.3e ", b1.fTM[iy * 6 + ix]);
    }
    std::cout << std::endl;
  }

  cudaDeviceSynchronize();
  bsympass4<<<NBlocks, NThreadsPerBlock>>>(particles, b1.fMa, b1.fMb, b1.fM1, b1.fM2, b1.fK1Ld, b1.fK1Lg, b1.fK2Ld, b1.fK2Lg, b1.fLg, b1.fLd,  b1.fdL, b1.fR, b1.fNKick, NParticles);

  cudaDeviceSynchronize();
  if (true)
  for (int i = 0; i != NParticles; ++i) {
    for (int j = 0; j != 6; ++j) {
      printf("%.6e ", particles[i*6 + j]);
    }
    std::cout << std::endl;
  }
  exit(0);


  Quad q1("q1",0.1,0.5,4,0,0, 0.0);
  std::cout << q1 << std::endl;
  //q1.PrintTM();
  for (int iy = 0; iy != 6; ++iy) {
    for (int ix = 0; ix != 6; ++ix) {
      printf("%+1.3e ", q1.fMb[iy * 6 + ix]);
    }
    std::cout << std::endl;
  }

  std::cout << (q1.fNKick) << std::endl;

  cudaDeviceSynchronize();
  qsympass4<<<NBlocks, NThreadsPerBlock>>>(particles, q1.fMa, q1.fMb, q1.fK1Lg, q1.fK1Ld, q1.fdL, q1.fNKick, NParticles);
  cudaDeviceSynchronize();

  std::cout << "hihihihi" << std::endl;
  if (true)
  for (int i = 0; i != NParticles; ++i) {
    for (int j = 0; j != 6; ++j) {
      printf("%.6e ", particles[i*6 + j]);
    }
    std::cout << std::endl;
  }


  exit(0);










  //q1.DoDxDyTilt();
//  for (int i = 0; i != NParticles; ++i) {
//    S[i] = 0;
//  }
//  for (int i = 0; i != q1.fNKick; ++i) {
//    xeqmdotx<<<10, 512>>>(q1.fMa, particles);
//    xupdate13<<<10, 512>>>(particles, q1.fK1Lg);
//    xeqmdotx<<<10, 512>>>(q1.fMb, particles);
//    xupdate13<<<10, 512>>>(particles, q1.fK1Ld);
//    xeqmdotx<<<10, 512>>>(q1.fMb, particles);
//    xupdate13<<<10, 512>>>(particles, q1.fK1Lg);
//    xeqmdotx<<<10, 512>>>(q1.fMa, particles);
//  }
  //q1.UnDoDxDyTilt();

  exit(0);

  Drif d0001("D0001", 1.234);
  d0001.PrintTM();
  d0001.PrintTwiss();

  double count = 0;

  if (false)
  for (int i = 0; i != NParticles; ++i) {
    for (int j = 0; j != 6; ++j) {
      printf("%.6e ", particles[i*6 + j]);
    }
    std::cout << std::endl;
  }

  cudaDeviceSynchronize();
  std::cout << "now calc" << std::endl;
  for (int i = 0; i != 1000; ++i) {
    drif<<<NBlocks, NThreadsPerBlock>>>(d0001.fTM, particles);
  }

  cudaDeviceSynchronize();
  if (false)
  for (int i = 0; i != NParticles; ++i) {
    for (int j = 0; j != 6; ++j) {
      printf("%.6e ", particles[i*6 + j]);
    }
    std::cout << std::endl;
  }

  //delete particles;
  cudaFree(particles);
  cudaFree(S);

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " [N]" << std::endl;
    return 1;
  }

  int const N = atoi(argv[1]);
  ParticleGen(N);

  return 0;
}
