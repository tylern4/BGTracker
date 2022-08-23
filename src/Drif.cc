#include "Drif.h"

#include <cuda.h>
#include <cuda_runtime_api.h>


Drif::Drif (std::string const Name, double const L)
{
  // Allocate managed memory
  cudaMallocManaged((void**) &fTM,    36*sizeof(double));
  cudaMallocManaged((void**) &fTX,     9*sizeof(double));
  cudaMallocManaged((void**) &fTY,     9*sizeof(double));

  fType = ElementType::Drif;

  fName = Name;
  this->L(L);
}

Drif::~Drif ()
{
  cudaFree(fTM);
  cudaFree(fTX);
  cudaFree(fTY);
}


void Drif::Update ()
{
  this->TransMatrix();
  this->TwissMatrix();
  return;
}


void Drif::TransMatrix ()
{
  // Update TM
  for (int iy = 0; iy != 6; ++iy) {
    for (int ix = 0; ix != 6; ++ix) {
      ix == iy ? fTM[iy*6 + ix] = 1 : fTM[iy*6 + ix] = 0;
    }
  }
  fTM[0*6 + 1] = fL;
  fTM[2*6 + 3] = fL;
  return;
}


void Drif::TwissMatrix ()
{
  fTX[0*3+0] = fTM[0*6+0]*fTM[0*6+0];
  fTX[0*3+1] = -2*fTM[0*6+0]*fTM[0*6+1];
  fTX[0*3+2] = fTM[0*6+1]*fTM[0*6+1];

  fTX[1*3+0] = -fTM[1*6+0]*fTM[0*6+0];
  fTX[1*3+1] = 1+2*fTM[0*6+1]*fTM[1*6+0];
  fTX[1*3+2] = -fTM[0*6+1]*fTM[1*6+1];

  fTX[2*3+0] = fTM[1*6+0]*fTM[1*6+0];
  fTX[2*3+1] = -2*fTM[1*6+1]*fTM[1*6+0];
  fTX[2*3+2] = fTM[1*6+1]*fTM[1*6+1];


  fTY[0*3+0] = fTM[2*6+2]*fTM[2*6+2];
  fTY[0*3+1] = -2*fTM[2*6+2]*fTM[2*6+3];
  fTY[0*3+2] = fTM[2*6+3]*fTM[2*6+3];

  fTY[1*3+0] = -fTM[3*6+2]*fTM[2*6+2];
  fTY[1*3+1] = 1+2*fTM[2*6+3]*fTM[3*6+2];
  fTY[1*3+2] = -fTM[2*6+3]*fTM[3*6+3];

  fTY[2*3+0] = fTM[3*6+2]*fTM[3*6+2];
  fTY[2*3+1] = -2*fTM[3*6+3]*fTM[3*6+2];
  fTY[2*3+2] = fTM[3*6+3]*fTM[3*6+3];

  return;
}

void Drif::PrintTM ()
{
  for (int iy = 0; iy != 6; ++iy) {
    for (int ix = 0; ix != 6; ++ix) {
      printf("%1.3e ", fTM[iy * 6 + ix]);
    }
    std::cout << std::endl;
  }
  return;
}


void Drif::PrintTwiss ()
{
  for (int iy = 0; iy != 3; ++iy) {
    for (int ix = 0; ix != 3; ++ix) {
      printf("%+1.3e ", fTX[iy * 3 + ix]);
    }
    std::cout << std::endl;
  }
  for (int iy = 0; iy != 3; ++iy) {
    for (int ix = 0; ix != 3; ++ix) {
      printf("%+1.3e ", fTY[iy * 3 + ix]);
    }
    std::cout << std::endl;
  }
  return;
}



std::string const& Drif::Name () const
{
  return fName;
}


double Drif::L () const {
  return fL;
}

void Drif::L (double const L)
{
  fL = L;
  if (L < 0) {
    std::cerr << "WARNING: " << fName << " has a negative length" << fL << std::endl;
  }
  this->Update();
  return;
}


