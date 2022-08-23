////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Aug 23 09:29:52 EDT 2022
//
////////////////////////////////////////////////////////////////////

#include "Sext.h"
#include "util.h"

#include <cuda_runtime_api.h>

Sext::Sext (
    std::string const& Name,
    double const L,
    double const K2,
    int    const NKick,
    double const Dx,
    double const Dy,
    double const Tilt
    )
{
  cudaMallocManaged((void**) &fTM,    36*sizeof(double));
  cudaMallocManaged((void**) &fMa,    36*sizeof(double));
  cudaMallocManaged((void**) &fMb,    36*sizeof(double));

  fType = ElementType::Sext;

  fName = Name;
  fL = L;
  fK2 = K2;
  fNKick = NKick;
  fDx = Dx;
  fDy = Dy;
  fTilt = Tilt;

  this->Update();
}


Sext::~Sext ()
{
  cudaFree(fTM);
  cudaFree(fMa);
  cudaFree(fMb);
}


void Sext::Update ()
{
  this->SetSympass();
  return;
}



void Sext::SetSympass ()
{
  double const a =  0.675603595979828664;
  double const b = -0.175603595979828664;
  double const g =  1.351207191959657328;
  double const d = -1.702414383919314656;

  fdL = fL / fNKick;

  eye(fMa, 6);
  eye(fMb, 6);

  fMa[0*6+1] = a * fdL;
  fMa[2*6+3] = fMa[0*6+1];
  fMb[0*6+1] = b * fdL;
  fMb[2*6+3] = fMb[0*6+1];

  fK2Lg = g * fK2 * fdL;
  fK2Ld = d * fK2 * fdL;
  return;
}
