////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Mon Aug 22 14:36:02 EDT 2022
//
////////////////////////////////////////////////////////////////////


#include "Bend.h"
#include "util.h"

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cmath>


Bend::Bend (
    std::string const& Name,
    double const L,
    double const Angle,
    double const E1,
    double const E2,
    double const K1,
    double const K2,
    double const NKick,
    double const HGap,
    double const FInt,
    double const Dx,
    double const Dy,
    double const Tilt
    )
{
  // Allocate managed memory
  cudaMallocManaged((void**) &fTM,    36*sizeof(double));
  cudaMallocManaged((void**) &fMa,    36*sizeof(double));
  cudaMallocManaged((void**) &fMb,    36*sizeof(double));
  cudaMallocManaged((void**) &fM1,    36*sizeof(double));
  cudaMallocManaged((void**) &fM2,    36*sizeof(double));

  fType = ElementType::Bend;

  fL      = L;
  fAngle  = Angle;
  fE1     = E1;
  fE2     = E2;
  fK1     = K1;
  fK2     = K2;
  fNKick  = NKick;
  fHGap   = HGap;
  fFInt   = FInt;
  fDx     = Dx;
  fDy     = Dy;
  fTilt   = Tilt;

  eye(fM1);
  eye(fM2);
  this->Update();

}


Bend::~Bend ()
{
  cudaFree(fTM);
  cudaFree(fMa);
  cudaFree(fMb);
  cudaFree(fM1);
  cudaFree(fM2);
}



void Bend::Update ()
{
  fR = fL / fAngle;

  this->TransMatrix();
  this->SetSympass();
  return;
}



void Bend::SetSympass ()
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

  fLg = g * fdL;
  fLd = d * fdL;
  fK1Lg = g * fK1 * fdL;
  fK1Ld = d * fK1 * fdL;
  fK2Lg = g * fK2 * fdL;
  fK2Ld = d * fK2 * fdL;
}



void Bend::TransMatrix ()
{
  double kx = fK1 + 1 / (fR * fR);

  eye(fTM, 6);

  double k;
  double p;
  if (kx > 0) {
    k = sqrt(kx);
    p = k * fL;

    fTM[0*6 + 0] = cos(p);
    fTM[0*6 + 1] = sin(p) / k;
    fTM[1*6 + 0] = -k * sin(p);
    fTM[1*6 + 1] = cos(p);

    fTM[0*6 + 5] = (1 - cos(p)) / (fR * kx);
    fTM[1*6 + 5] = sin(p) / (fR * k);

    fTM[4*6 + 0] = fTM[1*6 + 5];
    fTM[4*6 + 1] = fTM[0*6 + 5];
    fTM[4*6 + 5] = (k * fL - sin(p)) / fR / fR / k / kx;
  } else if (kx < 0) {
    k = sqrt(-kx);
    p = k * fL;

    fTM[0*6 + 0] = cosh(p);
    fTM[0*6 + 1] = sinh(p) / k;
    fTM[1*6 + 0] = k * sinh(p);
    fTM[1*6 + 1] = cosh(p);

    fTM[0*6 + 5] = (cos(p) - 1) / (-fR * kx);
    fTM[1*6 + 5] = sinh(p) / (fR * k);

    fTM[4*6 + 0] = fTM[1*6 + 5];
    fTM[4*6 + 1] = fTM[0*6 + 5];
    fTM[4*6 + 5] = (k * fL - sinh(p)) / fR / fR / k / kx;
  } else {
    fTM[0*6 + 1] = fL;
  }

  if (fK1 > 0) {
    k = sqrt(fK1);
    p = k * fL;

    fTM[2*6 + 2] = cosh(p);
    fTM[2*6 + 3] = sinh(p) / k;
    fTM[3*6 + 2] = k * sinh(p);
    fTM[3*6 + 3] = cosh(p);
  } else if (fK1 < 0) {
    k = sqrt(-fK1);
    p = k * fL;

    fTM[2*6 + 2] = cosh(p);
    fTM[2*6 + 3] = sin(p) / k;
    fTM[3*6 + 2] = -k * sin(p);
    fTM[3*6 + 3] = cos(p);
  } else {
    fTM[2*6 + 3] = fL;
  }

  // Entrance
  double vf = 0;
  if (fHGap != 0 && fFInt != 0) {
    vf = -1 / fR * fHGap * fFInt * (1 + pow(sin(fE1), 2)) / cos(fE1) * 2;
  }

  if (fE1 != 0 || vf != 0) {
    eye(fM1, 6);
    fM1[1*6 + 0] = tan(fE1) / fR;
    fM1[3*6 + 2] = -tan(fE1 + vf) / fR;
    MeqMdotN(fTM, fM1);
  }

  // Exit
  if (fHGap != 0 || fFInt != 0) {
    vf = -1 / fR * fHGap * fFInt * (1 + pow(sin(fE2), 2)) / cos(fE2) * 2;
  } else {
    vf = 0;
  }

  if (fE1 != 0 || vf != 0) {
    eye(fM2, 6);
    fM2[1*6 + 0] = tan(fE2) / fR;
    fM2[3*6 + 2] = -tan(fE2 + vf) / fR;
    MeqNdotM(fTM, fM2);
  }

  if (fTilt != 0) {
    double r1[36];
    double r0[36];
    
    rotmat(r1, -fTilt);
    rotmat(r0, fTilt);
    MeqMdotN(r1, fTM);
    LeqMdotN(fTM, r1, r0);

  }


  return;
}
