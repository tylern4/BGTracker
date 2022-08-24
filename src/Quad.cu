////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed Aug 17 17:51:09 EDT 2022
//
////////////////////////////////////////////////////////////////////

#include "Quad.h"
#include "util.h"

#include <cuda.h>
#include <cuda_runtime_api.h>
#include "BGCuda.h"

Quad::Quad(
    std::string const Name,
    double const L,
    double const K1,
    int const NKick,
    double const Dx,
    double const Dy,
    double const Tilt)
{
  // Allocate managed memory
  cudaMallocManaged((void **)&fTM, 36 * sizeof(double));
  cudaMallocManaged((void **)&fMa, 36 * sizeof(double));
  cudaMallocManaged((void **)&fMb, 36 * sizeof(double));

  fType = ElementType::Quad;

  fName = Name;
  fK1 = K1;
  fNKick = NKick;
  fDx = Dx;
  fDy = Dy;
  fTilt = Tilt;
  this->L(L);
}

Quad::~Quad()
{
  cudaFree(fTM);
  cudaFree(fMa);
  cudaFree(fMb);
}

void Quad::Update()
{
  this->TransMatrix();
  // this->TwissMatrix();
  this->SetSympass();
  return;
}

void Quad::run_qsympass4(int NParticles, int NTurns)
{

  double *particles;
  cudaMallocManaged((void **)&particles, NParticles * 6 * sizeof(double));

  // Can randomize this later..
  for (int i = 0; i != NParticles; ++i)
  {
    particles[i * 6 + 0] = 1e-8;
    particles[i * 6 + 1] = 0;
    particles[i * 6 + 2] = 0;
    particles[i * 6 + 3] = 0;
    particles[i * 6 + 4] = 0;
    particles[i * 6 + 5] = 0;
    // for (int j = 0; j != 6; ++j) {
    //   particles[i*6 + j] = (double)std::rand()/(float) RAND_MAX;
    // }
  }

  int const NThreadsPerBlock = 32;
  int const NBlocks = NParticles % NThreadsPerBlock == 0 ? ((int)NParticles / NThreadsPerBlock) : ((int)NParticles / NThreadsPerBlock) + 1;
  std::cout << "NBlocks: " << NBlocks << std::endl;
  std::cout << NParticles * 6 * sizeof(double) / 1024 / 1024 << " MB particle data" << std::endl;
  for (int iturn = 0; iturn != NTurns; ++iturn) {
    qsympass4<<<NBlocks, NThreadsPerBlock>>>(particles, this->fMa, this->fMb, this->fK1Lg, this->fK1Ld, this->fdL, this->fNKick, NParticles);
  }

}

void Quad::TransMatrix()
{
  // Update TM
  for (int iy = 0; iy != 6; ++iy)
  {
    for (int ix = 0; ix != 6; ++ix)
    {
      ix == iy ? fTM[iy * 6 + ix] = 1 : fTM[iy * 6 + ix] = 0;
    }
  }

  double const k = sqrt(abs(fK1));
  double const p = k * fL;
  if (fK1 > 0)
  {
    fTM[0 * 6 + 0] = cos(p);
    fTM[0 * 6 + 1] = sin(p) / k;
    fTM[1 * 6 + 0] = -k * sin(p);
    fTM[1 * 6 + 1] = cos(p);

    fTM[2 * 6 + 2] = cosh(p);
    fTM[2 * 6 + 3] = sinh(p) / k;
    fTM[3 * 6 + 2] = +k * sinh(p);
    fTM[3 * 6 + 3] = cosh(p);
  }
  else if (fK1 < 0)
  {
    fTM[0 * 6 + 0] = cosh(p);
    fTM[0 * 6 + 1] = sinh(p) / k;
    fTM[1 * 6 + 0] = +k * sinh(p);
    fTM[1 * 6 + 1] = cosh(p);

    fTM[2 * 6 + 2] = cos(p);
    fTM[2 * 6 + 3] = sin(p) / k;
    fTM[3 * 6 + 2] = -k * sin(p);
    fTM[3 * 6 + 3] = cos(p);
  }
  else
  {
    // super(quad,self)._transmatrix()
  }

  if (fTilt != 0)
  {
    rotmat(fTM, fTilt);
  }

  return;
}

void Quad::SetSympass()
{
  double const a = 0.675603595979828664;
  double const b = -0.175603595979828664;
  double const g = 1.351207191959657328;
  double const d = -1.702414383919314656;

  fdL = fL / fNKick;

  for (int iy = 0; iy != 6; ++iy)
  {
    for (int ix = 0; ix != 6; ++ix)
    {
      if (ix == iy)
      {
        fMa[iy * 6 + ix] = 1;
        fMb[iy * 6 + ix] = 1;
      }
      else
      {
        fMa[iy * 6 + ix] = 0;
        fMb[iy * 6 + ix] = 0;
      }
    }
  }
  fMa[0 * 6 + 1] = a * fdL;
  fMa[2 * 6 + 3] = fMa[0 * 6 + 1];
  fMb[0 * 6 + 1] = b * fdL;
  fMb[2 * 6 + 3] = fMb[0 * 6 + 1];

  fK1Lg = g * fK1 * fdL;
  fK1Ld = d * fK1 * fdL;
  return;
}

void Quad::DoDxDyTilt(double *particles, int const N)
{
  double r[36];
  if (fTilt != 0)
  {
    rotmat(r, fTilt);
  }
  for (int i = 0; i != N; ++i)
  {
    if (fDx != 0)
      particles[i * 6 + 0] -= fDx;
    if (fDy != 0)
      particles[i * 6 + 2] -= fDy;
    if (fTilt != 0)
    {
      VeqMdotV(r, &particles[i * 6]);
    }
  }
  return;
}

void Quad::UnDoDxDyTilt(double *particles, int const N)
{
  double r[36];
  if (fTilt != 0)
  {
    rotmat(r, -fTilt);
  }
  for (int i = 0; i != N; ++i)
  {
    if (fTilt != 0)
    {
      VeqMdotV(r, &particles[i * 6]);
    }
    if (fDy != 0)
      particles[i * 6 + 2] -= fDy;
    if (fDx != 0)
      particles[i * 6 + 0] -= fDx;
  }
  return;
}

void Quad::PrintTM()
{
  for (int iy = 0; iy != 6; ++iy)
  {
    for (int ix = 0; ix != 6; ++ix)
    {
      printf("%+1.3e ", fTM[iy * 6 + ix]);
    }
    std::cout << std::endl;
  }
  return;
}

std::string const &Quad::Name() const
{
  return fName;
}

double Quad::L() const
{
  return fL;
}

void Quad::L(double const L)
{
  fL = L;
  if (L < 0)
  {
    std::cerr << "WARNING: " << fName << " has a negative length" << fL << std::endl;
  }
  this->Update();
  return;
}

double Quad::K1() const
{
  return fK1;
}

void Quad::K1(double const K1)
{
  fK1 = K1;
  this->Update();
  return;
}

double Quad::Tilt() const
{
  return fTilt;
}

void Quad::Tilt(double const Tilt)
{
  fTilt = Tilt;
  this->Update();
  return;
}
