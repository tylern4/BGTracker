#ifndef __GUARD_Quad_h__
#define __GUARD_Quad_h__

#include "BGE.h"

#include <iostream>
#include <string>
#include <cmath>

void rotmat (double* m, double const angle=M_PI_4);
class Quad : public BGE
{
  public:
    Quad (
        std::string const Name = "Q01",
        double      const L = 0.25,
        double      const K1 = 1,
        int         const NKick = 4,
        double      const Dx = 0,
        double      const Dy = 0,
        double      const Tilt = 0);
    ~Quad ();

    void PrintTM ();
    void Update ();
    void TransMatrix ();
    void TwissMatrix ();
    void SetSympass ();

    void DoDxDyTilt (double*, int const);
    void UnDoDxDyTilt (double*, int const);
    
    std::string const& Name () const;
    double L () const;
    void   L (double const);
    double K1 () const;
    void   K1 (double const);
    double Tilt () const;
    void   Tilt (double const);

    std::string fName;
    double      fL;
    double      fK1;
    double      fNKick;
    double      fDx;
    double      fDy;
    double      fTilt;

    double      fdL;
    double      fK1Ld;
    double      fK1Lg;

    double *fTM;
    double *fMa;
    double *fMb;
    // tag?

};


inline std::ostream& operator << (std::ostream& os, Quad const& o)
{
  os << o.Name() << ": quad, L=" << o.L() << ",K1=" << o.K1() << ",tilt=" << o.Tilt() << std::endl;
  return os;
}

#endif
