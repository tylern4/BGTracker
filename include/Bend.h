#ifndef __GUARD__Bend_h
#define __GUARD__Bend_h

#include "BGE.h"

#include <iostream>
#include <string>

class Bend : public BGE
{
  public:
    Bend (
        std::string const& Name = "B01",
        double const L = 1,
        double const Angle = 1e-9,
        double const E1 = 0,
        double const E2 = 0,
        double const K1 = 0,
        double const K2 = 0,
        double const NKick = 10,
        double const HGap = 0,
        double const FInt = 0.5,
        double const Dx = 0,
        double const Dy = 0,
        double const Tilt = 0
        );
    ~Bend ();

    void Update ();
    void SetSympass ();
    void TransMatrix ();

    std::string fName;
    double fL;
    double fAngle;
    double fE1;
    double fE2;
    double fK1;
    double fK2;
    double fNKick;
    double fDx;
    double fDy;
    double fTilt;
    double fHGap;
    double fFInt;

    double fR;
    double fdL;
    double fLg;
    double fLd;
    double fK1Lg;
    double fK1Ld;
    double fK2Lg;
    double fK2Ld;
    double *fTM;
    double *fMa;
    double *fMb;
    double *fM1;
    double *fM2;

    // Tag ?

};




#endif
