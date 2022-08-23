#ifndef __GUARD__Sext_h
#define __GUARD__Sext_h

#include "BGE.h"

#include <string>

class Sext : public BGE
{
  public:
    Sext (
        std::string const& Name = "SEXT01",
        double const L = 0.25,
        double const K2 = 1,
        int    const NKick = 4,
        double const Dx = 0,
        double const Dy = 0,
        double const Tilt = 0
        );
    ~Sext ();


    void Update ();
    void SetSympass ();
    

    std::string fName;
    double fL;
    double fK2;
    int    fNKick;
    double fDx;
    double fDy;
    double fTilt;
    // tag?

    double fdL;
    double fK2Lg;
    double fK2Ld;
    double *fTM;
    double *fMa;
    double *fMb;


};







#endif
