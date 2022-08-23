#ifndef __GUARD_Drif_h__
#define __GUARD_Drif_h__

#include "BGE.h"

#include <iostream>
#include <string>

class Drif : public BGE {
  public:
    Drif (std::string const Name = "D01", double const L = 0);
    ~Drif ();

    void Update ();
    void TransMatrix ();
    void TwissMatrix ();
    void PrintTM ();
    void PrintTwiss ();

    std::string const& Name () const;
    double L () const;
    void   L (double const);

    double *fTM;
    double *fTX;
    double *fTY;

    std::string fName;
    double  fL;
    int     fNKick;
    double  fDx;
    double  fDy;
    double  fTilt;
    // tag ?
};



inline std::ostream& operator << (std::ostream& os, Drif const& o)
{
  os << o.Name() << ": drif, L=" << o.L() << std::endl;
  return os;
}

#endif
