#ifndef __GUARD__BGE_h
#define __GUARD__BGE_h


class BGE
{
  public:
    BGE () {}
    virtual ~BGE () {}


    enum ElementType {
      NONE,
      Drif,
      Quad,
      Bend,
      Sext
    };

    ElementType Type () const
    {
      return fType;
    }

  protected:
    ElementType fType;


};









#endif
