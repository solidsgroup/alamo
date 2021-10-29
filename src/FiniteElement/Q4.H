#ifndef FINITEELEMENT_Q4_
#define FINITEELEMENT_Q4_

#include "Set/Set.H"

namespace FiniteElement
{
class Q4
{
public:
  Q4() {} ;
  Q4(amrex::Real _dx, amrex::Real _dy): dx(_dx), dy(_dy) {} ;
  void Define(Set::Scalar a_dx, Set::Scalar a_dy)
  {
    dx = a_dx; dy = a_dy;
  }

  AMREX_FORCE_INLINE
  Set::Scalar Phi(int m, const Set::Vector x) const
  {
    if (m==0) return (dx-x(0))*(dy-x(1))/dx/dy;
    if (m==1) return (x(0)-0.0)*(dy-x(1))/dx/dy;
    if (m==2) return (dx-x(0))*(x(1)-0.0)/dx/dy;
    if (m==3) return (x(0)-0.0)*(x(1)-0.0)/dx/dy;
  }

  AMREX_FORCE_INLINE
  Set::Vector DPhi (int m, const Set::Vector x) const
  {
    if (m==0) return {AMREX_D_DECL(-(dy-x(1))/dx/dy,-(dx-x(0))/dx/dy,0)};
    if (m==1) return {AMREX_D_DECL((dy-x(1))/dx/dy,-x(0)/dx/dy,0)};
    if (m==2) return {AMREX_D_DECL(-x(1)/dx/dy,(dx-x(0))/dx/dy,0)};
    if (m==3) return {AMREX_D_DECL(x(1)/dx/dy, x(0)/dx/dy,0)};
  }

  Set::Vector Qpt (const int m) const { return qpts[m]; }
  Set::Scalar Qwt (const int m) const { return qwts[m]; }

private:  
  Set::Scalar dx, dy;
  const std::array<Set::Vector,4> qpts ={Set::Vector(0.5 - 0.5/sqrt(3.0),0.5 - 0.5/sqrt(3.0)),
                                        Set::Vector(0.5 - 0.5/sqrt(3.0),0.5 + 0.5/sqrt(3.0)),
                                        Set::Vector(0.5 + 0.5/sqrt(3.0),0.5 - 0.5/sqrt(3.0)),
                                        Set::Vector(0.5 + 0.5/sqrt(3.0),0.5 + 0.5/sqrt(3.0))};
  const std::array<Set::Scalar,4> qwts ={1.0, 1.0, 1.0, 1.0};

};
}
#endif
