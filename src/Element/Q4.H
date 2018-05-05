#ifndef ELEMENT_Q4_
#define ELEMENT_Q4_

#include "Set/Set.H"
#include "Element.H"

template<ShapeFunctionType sftype> struct Element{};

template <> struct Element<Q4>
{
public:
  Element<Q4>(amrex::Real _dx, amrex::Real _dy): dx(_dx), dy(_dy) {} ;

  template<int n>
  inline amrex::Real Phi(Set::Vector);

  template<int n>
  inline Set::Vector DPhi(Set::Vector);

  inline std::array<std::array<Set::Vector,4>,4> DPhis();

  inline std::array<amrex::Real,4> Ws();

  template<int n>
  inline Set::Vector QPoint();

  template<int n>
  inline amrex::Real QWeight();

private:
  amrex::Real dx, dy;
};

template <> inline amrex::Real Element<Q4>::Phi<1>(Set::Vector x)
{return (dx-x(0))*(dy-x(1))/dx/dy;}
 
template <> inline amrex::Real Element<Q4>::Phi<2>(Set::Vector x)
{return (x(0)-0.0)*(dy-x(1))/dx/dy;} 

template <> inline amrex::Real Element<Q4>::Phi<3>(Set::Vector x)
{return (dx-x(0))*(x(1)-0.0)/dx/dy;} 

template <> inline amrex::Real Element<Q4>::Phi<4>(Set::Vector x)
{return (x(0)-0.0)*(x(1)-0.0)/dx/dy;} 

template <> inline Set::Vector Element<Q4>::DPhi<1>(Set::Vector x)
{return {-(dy-x(1))/dx/dy,-(dx-x(0))/dx/dy};}

template <> inline Set::Vector Element<Q4>::DPhi<2>(Set::Vector x)
{return {(dy-x(1))/dx/dy,-x(0)/dx/dy};} 

template <> inline Set::Vector Element<Q4>::DPhi<3>(Set::Vector x)
{return {-x(1)/dx/dy,(dx-x(0))/dx/dy};} 

template <> inline Set::Vector Element<Q4>::DPhi<4>(Set::Vector x) 
{return {x(1)/dx/dy, x(0)/dx/dy};} 

template <> inline Set::Vector Element<Q4>::QPoint<1>()
{return {dx*(sqrt(3.0)-1.0)/2.0/sqrt(3.0), dy*(sqrt(3.0)-1.0)/2.0/sqrt(3.0)};}

template <> inline Set::Vector Element<Q4>::QPoint<2>()
{return {dx*(sqrt(3.0)+1.0)/2.0/sqrt(3.0), dy*(sqrt(3.0)-1.0)/2.0/sqrt(3.0)};}

template <> inline Set::Vector Element<Q4>::QPoint<3>()
{return {dx*(sqrt(3.0)-1.0)/2.0/sqrt(3.0), dy*(sqrt(3.0)+1.0)/2.0/sqrt(3.0)};}

template <> inline Set::Vector Element<Q4>::QPoint<4>()
{return {dx*(sqrt(3.0)+1.0)/2.0/sqrt(3.0), dy*(sqrt(3.0)+1.0)/2.0/sqrt(3.0)};}

template <> inline amrex::Real Element<Q4>::QWeight<1>() {return 0.25*dx*dy; }
template <> inline amrex::Real Element<Q4>::QWeight<2>() {return 0.25*dx*dy; }
template <> inline amrex::Real Element<Q4>::QWeight<3>() {return 0.25*dx*dy; }
template <> inline amrex::Real Element<Q4>::QWeight<4>() {return 0.25*dx*dy; }

inline std::array<std::array<Set::Vector,4>,4> Element<Q4>::DPhis()
{
  return
    {{{DPhi<1>(QPoint<1>()), DPhi<1>(QPoint<2>()), DPhi<1>(QPoint<3>()), DPhi<1>(QPoint<4>())},
	  {DPhi<2>(QPoint<1>()), DPhi<2>(QPoint<2>()), DPhi<2>(QPoint<3>()), DPhi<2>(QPoint<4>())},
	    {DPhi<3>(QPoint<1>()), DPhi<3>(QPoint<2>()), DPhi<3>(QPoint<3>()), DPhi<3>(QPoint<4>())},
	      {DPhi<4>(QPoint<1>()), DPhi<4>(QPoint<2>()), DPhi<4>(QPoint<3>()), DPhi<4>(QPoint<4>())}}};
}

inline std::array<amrex::Real,4> Element<Q4>::Ws()
{
  return {QWeight<1>(),QWeight<2>(),QWeight<3>(),QWeight<4>()};
}

#endif
