#ifndef PFBoundarySin_H
#define PFBoundarySin_H

#include "AMReX.H"
#include "PFBoundary.H"

class PFBoundarySin : public PFBoundary
{
public:
  PFBoundarySin(amrex::Real _theta0, amrex::Real _sigma0, amrex::Real _sigma1)
  {
    theta0 = _theta0;
    sigma0 = _sigma0;
    sigma1 = _sigma1;
  };
  amrex::Real W(amrex::Real theta)
  {
    
    return sin(2*(theta-theta0))*sin(2*(theta-theta0))*sigma1+sigma0;
  };
  amrex::Real DW(amrex::Real theta)
  {
    
    return 4*sin(2*(theta-theta0))*cos(2*(theta-theta0));;
  };
  
  //virtual amrex::Real operator()(amrex::Real theta)
  //{
  //amrex::Real d_gbenergy;
  //d_gbenergy= 
  //return d_gbenergy;    
  //};
  
private:
  amrex::Real theta0, sigma0, sigma1;
};


// class PFBoundaryTab : public PFBoundary
// {
// public:
//   PFBoundarySin(std::string filename)
//   {
//     // initialize thetas and energies with data from filename
//   };
//   virtual amrex::Real operator()(amrex::Real theta)
//   {
//     return 1.; // some interpolation function
//   };
// private:
//   std::vector<amrex::Real> thetas, energies;
// };

#endif
