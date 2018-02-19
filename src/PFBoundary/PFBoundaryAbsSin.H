#ifndef PFBoundaryAbsSin_H
#define PFBoundaryAbsSin_H

#include <iostream>
#include <fstream>

#include "AMReX.H"
#include "PFBoundary/PFBoundary.H"

#define PI 3.14159265 

class PFBoundaryAbsSin : public PFBoundary
{
public:
  PFBoundaryAbsSin(amrex::Real _theta0, amrex::Real _sigma0, amrex::Real _sigma1)
  {       
    theta0 = _theta0;
    sigma0 = _sigma0;
    sigma1 = _sigma1;
  };
  amrex::Real W(amrex::Real theta)
  {
    //sigma(theta)=sigma0+sigma1*fabs(sin(n*(theta-theta0)))
    //n=2:  
    return sigma0+sigma1*fabs(sin(2*(theta-theta0)));  
    
  };
  amrex::Real DW(amrex::Real theta)
  {
    //sigma'(theta)=(n*sigma1*sin(n*(theta-theta0))*cos(n*(theta-theta0)))/fabs(sin(n*(theta-theta0)))
    //n=2:
    return (2*sigma1*sin(2*(theta-theta0))*cos(2*(theta-theta0)))/fabs(sin(2*(theta-theta0)));
  };
  amrex::Real DDW(amrex::Real theta)
  {
    //sigma''(theta)=-n^2*sigma1*fabs(sin(n*(theta-theta0)))
    //n=2:
    return -4*sigma1*fabs(sin(2*(theta-theta0)));
  };
 
private:
  amrex::Real theta0, sigma0, sigma1;
};

#endif
