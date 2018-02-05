#ifndef PFBoundarySin_H
#define PFBoundarySin_H

#include <iostream>
#include <fstream>

#include "AMReX.H"
#include "PFBoundary.H"

#define PI 3.14159265 

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
    //theta = theta*PI/180;
    //theta0 = theta0*PI/180;
    //return sin(2.0*(theta-theta0))*sin(2.0*(theta-theta0))*sigma1+sigma0;
    return sigma0 + 0.5*sigma1*(1.0 - cos(4.0*(theta-theta0)));
  };
  amrex::Real DW(amrex::Real theta)
  {
    //theta = theta*PI/180;
    //theta0 = theta0*PI/180;
    //return 4.0*sigma1*sin(2.0*(theta-theta0))*cos(2.0*(theta-theta0));
    return 2.0*sigma1*sin(4.0*(theta-theta0));
  };
  amrex::Real DDW(amrex::Real theta)
  {
    //theta = theta*PI/180;
    //    theta0 = theta0*PI/180;
    return 8.0*sigma1*cos(4.0*(theta-theta0));
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

#endif
