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
    //return sin(2.0*(theta-theta0))*sin(2.0*(theta-theta0))*sigma1+sigma0;
    return sigma0 + 0.5*sigma1*(1.0 - cos(4.0*(theta-theta0)));
  };
  amrex::Real DW(amrex::Real theta)
  {
    //return 4.0*sigma1*sin(2.0*(theta-theta0))*cos(2.0*(theta-theta0));
    return 2.0*sigma1*sin(4.0*(theta-theta0));
  };
  amrex::Real DDW(amrex::Real theta)
  {
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


class PFBoundaryTab : public PFBoundary
{
public:
  PFBoundaryTab(std::string filename)
  {
    // initialize thetas and energies with data from filename
    // current directory ./alamo/Energy_Inputs/110_fcc_output.dat
    std::vector<amrex::Real>   thetas;
    std::vector<amrex::Real>   energies;
    //std::stringstream   element;
    thetas.clear();
    energies.clear();
    
    std::ifstream infile;
    infile.open(filename , std::ios_base::in);
    
    std::ifstream& operator >> (std::ifstream& infile, std::vector<amrex::Real> table);

    amrex::Real tmp_thetas, tmp_energies;
    std:: string line;
    //getline( infile , line );
    
    //std::stringstream ss( line );
    //double f;
    // while(ss >> f)
    // {
    //table.push_back(f);
    //}
    
    //    length_of_vector = table.size()/2;
    int i=0;
    while(getline(infile,line))
      {
	
	infile >> tmp_thetas >> tmp_energies;
	thetas[i]= tmp_thetas;
	energies[i]=tmp_energies;
	i++;
      }
    length_of_vector = thetas.size(); //or i, both should work



  };

  

  
  amrex::Real W(amrex::Real theta)
  {
    amrex::Real W;
 
    int i =(int)(2*theta);    
  
    W= (energies[(i)%length_of_vector]);
    
   // some interpolation function
};

 amrex::Real DW(amrex::Real theta)
 {
   amrex::Real d_W;
   int i =(int)(2*theta);    
   amrex::Real d_theta;
   d_theta= thetas[1]-thetas[0];
   
   d_W=(energies[(i+1)%length_of_vector]-energies[(i-1)%length_of_vector])/(2*d_theta);

   
   return d_W;
 };

private:
  std::vector<amrex::Real> table, thetas, energies;
  int length_of_vector;
  //std::vector<amrex::Real> thetas;
};

#endif
