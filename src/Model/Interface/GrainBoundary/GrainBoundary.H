#ifndef MODEL_INTERFACE_GRAINBOUNDARY_H
#define MODEL_INTERFACE_GRAINBOUNDARY_H

#include <AMReX.H>
#include <AMReX_AmrCore.H>

#include <iostream>
#include <fstream>

namespace Model
{
namespace Interface
{
namespace GrainBoundary
{
class GrainBoundary
{
public:
  GrainBoundary() {};
  virtual amrex::Real W(amrex::Real theta) = 0;
  virtual amrex::Real DW(amrex::Real theta) = 0;
  virtual amrex::Real DDW(amrex::Real theta) = 0;
  bool Test()
  {
    amrex::Real small = 1E-12;
    amrex::Real tolerance = 1E-3;
   
    for (int i = 0; i<20; i++)
      {
	amrex::Real theta = 2.0*pi*((amrex::Real)rand()/(amrex::Real)RAND_MAX);
	amrex::Real numerical_DW = (W(theta+small) - W(theta))/small;
	if (fabs(numerical_DW-DW(theta))/DW(theta) > tolerance)
	  {
	    std::cout << "DW Exact does not match DW Numerical" << std::endl;
	    std::cout << "       Theta:         " << theta << std::endl;
	    std::cout << "       DW Exact:      " << DW(theta) << std::endl;
	    std::cout << "       DW Numerical:  " << numerical_DW << std::endl;
	    std::cout << "       Difference:    " << fabs(numerical_DW-DW(theta)) << std::endl;
	    return false;
	  }
	amrex::Real numerical_DDW = (DW(theta+small) - DW(theta))/small;
	if (fabs(numerical_DDW-DDW(theta)) > tolerance)
	  {
	    std::cout << "DDW Exact does not match DDW Numerical" << std::endl;
	    std::cout << "       Theta:         " << theta << std::endl;
	    std::cout << "       DDW Exact:     " << DDW(theta) << std::endl;
	    std::cout << "       DDW Numerical: " << numerical_DDW << std::endl;
	    std::cout << "       Difference:    " << fabs(numerical_DDW-DDW(theta)) << std::endl;
	    return false;
	  }
      }
    return true;
  }

  void ExportToFile(std::string filename, amrex::Real dTheta)
  {
    std::ofstream outFile;
    outFile.open(filename);
    
    for(amrex::Real theta=0; theta<2*pi; theta=theta+dTheta)
      {
	outFile << theta << " " << W(theta) << std::endl;   
      }
    outFile.close();
    
  }

protected:
  static constexpr amrex::Real pi = 3.14159265359;
};
}
}
}

#endif
