#ifndef MODEL_SOLID_ELASTIC_H_
#define MODEL_SOLID_ELASTIC_H_

#include <AMReX.H>
#include "Model/Solid/Solid.H"

namespace Model
{
namespace Solid
{

class Elastic : public Solid
{
public:
  Elastic()
  {
    mu = 1.0; lambda=1.0;

    for (int i=0; i<2; i++)
      for (int j=0; j<2; j++)
	for (int k=0; k<2; k++)
	  for (int l=0; l<2; l++)
	    {
	      C[i][j][k][l] = 0.0;
	      if (i==k && j==l) 
		C[i][j][k][l] += mu;
	      if (i==l && j==k) 
		C[i][j][k][l] += mu;
	      if (i==j && k==l) 
		C[i][j][k][l] += lambda;
	    }
  };

  void 
  W(amrex::Real &w,
    const std::array<std::array<amrex::Real,2>,2> &gradu)
  {
    w = 0.0;
    for (int i=0; i<2; i++)
      for (int j=0; j<2; j++)
	for (int k=0; k<2; k++)
	  for (int l=0; l<2; l++)
	    w += 0.5*gradu[i][j]*C[i][j][k][l]*gradu[k][l];

  }

  
  void
  DW(std::array<std::array<amrex::Real,2>,2> &dw,
     const std::array<std::array<amrex::Real,2>,2> &gradu)
  {
    std::array<std::array<std::array<std::array<amrex::Real,2>,2>,2>,2> ddw;
    for (int i=0; i<2; i++)
      for (int j=0; j<2; j++)
	{
	  dw[i][j] = 0.0;
	  for (int k=0; k<2; k++)
	    for (int l=0; l<2; l++)
	      dw[i][j] += C[i][j][k][l]*gradu[k][l];
	}
  }
  
  void
  DDW(std::array<std::array<std::array<std::array<amrex::Real,2>,2>,2>,2> &ddw,
      const std::array<std::array<amrex::Real,2>,2> &gradu)
  {
    for (int i=0; i<2; i++)
      for (int j=0; j<2; j++)
	for (int k=0; k<2; k++)
	  for (int l=0; l<2; l++)
	    ddw[i][j][k][l] = C[i][j][k][l];
  }

private:
  amrex::Real mu, lambda;
  std::array<std::array<std::array<std::array<amrex::Real,2>,2>,2>,2> C;
};

}
}


#endif


