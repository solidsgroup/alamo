#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include <AMReX_AmrCore.H>

#include "Integrator/HeatConduction/HeatConduction.H"


int main (int argc, char* argv[])
{
  if (argc < 2)
    {
      std::cout << "Missing input file" << std::endl;
      exit(-1);
    }
  amrex::Initialize(argc,argv);

  {
    Integrator::HeatConduction myamr;
    myamr.InitData();
    myamr.Evolve();
  }

  amrex::Finalize();
  return 0;

} 
