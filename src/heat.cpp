#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include <AMReX_AmrCore.H>

#include "HeatConduction/Integrator.H"
#include "GeneralAMRIntegrator/GeneralAMRIntegrator.H"


int main (int argc, char* argv[])
{
  if (argc < 2)
    {
      std::cout << "Missing input file" << std::endl;
      exit(-1);
    }
  amrex::Initialize(argc,argv);

  {
    HeatConduction::Integrator myamr;
    myamr.InitData();
    myamr.Evolve();
  }

  amrex::Finalize();
  return 0;

} 
