#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include "AMReX_ParallelDescriptor.H"
#include "PFFlame.H"

#include "GeneralAMRIntegrator.H"
#include "GeneralAMRIntegratorBC.H"

int main (int argc, char* argv[])
{

  if (argc < 2)
    {
      std::cout << "Missing input file" << std::endl;
      exit(-1);
    }
  amrex::Initialize(argc,argv);

  PFFlame pfamr;
  pfamr.InitData();
  pfamr.Evolve();
  
  amrex::Finalize();
} 
