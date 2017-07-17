#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include "PFAmr.H"
//#include "PFBoundary.H"
#include "PFBoundarySin.H"

int main (int argc, char* argv[])
{
  if (argc < 2)
    {
      std::cout << "Missing input file" << std::endl;
      exit(-1);
    }
  amrex::Initialize(argc,argv);
  srand(1.0*amrex::ParallelDescriptor::MyProc());
  {
    PFAmr pfamr;
    pfamr.InitData();
    pfamr.Evolve();
  }
  
  amrex::Finalize();
} 
