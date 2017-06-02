#include <iostream>
#include <fstream>
#include <iomanip>

#include "PFAmr.H"

using namespace amrex;

int main (int argc, char* argv[])
{
  if (argc < 2)
    {
      std::cout << "Missing input file" << std::endl;
      exit(-1);
    }
  amrex::Initialize(argc,argv);
  srand(1.0*ParallelDescriptor::MyProc());
  PFAmr pfamr;
  pfamr.InitData();
  pfamr.Evolve();
  amrex::Finalize();
  return 0;
} 
