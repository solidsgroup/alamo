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
  srand(0);
  amrex::Initialize(argc,argv);
  PFAmr pfamr;
  pfamr.InitData();
  pfamr.Evolve();
  amrex::Finalize();
  return 0;
}
