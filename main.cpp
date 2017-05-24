#include <iostream>
#include <fstream>
#include <iomanip>

#include "PFAmr.H"

using namespace amrex;

int main (int argc, char* argv[])
{
  amrex::Initialize(argc,argv);
  PFAmr pfamr;
  pfamr.InitData();
  pfamr.Evolve();
  amrex::Finalize();
  return 0;
}
