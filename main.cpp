#include <iostream>
#include <fstream>
#include <iomanip>

#include "AmrAdv.H"

using namespace amrex;

int main (int argc, char* argv[])
{
  amrex::Initialize(argc,argv);
  AmrAdv amradv;
  amradv.InitData();
  //amradv.MakeNewGrids();
  amradv.Evolve();
  amrex::Finalize();
  return 0;
}
