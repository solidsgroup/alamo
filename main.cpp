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
  std::cout << __FILE__ << ":" << __LINE__ << std::endl;
  amradv.Evolve();
  std::cout << __FILE__ << ":" << __LINE__ << std::endl;
  amrex::Finalize();
  return 0;
}
