#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"

#if AMREX_SPACEDIM == 2

#include "PhaseFieldMicrostructure/PhaseFieldMicrostructure.H"

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
    PhaseFieldMicrostructure model;
    model.InitData();
    model.Evolve();
  }
  
  amrex::Finalize();
} 
#else
int main()
{
  std::cout << "This program works in 2D only, but AMREX_SPACEDIM=" << AMREX_SPACEDIM << std::endl;
}
#endif
