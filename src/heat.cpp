#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include <AMReX_AmrCore.H>

#include "HeatConduction/Integrator.H"
#include "GeneralAMRIntegrator/GeneralAMRIntegrator.H"


class Integrator : public GeneralAMRIntegrator
{
public:

  /// \brief Read in parameters and register field variables
  Integrator()  : GeneralAMRIntegrator(),   mybc(geom){
    RegisterNewFab(Temp,     mybc, number_of_components, number_of_ghost_cells, "Temp");
  };
  ~Integrator() {
    std::cout << "DEBUG " << __FILE__ << ":" << __LINE__ << std::endl;
  };

protected:

  void Initialize (int lev) {
  };

  void Advance (int lev, Real /*time*/, Real dt) {};


  void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/) {};

private:
  int number_of_components = 1;					///< Number of components
  int number_of_ghost_cells = 1;				///< Number of ghost cells
  amrex::Array<std::unique_ptr<amrex::MultiFab> > Temp;		///< Temperature field variable (current timestep)
  GeneralAMRIntegratorPhysBC mybc;				///< Stock generic boundary condition object
};



int main (int argc, char* argv[])
{
  if (argc < 2)
    {
      std::cout << "Missing input file" << std::endl;
      exit(-1);
    }
  amrex::Initialize(argc,argv);

  {
    Integrator myamr;
    myamr.InitData();
    myamr.Evolve();
    std::cout << "DEBUG " << __FILE__ << ":" << __LINE__ << std::endl;
  }
  std::cout << "DEBUG " << __FILE__ << ":" << __LINE__ << std::endl;

  amrex::Finalize();
  return 0;

  // srand(1.0*amrex::ParallelDescriptor::MyProc());
  // {
  //   HeatConduction::Integrator model;
  //   model.InitData();
  //   model.Evolve();
  //   std::cout << "Outside of Evolve"<< std::endl;
  // }
  
  // std::cout << "Done! Calling finalize now"<< std::endl;
  // amrex::Finalize();
  // std::cout << "Called finalize"<< std::endl;

  // return 0;
} 
