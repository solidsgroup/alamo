#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include "AMReX_ParallelDescriptor.H"
#include "PFFlame.H"

#include "GeneralAMRIntegrator.H"


class Example : public GeneralAMRIntegrator
{
  amrex::Array<std::unique_ptr<amrex::MultiFab> > solution;
  Example() : GeneralAMRIntegrator()
  {
    RegisterNewFab(solution,1,1);
  }

  void
  Advance (int lev, Real time, Real dt)
  {
    // TODO - replace this
    // std::swap(phi_old[0][lev], phi_new[0][lev]);
    const Real* dx = geom[lev].CellSize();

    amrex::Array<std::unique_ptr<amrex::MultiFab> > Sborder(number_of_fabs);

    // TODO - fix this
    for (int n=0; n<number_of_fabs; n++)
      Sborder[n].reset(new amrex::MultiFab(grids[lev], dmap[lev], number_of_grains+2,nghost)); 

    FillPatch(lev,t_old[lev],Sborder,0);

    // for ( MFIter mfi(*phi_new[0][lev],true); mfi.isValid(); ++mfi )
    //   {
    //     const Box& bx = mfi.tilebox();

    //     amrex::BaseFab<Real> &old_phi = (*Sborder[0])[mfi];
    //     amrex::BaseFab<Real> &new_phi = (*phi_new[0][lev])[mfi];

    //     amrex::Array<amrex::Real> Laplacian(number_of_fabs);
    //     for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
    // 	for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
    // 	  {
    // 	    new_phi(amrex::IntVect(i,j)) = i*j*dx[0]*dx[1]; // here's where we put the numerics
    // 	  }
    //   }
  }


};

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
    PFFlame pfamr;
    pfamr.InitData();
    pfamr.Evolve();
  }
  
  amrex::Finalize();
} 
