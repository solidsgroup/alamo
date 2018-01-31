#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include "AMReX_ParallelDescriptor.H"
#include "PFFlame.H"

#include "GeneralAMRIntegrator.H"
#include "GeneralAMRIntegratorBC.H"

class Example : public GeneralAMRIntegrator
{
public:
  Example() :
    GeneralAMRIntegrator(), 
    mybc(geom)
  {
    RegisterNewFab(Temp,     mybc, number_of_components, number_of_ghost_cells,"Temp");
    RegisterNewFab(Temp_old, mybc, number_of_components, number_of_ghost_cells,"Temp old");
  }

protected:
  void
  Advance (int lev, Real time, Real dt)
  {
    std::swap(*Temp[lev], *Temp_old[lev]);

    FillPatch(lev,t_old[lev],*Temp_old[lev],mybc,0);

    const Real* dx = geom[lev].CellSize();

    for ( amrex::MFIter mfi(*Temp[lev],true); mfi.isValid(); ++mfi )
         {
	   const amrex::Box& bx = mfi.tilebox();

	   amrex::BaseFab<Real> &Temp_old_box = (*Temp_old[lev])[mfi];
	   amrex::BaseFab<Real> &Temp_box = (*Temp[lev])[mfi];

	   for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
	     for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
	       {
		 Temp_box(amrex::IntVect(i,j))
		   = Temp_old_box(amrex::IntVect(i,j))
		   + dt * ((Temp_old_box(amrex::IntVect(i+1,j)) + Temp_old_box(amrex::IntVect(i-1,j)) - 2*Temp_old_box(amrex::IntVect(i,j))) / dx[0] / dx[0] +
			   (Temp_old_box(amrex::IntVect(i,j+1)) + Temp_old_box(amrex::IntVect(i,j+1)) - 2*Temp_old_box(amrex::IntVect(i,j))) / dx[1] / dx[1]  );
	       }
       }
  }

  void
  Initialize (int lev)
  {
    const amrex::Real width = geom[lev].ProbHi()[0] - geom[lev].ProbHi()[1];
    for (amrex::MFIter mfi(*Temp[lev],true); mfi.isValid(); ++mfi)
      {
	const amrex::Box& box = mfi.tilebox();
	amrex::BaseFab<Real> &Temp_box = (*Temp[lev])[mfi];
	amrex::BaseFab<Real> &Temp_old_box = (*Temp_old[lev])[mfi];
	for (int i = box.loVect()[0]-number_of_ghost_cells; i<=box.hiVect()[0]+number_of_ghost_cells; i++) 
	  for (int j = box.loVect()[1]-number_of_ghost_cells; j<=box.hiVect()[1]+number_of_ghost_cells; j++)
     	    {
	      Temp_box(amrex::IntVect(i,j),0) = width*i*j;
	      Temp_old_box(amrex::IntVect(i,j),0) = Temp_box(amrex::IntVect(i,j),0);
     	    }
      }
  }


  void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/)
  {

    const Real* dx      = geom[lev].CellSize();

    amrex::Array<int>  itags;
 	
    for (amrex::MFIter mfi(*Temp[lev],true); mfi.isValid(); ++mfi)
      {
	const amrex::Box&  bx  = mfi.tilebox();
	amrex::TagBox&     tag  = tags[mfi];
 	    
	amrex::BaseFab<Real> &Temp_box = (*Temp[lev])[mfi];

	for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
	  for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
	    {
	      //if (Temp_box(amrex::IntVect(i,j)) < 1.0) tag(amrex::IntVect(i,j)) = amrex::TagBox::SET;
	      if (i==0 && j==0) tag(amrex::IntVect(i,j)) = amrex::TagBox::SET;
	    }

      }

  }

private:
  int number_of_components = 1;
  int number_of_ghost_cells = 2;

  amrex::Array<std::unique_ptr<amrex::MultiFab> > Temp;
  amrex::Array<std::unique_ptr<amrex::MultiFab> > Temp_old;
  GeneralAMRIntegratorPhysBC mybc;
};

int main (int argc, char* argv[])
{

  if (argc < 2)
    {
      std::cout << "Missing input file" << std::endl;
      exit(-1);
    }
  amrex::Initialize(argc,argv);

  Example ex;
  ex.InitData();
  ex.Evolve();

  // srand(1.0*amrex::ParallelDescriptor::MyProc());
  // {
  //   PFFlame pfamr;
  //   pfamr.InitData();
  //   pfamr.Evolve();
  // }
  
  amrex::Finalize();
} 
