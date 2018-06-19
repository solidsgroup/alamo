#include "Diffusion.H"

#if BL_SPACEDIM == 2

Diffusion::Diffusion() :
  GeneralAMRIntegrator(), 
  mybc(geom)
{
  //
  // Read input parameters
  // Just reading diffusion parameters
  //
  amrex::ParmParse pp("df");
  pp.query("mu",mu);

  RegisterNewFab(conc,     mybc, number_of_components, number_of_ghost_cells, "conc");
  RegisterNewFab(conc_old, mybc, number_of_components, number_of_ghost_cells, "conc_old");
}

void
Diffusion::Advance (int lev, Real /*time*/, Real dt)
{
  std::swap(*conc[lev], *conc_old[lev]);

  const Real* dx = geom[lev].CellSize();

  for ( amrex::MFIter mfi(*conc[lev],true); mfi.isValid(); ++mfi )
    {
      const amrex::Box& bx = mfi.tilebox();

      amrex::BaseFab<Real> &conc_old_box = (*conc_old[lev])[mfi];
      amrex::BaseFab<Real> &conc_box = (*conc[lev])[mfi];

      for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
	for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
	  {
	    conc_box(amrex::IntVect(i,j))
	      = conc_old_box(amrex::IntVect(i,j))
	      + mu*dt*((conc_old_box(amrex::IntVect(i+1,j)) + conc_old_box(amrex::IntVect(i-1,j)) - 2*conc_old_box(amrex::IntVect(i,j))) / dx[0] / dx[0] +
		      (conc_old_box(amrex::IntVect(i,j+1)) + conc_old_box(amrex::IntVect(i,j-1)) - 2*conc_old_box(amrex::IntVect(i,j))) / dx[1] / dx[1]);
	  }
    }
}



void
Diffusion::Initialize (int lev)
{
  //const amrex::Real width = geom[lev].ProbHi()[0] - geom[lev].ProbHi()[1];
  for (amrex::MFIter mfi(*conc[lev],true); mfi.isValid(); ++mfi)
    {
      const amrex::Box& box = mfi.tilebox();
      amrex::BaseFab<Real> &conc_box = (*conc[lev])[mfi];
      amrex::BaseFab<Real> &conc_old_box = (*conc_old[lev])[mfi];
      for (int i = box.loVect()[0]-number_of_ghost_cells; i<=box.hiVect()[0]+number_of_ghost_cells; i++) 
	for (int j = box.loVect()[1]-number_of_ghost_cells; j<=box.hiVect()[1]+number_of_ghost_cells; j++)
	  {
	    conc_box(amrex::IntVect(i,j),0) = 0; 
	    conc_old_box(amrex::IntVect(i,j),0) = conc_box(amrex::IntVect(i,j),0);
	  }
    }
}


void
Diffusion::TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/)
{

  const Real* dx      = geom[lev].CellSize();

  amrex::Array<int>  itags;
 	
  for (amrex::MFIter mfi(*conc[lev],true); mfi.isValid(); ++mfi)
    {
      const amrex::Box&  bx  = mfi.tilebox();
      amrex::TagBox&     tag  = tags[mfi];
 	    
      amrex::BaseFab<Real> &conc_old_box = (*conc_old[lev])[mfi];

      for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
	for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
	  {
	    amrex::Real grad1 = (conc_old_box(amrex::IntVect(i+1,j)) - conc_old_box(amrex::IntVect(i-1,j)))/(2*dx[0]);
	    amrex::Real grad2 = (conc_old_box(amrex::IntVect(i,j+1)) - conc_old_box(amrex::IntVect(i,j-1)))/(2*dx[1]);

	    if ((grad1*grad1 + grad2*grad2)*dx[0]*dx[1] > 0.001) tag(amrex::IntVect(i,j)) = amrex::TagBox::SET;

	  }

    }

}

#endif
