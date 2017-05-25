
#include <AMReX_ParmParse.H>

#include <PFAmr.H>

using namespace amrex;

void
PFAmr::ErrorEst (int lev, TagBoxArray& tags, Real time, int /*ngrow*/)
{
  const Real* dx      = geom[lev].CellSize();
  const Real* prob_lo = geom[lev].ProbLo();

  const MultiFab& state = *phi_new[0][lev];

  {
    Array<int>  itags;
	
    for (MFIter mfi(state,true); mfi.isValid(); ++mfi)
      {
	const Box&  bx  = mfi.tilebox();

	TagBox&     tag  = tags[mfi];
	    
	amrex::BaseFab<Real> &new_phi_box = (*phi_new[0][lev])[mfi];

	for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
	  for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
	    {
	      amrex::Real gradx = (new_phi_box(amrex::IntVect(i+1,j)) - new_phi_box(amrex::IntVect(i,j)))/dx[0]/dx[0];
	      amrex::Real grady = (new_phi_box(amrex::IntVect(i,j+1)) - new_phi_box(amrex::IntVect(i,j)))/dx[1]/dx[1];
	      if (dx[0]*sqrt(gradx*gradx + grady*grady)>2.)
		{
		  tag(amrex::IntVect(i,j)) = TagBox::SET;
		}
	    }
      }
  }
}
