
#include <AMReX_ParmParse.H>

#include <PFFlame.H>

using namespace amrex;

void
PFFlame::ErrorEst (int lev, TagBoxArray& tags, Real time, int /*ngrow*/)
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
	      amrex::Real x = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];

	      // for (int n = 0; n < number_of_grains; n++)
	      // 	{
	      amrex::Real gradx = (new_phi_box(amrex::IntVect(i+1,j),0) - new_phi_box(amrex::IntVect(i-1,j),0))/(2.*dx[0]);
	      amrex::Real grady = (new_phi_box(amrex::IntVect(i,j+1),0) - new_phi_box(amrex::IntVect(i,j-1),0))/(2.*dx[1]);
	      //std::cout << x << " : " << new_phi_box(amrex::IntVect(i-1,j),0) << std::endl;
	      if (dx[0]*sqrt(gradx*gradx + grady*grady)>0.1)
		{
		  //std::cout << "tagging" << std::endl; 
		  tag(amrex::IntVect(i,j)) = TagBox::SET;
		}
		  //		}
	       // if (new_phi_box(amrex::IntVect(i,j)) > 0.1 && new_phi_box(amrex::IntVect(i,j)) < 0.9) tag(amrex::IntVect(i,j)) = TagBox::SET;
	    }
      }
    //exit(0);
  }
}
