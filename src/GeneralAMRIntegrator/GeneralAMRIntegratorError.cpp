
#include <AMReX_ParmParse.H>

#include "GeneralAMRIntegrator.H"

using namespace amrex;

void
GeneralAMRIntegrator::ErrorEst (int lev, TagBoxArray& tags, Real time, int /*ngrow*/)
{
  // const Real* dx      = geom[lev].CellSize();

  // const MultiFab& state = *phi_new[0][lev];
  // {
  //   Array<int>  itags;
	
  //   for (MFIter mfi(state,true); mfi.isValid(); ++mfi)
  //     {
  // 	const Box&  bx  = mfi.tilebox();

  // 	TagBox&     tag  = tags[mfi];
	    
  // 	amrex::BaseFab<Real> &new_phi_box = (*phi_new[0][lev])[mfi];

  // 	for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
  // 	  for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
  // 	    {
  // 		  if (i==0 && j==0)tag(amrex::IntVect(i,j)) = TagBox::SET;
  // 	    }
  // 	for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
  // 	  for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
  // 	    for (int k = bx.loVect()[2]; k<=bx.hiVect()[2]; k++)
  // 	      {
  // 		    if (i==0 && j==0 && k == 0 )tag(amrex::IntVect(i,j)) = TagBox::SET;
  // 	      }

  //     }
  // }
}
