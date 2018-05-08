#ifndef GeneralAMRIntegratorBC_H_
#define GeneralAMRIntegratorBC_H_

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCRec.H>
#include <AMReX_PhysBCFunct.H>

#if 0
Just added a variable for NEUMANN BC. 
The value is set to 5 because the it goes upto 4 in AMREX_BC_TYPES.H
Maybe we need to change it in the future. It is not being used anywhere
else for now.
#endif
#define NEUMANN 5


class GeneralAMRIntegratorPhysBC
  : public amrex::PhysBCFunctBase
{

public:
  GeneralAMRIntegratorPhysBC (amrex::Array<amrex::Geometry> &_geom)
    : geom(_geom)
  {
    {
      amrex::ParmParse pp("bc");
      amrex::Array<std::string> bc_hi_str(BL_SPACEDIM);
      amrex::Array<std::string> bc_lo_str(BL_SPACEDIM);
      pp.queryarr("hi",bc_hi_str,0,BL_SPACEDIM);
      pp.queryarr("lo",bc_lo_str,0,BL_SPACEDIM);

      for (int i=0;i<BL_SPACEDIM;i++)
	{
	  if (bc_hi_str[i] == "REFLECT_ODD"	) bc_hi[i] = REFLECT_ODD; 
	  if (bc_hi_str[i] == "INT_DIR"		) bc_hi[i] = INT_DIR;
	  if (bc_hi_str[i] == "REFLECT_EVEN"	) bc_hi[i] = REFLECT_EVEN;
	  if (bc_hi_str[i] == "FOEXTRAP"		) bc_hi[i] = FOEXTRAP;
	  if (bc_hi_str[i] == "EXT_DIR"		) bc_hi[i] = EXT_DIR;
	  if (bc_hi_str[i] == "HOEXTRAP"		) bc_hi[i] = HOEXTRAP;

	  if (bc_lo_str[i] == "REFLECT_ODD"	) bc_lo[i] = REFLECT_ODD;
	  if (bc_lo_str[i] == "INT_DIR"		) bc_lo[i] = INT_DIR;
	  if (bc_lo_str[i] == "REFLECT_EVEN"	) bc_lo[i] = REFLECT_EVEN;
	  if (bc_lo_str[i] == "FOEXTRAP"		) bc_lo[i] = FOEXTRAP;
	  if (bc_lo_str[i] == "EXT_DIR"		) bc_lo[i] = EXT_DIR;
	  if (bc_lo_str[i] == "HOEXTRAP"		) bc_lo[i] = HOEXTRAP;

    // Added for Neumann BC.
	  if (bc_hi_str[i] == "NEUMANN"		) bc_hi[i] = NEUMANN;
	  if (bc_lo_str[i] == "NEUMANN"		) bc_lo[i] = NEUMANN;
	}
    
      // todo -- add ability to specify Dirichlet/Neumann BC values
      
      if (bc_lo[0] == EXT_DIR) pp.getarr("lo_1",bc_lo_1); if(bc_hi[0] == EXT_DIR) pp.getarr("hi_1",bc_hi_1);
      if (bc_lo[1] == EXT_DIR) pp.getarr("lo_2",bc_lo_2); if(bc_hi[1] == EXT_DIR) pp.getarr("hi_2",bc_hi_2);

      // Added for Neumann BC
      // By default it is assumed that the flux values are along the normals. So for left boundary, a positive flux means
      // that things are leaving the domain. 
      if(bc_lo[0] == NEUMANN) pp.getarr("lo_1",bc_lo_1); if(bc_hi[0] == NEUMANN) pp.getarr("hi_1",bc_hi_1);
      if(bc_lo[1] == NEUMANN) pp.getarr("lo_2",bc_lo_2); if(bc_hi[1] == NEUMANN) pp.getarr("hi_2",bc_hi_2);

#if AMREX_SPACEDIM > 2
      if (bc_lo[2] == EXT_DIR) pp.getarr("lo_3",bc_lo_3); if(bc_hi[2] == EXT_DIR) pp.getarr("hi_3",bc_hi_3);
      // Added for Neumann BC
      if(bc_lo[2] == NEUMANN) pp.getarr("lo_3",bc_lo_3); if(bc_hi[2] == NEUMANN) pp.getarr("hi_3",bc_hi_3);
#endif
    }
  }

  virtual void FillBoundary (amrex::MultiFab& mf, int, int, amrex::Real /*time*/) override
  {
    amrex::Box domain(geom[lev].Domain());

    mf.FillBoundary(geom[lev].periodicity());

    // Added for Neumann BC
    const Real* dx = geom[lev].CellSize();

    for (amrex::MFIter mfi(mf,true); mfi.isValid(); ++mfi)
      {
	const amrex::Box& box = mfi.tilebox();

	amrex::BaseFab<amrex::Real> &phi_box = mf[mfi];

	for (int i = box.loVect()[0] - mf.nGrow(); i<=box.hiVect()[0] + mf.nGrow(); i++)
	  for (int j = box.loVect()[1] - mf.nGrow(); j<=box.hiVect()[1] + mf.nGrow(); j++)
#if BL_SPACEDIM>2
	  for (int k = box.loVect()[2] - mf.nGrow(); k<=box.hiVect()[2] + mf.nGrow(); k++)
#endif
	    for (int n = 0; n < mf.nComp(); n++)
	      {
		if (i < domain.loVect()[0]) // Left boundary

		{
		  if (bc_lo[0] == EXT_DIR)
                  {
		    phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = bc_lo_1[n];
		  }
		  else if(bc_lo[0] == NEUMANN)
		  {
		    phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = phi_box(amrex::IntVect(AMREX_D_DECL(i+1,j,k)),n) - bc_lo_1[n]*dx[0];
		  }
		}

		if (i > domain.hiVect()[0]) // Right boundary
		{
		  if (bc_hi[0] == EXT_DIR)
                  {
		    phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = bc_hi_1[n];
		  }
		  else if(bc_hi[0] == NEUMANN)
		  {
		    phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = phi_box(amrex::IntVect(AMREX_D_DECL(i-1,j,k)),n) - bc_hi_1[n]*dx[0];
		  }
		}

		if (j < domain.loVect()[1]) // Bottom boundary
		{
		  if (bc_lo[1] == EXT_DIR)
		  {
		    phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = bc_lo_2[n];
		  }
		  else if (bc_lo[1] == NEUMANN)
		  {
		    phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = phi_box(amrex::IntVect(AMREX_D_DECL(i,j+1,k)),n) - bc_lo_2[n]*dx[1];
		  }
		}

		if (j > domain.hiVect()[1]) // Top boundary
		{
		  if (bc_hi[1] == EXT_DIR)
		  {
		    phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = bc_hi_2[n];
		  }
		  else if (bc_hi[1] == NEUMANN)
		  {
		    phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = phi_box(amrex::IntVect(AMREX_D_DECL(i,j-1,k)),n) - bc_hi_2[n]*dx[1];
		  }
		}


#if BL_SPACEDIM>2
		if (k < domain.loVect()[2])
		{
		  if (bc_lo[2] == EXT_DIR)
		  {
		    phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = bc_lo_3[n];
		  }
		  else if (bc_lo[2] == NEUMANN)
		  {
		    phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k+1)),n) - bc_lo_3[n]*dx[2];
		  }
		}

		if (k > domain.hiVect()[2])
		{
		  if (bc_hi[2] == EXT_DIR)
		  {
		    phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = bc_hi_3[n];
		  }
		  else if(bc_hi[2] == NEUMANN)
		  {
		    phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k-1)),n) - bc_hi_3[n]*dx[2];
		  }
		}
#endif

	      }
      }
  }
  void SetLevel(int _lev) {lev=_lev;}
  amrex::BCRec GetBCRec() {return amrex::BCRec(bc_lo,bc_hi);}
protected:
  int bc_lo[BL_SPACEDIM];
  int bc_hi[BL_SPACEDIM];
  int bc_hi_flux[BL_SPACEDIM];
  int bc_lo_flux[BL_SPACEDIM];
  amrex::Vector<amrex::Real> bc_lo_1, bc_hi_1;
  amrex::Vector<amrex::Real> bc_lo_2, bc_hi_2;

#if BL_SPACEDIM > 2
  amrex::Vector<amrex::Real> bc_lo_3, bc_hi_3;
#endif

  std::vector<amrex::Geometry> &geom;
  int lev=0;
};

#endif
