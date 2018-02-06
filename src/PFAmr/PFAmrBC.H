#include <AMReX.H>
#if AMREX_SPACEDIM==2
#ifndef PFAmrBC_H_
#define PFAmrBC_H_

#include <AMReX_ParmParse.H>
#include <AMReX_BCRec.H>
#include <AMReX_PhysBCFunct.H>

// In this tutorial, we assume the boundaries are all periodic.  
// Thus there are no physical boundaries. 

class PFAmrPhysBC
  : public amrex::PhysBCFunctBase
{
public:
  PFAmrPhysBC (amrex::Array<amrex::Geometry> &_geom)
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
	  if (bc_hi_str[i] == "FOEXTRAP"	) bc_hi[i] = FOEXTRAP;
	  if (bc_hi_str[i] == "EXT_DIR"		) bc_hi[i] = EXT_DIR;
	  if (bc_hi_str[i] == "HOEXTRAP"	) bc_hi[i] = HOEXTRAP;

	  if (bc_lo_str[i] == "REFLECT_ODD"	) bc_lo[i] = REFLECT_ODD;
	  if (bc_lo_str[i] == "INT_DIR"		) bc_lo[i] = INT_DIR;
	  if (bc_lo_str[i] == "REFLECT_EVEN"	) bc_lo[i] = REFLECT_EVEN;
	  if (bc_lo_str[i] == "FOEXTRAP"	) bc_lo[i] = FOEXTRAP;
	  if (bc_lo_str[i] == "EXT_DIR"		) bc_lo[i] = EXT_DIR;
	  if (bc_lo_str[i] == "HOEXTRAP"	) bc_lo[i] = HOEXTRAP;
	}
    
      // todo -- add ability to specify Dirichlet/Neumann BC values
      
      if (bc_lo[0] == EXT_DIR) pp.getarr("lo_1",bc_lo_1); if(bc_hi[0] == EXT_DIR) pp.getarr("hi_1",bc_hi_1);
#if BL_SPACEDIM > 1
      if (bc_lo[1] == EXT_DIR) pp.getarr("lo_2",bc_lo_2); if(bc_hi[1] == EXT_DIR) pp.getarr("hi_2",bc_hi_2);
#endif
#if BL_SPACEDIM > 2
      if (BC_lo[2] == EXT_DIR) pp.getarr("lo_3",bc_lo_3); if(bc_hi[2] == EXT_DIR) pp.getarr("hi_3",bc_hi_3);
#endif


    }
  }

  // TODO - UPDATE THIS
  //PFAmrPhysBC (amrex::Array<amrex::Geometry> &_geom, int *bc_lo, int *bc_hi) : geom(_geom)  {}
  //virtual ~PFAmrPhysBC () {}

  virtual void FillBoundary (amrex::MultiFab& mf, int, int, amrex::Real time) override
  {
    amrex::Box domain(geom[lev].Domain());

    mf.FillBoundary(geom[lev].periodicity());

    for (amrex::MFIter mfi(mf,true); mfi.isValid(); ++mfi)
      {
	const amrex::Box& box = mfi.tilebox();

	amrex::BaseFab<amrex::Real> &phi_box = mf[mfi];

	for (int i = box.loVect()[0] - mf.nGrow(); i<=box.hiVect()[0] + mf.nGrow(); i++)
	  {
	    for (int j = box.loVect()[1] - mf.nGrow(); j<=box.hiVect()[1] + mf.nGrow(); j++)
	      for (int n = 0; n < mf.nComp(); n++)
		{
		  if (i < domain.loVect()[0]) // Left boundary
		    {
		      if (bc_lo[0] == EXT_DIR)
			phi_box(amrex::IntVect(i,j),n) = bc_lo_1[n];

		      // Nothing
		    }
		  if (i > domain.hiVect()[0]) // Right boundary
		    {
		      if (bc_hi[0] == EXT_DIR)
			phi_box(amrex::IntVect(i,j),n) = bc_hi_1[n];

		      // Nothing
		    }
		  if (j < domain.loVect()[1]) // Bottom boundary
		    {
		      if (bc_lo[1] == EXT_DIR)
			phi_box(amrex::IntVect(i,j),n) = bc_lo_2[n];
		      // if (n<2)
		      // 	{
		      // 	  if (j == domain.loVect()[1]-1) phi_box(amrex::IntVect(i,j),n) = (amrex::Real)(1-n);//phi_box(amrex::IntVect(i,j+1),n); //1.0;//
		      // 	  if (j == domain.loVect()[1]-2) phi_box(amrex::IntVect(i,j),n) = (amrex::Real)(1-n);//phi_box(amrex::IntVect(i,j+3),n); //1.0;//
		      // 	}
		    }
		  if (j > domain.hiVect()[1]) // Top boundary
		    {
		      if (bc_hi[1] == EXT_DIR)
			phi_box(amrex::IntVect(i,j),n) = bc_hi_2[n];

		      // if (n<2)
		      // 	{
		      // 	  if (j == domain.hiVect()[1]+1) phi_box(amrex::IntVect(i,j),n) = (amrex::Real)(n);//phi_box(amrex::IntVect(i,j-1),n); //0.0;//
		      // 	  if (j == domain.hiVect()[1]+2) phi_box(amrex::IntVect(i,j),n) = (amrex::Real)(n);//phi_box(amrex::IntVect(i,j-3),n); //0.0;//
		      // 	}
		    }
		}
	  }
      }
  }
  void SetLevel(int _lev) {lev=_lev;}
  amrex::BCRec GetBCRec() {return amrex::BCRec(bc_lo,bc_hi);}
private:
  int bc_lo[BL_SPACEDIM];
  int bc_hi[BL_SPACEDIM];
  amrex::Array<amrex::Real> bc_lo_1, bc_hi_1;
#if BL_SPACEDIM > 1
  amrex::Array<amrex::Real> bc_lo_2, bc_hi_2;
#endif
#if BL_SPACEDIM > 2
  amrex::Array<amrex::Real> bc_lo_3, bc_hi_3;
#endif

  amrex::Array<amrex::Geometry> &geom;
  int lev=0;
};

#endif
#endif
