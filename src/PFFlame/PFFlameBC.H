#ifndef PFFlameBC_H_
#define PFFlameBC_H_

#include <AMReX_PhysBCFunct.H>
#include "GeneralAMRIntegrator/GeneralAMRIntegratorBC.H"

class PFFlamePhysBC : public GeneralAMRIntegratorPhysBC
{
public:
  PFFlamePhysBC (amrex::Array<amrex::Geometry> &_geom, std::string name)
    : GeneralAMRIntegratorPhysBC(_geom) 
  {
    {
      amrex::ParmParse pp(name);
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
      
      if (bc_lo[0] == EXT_DIR) pp.getarr("lo_1",bc_lo_1); if(bc_hi[0] == EXT_DIR) pp.getarr("hi_1",bc_hi_1);
      if (bc_lo[1] == EXT_DIR) pp.getarr("lo_2",bc_lo_2); if(bc_hi[1] == EXT_DIR) pp.getarr("hi_2",bc_hi_2);
#if BL_SPACEDIM > 2
      if (bc_lo[2] == EXT_DIR) pp.getarr("lo_3",bc_lo_3); if(bc_hi[2] == EXT_DIR) pp.getarr("hi_3",bc_hi_3);
#endif
    }
  }


  virtual ~PFFlamePhysBC () {}

  virtual void FillBoundary (amrex::MultiFab& mf, int, int, amrex::Real /*time*/) override
  {
    amrex::Box domain(geom[lev].Domain());

    mf.FillBoundary(geom[lev].periodicity());

    for (amrex::MFIter mfi(mf,true); mfi.isValid(); ++mfi)
      {
   	const amrex::Box& box = mfi.tilebox();

   	amrex::BaseFab<amrex::Real> &mf_box = mf[mfi];

   	for (int i = box.loVect()[0] - mf.nGrow(); i<=box.hiVect()[0] + mf.nGrow(); i++)
   	  for (int j = box.loVect()[1] - mf.nGrow(); j<=box.hiVect()[1] + mf.nGrow(); j++)
#if BL_SPACEDIM > 2
	    for (int k = box.loVect()[2] - mf.nGrow(); k<=box.hiVect()[2] + mf.nGrow(); k++)
#endif
	      {
		if (bc_lo[0] == EXT_DIR)
		  if (i < domain.loVect()[0]) // Low 1 boundary
		    for (int n = 0 ; n < mf.nComp() ; n++) 
		      mf_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) =  bc_lo_1[n];

		if (bc_hi[0] == EXT_DIR)
		  if (i > domain.hiVect()[0]) // High 1 boundary
		    for (int n = 0 ; n < mf.nComp() ; n++) 
		      mf_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) =  bc_hi_1[n];

		if (bc_lo[1] == EXT_DIR)
		  if (j < domain.loVect()[1]) // Low 2 boundary
		    for (int n = 0 ; n < mf.nComp() ; n++) 
		      mf_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) =  bc_lo_2[n];

		if (bc_hi[1] == EXT_DIR)
		  if (j > domain.hiVect()[1]) // High 2 boundary
		    for (int n = 0 ; n < mf.nComp() ; n++) 
		      mf_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) =  bc_hi_2[n];

#if BL_SPACEDIM > 2
		if (bc_lo[2] == EXT_DIR)
		  if (k < domain.loVect()[2]) // Low 2 boundary
		    for (int n = 0 ; n < mf.nComp() ; n++) 
		      mf_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) =  bc_lo_3[n];

		if (bc_hi[2] == EXT_DIR)
		  if (k > domain.hiVect()[2]) // High 2 boundary
		    for (int n = 0 ; n < mf.nComp() ; n++) 
		      mf_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) =  bc_hi_3[n];
#endif
	      }
      }
  }


private:
  amrex::Array<amrex::Real> bc_lo_1, bc_hi_1;
  amrex::Array<amrex::Real> bc_lo_2, bc_hi_2;
#if BL_SPACEDIM > 2
  amrex::Array<amrex::Real> bc_lo_3, bc_hi_3;
#endif
  

};

#endif
