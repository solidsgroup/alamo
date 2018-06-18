#include "BC.H"

BC::BC::BC (amrex::Vector<amrex::Geometry> &_geom,
	    amrex::Vector<std::string> bc_hi_str,
	    amrex::Vector<std::string> bc_lo_str,
	    amrex::Vector<amrex::Real> _bc_lo_1,
	    amrex::Vector<amrex::Real> _bc_hi_1,
	    amrex::Vector<amrex::Real> _bc_lo_2,
	    amrex::Vector<amrex::Real> _bc_hi_2,
	    amrex::Vector<amrex::Real> _bc_lo_3,
	    amrex::Vector<amrex::Real> _bc_hi_3)
  : geom(_geom),
    bc_lo_1(_bc_lo_1), bc_hi_1(_bc_hi_1), 
    bc_lo_2(_bc_lo_2), bc_hi_2(_bc_hi_2), 
    bc_lo_3(_bc_lo_3), bc_hi_3(_bc_hi_3)
{
  for (int i=0;i<BL_SPACEDIM;i++)
    {
      if (bc_hi_str[i] == "REFLECT_ODD"	) bc_hi[i] = REFLECT_ODD; 
      if (bc_hi_str[i] == "INT_DIR"		) bc_hi[i] = INT_DIR;
      if (bc_hi_str[i] == "PERIODIC"		) bc_hi[i] = INT_DIR;
      if (bc_hi_str[i] == "REFLECT_EVEN"	) bc_hi[i] = REFLECT_EVEN;
      if (bc_hi_str[i] == "FOEXTRAP"		) bc_hi[i] = FOEXTRAP;
      if (bc_hi_str[i] == "EXT_DIR"		) bc_hi[i] = EXT_DIR;
      if (bc_hi_str[i] == "DIRICHLET"		) bc_hi[i] = EXT_DIR
      if (bc_hi_str[i] == "HOEXTRAP"		) bc_hi[i] = HOEXTRAP;

      if (bc_lo_str[i] == "REFLECT_ODD"	) bc_lo[i] = REFLECT_ODD;
      if (bc_lo_str[i] == "INT_DIR"		) bc_lo[i] = INT_DIR;
      if (bc_lo_str[i] == "PERIODIC"		) bc_lo[i] = INT_DIR;
      if (bc_lo_str[i] == "REFLECT_EVEN"	) bc_lo[i] = REFLECT_EVEN;
      if (bc_lo_str[i] == "FOEXTRAP"		) bc_lo[i] = FOEXTRAP;
      if (bc_lo_str[i] == "EXT_DIR"		) bc_lo[i] = EXT_DIR;
      if (bc_lo_str[i] == "DIRICHLET"		) bc_lo[i] = EXT_DIR
      if (bc_lo_str[i] == "HOEXTRAP"		) bc_lo[i] = HOEXTRAP;

      // Added for Neumann BC.
      if (bc_hi_str[i] == "NEUMANN"		) bc_hi[i] = NEUMANN;
      if (bc_lo_str[i] == "NEUMANN"		) bc_lo[i] = NEUMANN;
    }
}

void
BC::BC::FillBoundary (amrex::MultiFab& mf, int, int, amrex::Real /*time*/) 
{
  if ((bc_lo[0] == NEUMANN || bc_lo[0] == EXT_DIR) && bc_lo_1.size() < mf.nComp()) amrex::Abort("Not enough values specified for bc_lo_1");
  if ((bc_hi[0] == NEUMANN || bc_hi[0] == EXT_DIR) && bc_hi_1.size() < mf.nComp()) amrex::Abort("Not enough values specified for bc_hi_1");
  if ((bc_lo[1] == NEUMANN || bc_lo[1] == EXT_DIR) && bc_lo_2.size() < mf.nComp()) amrex::Abort("Not enough values specified for bc_lo_2");
  if ((bc_hi[1] == NEUMANN || bc_hi[1] == EXT_DIR) && bc_hi_2.size() < mf.nComp()) amrex::Abort("Not enough values specified for bc_hi_2");
#if AMREX_SPACEDIM>2
  if ((bc_lo[2] == NEUMANN || bc_lo[2] == EXT_DIR) && bc_lo_3.size() < mf.nComp()) amrex::Abort("Not enough values specified for bc_lo_3");
  if ((bc_hi[2] == NEUMANN || bc_hi[2] == EXT_DIR) && bc_hi_3.size() < mf.nComp()) amrex::Abort("Not enough values specified for bc_hi_3");
#endif

  amrex::Box domain(geom[lev].Domain());

  mf.FillBoundary(geom[lev].periodicity());

  // Added for Neumann BC
  const amrex::Real* dx = geom[lev].CellSize();

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
                    else if(bc_lo[0] == INT_DIR)
		      {
			phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = 
				phi_box(amrex::IntVect(AMREX_D_DECL(-i+box.loVecto()[0]+box.hiVect()[0],j,k)),n);
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
                    else if(bc_hi[0] == INT_DIR)
		      {
			phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = 
				phi_box(amrex::IntVect(AMREX_D_DECL(-i+box.loVecto()[0]+box.hiVect()[0],j,k)),n);
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
                    else if(bc_lo[1] == INT_DIR)
		      {
			phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = 
				phi_box(amrex::IntVect(AMREX_D_DECL(i,-j+box.loVecto()[1]+box.hiVect()[1],k)),n);
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
                    else if(bc_hi[1] == INT_DIR)
		      {
			phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = 
				phi_box(amrex::IntVect(AMREX_D_DECL(i,-j+box.loVecto()[1]+box.hiVect()[1],k)),n);
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
                    else if(bc_lo[2] == INT_DIR)
		      {
			phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = 
				phi_box(amrex::IntVect(AMREX_D_DECL(i,j,-k+box.loVecto()[2]+box.hiVect()[2])),n);
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
                    else if(bc_hi[2] == INT_DIR)
		      {
			phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = 
				phi_box(amrex::IntVect(AMREX_D_DECL(i,j,-k+box.loVecto()[2]+box.hiVect()[2])),n);
		      }
		  }
#endif

	      }
    }
}

void
BC::BC::SetLevel(int _lev) {lev=_lev;}

amrex::BCRec
BC::BC::GetBCRec() {return amrex::BCRec(bc_lo,bc_hi);}
