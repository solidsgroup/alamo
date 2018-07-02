#include "Constant.H"

namespace BC
{
Constant::Constant (amrex::Vector<amrex::Geometry> & _geom
	  ,amrex::Vector<std::string> bc_hi_str
	  ,amrex::Vector<std::string> bc_lo_str
	  ,amrex::Vector<amrex::Real> _bc_lo_1
	  ,amrex::Vector<amrex::Real> _bc_hi_1
	  ,amrex::Vector<amrex::Real> _bc_lo_2
	  ,amrex::Vector<amrex::Real> _bc_hi_2
#if AMREX_SPACEDIM > 2
	  ,amrex::Vector<amrex::Real> _bc_lo_3
	  ,amrex::Vector<amrex::Real> _bc_hi_3
#endif
	  )
	: BC(_geom)
	,bc_lo_1(_bc_lo_1), bc_hi_1(_bc_hi_1)
	,bc_lo_2(_bc_lo_2), bc_hi_2(_bc_hi_2)
#if AMREX_SPACEDIM > 2
	,bc_lo_3(_bc_lo_3), bc_hi_3(_bc_hi_3)
#endif
{
	for (int i=0;i<BL_SPACEDIM;i++)
	{
		bc_hi[i] = BCUtil::ReadString(bc_hi_str[i]);
		bc_lo[i] = BCUtil::ReadString(bc_lo_str[i]);

		if (BCUtil::IsPeriodic(bc_lo[i]) != BCUtil::IsPeriodic(bc_hi[i]))
			Util::Abort("Invalid BCs cannot be periodic on one side and not the other");
	}
}

void
Constant::FillBoundary (amrex::MultiFab& mf, int, int, amrex::Real /*time*/) 
{
	if ((BCUtil::IsNeumann(bc_lo[0]) || BCUtil::IsDirichlet(bc_lo[0])) && bc_lo_1.size() < mf.nComp()) Util::Abort("Not enough values specified for bc_lo_1");
	if ((BCUtil::IsNeumann(bc_hi[0]) || BCUtil::IsDirichlet(bc_hi[0])) && bc_hi_1.size() < mf.nComp()) Util::Abort("Not enough values specified for bc_hi_1");
	if ((BCUtil::IsNeumann(bc_lo[1]) || BCUtil::IsDirichlet(bc_lo[1])) && bc_lo_2.size() < mf.nComp()) Util::Abort("Not enough values specified for bc_lo_2");
	if ((BCUtil::IsNeumann(bc_hi[1]) || BCUtil::IsDirichlet(bc_hi[1])) && bc_hi_2.size() < mf.nComp()) Util::Abort("Not enough values specified for bc_hi_2");
#if AMREX_SPACEDIM>2
	if ((BCUtil::IsNeumann(bc_lo[2]) || BCUtil::IsDirichlet(bc_lo[2])) && bc_lo_3.size() < mf.nComp()) Util::Abort("Not enough values specified for bc_lo_3");
	if ((BCUtil::IsNeumann(bc_hi[2]) || BCUtil::IsDirichlet(bc_hi[2])) && bc_hi_3.size() < mf.nComp()) Util::Abort("Not enough values specified for bc_hi_3");
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
							if (BCUtil::IsDirichlet(bc_lo[0]))
							{
								phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = bc_lo_1[n];
							}
							else if(BCUtil::IsNeumann(bc_lo[0]))
							{
								phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = phi_box(amrex::IntVect(AMREX_D_DECL(i+1,j,k)),n) - bc_lo_1[n]*dx[0];
							}
							else if(BCUtil::IsPeriodic(bc_lo[0]))
							{
								phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = 
									phi_box(amrex::IntVect(AMREX_D_DECL(-i+box.loVect()[0]+box.hiVect()[0],j,k)),n);
							}
						}

						if (i > domain.hiVect()[0]) // Right boundary
						{
							if (BCUtil::IsDirichlet(bc_hi[0]))
							{
								phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = bc_hi_1[n];
							}
							else if(BCUtil::IsNeumann(bc_hi[0]))
							{
								phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = phi_box(amrex::IntVect(AMREX_D_DECL(i-1,j,k)),n) - bc_hi_1[n]*dx[0];
							}
							else if(BCUtil::IsPeriodic(bc_hi[0]))
							{
								phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = 
									phi_box(amrex::IntVect(AMREX_D_DECL(-i+box.loVect()[0]+box.hiVect()[0],j,k)),n);
							}
						}

						if (j < domain.loVect()[1]) // Bottom boundary
						{
							if (BCUtil::IsDirichlet(bc_lo[1]))
							{
								phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = bc_lo_2[n];
							}
							else if (BCUtil::IsNeumann(bc_lo[1]))
							{
								phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = phi_box(amrex::IntVect(AMREX_D_DECL(i,j+1,k)),n) - bc_lo_2[n]*dx[1];
							}
							else if(BCUtil::IsPeriodic(bc_lo[1]))
							{
								phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = 
									phi_box(amrex::IntVect(AMREX_D_DECL(i,-j+box.loVect()[1]+box.hiVect()[1],k)),n);
							}
						}

						if (j > domain.hiVect()[1]) // Top boundary
						{
							if (BCUtil::IsDirichlet(bc_hi[1]))
							{
								phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = bc_hi_2[n];
							}
							else if (BCUtil::IsNeumann(bc_hi[1]))
							{
								phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = phi_box(amrex::IntVect(AMREX_D_DECL(i,j-1,k)),n) - bc_hi_2[n]*dx[1];
							}
							else if(BCUtil::IsPeriodic(bc_hi[1]))
							{
								phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = 
									phi_box(amrex::IntVect(AMREX_D_DECL(i,-j+box.loVect()[1]+box.hiVect()[1],k)),n);
							}
						}


#if AMREX_SPACEDIM>2
						if (k < domain.loVect()[2])
						{
							if (BCUtil::IsDirichlet(bc_lo[2]))
							{
								phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = bc_lo_3[n];
							}
							else if (BCUtil::IsNeumann(bc_lo[2]))
							{
								phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k+1)),n) - bc_lo_3[n]*dx[2];
							}
							else if(BCUtil::IsPeriodic(bc_lo[2]))
							{
								phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = 
									phi_box(amrex::IntVect(AMREX_D_DECL(i,j,-k+box.loVect()[2]+box.hiVect()[2])),n);
							}
						}

						if (k > domain.hiVect()[2])
						{
							if (BCUtil::IsDirichlet(bc_hi[2]))
							{
								phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = bc_hi_3[n];
							}
							else if(BCUtil::IsNeumann(bc_hi[2]))
							{
								phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k-1)),n) - bc_hi_3[n]*dx[2];
							}
							else if(BCUtil::IsPeriodic(bc_hi[2]))
							{
								phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = 
									phi_box(amrex::IntVect(AMREX_D_DECL(i,j,-k+box.loVect()[2]+box.hiVect()[2])),n);
							}
						}
#endif

					}
	}
}

amrex::BCRec
Constant::GetBCRec() {return amrex::BCRec(bc_lo,bc_hi);}

amrex::Array<int,AMREX_SPACEDIM>
Constant::IsPeriodic()
{
	return {AMREX_D_DECL(BCUtil::IsPeriodic(bc_lo[0]),BCUtil::IsPeriodic(bc_lo[1]),BCUtil::IsPeriodic(bc_lo[2]))};
}
}


