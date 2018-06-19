#include "BC.H"
#include "BC/Util.H"

namespace BC
{
BC::BC (amrex::Vector<amrex::Geometry> &_geom
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
	: geom(_geom)
	,bc_lo_1(_bc_lo_1), bc_hi_1(_bc_hi_1)
	,bc_lo_2(_bc_lo_2), bc_hi_2(_bc_hi_2)
#if AMREX_SPACEDIM > 2
	,bc_lo_3(_bc_lo_3), bc_hi_3(_bc_hi_3)
#endif
{
	for (int i=0;i<BL_SPACEDIM;i++)
		{
			bc_hi[i] = Util::ReadString(bc_hi_str[i]);
			bc_lo[i] = Util::ReadString(bc_lo_str[i]);

			if (Util::IsPeriodic(bc_lo[i]) != Util::IsPeriodic(bc_hi[i]))
				amrex::Abort("Invalid BCs cannot be periodic on one side and not the other");
		}
}

void
BC::FillBoundary (amrex::MultiFab& mf, int, int, amrex::Real /*time*/) 
{
	if ((Util::IsNeumann(bc_lo[0]) || Util::IsDirichlet(bc_lo[0])) && bc_lo_1.size() < mf.nComp()) amrex::Abort("Not enough values specified for bc_lo_1");
	if ((Util::IsNeumann(bc_hi[0]) || Util::IsDirichlet(bc_hi[0])) && bc_hi_1.size() < mf.nComp()) amrex::Abort("Not enough values specified for bc_hi_1");
	if ((Util::IsNeumann(bc_lo[1]) || Util::IsDirichlet(bc_lo[1])) && bc_lo_2.size() < mf.nComp()) amrex::Abort("Not enough values specified for bc_lo_2");
	if ((Util::IsNeumann(bc_hi[1]) || Util::IsDirichlet(bc_hi[1])) && bc_hi_2.size() < mf.nComp()) amrex::Abort("Not enough values specified for bc_hi_2");
#if AMREX_SPACEDIM>2
	if ((Util::IsNeumann(bc_lo[2]) || Util::IsDirichlet(bc_lo[2])) && bc_lo_3.size() < mf.nComp()) amrex::Abort("Not enough values specified for bc_lo_3");
	if ((Util::IsNeumann(bc_hi[2]) || Util::IsDirichlet(bc_hi[2])) && bc_hi_3.size() < mf.nComp()) amrex::Abort("Not enough values specified for bc_hi_3");
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
										if (Util::IsDirichlet(bc_lo[0]))
											{
												phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = bc_lo_1[n];
											}
										else if(Util::IsNeumann(bc_lo[0]))
											{
												phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = phi_box(amrex::IntVect(AMREX_D_DECL(i+1,j,k)),n) - bc_lo_1[n]*dx[0];
											}
										else if(Util::IsPeriodic(bc_lo[0]))
											{
												phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = 
													phi_box(amrex::IntVect(AMREX_D_DECL(-i+box.loVect()[0]+box.hiVect()[0],j,k)),n);
											}
									}

								if (i > domain.hiVect()[0]) // Right boundary
									{
										if (Util::IsDirichlet(bc_hi[0]))
											{
												phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = bc_hi_1[n];
											}
										else if(Util::IsNeumann(bc_hi[0]))
											{
												phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = phi_box(amrex::IntVect(AMREX_D_DECL(i-1,j,k)),n) - bc_hi_1[n]*dx[0];
											}
										else if(Util::IsPeriodic(bc_hi[0]))
											{
												phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = 
													phi_box(amrex::IntVect(AMREX_D_DECL(-i+box.loVect()[0]+box.hiVect()[0],j,k)),n);
											}
									}

								if (j < domain.loVect()[1]) // Bottom boundary
									{
										if (Util::IsDirichlet(bc_lo[1]))
											{
												phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = bc_lo_2[n];
											}
										else if (Util::IsNeumann(bc_lo[1]))
											{
												phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = phi_box(amrex::IntVect(AMREX_D_DECL(i,j+1,k)),n) - bc_lo_2[n]*dx[1];
											}
										else if(Util::IsPeriodic(bc_lo[1]))
											{
												phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = 
													phi_box(amrex::IntVect(AMREX_D_DECL(i,-j+box.loVect()[1]+box.hiVect()[1],k)),n);
											}
									}

								if (j > domain.hiVect()[1]) // Top boundary
									{
										if (Util::IsDirichlet(bc_hi[1]))
											{
												phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = bc_hi_2[n];
											}
										else if (Util::IsNeumann(bc_hi[1]))
											{
												phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = phi_box(amrex::IntVect(AMREX_D_DECL(i,j-1,k)),n) - bc_hi_2[n]*dx[1];
											}
										else if(Util::IsPeriodic(bc_hi[1]))
											{
												phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = 
													phi_box(amrex::IntVect(AMREX_D_DECL(i,-j+box.loVect()[1]+box.hiVect()[1],k)),n);
											}
									}


#if AMREX_SPACEDIM>2
								if (k < domain.loVect()[2])
									{
										if (Util::IsDirichlet(bc_lo[2]))
											{
												phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = bc_lo_3[n];
											}
										else if (Util::IsNeumann(bc_lo[2]))
											{
												phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k+1)),n) - bc_lo_3[n]*dx[2];
											}
										else if(Util::IsPeriodic(bc_lo[2]))
											{
												phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = 
													phi_box(amrex::IntVect(AMREX_D_DECL(i,j,-k+box.loVect()[2]+box.hiVect()[2])),n);
											}
									}

								if (k > domain.hiVect()[2])
									{
										if (Util::IsDirichlet(bc_hi[2]))
											{
												phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = bc_hi_3[n];
											}
										else if(Util::IsNeumann(bc_hi[2]))
											{
												phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k-1)),n) - bc_hi_3[n]*dx[2];
											}
										else if(Util::IsPeriodic(bc_hi[2]))
											{
												phi_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n) = 
													phi_box(amrex::IntVect(AMREX_D_DECL(i,j,-k+box.loVect()[2]+box.hiVect()[2])),n);
											}
									}
#endif

							}
		}
}

void
BC::SetLevel(int _lev) {lev=_lev;}

amrex::BCRec
BC::GetBCRec() {return amrex::BCRec(bc_lo,bc_hi);}

amrex::Array<int,AMREX_SPACEDIM>
BC::IsPeriodic()
{
	return {AMREX_D_DECL(Util::IsPeriodic(bc_lo[0]),Util::IsPeriodic(bc_lo[1]),Util::IsPeriodic(bc_lo[2]))};
}


}
