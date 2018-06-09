#include <AMReX_MultiFabUtil.H>
#include "eigen3/Eigen/Core"
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>
#include "Set/Set.H"

#include <AMReX_ArrayLim.H>

#include "Elastic.H"

/// \fn Operator::Elastic::MLPFStiffnessMatrix
///
/// Relay to the define function
/// Also define elastic constants here.
Operator::Elastic::Elastic::Elastic (const Vector<Geometry>& a_geom,
									 const Vector<BoxArray>& a_grids,
									 const Vector<DistributionMapping>& a_dmap,
									 const LPInfo& a_info)
{
	define(a_geom, a_grids, a_dmap, a_info);
}

Operator::Elastic::Elastic::~Elastic ()
{}

void
Operator::Elastic::Elastic::Fapply (int amrlev, ///<[in] AMR Level
									int mglev,  ///<[in]
									MultiFab& f,///<[out] The force vector
									const MultiFab& u ///<[in] The displacements vector
									) const
{
	BL_PROFILE("Operator::Elastic::Elastic::Fapply()");

	const Real* DX = m_geom[amrlev][mglev].CellSize();
  
#if AMREX_SPACEDIM == 1
	static amrex::IntVect dx(1);
#elif AMREX_SPACEDIM == 2
	static amrex::IntVect dx(1,0), dy(0,1);
#elif AMREX_SPACEDIM == 3
	static amrex::IntVect dx(1,0,0), dy(0,1,0), dz(0,0,1);
#endif 

	for (MFIter mfi(f, true); mfi.isValid(); ++mfi)
		{
			const Box& bx = mfi.tilebox();
			const amrex::BaseFab<amrex::Real> &ufab  = u[mfi];
			amrex::BaseFab<amrex::Real>       &ffab  = f[mfi];

			for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++)
				for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++)
#if AMREX_SPACEDIM > 2
					for (int m3 = bx.loVect()[2]; m3<=bx.hiVect()[2]; m3++)
#endif
						{
							amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));
							for (int i=0; i<AMREX_SPACEDIM; i++)
								{
									ffab(amrex::IntVect(AMREX_D_DECL(m1,m2,m3)),i) = 0.0;
									for (int k=0; k<AMREX_SPACEDIM; k++)
										{
											// C_{ijkl} u_{k,jl}
						
											ffab(m,i) -=
												C(i,0,k,0,m,amrlev,mglev,mfi) * (ufab(m+dx,k) - 2.0*ufab(m,k) + ufab(m-dx,k))/DX[0]/DX[0]
												+
												(C(i,0,k,1,m,amrlev,mglev,mfi) + C(i,1,k,0,m,amrlev,mglev,mfi)) * (ufab(m+dx+dy,k) + ufab(m-dx-dy,k) - ufab(m+dx-dy,k) - ufab(m-dx+dy,k))/(2.0*DX[0])/(2.0*DX[1])
												+
												C(i,1,k,1,m,amrlev,mglev,mfi) *(ufab(m+dy,k) - 2.0*ufab(m,k) + ufab(m-dy,k))/DX[1]/DX[1]
#if AMREX_SPACEDIM > 2
												+
												(C(i,0,k,2,m,amrlev,mglev,mfi) + C(i,2,k,0,m,amrlev,mglev,mfi)) * (ufab(m+dx+dz,k) + ufab(m-dx-dz,k) - ufab(m+dx-dz,k) - ufab(m-dx+dz,k))/(2.0*DX[0])/(2.0*DX[2])
												+
												(C(i,1,k,2,m,amrlev,mglev,mfi) + C(i,2,k,1,m,amrlev,mglev,mfi)) * (ufab(m+dy+dz,k) + ufab(m-dy-dz,k) - ufab(m+dy-dz,k) - ufab(m-dy+dz,k))/(2.0*DX[1])/(2.0*DX[2])
												+
												C(i,2,k,2,m,amrlev,mglev,mfi) *(ufab(m+dz,k) - 2.0*ufab(m,k) + ufab(m-dz,k))/DX[2]/DX[2]
#endif
												;

											// C_{ijkl,j} u_{k,l}
#if AMREX_SPACEDIM == 2
											ffab(m,i) -=
												((C(i,0,k,0,m+dx,amrlev,mglev,mfi) - C(i,0,k,0,m-dx,amrlev,mglev,mfi))/(2.0*DX[0]) +
												 (C(i,1,k,0,m+dy,amrlev,mglev,mfi) - C(i,1,k,0,m-dy,amrlev,mglev,mfi))/(2.0*DX[1])) *
												((ufab(m+dx,k) - ufab(m-dx,k))/(2.0*DX[0]))
												+
												((C(i,0,k,1,m+dx,amrlev,mglev,mfi) - C(i,0,k,1,m-dx,amrlev,mglev,mfi))/(2.0*DX[0]) +
												 (C(i,1,k,1,m+dy,amrlev,mglev,mfi) - C(i,1,k,1,m-dy,amrlev,mglev,mfi))/(2.0*DX[1])) *
												((ufab(m+dy,k) - ufab(m-dy,k))/(2.0*DX[1]));
#elif AMREX_SPACEIM == 3
											ffab(m,i) -=
												((C(i,0,k,0,m+dx,amrlev,mglev,mfi) - C(i,0,k,0,m-dx,amrlev,mglev,mfi))/(2.0*DX[0]) +
												 (C(i,1,k,0,m+dy,amrlev,mglev,mfi) - C(i,1,k,0,m-dy,amrlev,mglev,mfi))/(2.0*DX[1]) +
												 (C(i,2,k,0,m+dz,amrlev,mglev,mfi) - C(i,2,k,0,m-dz,amrlev,mglev,mfi))/(2.0*DX[2])) *
												((ufab(m+dx,k) - ufab(m-dx,k))/(2.0*DX[0]))
												+
												((C(i,0,k,1,m+dx,amrlev,mglev,mfi) - C(i,0,k,1,m-dx,amrlev,mglev,mfi))/(2.0*DX[0]) +
												 (C(i,1,k,1,m+dy,amrlev,mglev,mfi) - C(i,1,k,1,m-dy,amrlev,mglev,mfi))/(2.0*DX[1]) +
												 (C(i,2,k,1,m+dz,amrlev,mglev,mfi) - C(i,2,k,1,m-dz,amrlev,mglev,mfi))/(2.0*DX[2])) *
												((ufab(m+dy,k) - ufab(m-dy,k))/(2.0*DX[1]))
												+
												((C(i,0,k,2,m+dx,amrlev,mglev,mfi) - C(i,0,k,2,m-dx,amrlev,mglev,mfi))/(2.0*DX[0]) +
												 (C(i,1,k,2,m+dy,amrlev,mglev,mfi) - C(i,1,k,2,m-dy,amrlev,mglev,mfi))/(2.0*DX[1]) +
												 (C(i,2,k,2,m+dz,amrlev,mglev,mfi) - C(i,2,k,2,m-dz,amrlev,mglev,mfi))/(2.0*DX[2])) *
												((ufab(m+dz,k) - ufab(m-dz,k))/(2.0*DX[2]));
#endif
										}
								}
						}
		}
}


/// \fn Operator::Elastic::Fsmooth
///
/// Perform one half Gauss-Seidel iteration corresponding to the operator specified
/// in Operator::Elastic::Fapply.
/// The variable redblack corresponds to whether to smooth "red" nodes or "black"
/// nodes, where red and black nodes are distributed in a checkerboard pattern.
///
/// \todo Extend to 3D
///
void
Operator::Elastic::Elastic::Fsmooth (int amrlev,          ///<[in] AMR level
									 int mglev,           ///<[in]
									 MultiFab& u,       ///<[inout] Solution (displacement field)
									 const MultiFab& rhs, ///<[in] Body force vectors (rhs=right hand side)
									 int redblack         ///<[in] Smooth even vs. odd modes
									 ) const
{
	BL_PROFILE("Operator::Elastic::Elastic::Fsmooth()");

	const Real* DX = m_geom[amrlev][mglev].CellSize();

#if AMREX_SPACEDIM == 1
	static amrex::IntVect dx(1);
#elif AMREX_SPACEDIM == 2
	static amrex::IntVect dx(1,0), dy(0,1);
#elif AMREX_SPACEDIM == 3
	static amrex::IntVect dx(1,0,0), dy(0,1,0), dz(0,0,1);
#endif 

	for (MFIter mfi(u,MFItInfo().EnableTiling().SetDynamic(true));
		 mfi.isValid(); ++mfi)
		{
			const Box&       bx     = mfi.tilebox();
			FArrayBox&       ufab    = u[mfi];
			const FArrayBox& rhsfab  = rhs[mfi];

			for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++)
				for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++)
#if AMREX_SPACEDIM > 2
					for (int m3 = bx.loVect()[2]; m3<=bx.hiVect()[2]; m3++)
#endif
						{
							if ((AMREX_D_TERM(m1, + m2, + m3))%2 == redblack) continue;
							amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));

							for (int i=0; i<AMREX_SPACEDIM; i++)
								{
									amrex::Real rho = 0.0, aa = 0.0;
									for (int k=0; k<AMREX_SPACEDIM; k++)
										{
											// C_{ijkl} u_{k,jl}
											rho -=
												C(i,0,k,0,m,amrlev,mglev,mfi)
												*(ufab(m+dx,k) - (i==k ? 0.0 : 2.0*ufab(m,k)) + ufab(m-dx,k))/DX[0]/DX[0]
												+
												(C(i,0,k,1,m,amrlev,mglev,mfi) + C(i,1,k,0,m,amrlev,mglev,mfi))
												*(ufab(m+dx+dy,k) + ufab(m-dx-dy,k) - ufab(m+dx-dy,k) - ufab(m-dx+dy,k))/(2.0*DX[0])/(2.0*DX[1])
												+
												C(i,1,k,1,m,amrlev,mglev,mfi)
												*(ufab(m+dy,k) - (i==k ? 0.0 : 2.0*ufab(m,k)) + ufab(m-dy,k))/DX[1]/DX[1]
												;
#if AMREX_SPACEDIM == 3
											rho -=
												(C(i,0,k,2,m,amrlev,mglev,mfi) + C(i,2,k,0,m,amrlev,mglev,mfi))
												*(ufab(m+dx+dz,k) + ufab(m-dx-dz,k) - ufab(m+dx-dz,k) - ufab(m-dx+dz,k))/(2.0*DX[0])/(2.0*DX[2])
												+
												(C(i,1,k,2,m,amrlev,mglev,mfi) + C(i,2,k,1,m,amrlev,mglev,mfi))
												*(ufab(m+dy+dz,k) + ufab(m-dy-dz,k) - ufab(m+dy-dz,k) - ufab(m-dy+dz,k))/(2.0*DX[1])/(2.0*DX[2])
												+
												C(i,2,k,2,m,amrlev,mglev,mfi)
												*(ufab(m+dz,k) - (i==k ? 0.0 : 2.0*ufab(m,k)) + ufab(m-dz,k))/DX[2]/DX[2];
#endif
											


											// C_{ijkl,j} u_{k,l}
#if AMREX_SPACEDIM == 2
											rho -=
												((C(i,0,k,0,m+dx,amrlev,mglev,mfi) - C(i,0,k,0,m-dx,amrlev,mglev,mfi))/(2.0*DX[0]) +
												 (C(i,1,k,0,m+dy,amrlev,mglev,mfi) - C(i,1,k,0,m-dy,amrlev,mglev,mfi))/(2.0*DX[1])) *
												((ufab(m+dx,k) - ufab(m-dx,k))/(2.0*DX[0]))
												+
												((C(i,0,k,1,m+dx,amrlev,mglev,mfi) - C(i,0,k,1,m-dx,amrlev,mglev,mfi))/(2.0*DX[0]) +
												 (C(i,1,k,1,m+dy,amrlev,mglev,mfi) - C(i,1,k,1,m-dy,amrlev,mglev,mfi))/(2.0*DX[1])) *
												((ufab(m+dy,k) - ufab(m-dy,k))/(2.0*DX[1]));
#elif AMREX_SPACEDIM == 3
											rho -=
												((C(i,0,k,0,m+dx,amrlev,mglev,mfi) - C(i,0,k,0,m-dx,amrlev,mglev,mfi))/(2.0*DX[0]) +
												 (C(i,1,k,0,m+dy,amrlev,mglev,mfi) - C(i,1,k,0,m-dy,amrlev,mglev,mfi))/(2.0*DX[1]) +
												 (C(i,2,k,0,m+dz,amrlev,mglev,mfi) - C(i,2,k,0,m-dz,amrlev,mglev,mfi))/(2.0*DX[2])) *
												((ufab(m+dx,k) - ufab(m-dx,k))/(2.0*DX[0]))
												+
												((C(i,0,k,1,m+dx,amrlev,mglev,mfi) - C(i,0,k,1,m-dx,amrlev,mglev,mfi))/(2.0*DX[0]) +
												 (C(i,1,k,1,m+dy,amrlev,mglev,mfi) - C(i,1,k,1,m-dy,amrlev,mglev,mfi))/(2.0*DX[1]) +
												 (C(i,2,k,1,m+dz,amrlev,mglev,mfi) - C(i,2,k,1,m-dz,amrlev,mglev,mfi))/(2.0*DX[2])) *
												((ufab(m+dy,k) - ufab(m-dy,k))/(2.0*DX[1]))
												+
												((C(i,0,k,2,m+dx,amrlev,mglev,mfi) - C(i,0,k,2,m-dx,amrlev,mglev,mfi))/(2.0*DX[0]) +
												 (C(i,1,k,2,m+dy,amrlev,mglev,mfi) - C(i,1,k,2,m-dy,amrlev,mglev,mfi))/(2.0*DX[1]) +
												 (C(i,2,k,2,m+dz,amrlev,mglev,mfi) - C(i,2,k,2,m-dz,amrlev,mglev,mfi))/(2.0*DX[2])) *
												((ufab(m+dz,k) - ufab(m-dz,k))/(2.0*DX[2]));
#endif

										}

									aa -=
										-2.0*C(i,0,i,0,m,amrlev,mglev,mfi)/DX[0]/DX[0]
										-2.0*C(i,1,i,1,m,amrlev,mglev,mfi)/DX[1]/DX[1]
#if AMREX_SPACEDIM > 2
										-2.0*C(i,2,i,2,m,amrlev,mglev,mfi)/DX[2]/DX[2]
#endif
										;

									//std::cout << "nans not detetected, rho=" << rho << ", aa=" << aa << std::endl;
									if (rho != rho) std::cout << "nans detetected, rho=" << rho << ", aa=" << aa << std::endl;
									if (rho != rho) amrex::Abort("nans detected");

									ufab(m,i) = (rhsfab(m,i) - rho) / aa;
								}
						}
		}
}

/// \fn Operator::Elastic::FFlux
///
/// Compute the "flux" corresponding to the operator in Operator::Elastic::Fapply.
/// Because the operator is self-adjoint and positive-definite, the flux is not
/// required for adequate convergence (?)
/// Therefore, the fluxes are simply set to zero and returned.
///
/// \todo Extend to 3D
///
void
Operator::Elastic::Elastic::FFlux (int /*amrlev*/, const MFIter& /*mfi*/,
								   const std::array<FArrayBox*,AMREX_SPACEDIM>& sigmafab,
								   const FArrayBox& /*ufab*/, const int /*face_only*/) const
{

	amrex::BaseFab<amrex::Real> &fxfab = *sigmafab[0];
	amrex::BaseFab<amrex::Real> &fyfab = *sigmafab[1];
	fxfab.setVal(0.0);
	fyfab.setVal(0.0);
}


void
Operator::Elastic::Elastic::Stress (FArrayBox& sigmafab,
									const FArrayBox& ufab,
									int amrlev, const MFIter& mfi) const
{
	/// \todo add assert for sigmafab.ncomp=3 (for SPACEDIM=2) and ncomp=6 (for SPACEDIM=3)

	sigmafab.setVal(0.0);

	const amrex::Real* DX = m_geom[amrlev][0].CellSize();

#if AMREX_SPACEDIM==2
	amrex::IntVect dx(1,0);
	amrex::IntVect dy(0,1);
#elif AMREX_SPACEDIM==3
	amrex::IntVect dx(1,0,0);
	amrex::IntVect dy(0,1,0);
	amrex::IntVect dz(0,0,1);
#endif
	const Box& bx = mfi.tilebox();
	for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++)
		for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++)
#if AMREX_SPACEDIM==3
			for(int m3 = bx.loVect()[2]; m3<=bx.hiVect()[2]; m3++)
#endif
		{
				amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));

	 			amrex::Real du1_dx1 = (ufab(m+dx,0) - ufab(m-dx,0))/(2.0*DX[0]);
	 			amrex::Real du1_dx2 = (ufab(m+dy,0) - ufab(m-dy,0))/(2.0*DX[1]);
#if AMREX_SPACEDIM==3
				amrex::Real du1_dx3 = (ufab(m+dy,0) - ufab(m-dy,0))/(2.0*DX[2]);
#endif
	 			amrex::Real du2_dx1 = (ufab(m+dx,1) - ufab(m-dx,1))/(2.0*DX[0]);
	 			amrex::Real du2_dx2 = (ufab(m+dy,1) - ufab(m-dy,1))/(2.0*DX[1]);
#if AMREX_SPACEDIM==3
				amrex::Real du2_dx3 = (ufab(m+dz,1) - ufab(m-dz,1))/(2.0*DX[2]);
				amrex::Real du3_dx1 = (ufab(m+dx,2) - ufab(m-dx,2))/(2.0*DX[0]);
				amrex::Real du3_dx2 = (ufab(m+dy,2) - ufab(m-dy,2))/(2.0*DX[1]);
				amrex::Real du3_dx3 = (ufab(m+dz,2) - ufab(m-dz,2))/(2.0*DX[2]);
#endif

	 			for (int i=0; i<AMREX_SPACEDIM; i++)
	 				for (int j=0; j<AMREX_SPACEDIM; j++)
					{
						int voigt = (i+1)*(i==j ? 1:0) + (1- (i==j ? 1:0))*
							(AMREX_SPACEDIM*AMREX_SPACEDIM - (i+1) - (j+1)) - 1;
						sigmafab(m,voigt) =
							C(i,j,0,0,m,amrlev,0,mfi)*du1_dx1 +
							C(i,j,0,1,m,amrlev,0,mfi)*du1_dx2 +
							C(i,j,1,0,m,amrlev,0,mfi)*du2_dx1 +
							C(i,j,1,1,m,amrlev,0,mfi)*du2_dx2;
#if AMREX_SPACEDIM == 3
						sigmafab(m,voigt) += 
							C(i,j,0,2,m,amrlev,0,mfi)*du1_dx3 + 
							C(i,j,1,2,m,amrlev,0,mfi)*du2_dx3 + 
							C(i,j,2,0,m,amrlev,0,mfi)*du3_dx1 + 
							C(i,j,2,1,m,amrlev,0,mfi)*du3_dx2 + 
							C(i,j,2,2,m,amrlev,0,mfi)*du3_dx3;
							
#endif
					}
		}

}

void
Operator::Elastic::Elastic::Energy (FArrayBox& energyfab,
									const FArrayBox& ufab,
									int amrlev, const MFIter& mfi) const
{
	/// \todo RE-IMPLEMENT in 2D and 3D

	const amrex::Real* DX = m_geom[amrlev][0].CellSize();

	energyfab.setVal(0.0);
	
#if AMREX_SPACEDIM==2
	amrex::IntVect dx(1,0);
	amrex::IntVect dy(0,1);
#elif AMREX_SPACEDIM==3
	amrex::IntVect dx(1,0,0);
	amrex::IntVect dy(0,1,0);
	amrex::IntVect dz(0,0,1);
#endif

	const Box& bx = mfi.tilebox();
	for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++)
#if AMREX_SPACEDIM>1
		for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++)
#if AMREX_SPACEDIM==3
			for (int m3 = bx.loVect()[2]; m3<=bx.hiVect()[2]; m3++)
#endif
#endif
		{
	 			amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));
				Set::Matrix gradu;
#if AMREX_SPACEDIM==2
				gradu <<
					(ufab(m+dx,0) - ufab(m-dx,0))/(2.0*DX[0]), (ufab(m+dy,0) - ufab(m-dy,0))/(2.0*DX[1]),
					(ufab(m+dx,1) - ufab(m-dx,1))/(2.0*DX[0]), (ufab(m+dy,1) - ufab(m-dy,1))/(2.0*DX[1]);
#elif AMREX_SPACEDIM==3
				gradu << 
					(ufab(m+dx,0) - ufab(m-dx,0))/(2.0*DX[0]), (ufab(m+dy,0) - ufab(m-dy,0))/(2.0*DX[1]), (ufab(m+dz,0) - ufab(m-dz,0))/(2.0*DX[2]),
					(ufab(m+dx,1) - ufab(m-dx,1))/(2.0*DX[0]), (ufab(m+dy,1) - ufab(m-dy,1))/(2.0*DX[1]), (ufab(m+dz,1) - ufab(m-dz,1))/(2.0*DX[2]),
					(ufab(m+dx,2) - ufab(m-dx,2))/(2.0*DX[0]), (ufab(m+dy,2) - ufab(m-dy,2))/(2.0*DX[1]), (ufab(m+dz,2) - ufab(m-dz,2))/(2.0*DX[2]);
#endif

	 			energyfab(m) = 0.0;

	 			for (int i=0; i<AMREX_SPACEDIM; i++)
	 				for (int j=0; j<AMREX_SPACEDIM; j++)
	 					for (int k=0; k<AMREX_SPACEDIM; k++)
	 						for (int l=0; l<AMREX_SPACEDIM; l++)
								energyfab(m) += gradu(i,j) * C(i,j,k,l,m,amrlev,0,mfi) * gradu(k,l);
				
	 			energyfab(m) *= 0.5;
		}

}
