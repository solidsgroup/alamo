#include <AMReX_MultiFabUtil.H>
#include "eigen3/Eigen/Core"
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>

#include <AMReX_ArrayLim.H>

#include "FiniteElastic.H"

/// \fn Operator::FiniteElastic::MLPFStiffnessMatrix
///
/// Relay to the define function
/// Also define elastic constants here.
Operator::FiniteElastic::FiniteElastic::FiniteElastic (const Vector<Geometry>& a_geom,
					  const Vector<BoxArray>& a_grids,
					  const Vector<DistributionMapping>& a_dmap,
					  const LPInfo& a_info)
{
  define(a_geom, a_grids, a_dmap, a_info);
}

Operator::FiniteElastic::FiniteElastic::~FiniteElastic ()
{}

void
Operator::FiniteElastic::FiniteElastic::Fapply (int amrlev, ///<[in] AMR Level
			     int mglev,  ///<[in]
			     MultiFab& f,///<[out] The force vector
			     const MultiFab& u ///<[in] The displacements vector
			     ) const
{
  const Real* DX = m_geom[amrlev][mglev].CellSize();
  
  static amrex::IntVect dx(1,0), dy(0,1);
  for (MFIter mfi(f, true); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();
		const amrex::BaseFab<amrex::Real> &ufab  = u[mfi];
		amrex::BaseFab<amrex::Real>       &ffab  = f[mfi];

		for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++)
			for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++)
			{
				amrex::IntVect m(m1,m2);
				for (int i=0; i<AMREX_SPACEDIM; i++)
				{
					ffab(amrex::IntVect(m1,m2),i) = 0.0;
					for (int k=0; k<AMREX_SPACEDIM; k++)
					{
						// C_{ijkl} u_{k,jl}
						
						ffab(m,i) -=
							C(i,0,k,0,m,amrlev,mglev,mfi)
							*(ufab(m+dx,k) - 2.0*ufab(m,k) + ufab(m-dx,k))/DX[0]/DX[0]
							+
							(C(i,0,k,1,m,amrlev,mglev,mfi) + C(i,1,k,0,m,amrlev,mglev,mfi))
							*(ufab(m+dx+dy,k) + ufab(m-dx-dy,k) - ufab(m+dx-dy,k) - ufab(m-dx+dy,k))/(2.0*DX[0])/(2.0*DX[1])
							+
							C(i,1,k,1,m,amrlev,mglev,mfi)
							*(ufab(m+dy,k) - 2.0*ufab(m,k) + ufab(m-dy,k))/DX[1]/DX[1];

						// // C_{ijkl,j} u_{k,l}

						ffab(m,i) -=
							((C(i,0,k,0,m+dx,amrlev,mglev,mfi) - C(i,0,k,0,m-dx,amrlev,mglev,mfi))/(2.0*DX[0]) +
							 (C(i,1,k,0,m+dy,amrlev,mglev,mfi) - C(i,1,k,0,m-dy,amrlev,mglev,mfi))/(2.0*DX[1])) *
							((ufab(m+dx,k) - ufab(m-dx,k))/(2.0*DX[0]))
							+
							((C(i,0,k,1,m+dx,amrlev,mglev,mfi) - C(i,0,k,1,m-dx,amrlev,mglev,mfi))/(2.0*DX[0]) +
							 (C(i,1,k,1,m+dy,amrlev,mglev,mfi) - C(i,1,k,1,m-dy,amrlev,mglev,mfi))/(2.0*DX[1])) *
							((ufab(m+dy,k) - ufab(m-dy,k))/(2.0*DX[1]));

					}
				}
			}
	}
}


/// \fn Operator::FiniteElastic::Fsmooth
///
/// Perform one half Gauss-Seidel iteration corresponding to the operator specified
/// in Operator::FiniteElastic::Fapply.
/// The variable redblack corresponds to whether to smooth "red" nodes or "black"
/// nodes, where red and black nodes are distributed in a checkerboard pattern.
///
/// \todo Extend to 3D
///
void
Operator::FiniteElastic::FiniteElastic::Fsmooth (int amrlev,          ///<[in] AMR level
			      int mglev,           ///<[in]
			      MultiFab& u,       ///<[inout] Solution (displacement field)
			      const MultiFab& rhs, ///<[in] Body force vectors (rhs=right hand side)
			      int redblack         ///<[in] Smooth even vs. odd modes
			      ) const
{
	BL_PROFILE("Operator::FiniteElastic::Fsmooth()");

	const Real* dx = m_geom[amrlev][mglev].CellSize();

  for (MFIter mfi(u,MFItInfo().EnableTiling().SetDynamic(true));
       mfi.isValid(); ++mfi)
	{
		const Box&       tbx     = mfi.tilebox();
		FArrayBox&       ufab    = u[mfi];
		const FArrayBox& rhsfab  = rhs[mfi];

		for (int n = tbx.loVect()[1]; n<=tbx.hiVect()[1]; n++)
		{
			int noffset = (tbx.loVect()[0] + n + redblack)%2;
			for (int m = tbx.loVect()[0] + noffset; m <= tbx.hiVect()[0]; m+= 2)
			{
				for (int i=0; i<AMREX_SPACEDIM; i++)
				{
					amrex::Real rho = 0.0, aa = 0.0;
					for (int k=0; k<AMREX_SPACEDIM; k++)
					{
						rho -=
							C(i,0,k,0,amrex::IntVect(m,n),amrlev,mglev,mfi)
							*(ufab(amrex::IntVect(m+1,n),k) - (i==k ? 0.0 : 2.0*ufab(amrex::IntVect(m,n),k)) + ufab(amrex::IntVect(m-1,n),k))/dx[0]/dx[0]
							+
							(C(i,0,k,1,amrex::IntVect(m,n),amrlev,mglev,mfi) + C(i,1,k,0,amrex::IntVect(m,n),amrlev,mglev,mfi))
							*(ufab(amrex::IntVect(m+1,n+1),k) + ufab(amrex::IntVect(m-1,n-1),k) - ufab(amrex::IntVect(m+1,n-1),k) - ufab(amrex::IntVect(m-1,n+1),k))/(2.0*dx[0])/(2.0*dx[1])
							+
							C(i,1,k,1,amrex::IntVect(m,n),amrlev,mglev,mfi)
							*(ufab(amrex::IntVect(m,n+1),k) - (i==k ? 0.0 : 2.0*ufab(amrex::IntVect(m,n),k)) + ufab(amrex::IntVect(m,n-1),k))/dx[1]/dx[1];

						// C_{ijkl,j} u_{k,l}

						rho -=
							((C(i,0,k,0,amrex::IntVect(m+1,n),amrlev,mglev,mfi) - C(i,0,k,0,amrex::IntVect(m-1,n),amrlev,mglev,mfi))/(2.0*dx[0]) +
							(C(i,1,k,0,amrex::IntVect(m,n+1),amrlev,mglev,mfi) - C(i,1,k,0,amrex::IntVect(m,n-1),amrlev,mglev,mfi))/(2.0*dx[1])) *
							((ufab(amrex::IntVect(m+1,n),k) - ufab(amrex::IntVect(m-1,n),k))/(2.0*dx[0]))
							+
							((C(i,0,k,1,amrex::IntVect(m+1,n),amrlev,mglev,mfi) - C(i,0,k,1,amrex::IntVect(m-1,n),amrlev,mglev,mfi))/(2.0*dx[0]) +
							(C(i,1,k,1,amrex::IntVect(m,n+1),amrlev,mglev,mfi) - C(i,1,k,1,amrex::IntVect(m,n-1),amrlev,mglev,mfi))/(2.0*dx[1])) *
							((ufab(amrex::IntVect(m,n+1),k) - ufab(amrex::IntVect(m,n-1),k))/(2.0*dx[1]));
					}

					aa -=
						-2.0*C(i,0,i,0,amrex::IntVect(m,n),amrlev,mglev,mfi)/dx[0]/dx[0]
						-2.0*C(i,1,i,1,amrex::IntVect(m,n),amrlev,mglev,mfi)/dx[1]/dx[1];

					//std::cout << "nans not detetected, rho=" << rho << ", aa=" << aa << std::endl;
					if (rho != rho) std::cout << "nans detetected, rho=" << rho << ", aa=" << aa << std::endl;
					if (rho != rho) amrex::Abort("nans detected");

					ufab(amrex::IntVect(m,n),i) = (rhsfab(amrex::IntVect(m,n),i) - rho) / aa;
				}
			}
		}
	}
}

/// \fn Operator::FiniteElastic::FFlux
///
/// Compute the "flux" corresponding to the operator in Operator::FiniteElastic::Fapply.
/// Because the operator is self-adjoint and positive-definite, the flux is not
/// required for adequate convergence (?)
/// Therefore, the fluxes are simply set to zero and returned.
///
/// \todo Extend to 3D
///
void
Operator::FiniteElastic::FiniteElastic::FFlux (int amrlev, const MFIter& mfi,
					       const std::array<FArrayBox*,AMREX_SPACEDIM>& sigmafab,
					       const FArrayBox& ufab, const int /*face_only*/) const
{
  // THIS ONLY HAPPENS WHEN MULTIPLE AMR LEVELS ARE USED...

  amrex::BaseFab<amrex::Real> &fxfab = *sigmafab[0];
  amrex::BaseFab<amrex::Real> &fyfab = *sigmafab[1];
  fxfab.setVal(0.0);
  fyfab.setVal(0.0);

  const amrex::Real* DX = m_geom[amrlev][0].CellSize();

  amrex::IntVect dx(1,0);
  amrex::IntVect dy(0,1);

  const Box& bx = mfi.tilebox();
  for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++)
    for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++)
		{
			amrex::IntVect m(m1,m2);

			amrex::Real du1_dx1 = (ufab(m+dx,0) - ufab(m-dx,0))/(2.0*DX[0])  - 1.0;
			amrex::Real du1_dx2 = (ufab(m+dy,0) - ufab(m-dy,0))/(2.0*DX[1]);
			amrex::Real du2_dx1 = (ufab(m+dx,1) - ufab(m-dx,1))/(2.0*DX[0]);
			amrex::Real du2_dx2 = (ufab(m+dy,1) - ufab(m-dy,1))/(2.0*DX[1])  - 1.0;


			
			int i,j;

			i=0; j=0;
			fxfab(m,j) = 
			  C(i,j,0,0,m,amrlev,0,mfi)*du1_dx1 + C(i,j,0,1,m,amrlev,0,mfi)*du1_dx2 + C(i,j,1,0,m,amrlev,0,mfi)*du2_dx1 + C(i,j,1,1,m,amrlev,0,mfi)*du2_dx2;
			
			i=0; j=1;
			fxfab(m,j) = 
			  C(i,j,0,0,m,amrlev,0,mfi)*du1_dx1 + C(i,j,0,1,m,amrlev,0,mfi)*du1_dx2 + C(i,j,1,0,m,amrlev,0,mfi)*du2_dx1 + C(i,j,1,1,m,amrlev,0,mfi)*du2_dx2;
			fyfab(m,0) = fxfab(m,1);

			i=1; j=1;
			fyfab(m,j) = 
			  C(i,j,0,0,m,amrlev,0,mfi)*du1_dx1 + C(i,j,0,1,m,amrlev,0,mfi)*du1_dx2 + C(i,j,1,0,m,amrlev,0,mfi)*du2_dx1 + C(i,j,1,1,m,amrlev,0,mfi)*du2_dx2;


		}  
}


void
Operator::FiniteElastic::FiniteElastic::Stress (FArrayBox& sigmafab,
				    const FArrayBox& ufab,
				    int amrlev, const MFIter& mfi) const
{
  /// \todo add assert for sigmafab.ncomp=3 (for SPACEDIM=2) and ncomp=6 (for SPACEDIM=3)

  const amrex::Real* DX = m_geom[amrlev][0].CellSize();

  amrex::IntVect dx(1,0);
  amrex::IntVect dy(0,1);

  const Box& bx = mfi.tilebox();
  for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++)
    for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++)
		{
			amrex::IntVect m(m1,m2);

			amrex::Real du1_dx1 = (ufab(m+dx,0) - ufab(m-dx,0))/(2.0*DX[0]);
			amrex::Real du1_dx2 = (ufab(m+dy,0) - ufab(m-dy,0))/(2.0*DX[1]);
			amrex::Real du2_dx1 = (ufab(m+dx,1) - ufab(m-dx,1))/(2.0*DX[0]);
			amrex::Real du2_dx2 = (ufab(m+dy,1) - ufab(m-dy,1))/(2.0*DX[1]);

			for (int i=0; i<AMREX_SPACEDIM; i++)
				for (int j=0; j<AMREX_SPACEDIM; j++)
				{
					int voigt;
					if (i==0 && j==0) voigt = 0;
					if (i==1 && j==1) voigt = 1;
					if (i==0 && j==1) voigt = 2;
					if (i==1 && j==0) continue;

					sigmafab(m,voigt) =
						C(i,j,0,0,m,amrlev,0,mfi)*du1_dx1 +
						C(i,j,0,1,m,amrlev,0,mfi)*du1_dx2 +
						C(i,j,1,0,m,amrlev,0,mfi)*du2_dx1 +
						C(i,j,1,1,m,amrlev,0,mfi)*du2_dx2;
				}
			}

}

void
Operator::FiniteElastic::FiniteElastic::Energy (FArrayBox& energyfab,
				    const FArrayBox& ufab,
				    int amrlev, const MFIter& mfi) const
{
  const amrex::Real* DX = m_geom[amrlev][0].CellSize();

  amrex::IntVect dx(1,0);
  amrex::IntVect dy(0,1);

  const Box& bx = mfi.tilebox();
  for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++)
    for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++)
  	{
			amrex::IntVect m(m1,m2);
			amrex::Real gradu[2][2] = {{(ufab(m+dx,0) - ufab(m-dx,0))/(2.0*DX[0]),
				(ufab(m+dy,0) - ufab(m-dy,0))/(2.0*DX[1])},
				{(ufab(m+dx,1) - ufab(m-dx,1))/(2.0*DX[0]),
				(ufab(m+dy,1) - ufab(m-dy,1))/(2.0*DX[1])}};

			energyfab(m) = 0.0;

			for (int i=0; i<AMREX_SPACEDIM; i++)
				for (int j=0; j<AMREX_SPACEDIM; j++)
					for (int k=0; k<AMREX_SPACEDIM; k++)
						for (int l=0; l<AMREX_SPACEDIM; l++)
							energyfab(m) += gradu[i][j] * C(i,j,k,l,m,amrlev,0,mfi) * gradu[k][l];
			energyfab(m) *= 0.5;
    }

}
