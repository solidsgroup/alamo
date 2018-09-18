#include <AMReX_MultiFabUtil.H>
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>
#include "Set/Set.H"

#include <AMReX_ArrayLim.H>

#include "Set/Set.H"
#include "Elastic.H"



namespace Operator
{
namespace Elastic
{
/// \fn Operator::Elastic::MLPFStiffnessMatrix
///
/// Relay to the define function
/// Also define elastic constants here.
Elastic::Elastic (const Vector<Geometry>& a_geom,
		  const Vector<BoxArray>& a_grids,
		  const Vector<DistributionMapping>& a_dmap,
		  const LPInfo& a_info)
{
	define(a_geom, a_grids, a_dmap// , a_bc, a_info
	       );
}

Elastic::~Elastic ()
{}

// void
// Elastic::SetEigenstrain(amrex::Vector<std::unique_ptr<amrex::MultiFab> > &eigenstrain)
// {
// 	usingEigenstrain = true;
// 	AMREX_ASSERT(eigenstrain[0]->nComp() == AMREX_SPACEDIM*AMREX_SPACEDIM);
// 	RegisterNewFab(eigenstrain);
// }

// void
// Elastic::SetEigenstrain(amrex::Vector<amrex::MultiFab> &eigenstrain)
// {
// 	usingEigenstrain = true;
// 	AMREX_ASSERT(eigenstrain[0].nComp() == AMREX_SPACEDIM*AMREX_SPACEDIM);
// 	RegisterNewFab(eigenstrain);
// }

void
Elastic::Fapply (int amrlev, ///<[in] AMR Level
		 int mglev,  ///<[in]
		 MultiFab& f,///<[out] The force vector
		 const MultiFab& u ///<[in] The displacements vector
		 ) const
{
	BL_PROFILE("Operator::Elastic::Elastic::Fapply()");

	amrex::Box domain(m_geom[amrlev][mglev].Domain());
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
		const amrex::FArrayBox &ufab    = u[mfi];
		amrex::FArrayBox       &ffab    = f[mfi];
		// if(ffab.contains_inf()) Util::Abort("Inf in ffab [before update]");
		// if(ffab.contains_nan()) Util::Abort("Nan in ffab [before update]");
		if(ufab.contains_inf()) Util::Abort("Inf in ufab [before update]");
		if(ufab.contains_nan()) Util::Abort("Nan in ufab [before update]");
	}

	for (MFIter mfi(f, true); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();
		const amrex::FArrayBox &ufab    = u[mfi];
		amrex::FArrayBox       &ffab    = f[mfi];

		AMREX_D_TERM(for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++),
			     for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++),
			     for (int m3 = bx.loVect()[2]; m3<=bx.hiVect()[2]; m3++))
		{
			amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));
			for (int i=0; i<AMREX_SPACEDIM; i++)
			{
				ffab(amrex::IntVect(AMREX_D_DECL(m1,m2,m3)),i) = 0.0;

				for (int k=0; k<AMREX_SPACEDIM; k++)
				{

					if (m1 == domain.loVect()[0] || m1 == domain.hiVect()[0]+1)// ||
						//     //m2 == domain.loVect()[1] || m2 == domain.hiVect()[1] +1
					//     //m3 == domain.loVect()[2]
					//     //m2 == domain.loVect()[1]
					//     m3 == domain.loVect()[2]
					//     )
					{
					 	ffab(m,k) = ufab(m,k);
					 	continue;
					}
					// else
					// {
					// 	ffab(m,k) = 10.0*ufab(m,k);
					// 	continue;
					// }

					Set::Vector gradu_k; // gradu_k(l) = u_{k,l}
					AMREX_D_TERM(gradu_k(0) = (ufab(m+dx,k) - ufab(m-dx,k))/(2.0*DX[0]);,
						     gradu_k(1) = (ufab(m+dy,k) - ufab(m-dy,k))/(2.0*DX[1]);,
						     gradu_k(2) = (ufab(m+dz,k) - ufab(m-dz,k))/(2.0*DX[2]);)

						Set::Matrix gradgradu_k; // gradgradu_k(l,j) = u_{k,lj}
					AMREX_D_TERM(gradgradu_k(0,0) = (ufab(m+dx,k) - 2.0*ufab(m,k) + ufab(m-dx,k))/DX[0]/DX[0];
						     ,// 2D
						     gradgradu_k(0,1) = (ufab(m+dx+dy,k) + ufab(m-dx-dy,k) - ufab(m+dx-dy,k) - ufab(m-dx+dy,k))/(2.0*DX[0])/(2.0*DX[1]);
						     gradgradu_k(1,0) = gradgradu_k(0,1);
						     gradgradu_k(1,1) = (ufab(m+dy,k) - 2.0*ufab(m,k) + ufab(m-dy,k))/DX[1]/DX[1];
						     ,// 3D
						     gradgradu_k(0,2) = (ufab(m+dx+dz,k) + ufab(m-dx-dz,k) - ufab(m+dx-dz,k) - ufab(m-dx+dz,k))/(2.0*DX[0])/(2.0*DX[2]);
						     gradgradu_k(1,2) = (ufab(m+dy+dz,k) + ufab(m-dy-dz,k) - ufab(m+dy-dz,k) - ufab(m-dy+dz,k))/(2.0*DX[0])/(2.0*DX[2]);
						     gradgradu_k(2,0) = gradgradu_k(0,2);
						     gradgradu_k(2,1) = gradgradu_k(1,2);
						     gradgradu_k(2,2) = (ufab(m+dz,k) - 2.0*ufab(m,k) + ufab(m-dz,k))/DX[2]/DX[2];);

					Set::Vector C_ik; // C_ik(l) = C_{ijkl,j}
					AMREX_D_TERM(C_ik(0) = AMREX_D_TERM(+ (C(i,0,k,0,m+dx,amrlev,mglev,mfi) - C(i,0,k,0,m-dx,amrlev,mglev,mfi))/(2.0*DX[0]),
									    + (C(i,1,k,0,m+dy,amrlev,mglev,mfi) - C(i,1,k,0,m-dy,amrlev,mglev,mfi))/(2.0*DX[1]),
									    + (C(i,2,k,0,m+dz,amrlev,mglev,mfi) - C(i,2,k,0,m-dz,amrlev,mglev,mfi))/(2.0*DX[2]));,
						     C_ik(1) = AMREX_D_TERM(+ (C(i,0,k,1,m+dx,amrlev,mglev,mfi) - C(i,0,k,1,m-dx,amrlev,mglev,mfi))/(2.0*DX[0]),
									    + (C(i,1,k,1,m+dy,amrlev,mglev,mfi) - C(i,1,k,1,m-dy,amrlev,mglev,mfi))/(2.0*DX[1]),
									    + (C(i,2,k,1,m+dz,amrlev,mglev,mfi) - C(i,2,k,1,m-dz,amrlev,mglev,mfi))/(2.0*DX[2]));,
						     C_ik(2) = AMREX_D_TERM(+ (C(i,0,k,2,m+dx,amrlev,mglev,mfi) - C(i,0,k,2,m-dx,amrlev,mglev,mfi))/(2.0*DX[0]),
									    + (C(i,1,k,2,m+dy,amrlev,mglev,mfi) - C(i,1,k,2,m-dy,amrlev,mglev,mfi))/(2.0*DX[1]),
									    + (C(i,2,k,2,m+dz,amrlev,mglev,mfi) - C(i,2,k,2,m-dz,amrlev,mglev,mfi))/(2.0*DX[2])););
							
					for (int l=0; l<AMREX_SPACEDIM; l++)
					{
						// f_i -= C_{ijkl} (u_{k,lj})
						for (int j=0; j<AMREX_SPACEDIM; j++)
							ffab(m,i) -= C(i,j,k,l,m,amrlev,mglev,mfi) * (gradgradu_k(j,l));

						// f_i -= C_{ijkl,j} u_{k,l}
						ffab(m,i) -= C_ik(l) * gradu_k(l);
					}
					if(std::isinf(ffab(m,i))) Util::Abort("Inf in ffab(m,k) [after update]");
					if(std::isnan(ffab(m,i))) Util::Abort("Nan in ffab(m,k) [after update]");
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
Elastic::Fsmooth (int amrlev,          ///<[in] AMR level
		  int mglev,           ///<[in]
		  MultiFab& u,       ///<[inout] Solution (displacement field)
		  const MultiFab& rhs ///<[in] Body force vectors (rhs=right hand side)
		  ) const
{
	for (int redblack = 0; redblack < 2; redblack++)
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

			AMREX_D_TERM(for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++),
				     for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++),
				     for (int m3 = bx.loVect()[2]; m3<=bx.hiVect()[2]; m3++))
			{
				if ((AMREX_D_TERM(m1, + m2, + m3))%2 == redblack) continue;
				amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));

				for (int i=0; i<AMREX_SPACEDIM; i++)
				{
					amrex::Real rho = 0.0, aa = 0.0;
					for (int k=0; k<AMREX_SPACEDIM; k++)
					{
						AMREX_D_TERM(if(std::isnan(ufab(m,k))) std::cout << "Nan in ufab(m,k)" << std::endl;
							     ,
							     if(std::isnan(ufab(m+dx,k))) std::cout << "Nan in ufab(m+dx,k)" << std::endl;
							     if(std::isnan(ufab(m-dx,k))) std::cout << "Nan in ufab(m-dx,k)" << std::endl;
							     if(std::isnan(ufab(m+dy,k))) std::cout << "Nan in ufab(m+dy,k)" << std::endl;
							     if(std::isnan(ufab(m-dy,k))) std::cout << "Nan in ufab(m-dy,k)" << std::endl;
							     if(std::isnan(ufab(m+dx+dy,k))) std::cout << "Nan in ufab(m+dx+dy,k)" << std::endl;
							     if(std::isnan(ufab(m-dx+dy,k))) std::cout << "Nan in ufab(m-dx+dy,k)" << std::endl;
							     ,
							     if(std::isnan(ufab(m+dz,k))) std::cout << "Nan in ufab(m+dz,k)" << std::endl;
							     if(std::isnan(ufab(m-dz,k))) std::cout << "Nan in ufab(m-dz,k)" << std::endl;
							     if(std::isnan(ufab(m+dx+dz,k))) std::cout << "Nan in ufab(m+dx+dz,k)" << std::endl;
							     if(std::isnan(ufab(m+dx-dy,k))) std::cout << "Nan in ufab(m+dx-dz,k)" << std::endl;
							     if(std::isnan(ufab(m-dx+dy,k))) std::cout << "Nan in ufab(m-dx+dy,k)" << std::endl;);


						Set::Vector gradu_k; // gradu_k(l) = u_{k,l}
						AMREX_D_TERM(gradu_k(0) = (ufab(m+dx,k) - ufab(m-dx,k))/(2.0*DX[0]);,
							     gradu_k(1) = (ufab(m+dy,k) - ufab(m-dy,k))/(2.0*DX[1]);,
							     gradu_k(2) = (ufab(m+dz,k) - ufab(m-dz,k))/(2.0*DX[1]););

						Set::Matrix gradgradu_k; // gradgradu_k(l,j) = u_{k,lj}
						AMREX_D_TERM(gradgradu_k(0,0) = (ufab(m+dx,k) - (i==k ? 0.0 : 2.0*ufab(m,k)) + ufab(m-dx,k))/DX[0]/DX[0];
							     ,// 2D
							     gradgradu_k(1,1) = (ufab(m+dy,k) - (i==k ? 0.0 : 2.0*ufab(m,k)) + ufab(m-dy,k))/DX[1]/DX[1];
							     gradgradu_k(0,1) = (ufab(m+dx+dy,k) + ufab(m-dx-dy,k) - ufab(m+dx-dy,k) - ufab(m-dx+dy,k))/(2.0*DX[0])/(2.0*DX[1]);
							     gradgradu_k(1,0) = gradgradu_k(0,1);
							     ,// 3D
							     gradgradu_k(2,2) = (ufab(m+dz,k) - (i==k ? 0.0 : 2.0*ufab(m,k)) + ufab(m-dz,k))/DX[2]/DX[2];
							     gradgradu_k(0,2) = (ufab(m+dx+dz,k) + ufab(m-dx-dz,k) - ufab(m+dx-dz,k) - ufab(m-dx+dz,k))/(2.0*DX[0])/(2.0*DX[2]);
							     gradgradu_k(1,2) = (ufab(m+dy+dz,k) + ufab(m-dy-dz,k) - ufab(m+dy-dz,k) - ufab(m-dy+dz,k))/(2.0*DX[0])/(2.0*DX[2]);
							     gradgradu_k(2,0) = gradgradu_k(0,2);
							     gradgradu_k(2,1) = gradgradu_k(1,2););

						Set::Vector C_ik; // C_ik(l) = C_{ijkl,j}
						AMREX_D_TERM(C_ik(0) = AMREX_D_TERM(+ (C(i,0,k,0,m+dx,amrlev,mglev,mfi) - C(i,0,k,0,m-dx,amrlev,mglev,mfi))/(2.0*DX[0]),
										    + (C(i,1,k,0,m+dy,amrlev,mglev,mfi) - C(i,1,k,0,m-dy,amrlev,mglev,mfi))/(2.0*DX[1]),
										    + (C(i,2,k,0,m+dz,amrlev,mglev,mfi) - C(i,2,k,0,m-dz,amrlev,mglev,mfi))/(2.0*DX[2]));,
							     C_ik(1) = AMREX_D_TERM(+ (C(i,0,k,1,m+dx,amrlev,mglev,mfi) - C(i,0,k,1,m-dx,amrlev,mglev,mfi))/(2.0*DX[0]),
										    + (C(i,1,k,1,m+dy,amrlev,mglev,mfi) - C(i,1,k,1,m-dy,amrlev,mglev,mfi))/(2.0*DX[1]),
										    + (C(i,2,k,1,m+dz,amrlev,mglev,mfi) - C(i,2,k,1,m-dz,amrlev,mglev,mfi))/(2.0*DX[2]));,
							     C_ik(2) = AMREX_D_TERM(+ (C(i,0,k,2,m+dx,amrlev,mglev,mfi) - C(i,0,k,2,m-dx,amrlev,mglev,mfi))/(2.0*DX[0]),
										    + (C(i,1,k,2,m+dy,amrlev,mglev,mfi) - C(i,1,k,2,m-dy,amrlev,mglev,mfi))/(2.0*DX[1]),
										    + (C(i,2,k,2,m+dz,amrlev,mglev,mfi) - C(i,2,k,2,m-dz,amrlev,mglev,mfi))/(2.0*DX[2])););

						for (int l=0; l<AMREX_SPACEDIM; l++)
						{
							// rho -= C_{ijkl} (u_{k,lj})
							for (int j=0; j<AMREX_SPACEDIM; j++)
								rho -= C(i,j,k,l,m,amrlev,mglev,mfi) * (gradgradu_k(j,l));

							// rho -= C_{ijkl,j} u_{k,l}
							rho -= C_ik(l) * gradu_k(l);
						}
					}

					aa -= AMREX_D_TERM( - 2.0*C(i,0,i,0,m,amrlev,mglev,mfi)/DX[0]/DX[0],
							    - 2.0*C(i,1,i,1,m,amrlev,mglev,mfi)/DX[1]/DX[1], 
							    - 2.0*C(i,2,i,2,m,amrlev,mglev,mfi)/DX[2]/DX[2]);

					if (std::isnan(rho))
					{
						std::cout << "WARNING: nans detetected, rho=" << rho << ", aa=" << aa << std::endl;
						amrex::Abort("nans detected");
					}
					else
					{
						ufab(m,i) = (rhsfab(m,i) - rho) / aa;
					}
				}
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
Elastic::FFlux (int /*amrlev*/, const MFIter& /*mfi*/,
		const std::array<FArrayBox*,AMREX_SPACEDIM>& sigmafab,
		const FArrayBox& /*ufab*/, const int /*face_only*/) const
{
	// amrex::BaseFab<amrex::Real> &fxfab = *sigmafab[0];
	// amrex::BaseFab<amrex::Real> &fyfab = *sigmafab[1];
	// fxfab.setVal(0.0);
	// fyfab.setVal(0.0);

	amrex::BaseFab<amrex::Real> AMREX_D_DECL( &fxfab = *sigmafab[0],
	 					  &fyfab = *sigmafab[1],
	 					  &fzfab = *sigmafab[2] ) ;
	AMREX_D_TERM(fxfab.setVal(0.0);,
	 	     fyfab.setVal(0.0);,
	 	     fzfab.setVal(0.0););

}

// void Elastic::AddEigenstrainToRHS (amrex::Vector<amrex::MultiFab>& rhsfab) const
// {
// 	for (int amrlev=0; amrlev < rhsfab.size(); amrlev ++)
// 	{
// 		for (amrex::MFIter mfi(rhsfab[amrlev], true); mfi.isValid(); ++mfi)
// 		{
// 			amrex::FArrayBox       &rhs    = rhsfab[amrlev][mfi];
// 			AddEigenstrainToRHS(rhs,amrlev,mfi);
// 		}
// 	}
// }

// void Elastic::AddEigenstrainToRHS (FArrayBox& rhsfab,
// 				   int amrlev, const MFIter& mfi) const
// {
// 	BL_PROFILE("Operator::Elastic::AddEigenstrainToRHS()");

// 	if (!usingEigenstrain) return;

// 	const Real* DX = m_geom[amrlev][0].CellSize();
  
// #if AMREX_SPACEDIM == 1
// 	static amrex::IntVect dx(1);
// #elif AMREX_SPACEDIM == 2
// 	static amrex::IntVect dx(1,0), dy(0,1);
// #elif AMREX_SPACEDIM == 3
// 	static amrex::IntVect dx(1,0,0), dy(0,1,0), dz(0,0,1);
// #endif 

// 	const Box& bx = mfi.tilebox();
// 	const amrex::FArrayBox &eps0fab = GetFab(0,amrlev,0,mfi);

// 	for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++)
// 		for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++)
// #if AMREX_SPACEDIM > 2
// 			for (int m3 = bx.loVect()[2]; m3<=bx.hiVect()[2]; m3++)
// #endif
// 			{
// 				amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));
// 				for (int i=0; i<AMREX_SPACEDIM; i++)
// 				{
// 					for (int k=0; k<AMREX_SPACEDIM; k++)
// 					{

// 						Set::Matrix gradeps0_k; //  gradeps0_k(l,j) == eps0_{kl,j}
// 						gradeps0_k(0,0) = (eps0fab(m+dx,k*AMREX_SPACEDIM + 0) - eps0fab(m-dx,k*AMREX_SPACEDIM + 0)) / (2.0*DX[0]); // eps0_{k0,0}
// 						gradeps0_k(0,1) = (eps0fab(m+dy,k*AMREX_SPACEDIM + 0) - eps0fab(m-dy,k*AMREX_SPACEDIM + 0)) / (2.0*DX[1]); // eps0_{k0,1}
// 						gradeps0_k(1,0) = (eps0fab(m+dx,k*AMREX_SPACEDIM + 1) - eps0fab(m-dx,k*AMREX_SPACEDIM + 1)) / (2.0*DX[0]); // eps0_{k1,0}
// 						gradeps0_k(1,1) = (eps0fab(m+dy,k*AMREX_SPACEDIM + 1) - eps0fab(m-dy,k*AMREX_SPACEDIM + 1)) / (2.0*DX[1]); // eps0_{k1,1}
// #if AMREX_SPACEDIM > 2
// 						gradeps0_k(0,2) = (eps0fab(m+dz,k*AMREX_SPACEDIM + 0) - eps0fab(m-dz,k*AMREX_SPACEDIM + 0)) / (2.0*DX[2]); // eps0_{k0,1}
// 						gradeps0_k(1,2) = (eps0fab(m+dz,k*AMREX_SPACEDIM + 1) - eps0fab(m-dz,k*AMREX_SPACEDIM + 1)) / (2.0*DX[2]); // eps0_{k0,1}
// 						gradeps0_k(2,0) = (eps0fab(m+dx,k*AMREX_SPACEDIM + 2) - eps0fab(m-dx,k*AMREX_SPACEDIM + 2)) / (2.0*DX[0]); // eps0_{k0,0}
// 						gradeps0_k(2,1) = (eps0fab(m+dy,k*AMREX_SPACEDIM + 2) - eps0fab(m-dy,k*AMREX_SPACEDIM + 2)) / (2.0*DX[1]); // eps0_{k0,1}
// 						gradeps0_k(2,2) = (eps0fab(m+dz,k*AMREX_SPACEDIM + 2) - eps0fab(m-dz,k*AMREX_SPACEDIM + 2)) / (2.0*DX[2]); // eps0_{k0,1}
// #endif

// 						Set::Vector eps0_k; // eps0_k(l) = eps0_{kl}
// 						eps0_k(0) = eps0fab(m,k*AMREX_SPACEDIM + 0);
// 						eps0_k(1) = eps0fab(m,k*AMREX_SPACEDIM + 1);
// #if AMREX_SPACEDIM > 2
// 						eps0_k(2) = eps0fab(m,k*AMREX_SPACEDIM + 2);
// #endif


// 						// C_{ijkl} (u_{k,lj} - eps0_{kl,j})

// 						for (int j=0; j<AMREX_SPACEDIM; j++)
// 							for (int l=0; l<AMREX_SPACEDIM; l++)
// 								rhsfab(m,i) -= C(i,j,k,l,m,amrlev,0,mfi) * (gradeps0_k(j,l));


// 						// C_{ijkl,j} (u_{k,l} - eps0_{k,l}
// #if AMREX_SPACEDIM == 2
// 						rhsfab(m,i) -=
// 							((C(i,0,k,0,m+dx,amrlev,0,mfi) - C(i,0,k,0,m-dx,amrlev,0,mfi))/(2.0*DX[0]) +
// 							 (C(i,1,k,0,m+dy,amrlev,0,mfi) - C(i,1,k,0,m-dy,amrlev,0,mfi))/(2.0*DX[1])) *
// 							(eps0_k(0))
// 							+
// 							((C(i,0,k,1,m+dx,amrlev,0,mfi) - C(i,0,k,1,m-dx,amrlev,0,mfi))/(2.0*DX[0]) +
// 							 (C(i,1,k,1,m+dy,amrlev,0,mfi) - C(i,1,k,1,m-dy,amrlev,0,mfi))/(2.0*DX[1])) *
// 							(eps0_k(1));
// #elif AMREX_SPACEIM == 3
// 						rhsfab(m,i) -=
// 							((C(i,0,k,0,m+dx,amrlev,0,mfi) - C(i,0,k,0,m-dx,amrlev,0,mfi))/(2.0*DX[0]) +
// 							 (C(i,1,k,0,m+dy,amrlev,0,mfi) - C(i,1,k,0,m-dy,amrlev,0,mfi))/(2.0*DX[1]) +
// 							 (C(i,2,k,0,m+dz,amrlev,0,mfi) - C(i,2,k,0,m-dz,amrlev,0,mfi))/(2.0*DX[2])) *
// 							(eps0_k(0))
// 							+
// 							((C(i,0,k,1,m+dx,amrlev,0,mfi) - C(i,0,k,1,m-dx,amrlev,0,mfi))/(2.0*DX[0]) +
// 							 (C(i,1,k,1,m+dy,amrlev,0,mfi) - C(i,1,k,1,m-dy,amrlev,0,mfi))/(2.0*DX[1]) +
// 							 (C(i,2,k,1,m+dz,amrlev,0,mfi) - C(i,2,k,1,m-dz,amrlev,0,mfi))/(2.0*DX[2])) *
// 							(eps0_k(1))
// 							+
// 							((C(i,0,k,2,m+dx,amrlev,0,mfi) - C(i,0,k,2,m-dx,amrlev,0,mfi))/(2.0*DX[0]) +
// 							 (C(i,1,k,2,m+dy,amrlev,0,mfi) - C(i,1,k,2,m-dy,amrlev,0,mfi))/(2.0*DX[1]) +
// 							 (C(i,2,k,2,m+dz,amrlev,0,mfi) - C(i,2,k,2,m-dz,amrlev,0,mfi))/(2.0*DX[2])) *
// 							(eps0_k(2));
// #endif
// 					}
// 				}
// 			}
// }

void
Elastic::Stress (FArrayBox& sigmafab,
		 const FArrayBox& ufab,
		 int amrlev, const MFIter& mfi,
		 bool voigt) const
{
	if (voigt)
		AMREX_ASSERT(sigmafab.nComp() == (AMREX_SPACEDIM*(AMREX_SPACEDIM-1)/2));
	else
		AMREX_ASSERT(sigmafab.nComp() == AMREX_SPACEDIM*AMREX_SPACEDIM);

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

				Set::Matrix eps0 = Set::Matrix::Zero();
// 				if (usingEigenstrain)
// 				{
// 					const FArrayBox &eps0fab = GetFab(0,amrlev,0,mfi);
// #if AMREX_SPACEDIM ==2
// 					eps0 << eps0fab(m,0), eps0fab(m,1), eps0fab(m,2), eps0fab(m,3);
// #elif AMREX_SPACEDIM ==3
// 					eps0 << eps0fab(m,0), eps0fab(m,1), eps0fab(m,2), eps0fab(m,3), eps0fab(m,4), eps0fab(m,5), eps0fab(m,6), eps0fab(m,7), eps0fab(m,8);
// #endif
// 				}

				for (int i=0; i<AMREX_SPACEDIM; i++)
					for (int j=0; j<AMREX_SPACEDIM; j++)
					{
						for (int k=0; k<AMREX_SPACEDIM; k++)
							for (int l=0; l<AMREX_SPACEDIM; l++)
								sigmafab(m,i*AMREX_SPACEDIM + j) += C(i,j,k,l,m,amrlev,0,mfi) * (gradu(k,l) - eps0(k,l));
					}
			}
}

void
Elastic::Energy (FArrayBox& energyfab,
		 const FArrayBox& ufab,
		 int amrlev, const MFIter& mfi) const
{
	/// \todo RE-IMPLEMENT in 2D and 3D

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

				Set::Matrix eps0 = Set::Matrix::Zero();
// 				if (usingEigenstrain)
// 				{
// 					const FArrayBox &eps0fab = GetFab(0,amrlev,0,mfi);
// #if AMREX_SPACEDIM ==2
// 					eps0 << eps0fab(m,0), eps0fab(m,1), eps0fab(m,2), eps0fab(m,3);
// #elif AMREX_SPACEDIM ==3
// 					eps0 << eps0fab(m,0), eps0fab(m,1), eps0fab(m,2), eps0fab(m,3), eps0fab(m,4), eps0fab(m,5), eps0fab(m,6), eps0fab(m,7), eps0fab(m,8);
// #endif
// 				}

				energyfab(m) = 0.0;

				for (int i=0; i<AMREX_SPACEDIM; i++)
					for (int j=0; j<AMREX_SPACEDIM; j++)
						for (int k=0; k<AMREX_SPACEDIM; k++)
							for (int l=0; l<AMREX_SPACEDIM; l++)
								energyfab(m) += 0.5 * (gradu(i,j) - eps0(i,j)) * C(i,j,k,l,m,amrlev,0,mfi) * (gradu(k,l) - eps0(k,l));
			}
}






}
}
