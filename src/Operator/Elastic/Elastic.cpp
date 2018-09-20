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



void
Elastic::Fapply (int amrlev,
		 int mglev,
		 const amrex::Box& bx,
		 amrex::FArrayBox& ffab,
		 const amrex::FArrayBox& ufab
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

				// DIRICHLET boundary conditions
				if (false
				    || m1 == domain.loVect()[0]    
				    || m1 == domain.hiVect()[0]+1  
				    || m2 == domain.loVect()[1]    
				    || m2 == domain.hiVect()[1]+1  
				    || m3 == domain.loVect()[2]    
				    || m3 == domain.hiVect()[2]+1  
				    )
				{
					ffab(m,k) = ufab(m,k);
					continue;
				}
				// IGNORE boundaries if not part of domain boundary
				if ( m1 == bx.loVect()[0] ||
				     m1 == bx.hiVect()[0] ||
				     m2 == bx.loVect()[1] ||
				     m2 == bx.hiVect()[1] ||
				     m3 == bx.loVect()[2] ||
				     m3 == bx.hiVect()[2]) continue;

				// ELSE perform operator

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
				AMREX_D_TERM(C_ik(0) = AMREX_D_TERM(+ (C(i,0,k,0,m+dx) - C(i,0,k,0,m-dx))/(2.0*DX[0]),
								    + (C(i,1,k,0,m+dy) - C(i,1,k,0,m-dy))/(2.0*DX[1]),
								    + (C(i,2,k,0,m+dz) - C(i,2,k,0,m-dz))/(2.0*DX[2]));,
					     C_ik(1) = AMREX_D_TERM(+ (C(i,0,k,1,m+dx) - C(i,0,k,1,m-dx))/(2.0*DX[0]),
								    + (C(i,1,k,1,m+dy) - C(i,1,k,1,m-dy))/(2.0*DX[1]),
								    + (C(i,2,k,1,m+dz) - C(i,2,k,1,m-dz))/(2.0*DX[2]));,
					     C_ik(2) = AMREX_D_TERM(+ (C(i,0,k,2,m+dx) - C(i,0,k,2,m-dx))/(2.0*DX[0]),
								    + (C(i,1,k,2,m+dy) - C(i,1,k,2,m-dy))/(2.0*DX[1]),
								    + (C(i,2,k,2,m+dz) - C(i,2,k,2,m-dz))/(2.0*DX[2])););
							
				for (int l=0; l<AMREX_SPACEDIM; l++)
				{
					// f_i -= C_{ijkl} (u_{k,lj})
					for (int j=0; j<AMREX_SPACEDIM; j++)
						ffab(m,i) -= C(i,j,k,l,m) * (gradgradu_k(j,l));

					// f_i -= C_{ijkl,j} u_{k,l}
					ffab(m,i) -= C_ik(l) * gradu_k(l);
				}
			}
		}
	}
}



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

	for (MFIter mfi(f, true); mfi.isValid(); ++mfi)
		if(u[mfi].contains_inf()) Util::Abort("Inf in ufab [before update]");

	for (MFIter mfi(f, true); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();
		const amrex::FArrayBox &ufab    = u[mfi];
		amrex::FArrayBox       &ffab    = f[mfi];
		Fapply(amrlev, mglev, bx, ffab, ufab);
	}

	for (MFIter mfi(f, true); mfi.isValid(); ++mfi)
		if(f[mfi].contains_nan()) Util::Abort("Nan in ffab [after update]");
}


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

		amrex::Box domain(m_geom[amrlev][mglev].Domain());
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
						// DIRICHLET BCs on operator
						if (true
						    || m1 == domain.loVect()[0]     
						    || m1 == domain.hiVect()[0]+1   
						    || m2 == domain.loVect()[1]     
						    || m2 == domain.hiVect()[1]+1   
						    || m3 == domain.loVect()[2]     
						    || m3 == domain.hiVect()[2]+1   
						    )
						{
							ufab(m,k) = rhsfab(m,k);
							continue;
						}
						// DO NOTHING on interior boundaries
						else if ( m1 == bx.loVect()[0] ||
							  m1 == bx.hiVect()[0] ||
							  m2 == bx.loVect()[1] ||
							  m2 == bx.hiVect()[1] ||
							  m3 == bx.loVect()[2] ||
							  m3 == bx.hiVect()[2]) continue;



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
							     gradu_k(2) = (ufab(m+dz,k) - ufab(m-dz,k))/(2.0*DX[2]););

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
						AMREX_D_TERM(C_ik(0) = AMREX_D_TERM(+ (C(i,0,k,0,m+dx) - C(i,0,k,0,m-dx))/(2.0*DX[0]),
										    + (C(i,1,k,0,m+dy) - C(i,1,k,0,m-dy))/(2.0*DX[1]),
										    + (C(i,2,k,0,m+dz) - C(i,2,k,0,m-dz))/(2.0*DX[2]));,
							     C_ik(1) = AMREX_D_TERM(+ (C(i,0,k,1,m+dx) - C(i,0,k,1,m-dx))/(2.0*DX[0]),
										    + (C(i,1,k,1,m+dy) - C(i,1,k,1,m-dy))/(2.0*DX[1]),
										    + (C(i,2,k,1,m+dz) - C(i,2,k,1,m-dz))/(2.0*DX[2]));,
							     C_ik(2) = AMREX_D_TERM(+ (C(i,0,k,2,m+dx) - C(i,0,k,2,m-dx))/(2.0*DX[0]),
										    + (C(i,1,k,2,m+dy) - C(i,1,k,2,m-dy))/(2.0*DX[1]),
										    + (C(i,2,k,2,m+dz) - C(i,2,k,2,m-dz))/(2.0*DX[2])););

						for (int l=0; l<AMREX_SPACEDIM; l++)
						{
							// rho -= C_{ijkl} (u_{k,lj})
							for (int j=0; j<AMREX_SPACEDIM; j++)
								rho -= C(i,j,k,l,m) * (gradgradu_k(j,l));

							// rho -= C_{ijkl,j} u_{k,l}
							rho -= C_ik(l) * gradu_k(l);
						}
					}

					aa -= AMREX_D_TERM( - 2.0*C(i,0,i,0,m)/DX[0]/DX[0],
							    - 2.0*C(i,1,i,1,m)/DX[1]/DX[1], 
							    - 2.0*C(i,2,i,2,m)/DX[2]/DX[2]);

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


//
// APPLICATION-SPECIFIC
//

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


void
Elastic::FineResidualContribution (int amrlev, int mglev,
				   const amrex::Box &c,
				   const amrex::Box &cg,
				   amrex::FArrayBox &f,
				   const amrex::FArrayBox &phi,
				   const amrex::FArrayBox &Ax, 
				   const amrex::IArrayBox& dmsk,
				   const amrex::Real *dxinv) const
{
	

//    integer, dimension(2), intent(in) :: clo, chi, cglo, cghi, flo, fhi, xlo, xhi, &
//         slo, shi, alo, ahi, mlo, mhi
//    real(amrex_real), intent(inout) :: f  (flo(1):fhi(1),flo(2):fhi(2))
//    real(amrex_real), intent(in   ) :: x  (xlo(1):xhi(1),xlo(2):xhi(2))
//    real(amrex_real), intent(in   ) :: sig(slo(1):shi(1),slo(2):shi(2))
//    real(amrex_real), intent(inout) :: Ax (alo(1):ahi(1),alo(2):ahi(2))
//    integer         , intent(in   ) :: msk(mlo(1):mhi(1),mlo(2):mhi(2))
//    real(amrex_real), intent(in) :: dxinv(2)
//
//    integer, dimension(2) :: lo, hi, glo, ghi, gtlo, gthi
//    integer :: i, j, ii, jj, step
//    real(amrex_real) :: facx, facy, fxy, f2xmy, fmx2y, fm, fp
//    real(amrex_real), parameter :: rfd = 0.25d0
//    real(amrex_real), parameter :: chip = 0.5d0
//    real(amrex_real), parameter :: chip2 = 0.25d0
//
//    facx = (1.d0/6.d0)*dxinv(1)*dxinv(1)
//    facy = (1.d0/6.d0)*dxinv(2)*dxinv(2)
//    fxy = facx + facy
//    f2xmy = 2.d0*facx - facy
//    fmx2y = 2.d0*facy - facx
//
//    lo = 2*clo
//    hi = 2*chi
//    glo = 2*cglo
//    ghi = 2*cghi
//
//    gtlo = max(lo-1,glo)
//    gthi = min(hi+1,ghi)
//


//    do jj = gtlo(2), gthi(2)
//       if (jj .eq. glo(2) .or. jj .eq. ghi(2)) then
//          step = 1
//       else
//          step = gthi(1)-gtlo(1)
//       end if
//       do ii = gtlo(1), gthi(1), step
//          if (ii.eq.glo(1) .or. ii.eq.ghi(1) .or. step .eq. 1) then
//             Ax(ii,jj) = x(ii-1,jj-1)*fxy*sig(ii-1,jj-1) &
//                  +      x(ii+1,jj-1)*fxy*sig(ii  ,jj-1) &
//                  +      x(ii-1,jj+1)*fxy*sig(ii-1,jj  ) &
//                  +      x(ii+1,jj+1)*fxy*sig(ii  ,jj  ) &
//                  +      x(ii-1,jj)*f2xmy*(sig(ii-1,jj-1)+sig(ii-1,jj  )) &
//                  +      x(ii+1,jj)*f2xmy*(sig(ii  ,jj-1)+sig(ii  ,jj  )) &
//                  +      x(ii,jj-1)*fmx2y*(sig(ii-1,jj-1)+sig(ii  ,jj-1)) &
//                  +      x(ii,jj+1)*fmx2y*(sig(ii-1,jj  )+sig(ii  ,jj  )) &
//                  +      x(ii,jj)*(-2.d0)*fxy*(sig(ii-1,jj-1)+sig(ii,jj-1)+sig(ii-1,jj)+sig(ii,jj))
//
//             if (is_rz) then
//                fp = facy / (2*ii+1)
//                fm = facy / (2*ii-1)
//                Ax(ii,jj) = Ax(ii,jj) + (fm*sig(ii-1,jj  )-fp*sig(ii,jj  ))*(x(ii,jj+1)-x(ii,jj)) &
//                     &                + (fm*sig(ii-1,jj-1)-fp*sig(ii,jj-1))*(x(ii,jj-1)-x(ii,jj))
//             end if
//          end if
//       end do
//    end do
//
//    do j = clo(2), chi(2)
//       jj = 2*j
//       if (j .eq. cglo(2) .or. j .eq. cghi(2)) then
//          step = 1
//       else
//          step = chi(1)-clo(1)
//       end if
//       do i = clo(1), chi(1), step
//          ii = 2*i
//          if (msk(ii,jj) .eq. dirichlet) then
//             f(i,j) = f(i,j) + rfd*(Ax(ii,jj) &
//                  + chip*(Ax(ii-1,jj)+Ax(ii+1,jj)+Ax(ii,jj-1)+Ax(ii,jj+1)) &
//                  + chip2*(Ax(ii-1,jj-1)+Ax(ii+1,jj-1)+Ax(ii-1,jj+1)+Ax(ii+1,jj+1)))
//          end if
//       end do
//    end do



}



}
}
