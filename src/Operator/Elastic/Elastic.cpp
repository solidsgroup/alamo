#include <AMReX_MultiFabUtil.H>
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>
#include "Set/Set.H"
#include "Util/Color.H"
#include <AMReX_ArrayLim.H>

#include "Set/Set.H"
#include "Elastic.H"

#define TRACER	std::cout << Color::FG::Yellow << __FILE__ << ":" << __LINE__ << Color::FG::Default << " " << __func__ << std::endl;

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

void
Elastic::define (const Vector<Geometry>& a_geom,
		 const Vector<BoxArray>& a_grids,
		 const Vector<DistributionMapping>& a_dmap,
		 const LPInfo& a_info,
		 const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
	Operator::define(a_geom,a_grids,a_dmap,a_info,a_factory);

	model = new Model::Solid::Elastic::Elastic(2.6, 6.0);

	// DEFINE Cijkl ( to be replaced with RegisterNewFab  eventually)

	m_a_coeffs.resize(m_a_coeffs.size() + 1);
	m_a_coeffs[m_num_a_fabs].resize(m_num_amr_levels);
	for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
	{
		m_a_coeffs[m_num_a_fabs][amrlev].resize(m_num_mg_levels[amrlev]);
		for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
		{
			m_a_coeffs[m_num_a_fabs][amrlev][mglev].define(m_grids[amrlev][mglev],
								       m_dmap[amrlev][mglev],
								       model->NComp(),
								       0);
			m_a_coeffs[m_num_a_fabs][amrlev][mglev].setVal(0.0);

			for (MFIter mfi(m_a_coeffs[m_num_a_fabs][amrlev][mglev], true); mfi.isValid(); ++mfi)
			{
				const Box& bx = mfi.tilebox();
				amrex::FArrayBox       &C    = m_a_coeffs[m_num_a_fabs][amrlev][mglev][mfi];

				AMREX_D_TERM(for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++),
					     for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++),
					     for (int m3 = bx.loVect()[2]; m3<=bx.hiVect()[2]; m3++))
				{
					amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));
					for (int i=0; i<AMREX_SPACEDIM; i++)
						for (int j=0; j<AMREX_SPACEDIM; j++)
							for (int k=0; k<AMREX_SPACEDIM; k++)
								for (int l=0; l<AMREX_SPACEDIM; l++)
									if ((*model)(i,j,k,l) >= 0)
										C(m,(*model)(i,j,k,l)) = (*model).C(i,j,k,l);
				}
			}
			
		}
		// amrex::MultiFab::Copy(m_a_coeffs[m_num_a_fabs][amrlev][0],
		// 		      input[amrlev], 0, 0,
		// 		      input[amrlev].nComp(),
		// 		      input[amrlev].nGrow());
	}
	m_num_a_fabs++;


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
			bool	xmin = (m1 == domain.loVect()[0]),
				xmax = (m1 == domain.hiVect()[0] + 1),
				ymin = (m2 == domain.loVect()[1]),
				ymax = (m2 == domain.hiVect()[1] + 1),
				zmin = (m3 == domain.loVect()[2]),
				zmax = (m3 == domain.hiVect()[2] + 1);

			for (int i=0; i<AMREX_SPACEDIM; i++)
			{
				ffab(amrex::IntVect(AMREX_D_DECL(m1,m2,m3)),i) = 0.0;

				for (int k=0; k<AMREX_SPACEDIM; k++)
				{
					Set::Vector gradu_k; // gradu_k(l) = u_{k,l}
					AMREX_D_TERM(gradu_k(0) = ((!xmax ? ufab(m+dx,k) : ufab(m,k)) - (!xmin ? ufab(m-dx,k) : ufab(m,k)))/((xmin || xmax ? 1.0 : 2.0)*DX[0]);,
						     gradu_k(1) = ((!ymax ? ufab(m+dy,k) : ufab(m,k)) - (!ymin ? ufab(m-dy,k) : ufab(m,k)))/((ymin || ymax ? 1.0 : 2.0)*DX[1]);,
						     gradu_k(2) = ((!zmax ? ufab(m+dz,k) : ufab(m,k)) - (!zmin ? ufab(m-dz,k) : ufab(m,k)))/((zmin || zmax ? 1.0 : 2.0)*DX[2]););

					if (xmin)
					{
						if (m_bc_xlo[k] == BC::Displacement)
							ffab(m,k) = ufab(m,k);
						else if (m_bc_xlo[k] == BC::Traction) 
							for (int l=0; l<AMREX_SPACEDIM; l++)
								ffab(m,i) -= C(i,0,k,l,m,amrlev,mglev,mfi) * gradu_k(l);
						else Util::Abort("Invalid BC");
					}
					if (xmax)
					{
						if (m_bc_xhi[k] == BC::Displacement)
							ffab(m,k) = ufab(m,k);
						else if (m_bc_xlo[k] == BC::Traction) 
							for (int l=0; l<AMREX_SPACEDIM; l++)
								ffab(m,i) += C(i,0,k,l,m,amrlev,mglev,mfi) * gradu_k(l);
						else Util::Abort("Invalid BC");
					}

					if (ymin)
					{
						if (m_bc_ylo[k] == BC::Displacement)
							ffab(m,k) = ufab(m,k);
						else if (m_bc_ylo[k] == BC::Traction) 
							for (int l=0; l<AMREX_SPACEDIM; l++)
								ffab(m,i) -= C(i,1,k,l,m,amrlev,mglev,mfi) * gradu_k(l);
						else Util::Abort("Invalid BC");
					}
					if (ymax)
					{
						if (m_bc_yhi[k] == BC::Displacement)
							ffab(m,k) = ufab(m,k);
						else if (m_bc_ylo[k] == BC::Traction) 
							for (int l=0; l<AMREX_SPACEDIM; l++)
								ffab(m,i) += C(i,1,k,l,m,amrlev,mglev,mfi) * gradu_k(l);
						else Util::Abort("Invalid BC");
					}

					if (zmin)
					{
						if (m_bc_zlo[k] == BC::Displacement)
							ffab(m,k) = ufab(m,k);
						else if (m_bc_zlo[k] == BC::Traction) 
							for (int l=0; l<AMREX_SPACEDIM; l++)
								ffab(m,i) -= C(i,2,k,l,m,amrlev,mglev,mfi) * gradu_k(l);
						else Util::Abort("Invalid BC");
					}
					if (zmax)
					{
						if (m_bc_zhi[k] == BC::Displacement)
							ffab(m,k) = ufab(m,k);
						else if (m_bc_zlo[k] == BC::Traction) 
							for (int l=0; l<AMREX_SPACEDIM; l++)
								ffab(m,i) += C(i,2,k,l,m,amrlev,mglev,mfi) * gradu_k(l);
						else Util::Abort("Invalid BC");
					}
					if (xmin || xmax || ymin || ymax || zmin || zmax) continue;
					





					Set::Matrix gradgradu_k; // gradgradu_k(l,j) = u_{k,lj}
					AMREX_D_TERM(gradgradu_k(0,0) = (ufab(m+dx,k) - 2.0*ufab(m,k) + ufab(m-dx,k))/DX[0]/DX[0];
						     ,// 2D
						     gradgradu_k(0,1) = (ufab(m+dx+dy,k) + ufab(m-dx-dy,k) - ufab(m+dx-dy,k) - ufab(m-dx+dy,k))/(2.0*DX[0])/(2.0*DX[1]);
						     gradgradu_k(1,0) = gradgradu_k(0,1);
						     gradgradu_k(1,1) = (ufab(m+dy,k) - 2.0*ufab(m,k) + ufab(m-dy,k))/DX[1]/DX[1];
						     ,// 3D
						     gradgradu_k(0,2) = (ufab(m+dx+dz,k) + ufab(m-dx-dz,k) - ufab(m+dx-dz,k) - ufab(m-dx+dz,k))/(2.0*DX[0])/(2.0*DX[2]);
						     gradgradu_k(1,2) = (ufab(m+dy+dz,k) + ufab(m-dy-dz,k) - ufab(m+dy-dz,k) - ufab(m-dy+dz,k))/(2.0*DX[1])/(2.0*DX[2]);
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
							
					Set::Scalar diag = 0.0;

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


void
Elastic::Fsmooth (int amrlev,
		  int mglev,
		  MultiFab& u,
		  const MultiFab& rhs
		  ) const
{
	BL_PROFILE("Operator::Elastic::Elastic::Fsmooth()");

	for (int redblack = 0; redblack < 2; redblack++)
	{
		amrex::Box domain(m_geom[amrlev][mglev].Domain());
		const Real* DX = m_geom[amrlev][mglev].CellSize();

		static amrex::IntVect AMREX_D_DECL(dx(AMREX_D_DECL(1,0,0)),
						   dy(AMREX_D_DECL(0,1,0)),
						   dz(AMREX_D_DECL(0,0,1)));

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

				bool	xmin = (m1 == domain.loVect()[0]),
					xmax = (m1 == domain.hiVect()[0] + 1),
					ymin = (m2 == domain.loVect()[1]),
					ymax = (m2 == domain.hiVect()[1] + 1),
					zmin = (m3 == domain.loVect()[2]),
					zmax = (m3 == domain.hiVect()[2] + 1);


				for (int i=0; i<AMREX_SPACEDIM; i++)
				{
					amrex::Real rho = 0.0, aa = 0.0;
					for (int k=0; k<AMREX_SPACEDIM; k++)
					{
						Set::Vector OffDiag_gradu_k; // gradu_k(l) = u_{k,l}
						AMREX_D_TERM(OffDiag_gradu_k(0) = ((!xmax ? ufab(m+dx,k) : (i==k? 0.0 : ufab(m,k))) - (!xmin ? ufab(m-dx,k) : (i==k ? 0.0 : ufab(m,k))))/((xmin || xmax ? 1.0 : 2.0)*DX[0]);,
							     OffDiag_gradu_k(1) = ((!ymax ? ufab(m+dy,k) : (i==k? 0.0 : ufab(m,k))) - (!ymin ? ufab(m-dy,k) : (i==k ? 0.0 : ufab(m,k))))/((ymin || ymax ? 1.0 : 2.0)*DX[1]);,
							     OffDiag_gradu_k(2) = ((!zmax ? ufab(m+dz,k) : (i==k? 0.0 : ufab(m,k))) - (!zmin ? ufab(m-dz,k) : (i==k ? 0.0 : ufab(m,k))))/((zmin || zmax ? 1.0 : 2.0)*DX[2]););
						Set::Vector Diag_gradu_k; // gradu_k(l) = u_{k,l}
						AMREX_D_TERM(Diag_gradu_k(0) = ((!xmax ? 0.0 : (i==k? 1.0 : 0.0)) - (!xmin ? 0.0 : (i==k ? 1.0 : 0.0)))/((xmin || xmax ? 1.0 : 2.0)*DX[0]);,
							     Diag_gradu_k(1) = ((!ymax ? 0.0 : (i==k? 1.0 : 0.0)) - (!ymin ? 0.0 : (i==k ? 1.0 : 0.0)))/((ymin || ymax ? 1.0 : 2.0)*DX[1]);,
							     Diag_gradu_k(2) = ((!zmax ? 0.0 : (i==k? 1.0 : 0.0)) - (!zmin ? 0.0 : (i==k ? 1.0 : 0.0)))/((zmin || zmax ? 1.0 : 2.0)*DX[2]););

						if (xmin)
						{
							if (m_bc_xlo[k] == BC::Displacement)
							{rho = 0.0; aa = 1.0;}
							else if (m_bc_xlo[k] == BC::Traction) 
								for (int l=0; l<AMREX_SPACEDIM; l++)
								{
									rho -= C(i,0,k,l,m,amrlev,mglev,mfi) * OffDiag_gradu_k(l);
									aa  -= C(i,0,k,l,m,amrlev,mglev,mfi) * Diag_gradu_k(l);
								}
							else Util::Abort("Invalid BC");
						}
						if (xmax)
						{
							if (m_bc_xhi[k] == BC::Displacement)
							{rho = 0.0; aa = 1.0;}
							else if (m_bc_xlo[k] == BC::Traction) 
								for (int l=0; l<AMREX_SPACEDIM; l++)
								{
									rho += C(i,0,k,l,m,amrlev,mglev,mfi) * OffDiag_gradu_k(l);
									aa  += C(i,0,k,l,m,amrlev,mglev,mfi) * Diag_gradu_k(l);
								}
							else Util::Abort("Invalid BC");
						}

						if (ymin)
						{
							if (m_bc_ylo[k] == BC::Displacement)
							{rho = 0.0; aa = 1.0;}
							else if (m_bc_ylo[k] == BC::Traction) 
								for (int l=0; l<AMREX_SPACEDIM; l++)
								{
									rho -= C(i,1,k,l,m,amrlev,mglev,mfi) * OffDiag_gradu_k(l);
									aa  -= C(i,1,k,l,m,amrlev,mglev,mfi) * Diag_gradu_k(l);
								}
							else Util::Abort("Invalid BC");
						}
						if (ymax)
						{
							if (m_bc_yhi[k] == BC::Displacement)
							{rho = 0.0; aa = 1.0;}
							else if (m_bc_ylo[k] == BC::Traction) 
								for (int l=0; l<AMREX_SPACEDIM; l++)
								{
									rho += C(i,1,k,l,m,amrlev,mglev,mfi) * OffDiag_gradu_k(l);
									aa  += C(i,1,k,l,m,amrlev,mglev,mfi) * Diag_gradu_k(l);
								}
							else Util::Abort("Invalid BC");
						}

						if (zmin)
						{
							if (m_bc_zlo[k] == BC::Displacement)
							{rho = 0.0; aa = 1.0;}
							else if (m_bc_zlo[k] == BC::Traction) 
								for (int l=0; l<AMREX_SPACEDIM; l++)
								{
									rho -= C(i,2,k,l,m,amrlev,mglev,mfi) * OffDiag_gradu_k(l);
									aa  -= C(i,2,k,l,m,amrlev,mglev,mfi) * Diag_gradu_k(l);
								}
							else Util::Abort("Invalid BC");
						}
						if (zmax)
						{
							if (m_bc_zhi[k] == BC::Displacement)
							{rho = 0.0; aa = 1.0;}
							else if (m_bc_zlo[k] == BC::Traction) 
								for (int l=0; l<AMREX_SPACEDIM; l++)
								{
									rho += C(i,2,k,l,m,amrlev,mglev,mfi) * OffDiag_gradu_k(l);
									aa  += C(i,2,k,l,m,amrlev,mglev,mfi) * Diag_gradu_k(l);
								}
							else Util::Abort("Invalid BC");
						}
						if (xmin || xmax || ymin || ymax || zmin || zmax) continue;



					
						Set::Matrix OffDiag_gradgradu_k; // gradgradu_k(l,j) = u_{k,lj}
						AMREX_D_TERM(OffDiag_gradgradu_k(0,0) = (ufab(m+dx,k) - (i==k ? 0.0 : 2.0*ufab(m,k)) + ufab(m-dx,k))/DX[0]/DX[0];
							     ,// 2D
							     OffDiag_gradgradu_k(0,1) = (ufab(m+dx+dy,k) + ufab(m-dx-dy,k) - ufab(m+dx-dy,k) - ufab(m-dx+dy,k))/(2.0*DX[0])/(2.0*DX[1]);
							     OffDiag_gradgradu_k(1,0) = OffDiag_gradgradu_k(0,1);
							     OffDiag_gradgradu_k(1,1) = (ufab(m+dy,k) - (i==k ? 0.0 : 2.0*ufab(m,k)) + ufab(m-dy,k))/DX[1]/DX[1];
							     ,// 3D
							     OffDiag_gradgradu_k(0,2) = (ufab(m+dx+dz,k) + ufab(m-dx-dz,k) - ufab(m+dx-dz,k) - ufab(m-dx+dz,k))/(2.0*DX[0])/(2.0*DX[2]);
							     OffDiag_gradgradu_k(1,2) = (ufab(m+dy+dz,k) + ufab(m-dy-dz,k) - ufab(m+dy-dz,k) - ufab(m-dy+dz,k))/(2.0*DX[1])/(2.0*DX[2]);
							     OffDiag_gradgradu_k(2,0) = OffDiag_gradgradu_k(0,2);
							     OffDiag_gradgradu_k(2,1) = OffDiag_gradgradu_k(1,2);
							     OffDiag_gradgradu_k(2,2) = (ufab(m+dz,k) - (i==k ? 0.0 : 2.0*ufab(m,k)) + ufab(m-dz,k))/DX[2]/DX[2];);

						Set::Matrix Diag_gradgradu_k; // gradgradu_k(l,j) = u_{k,lj}
						AMREX_D_TERM(Diag_gradgradu_k(0,0) = (0.0 - (i==k ? 2.0 : 0.0) + 0.0)/DX[0]/DX[0];
							     ,// 2D
							     Diag_gradgradu_k(0,1) = 0.0;
							     Diag_gradgradu_k(1,0) = 0.0;
							     Diag_gradgradu_k(1,1) = (0.0 - (i==k ? 2.0 : 0.0) + 0.0)/DX[1]/DX[1];
							     ,// 3D
							     Diag_gradgradu_k(0,2) = 0.0;
							     Diag_gradgradu_k(1,2) = 0.0;
							     Diag_gradgradu_k(2,0) = 0.0;
							     Diag_gradgradu_k(2,1) = 0.0;
							     Diag_gradgradu_k(2,2) = (0.0 - (i==k ? 2.0 : 0.0) + 0.0)/DX[2]/DX[2];);

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
								rho -= C(i,j,k,l,m,amrlev,mglev,mfi) * (OffDiag_gradgradu_k(j,l));
							for (int j=0; j<AMREX_SPACEDIM; j++)
								aa  -= C(i,j,k,l,m,amrlev,mglev,mfi) * (Diag_gradgradu_k(j,l));

							// f_i -= C_{ijkl,j} u_{k,l}
							rho -= C_ik(l) * OffDiag_gradu_k(l);
							aa  -= C_ik(l) * Diag_gradu_k(l);
						}


						// ufab(m,i) = (rhsfab(m,i) - rho) / aa;
					}


					// std::cout << "u = " << ufab(m,i) << "  rho = " << rho << " aa = " << aa << std::endl;
					// if (fabs(aa) < 1E-8)
					// 	Util::Abort("Defective aa in fsmooth");


					ufab(m,i) = (rhsfab(m,i) - rho) / aa;
				}
			}
		}
	}
}

void
Elastic::normalize (int amrlev, int mglev, MultiFab& mf) const
{
	//	return;
	BL_PROFILE("Operator::Elastic::Elastic::Fapply()");
	amrex::Box domain(m_geom[amrlev][mglev].Domain());
	const Real* DX = m_geom[amrlev][mglev].CellSize();
	
	static amrex::IntVect AMREX_D_DECL(dx(AMREX_D_DECL(1,0,0)),
					   dy(AMREX_D_DECL(0,1,0)),
					   dz(AMREX_D_DECL(0,0,1)));

	for (MFIter mfi(mf, true); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();
		amrex::FArrayBox &mffab    = mf[mfi];

		AMREX_D_TERM(for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++),
			     for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++),
			     for (int m3 = bx.loVect()[2]; m3<=bx.hiVect()[2]; m3++))
		{
			amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));
			bool	xmin = (m1 == domain.loVect()[0]),
				xmax = (m1 == domain.hiVect()[0] + 1),
				ymin = (m2 == domain.loVect()[1]),
				ymax = (m2 == domain.hiVect()[1] + 1),
				zmin = (m3 == domain.loVect()[2]),
				zmax = (m3 == domain.hiVect()[2] + 1);

			for (int i=0; i<AMREX_SPACEDIM; i++)
			{
				Set::Scalar aa = 0.0;

				for (int k=0; k<AMREX_SPACEDIM; k++)
				{
					Set::Vector gradu_k; // gradu_k(l) = u_{k,l}
					AMREX_D_TERM(gradu_k(0) = ((!xmax ? 0.0/*ufab(m+dx,k)*/ : (i==k ? 1.0 : 0.0)/*ufab(m,k)*/) - (!xmin ? 0.0/*ufab(m-dx,k)*/ : (i==k ? 1.0 : 0.0)/*ufab(m,k)*/))/((xmin || xmax ? 1.0 : 2.0)*DX[0]);,
						     gradu_k(1) = ((!ymax ? 0.0/*ufab(m+dy,k)*/ : (i==k ? 1.0 : 0.0)/*ufab(m,k)*/) - (!ymin ? 0.0/*ufab(m-dy,k)*/ : (i==k ? 1.0 : 0.0)/*ufab(m,k)*/))/((ymin || ymax ? 1.0 : 2.0)*DX[1]);,
						     gradu_k(2) = ((!zmax ? 0.0/*ufab(m+dz,k)*/ : (i==k ? 1.0 : 0.0)/*ufab(m,k)*/) - (!zmin ? 0.0/*ufab(m-dz,k)*/ : (i==k ? 1.0 : 0.0)/*ufab(m,k)*/))/((zmin || zmax ? 1.0 : 2.0)*DX[2]););


					if (xmin)
					{
						if (m_bc_xlo[k] == BC::Displacement)
							aa = 1.0;
						else if (m_bc_xlo[k] == BC::Traction) 
							for (int l=0; l<AMREX_SPACEDIM; l++)
								aa  -= C(i,0,k,l,m,amrlev,mglev,mfi) * gradu_k(l);
					}
					if (xmax)
					{
						if (m_bc_xhi[k] == BC::Displacement)
							aa = 1.0;
						else if (m_bc_xlo[k] == BC::Traction) 
							for (int l=0; l<AMREX_SPACEDIM; l++)
								aa  += C(i,0,k,l,m,amrlev,mglev,mfi) * gradu_k(l);

					}

					if (ymin)
					{
						if (m_bc_ylo[k] == BC::Displacement)
							aa = 1.0;
						else if (m_bc_ylo[k] == BC::Traction) 
							for (int l=0; l<AMREX_SPACEDIM; l++)
								aa  -= C(i,1,k,l,m,amrlev,mglev,mfi) * gradu_k(l);
					}
					if (ymax)
					{
						if (m_bc_yhi[k] == BC::Displacement)
							aa = 1.0;
						else if (m_bc_ylo[k] == BC::Traction) 
							for (int l=0; l<AMREX_SPACEDIM; l++)
								aa  += C(i,1,k,l,m,amrlev,mglev,mfi) * gradu_k(l);

					}

					if (zmin)
					{
						if (m_bc_zlo[k] == BC::Displacement)
							aa = 1.0;
						else if (m_bc_zlo[k] == BC::Traction) 
							for (int l=0; l<AMREX_SPACEDIM; l++)
								aa  -= C(i,2,k,l,m,amrlev,mglev,mfi) * gradu_k(l);

					}
					if (zmax)
					{
						if (m_bc_zhi[k] == BC::Displacement)
							aa = 1.0;
						else if (m_bc_zlo[k] == BC::Traction) 
							for (int l=0; l<AMREX_SPACEDIM; l++)
								aa  += C(i,2,k,l,m,amrlev,mglev,mfi) * gradu_k(l);
					}


					if (xmin || xmax || ymin || ymax || zmin || zmax) continue;


					Set::Matrix gradgradu_k; // gradgradu_k(l,j) = u_{k,lj}
					AMREX_D_TERM(gradgradu_k(0,0) = (/*ufab(m+dx,k)*/ - (i==k ? 2.0/**ufab(m,k)*/ : 0) /* + ufab(m-dx,k)*/)/DX[0]/DX[0];
						     ,// 2D
						     gradgradu_k(0,1) = 0.0 /*(ufab(m+dx+dy,k) + ufab(m-dx-dy,k) - ufab(m+dx-dy,k) - ufab(m-dx+dy,k))/(2.0*DX[0])/(2.0*DX[1])*/;
						     gradgradu_k(1,0) = 0.0 /*gradgradu_k(0,1)*/;
						     gradgradu_k(1,1) = (/*ufab(m+dy,k)*/ - (i==k ? 2.0/**ufab(m,k)*/ : 0) /*+ ufab(m-dy,k)*/)/DX[1]/DX[1];
						     ,// 3D
						     gradgradu_k(0,2) = 0.0 /*(ufab(m+dx+dz,k) + ufab(m-dx-dz,k) - ufab(m+dx-dz,k) - ufab(m-dx+dz,k))/(2.0*DX[0])/(2.0*DX[2])*/;
						     gradgradu_k(1,2) = 0.0 /*(ufab(m+dy+dz,k) + ufab(m-dy-dz,k) - ufab(m+dy-dz,k) - ufab(m-dy+dz,k))/(2.0*DX[1])/(2.0*DX[2])*/;
						     gradgradu_k(2,0) = 0.0 /*gradgradu_k(0,2)*/;
						     gradgradu_k(2,1) = 0.0 /*gradgradu_k(1,2)*/;
						     gradgradu_k(2,2) = (/*ufab(m+dz,k)*/ - (i==k ? 2.0/**ufab(m,k)*/ : 0) /*+ ufab(m-dz,k)*/)/DX[2]/DX[2];);

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
							aa -= C(i,j,k,l,m,amrlev,mglev,mfi) * (gradgradu_k(j,l));

						// f_i -= C_{ijkl,j} u_{k,l}
						if (i==l)
							aa -= C_ik(l) * gradu_k(l);
					}
				}
				mffab(m,i) = mffab(m,i) / aa;
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


void
Elastic::Stress (FArrayBox& sigmafab,
		 const FArrayBox& ufab,
		 int amrlev, const MFIter& mfi,
		 bool voigt) const
{
	amrex::Box domain(m_geom[amrlev][0].Domain());
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

				bool	xmin = (m1 == domain.loVect()[0]),
					xmax = (m1 == domain.hiVect()[0] + 1),
					ymin = (m2 == domain.loVect()[1]),
					ymax = (m2 == domain.hiVect()[1] + 1),
					zmin = (m3 == domain.loVect()[2]),
					zmax = (m3 == domain.hiVect()[2] + 1);



				Set::Matrix gradu;

				Set::Vector gradu_k; // gradu_k(l) = u_{k,l}
				AMREX_D_TERM(gradu(0,0) = ((!xmax ? ufab(m+dx,0) : ufab(m,0)) - (!xmin ? ufab(m-dx,0) : ufab(m,0)))/((xmin || xmax ? 1.0 : 2.0)*DX[0]);
					     ,
					     gradu(0,1) = ((!ymax ? ufab(m+dy,0) : ufab(m,0)) - (!ymin ? ufab(m-dy,0) : ufab(m,0)))/((ymin || ymax ? 1.0 : 2.0)*DX[1]);
					     gradu(1,0) = ((!xmax ? ufab(m+dx,1) : ufab(m,1)) - (!xmin ? ufab(m-dx,1) : ufab(m,1)))/((xmin || xmax ? 1.0 : 2.0)*DX[0]);
					     gradu(1,1) = ((!ymax ? ufab(m+dy,1) : ufab(m,1)) - (!ymin ? ufab(m-dy,1) : ufab(m,1)))/((ymin || ymax ? 1.0 : 2.0)*DX[1]);
					     ,
					     gradu(0,2) = ((!zmax ? ufab(m+dz,0) : ufab(m,0)) - (!zmin ? ufab(m-dz,0) : ufab(m,0)))/((zmin || zmax ? 1.0 : 2.0)*DX[2]);
					     gradu(1,2) = ((!zmax ? ufab(m+dz,1) : ufab(m,1)) - (!zmin ? ufab(m-dz,1) : ufab(m,1)))/((zmin || zmax ? 1.0 : 2.0)*DX[2]);
					     gradu(2,0) = ((!xmax ? ufab(m+dx,2) : ufab(m,2)) - (!xmin ? ufab(m-dx,2) : ufab(m,2)))/((xmin || xmax ? 1.0 : 2.0)*DX[0]);
					     gradu(2,1) = ((!ymax ? ufab(m+dy,2) : ufab(m,2)) - (!ymin ? ufab(m-dy,2) : ufab(m,2)))/((ymin || ymax ? 1.0 : 2.0)*DX[1]);
					     gradu(2,2) = ((!zmax ? ufab(m+dz,2) : ufab(m,2)) - (!zmin ? ufab(m-dz,2) : ufab(m,2)))/((zmin || zmax ? 1.0 : 2.0)*DX[2]);
					     );

				Set::Matrix eps0 = Set::Matrix::Zero();

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
