#include "Model/Solid/Elastic/Isotropic/Isotropic.H"
#include "Elastic.H"

namespace Operator
{
namespace Elastic
{

template<class T>
Elastic<T>::Elastic (const Vector<Geometry>& a_geom,
		     const Vector<BoxArray>& a_grids,
		     const Vector<DistributionMapping>& a_dmap,
		     const LPInfo& a_info)
{
	define(a_geom, a_grids, a_dmap, a_info);
}

template<class T>
Elastic<T>::~Elastic ()
{}

template<class T>
void
Elastic<T>::define (const Vector<Geometry>& a_geom,
		 const Vector<BoxArray>& a_grids,
		 const Vector<DistributionMapping>& a_dmap,
		 const LPInfo& a_info,
		 const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
	Operator::define(a_geom,a_grids,a_dmap,a_info,a_factory);

	model.resize(m_num_amr_levels);
	for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
	{
		model[amrlev].resize(m_num_mg_levels[amrlev]);
		for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
		{
			//\todo replace number of ghost cells with 0
			model[amrlev][mglev].reset(new amrex::FabArray<amrex::BaseFab<T> >(m_grids[amrlev][mglev], m_dmap[amrlev][mglev], 1, 1));
		}
	}
}

template <class T>
void
Elastic<T>::SetModel (int amrlev, const amrex::FabArray<amrex::BaseFab<T> >& a_model)
{
	for (MFIter mfi(a_model, true); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();
		amrex::BaseFab<T> &modelfab = (*(model[amrlev][0]))[mfi];
		const amrex::BaseFab<T> &a_modelfab = a_model[mfi];

		AMREX_D_TERM(for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++),
			     for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++),
			     for (int m3 = bx.loVect()[2]; m3<=bx.hiVect()[2]; m3++))
		{
			amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));
			modelfab(m) = a_modelfab(m);
		}
	}
}

template<class T>
void
Elastic<T>::Fapply (int amrlev, int mglev, MultiFab& f, const MultiFab& u) const
{
	BL_PROFILE("Operator::Elastic::Elastic::Fapply()");

	amrex::Box domain(m_geom[amrlev][mglev].Domain());
	const Real* DX = m_geom[amrlev][mglev].CellSize();
	
	static std::array<amrex::IntVect,AMREX_SPACEDIM> dx = {AMREX_D_DECL(amrex::IntVect(AMREX_D_DECL(1,0,0)),
									    amrex::IntVect(AMREX_D_DECL(0,1,0)),
									    amrex::IntVect(AMREX_D_DECL(0,0,1)))};
	// DEBUG: Check to see if Fapply is getting passed bad values.
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
		amrex::BaseFab<T> &C = (*(model[amrlev][mglev]))[mfi];
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

			// The displacement gradient tensor
			Set::Matrix gradu; // gradu(i,j) = u_{i,j)

			// The gradient of the displacement gradient tensor
			std::array<Set::Matrix,AMREX_SPACEDIM> gradgradu; // gradgradu[k](l,j) = u_{k,lj}

			// Fill gradu and gradgradu
			for (int i = 0; i < AMREX_SPACEDIM; i++)
			{
				// Note: the (?:) block modifies the stencil if the node is on a boundary
				AMREX_D_TERM(gradu(i,0) = ((!xmax ? ufab(m+dx[0],i) : ufab(m,i)) - (!xmin ? ufab(m-dx[0],i) : ufab(m,i)))/((xmin || xmax ? 1.0 : 2.0)*DX[0]);,
					     gradu(i,1) = ((!ymax ? ufab(m+dx[1],i) : ufab(m,i)) - (!ymin ? ufab(m-dx[1],i) : ufab(m,i)))/((ymin || ymax ? 1.0 : 2.0)*DX[1]);,
					     gradu(i,2) = ((!zmax ? ufab(m+dx[2],i) : ufab(m,i)) - (!zmin ? ufab(m-dx[2],i) : ufab(m,i)))/((zmin || zmax ? 1.0 : 2.0)*DX[2]););
			
				AMREX_D_TERM(gradgradu[i](0,0) = (ufab(m+dx[0],i) - 2.0*ufab(m,i) + ufab(m-dx[0],i))/DX[0]/DX[0];
					     ,// 2D
					     gradgradu[i](0,1) = (ufab(m+dx[0]+dx[1],i) + ufab(m-dx[0]-dx[1],i) - ufab(m+dx[0]-dx[1],i) - ufab(m-dx[0]+dx[1],i))/(2.0*DX[0])/(2.0*DX[1]);
					     gradgradu[i](1,0) = gradgradu[i](0,1);
					     gradgradu[i](1,1) = (ufab(m+dx[1],i) - 2.0*ufab(m,i) + ufab(m-dx[1],i))/DX[1]/DX[1];
					     ,// 3D
					     gradgradu[i](0,2) = (ufab(m+dx[0]+dx[2],i) + ufab(m-dx[0]-dx[2],i) - ufab(m+dx[0]-dx[2],i) - ufab(m-dx[0]+dx[2],i))/(2.0*DX[0])/(2.0*DX[2]);
					     gradgradu[i](1,2) = (ufab(m+dx[1]+dx[2],i) + ufab(m-dx[1]-dx[2],i) - ufab(m+dx[1]-dx[2],i) - ufab(m-dx[1]+dx[2],i))/(2.0*DX[1])/(2.0*DX[2]);
					     gradgradu[i](2,0) = gradgradu[i](0,2);
					     gradgradu[i](2,1) = gradgradu[i](1,2);
					     gradgradu[i](2,2) = (ufab(m+dx[2],i) - 2.0*ufab(m,i) + ufab(m-dx[2],i))/DX[2]/DX[2];);
			}

			// Strain tensor
			Set::Matrix eps = 0.5*(gradu + gradu.transpose());
			// Stress tensor computed using the model fab
			Set::Matrix sig = C(m)(eps);

			//
			// Boundary conditions
			//
			// BCs are implemented as boundary operators.
			//
			// ┌                      ┐ ┌                     ┐   ┌              ┐
			// │              |       │ │ interior		  │   │ body	     │
			// │  Div C Grad  │       │ │ displacements	  │   │ forces	     │
			// │              │       │ │			  │ = │		     │
			// │ ─────────────┼────── │ │ ──────────────────  │   │ ──────────── │
			// │              │ Bndry │ │ bndry displacements │   │ bndry values │
			// └                      ┘ └			  ┘   └		     ┘
			//
			// For displacement:
			//   (Bndry)(u) = u
			// For traction:
			//   (Bndry)(u) = C Grad (u) n   (n is surface normal)
			//
			// The displacement values or traction values are set as the boundary
			// values of the rhs fab.
			//
			for (int i = 0; i < AMREX_SPACEDIM; i++) // iterate over DIMENSIONS
			{
				for (int j = 0; j < AMREX_SPACEDIM; j++) // iterate over FACES
				{
					if (m[j] == domain.loVect()[j])
						if (m_bc_lo[j][i] == BC::Displacement)
							ffab(m,i) = ufab(m,i);
						else if (m_bc_lo[j][i] == BC::Traction) 
							ffab(m,i) = -sig(i,j);
						else Util::Abort("Invalid BC");
					if (m[j] == domain.hiVect()[j] + 1)
						if (m_bc_hi[j][i] == BC::Displacement)
							ffab(m,i) = ufab(m,i);
						else if (m_bc_hi[j][i] == BC::Traction) 
							ffab(m,i) = +sig(i,j);
						else Util::Abort("Invalid BC");
				}
			}
			if (xmax || xmin || ymax || ymin || zmax || zmin) continue;
			
			//
			// Operator
			//
			// The return value is
			//    f = C(grad grad u) + grad(C)*grad(u)
			// In index notation
			//    f_i = C_{ijkl,j} u_{k,l}  +  C_{ijkl}u_{k,lj}
			//
			Set::Vector f =
				C(m)(gradgradu) + 
				AMREX_D_TERM(((C(m+dx[0]) - C(m-dx[0]))/2.0/DX[0])(gradu).col(0),
				 	     + ((C(m+dx[1]) - C(m-dx[1]))/2.0/DX[1])(gradu).col(1),
				  	     + ((C(m+dx[2]) - C(m-dx[2]))/2.0/DX[2])(gradu).col(2));
			for (int i = 0; i < AMREX_SPACEDIM; i++)
				ffab(m,i) = f(i);
		}
	}
}


template<class T>
void
Elastic<T>::Fsmooth (int amrlev,
		  int mglev,
		  MultiFab& u,
		  const MultiFab& rhs
		  ) const
{
	Util::Abort("Fsmooth is under construction and not fully implemented");
	/*
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

				bool min[AMREX_SPACEDIM], max[AMREX_SPACEDIM];

				min[0] = (m1 == domain.loVect()[0]);
				max[0] = (m1 == domain.hiVect()[0] + 1);
				min[1] = (m2 == domain.loVect()[1]);
				max[1] = (m2 == domain.hiVect()[1] + 1);
				min[2] = (m3 == domain.loVect()[2]);
				max[2] = (m3 == domain.hiVect()[2] + 1);


				Set::Matrix Dgradu, ODgradu;// gradu(i,j) = u_{i,j)
				std::array<Set::Matrix,AMREX_SPACEDIM> Dgradgradu, ODgradgradu; // gradgradu[k](l,j) = u_{k,lj}
				for (int i = 0; i < AMREX_SPACEDIM; i++)
				{
					for (int j = 0; j < AMREX_SPACEDIM; j++)
					{
						ODgradu(i,j) = ((!max[j] ? ufab(m+dx[0],i) : (i==j ? 0.0 : ufab(m,i))) - (!min[j] ? ufab(m-dx[0],i) : (i==j ? 0.0 : ufab(m,i))))
							/((min[j] || max[j] ? 1.0 : 2.0)*DX[0]);

						Dgradu(i,j) = ((!max[j] ? 0.0 : (i==j ? 1.0 : 0.0)) - (!min[j] ? 0.0 : (i==j ? 1.0 : 0.0)))
							/((min[j] || max[j] ? 1.0 : 2.0)*DX[0]);
					}
					
					AMREX_D_TERM(ODgradgradu[k](0,0) = (ufab(m+dx,k) - (i==k ? 0.0 : 2.0*ufab(m,k)) + ufab(m-dx,k))/DX[0]/DX[0];
						     ,// 2D
						     ODgradgradu[k](0,1) = (ufab(m+dx+dy,k) + ufab(m-dx-dy,k) - ufab(m+dx-dy,k) - ufab(m-dx+dy,k))/(2.0*DX[0])/(2.0*DX[1]);
						     ODgradgradu[k](1,0) = OffDiag_gradgradu[k](0,1);
						     ODgradgradu[k](1,1) = (ufab(m+dy,k) - (i==k ? 0.0 : 2.0*ufab(m,k)) + ufab(m-dy,k))/DX[1]/DX[1];
						     ,// 3D
						     ODgradgradu[k](0,2) = (ufab(m+dx+dz,k) + ufab(m-dx-dz,k) - ufab(m+dx-dz,k) - ufab(m-dx+dz,k))/(2.0*DX[0])/(2.0*DX[2]);
						     ODgradgradu[k](1,2) = (ufab(m+dy+dz,k) + ufab(m-dy-dz,k) - ufab(m+dy-dz,k) - ufab(m-dy+dz,k))/(2.0*DX[1])/(2.0*DX[2]);
						     ODgradgradu[k](2,0) = OffDiag_gradgradu[k](0,2);
						     ODgradgradu[k](2,1) = OffDiag_gradgradu[k](1,2);
						     ODgradgradu[k](2,2) = (ufab(m+dz,k) - (i==k ? 0.0 : 2.0*ufab(m,k)) + ufab(m-dz,k))/DX[2]/DX[2];);
					

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


				}


				for (int i=0; i<AMREX_SPACEDIM; i++)
				{



					amrex::Real rho = 0.0, aa = 0.0;
					for (int k=0; k<AMREX_SPACEDIM; k++)
					{
						Set::Vector OffDiag_gradu_k; // gradu_k(l) = u_{k,l}
						AMREX_D_TERM(OffDiag_gradu_k(0) = ((!max[0] ? ufab(m+dx,k) : (i==k? 0.0 : ufab(m,k))) - (!min[0] ? ufab(m-dx,k) : (i==k ? 0.0 : ufab(m,k))))/((min[0] || max[0] ? 1.0 : 2.0)*DX[0]);,
							     OffDiag_gradu_k(1) = ((!max[1] ? ufab(m+dy,k) : (i==k? 0.0 : ufab(m,k))) - (!min[1] ? ufab(m-dy,k) : (i==k ? 0.0 : ufab(m,k))))/((min[1] || max[1] ? 1.0 : 2.0)*DX[1]);,
							     OffDiag_gradu_k(2) = ((!max[2] ? ufab(m+dz,k) : (i==k? 0.0 : ufab(m,k))) - (!min[2] ? ufab(m-dz,k) : (i==k ? 0.0 : ufab(m,k))))/((min[2] || max[2] ? 1.0 : 2.0)*DX[2]););
						Set::Vector Diag_gradu_k; // gradu_k(l) = u_{k,l}
						AMREX_D_TERM(Diag_gradu_k(0) = ((!max[0] ? 0.0 : (i==k? 1.0 : 0.0)) - (!min[0] ? 0.0 : (i==k ? 1.0 : 0.0)))/((min[0] || max[0] ? 1.0 : 2.0)*DX[0]);,
							     Diag_gradu_k(1) = ((!max[1] ? 0.0 : (i==k? 1.0 : 0.0)) - (!min[1] ? 0.0 : (i==k ? 1.0 : 0.0)))/((min[1] || max[1] ? 1.0 : 2.0)*DX[1]);,
							     Diag_gradu_k(2) = ((!max[2] ? 0.0 : (i==k? 1.0 : 0.0)) - (!min[2] ? 0.0 : (i==k ? 1.0 : 0.0)))/((min[2] || max[2] ? 1.0 : 2.0)*DX[2]););

						if (min[0])
						{
							if (m_bc_lo[0][k] == BC::Displacement)
							{rho = 0.0; aa = 1.0;}
							else if (m_bc_lo[0][k] == BC::Traction) 
								for (int l=0; l<AMREX_SPACEDIM; l++)
								{
									rho -= C(i,0,k,l,m,amrlev,mglev,mfi) * OffDiag_gradu_k(l);
									aa  -= C(i,0,k,l,m,amrlev,mglev,mfi) * Diag_gradu_k(l);
								}
							else Util::Abort("Invalid BC");
						}
						if (max[0])
						{
							if (m_bc_hi[0][k] == BC::Displacement)
							{rho = 0.0; aa = 1.0;}
							else if (m_bc_hi[0][k] == BC::Traction) 
								for (int l=0; l<AMREX_SPACEDIM; l++)
								{
									rho += C(i,0,k,l,m,amrlev,mglev,mfi) * OffDiag_gradu_k(l);
									aa  += C(i,0,k,l,m,amrlev,mglev,mfi) * Diag_gradu_k(l);
								}
							else Util::Abort("Invalid BC");
						}

						if (min[1])
						{
							if (m_bc_lo[1][k] == BC::Displacement)
							{rho = 0.0; aa = 1.0;}
							else if (m_bc_lo[1][k] == BC::Traction) 
								for (int l=0; l<AMREX_SPACEDIM; l++)
								{
									rho -= C(i,1,k,l,m,amrlev,mglev,mfi) * OffDiag_gradu_k(l);
									aa  -= C(i,1,k,l,m,amrlev,mglev,mfi) * Diag_gradu_k(l);
								}
							else Util::Abort("Invalid BC");
						}
						if (max[1])
						{
							if (m_bc_hi[1][k] == BC::Displacement)
							{rho = 0.0; aa = 1.0;}
							else if (m_bc_hi[1][k] == BC::Traction) 
								for (int l=0; l<AMREX_SPACEDIM; l++)
								{
									rho += C(i,1,k,l,m,amrlev,mglev,mfi) * OffDiag_gradu_k(l);
									aa  += C(i,1,k,l,m,amrlev,mglev,mfi) * Diag_gradu_k(l);
								}
							else Util::Abort("Invalid BC");
						}

						if (min[2])
						{
							if (m_bc_lo[2][k] == BC::Displacement)
							{rho = 0.0; aa = 1.0;}
							else if (m_bc_lo[2][k] == BC::Traction) 
								for (int l=0; l<AMREX_SPACEDIM; l++)
								{
									rho -= C(i,2,k,l,m,amrlev,mglev,mfi) * OffDiag_gradu_k(l);
									aa  -= C(i,2,k,l,m,amrlev,mglev,mfi) * Diag_gradu_k(l);
								}
							else Util::Abort("Invalid BC");
						}
						if (max[2])
						{
							if (m_bc_hi[2][k] == BC::Displacement)
							{rho = 0.0; aa = 1.0;}
							else if (m_bc_lo[2][k] == BC::Traction) 
								for (int l=0; l<AMREX_SPACEDIM; l++)
								{
									rho += C(i,2,k,l,m,amrlev,mglev,mfi) * OffDiag_gradu_k(l);
									aa  += C(i,2,k,l,m,amrlev,mglev,mfi) * Diag_gradu_k(l);
								}
							else Util::Abort("Invalid BC");
						}
						if (min[0] || max[0] || min[1] || max[1] || min[2] || max[2]) continue;



					
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
		}*/
}

template<class T>
void
Elastic<T>::normalize (int amrlev, int mglev, MultiFab& mf) const
{
	// See FApply for documentation. 

	amrex::Box domain(m_geom[amrlev][mglev].Domain());
	const Real* DX = m_geom[amrlev][mglev].CellSize();
	
	static std::array<amrex::IntVect,AMREX_SPACEDIM> dx = {AMREX_D_DECL(amrex::IntVect(AMREX_D_DECL(1,0,0)),
									    amrex::IntVect(AMREX_D_DECL(0,1,0)),
									    amrex::IntVect(AMREX_D_DECL(0,0,1)))};
	for (MFIter mfi(mf, true); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();
		amrex::BaseFab<T> &C = (*(model[amrlev][mglev]))[mfi];
		amrex::FArrayBox       &mffab    = mf[mfi];

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

			Set::Matrix gradu; // gradu(i,j) = u_{i,j)
			std::array<Set::Matrix,AMREX_SPACEDIM> gradgradu; // gradgradu[k](l,j) = u_{k,lj}

			for (int i = 0; i < AMREX_SPACEDIM; i++)
			{
				Set::Scalar aa = 0.0;
				for (int k = 0; k < AMREX_SPACEDIM; k++)
				{
					AMREX_D_TERM(gradu(k,0) = ((!xmax ? 0.0 : (i==k ? 1.0 : 0.0)) - (!xmin ? 0.0 : (i==k ? 1.0 : 0.0)))/((xmin || xmax ? 1.0 : 2.0)*DX[0]);,
						     gradu(k,1) = ((!ymax ? 0.0 : (i==k ? 1.0 : 0.0)) - (!ymin ? 0.0 : (i==k ? 1.0 : 0.0)))/((ymin || ymax ? 1.0 : 2.0)*DX[1]);,
						     gradu(k,2) = ((!zmax ? 0.0 : (i==k ? 1.0 : 0.0)) - (!zmin ? 0.0 : (i==k ? 1.0 : 0.0)))/((zmin || zmax ? 1.0 : 2.0)*DX[2]););
			
					AMREX_D_TERM(gradgradu[k](0,0) = (i==k ? -2.0 : 0.0)/DX[0]/DX[0];
						     ,// 2D
						     gradgradu[k](0,1) = 0.0;
						     gradgradu[k](1,0) = 0.0;
						     gradgradu[k](1,1) = (i==k ? -2.0 : 0.0)/DX[1]/DX[1];
						     ,// 3D
						     gradgradu[k](0,2) = 0.0;
						     gradgradu[k](1,2) = 0.0;
						     gradgradu[k](2,0) = 0.0;
						     gradgradu[k](2,1) = 0.0;
						     gradgradu[k](2,2) = (i==k ? -2.0 : 0.0))/DX[2]/DX[2];
				}

				Set::Matrix eps = 0.5*(gradu + gradu.transpose());
				Set::Matrix sig = C(m)(eps);

				if (xmax || xmin || ymax || ymin || zmax || zmin) 
				{
					for (int k = 0; k < AMREX_SPACEDIM; k++) // iterate over DIMENSIONS
					{
						for (int j = 0; j < AMREX_SPACEDIM; j++) // iterate over FACES
						{
							if (m[j] == domain.loVect()[j])
								if (m_bc_lo[j][k] == BC::Displacement)
									aa += 1.0;
								else if (m_bc_lo[j][k] == BC::Traction) 
									aa -= sig(k,j);
								else Util::Abort("Invalid BC");
							if (m[j] == domain.hiVect()[j] + 1)
								if (m_bc_hi[j][k] == BC::Displacement)
									aa += 1.0;
								else if (m_bc_hi[j][k] == BC::Traction) 
									aa += sig(k,j);
								else Util::Abort("Invalid BC");
						}
					}
					if (fabs(aa) < 1E-10) Util::Abort("Singular boundary operator: diagonal = 0");
				}
				else
				{
					Set::Vector f =
						C(m)(gradgradu) + 
						AMREX_D_TERM(((C(m+dx[0]) - C(m-dx[0]))/2.0/DX[0])(eps).col(0),
							     + ((C(m+dx[1]) - C(m-dx[1]))/2.0/DX[1])(eps).col(1),
							     + ((C(m+dx[2]) - C(m-dx[2]))/2.0/DX[2])(eps).col(2));
					aa += f(i);
				}

				if (fabs(aa) < 1E-10) Util::Abort("Singular operator: diagonal = 0");

				mffab(m,i) /= aa;
			}
		}
	}
}


template<class T>
void
Elastic<T>::FFlux (int /*amrlev*/, const MFIter& /*mfi*/,
		const std::array<FArrayBox*,AMREX_SPACEDIM>& sigmafab,
		const FArrayBox& /*ufab*/, const int /*face_only*/) const
{
	amrex::BaseFab<amrex::Real> AMREX_D_DECL( &fxfab = *sigmafab[0],
	 					  &fyfab = *sigmafab[1],
	 					  &fzfab = *sigmafab[2] ) ;
	AMREX_D_TERM(fxfab.setVal(0.0);,
	 	     fyfab.setVal(0.0);,
	 	     fzfab.setVal(0.0););

}


template<class T>
void
Elastic<T>::Stress (FArrayBox& sigmafab,
		 const FArrayBox& ufab,
		 int amrlev, const MFIter& mfi,
		 bool voigt) const
{
	Util::Abort("Stress() is not yet implemented - do not use");
	/*
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
	*/
}

template<class T>
void
Elastic<T>::Energy (FArrayBox& energyfab,
		 const FArrayBox& ufab,
		 int amrlev, const MFIter& mfi) const
{
	Util::Abort("Energy() is not yet implemented - do not use");
	/*
	const amrex::Real* DX = m_geom[amrlev][0].CellSize();

	static std::array<amrex::IntVect,AMREX_SPACEDIM> dx = {AMREX_D_DECL(amrex::IntVect(AMREX_D_DECL(1,0,0)),
									    amrex::IntVect(AMREX_D_DECL(0,1,0)),
									    amrex::IntVect(AMREX_D_DECL(0,0,1)))};

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
	*/
}



template class Elastic<Model::Solid::Elastic::Isotropic::Isotropic>;


}
}
