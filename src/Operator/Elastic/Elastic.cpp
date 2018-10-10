#include "Model/Solid/Elastic/Isotropic/Isotropic.H"
#include "Model/Solid/Elastic/Cubic/Cubic.H"
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
	
	// DEBUG: Check to see if Fapply is getting passed bad values.
	for (MFIter mfi(f, true); mfi.isValid(); ++mfi)
	{
		const amrex::FArrayBox &ufab    = u[mfi];
		if(ufab.contains_inf()) Util::Abort(INFO, "Inf in ufab [before update]");
		if(ufab.contains_nan()) Util::Abort(INFO, "Nan in ufab [before update]");
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
			bool    AMREX_D_DECL(xmin = (m1 == domain.loVect()[0]),
					     ymin = (m2 == domain.loVect()[1]),
					     zmin = (m3 == domain.loVect()[2])),
				AMREX_D_DECL(xmax = (m1 == domain.hiVect()[0] + 1),
					     ymax = (m2 == domain.hiVect()[1] + 1),
					     zmax = (m3 == domain.hiVect()[2] + 1));

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
					{
						if (m_bc_lo[j][i] == BC::Displacement)
							ffab(m,i) = ufab(m,i);
						else if (m_bc_lo[j][i] == BC::Traction) 
							ffab(m,i) = -sig(i,j);
						else Util::Abort(INFO, "Invalid BC");
					}
					if (m[j] == domain.hiVect()[j] + 1)
					{
						if (m_bc_hi[j][i] == BC::Displacement)
							ffab(m,i) = ufab(m,i);
						else if (m_bc_hi[j][i] == BC::Traction) 
							ffab(m,i) = +sig(i,j);
						else Util::Abort(INFO, "Invalid BC");
					}
				}
			}
			if (AMREX_D_TERM(xmax || xmin,
					 || ymax || ymin,
					 || zmax || zmin)) continue;
			
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
	amrex::Box domain(m_geom[amrlev][mglev].Domain());
	const Real* DX = m_geom[amrlev][mglev].CellSize();
	
	for (int redblack = 0; redblack < 2; redblack++)
	{
		for (MFIter mfi(rhs, true); mfi.isValid(); ++mfi)
		{
			const Box           &bx     = mfi.tilebox();
			amrex::BaseFab<T>   &C      = (*(model[amrlev][mglev]))[mfi];
			amrex::FArrayBox    &ufab   = u[mfi];
			const amrex::FArrayBox    &rhsfab   = rhs[mfi];

			AMREX_D_TERM(for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++),
				     for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++),
				     for (int m3 = bx.loVect()[2]; m3<=bx.hiVect()[2]; m3++))
			{
				if ((AMREX_D_TERM(m1, + m2, + m3))%2 == redblack) continue;

				amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));
				bool    AMREX_D_DECL(xmin = (m1 == domain.loVect()[0]),
						     ymin = (m2 == domain.loVect()[1]),
						     zmin = (m3 == domain.loVect()[2])),
					AMREX_D_DECL(xmax = (m1 == domain.hiVect()[0] + 1),
						     ymax = (m2 == domain.hiVect()[1] + 1),
						     zmax = (m3 == domain.hiVect()[2] + 1));

				Set::Matrix graduD, graduU; // gradu(i,j) = u_{i,j)
				std::array<Set::Matrix,AMREX_SPACEDIM> gradgraduD,gradgraduU; // gradgradu[k](l,j) = u_{k,lj}

				for (int i = 0; i < AMREX_SPACEDIM; i++)
				{
					Set::Scalar aa = 0.0, rho = 0.0;
					for (int k = 0; k < AMREX_SPACEDIM; k++)
					{
						AMREX_D_TERM(graduD(k,0) = ((!xmax ? 0.0 : (i==k ? 1.0 : 0.0)) - (!xmin ? 0.0 : (i==k ? 1.0 : 0.0)))/((xmin || xmax ? 1.0 : 2.0)*DX[0]);,
							     graduD(k,1) = ((!ymax ? 0.0 : (i==k ? 1.0 : 0.0)) - (!ymin ? 0.0 : (i==k ? 1.0 : 0.0)))/((ymin || ymax ? 1.0 : 2.0)*DX[1]);,
							     graduD(k,2) = ((!zmax ? 0.0 : (i==k ? 1.0 : 0.0)) - (!zmin ? 0.0 : (i==k ? 1.0 : 0.0)))/((zmin || zmax ? 1.0 : 2.0)*DX[2]););

						AMREX_D_TERM(graduU(i,0) = ((!xmax ? ufab(m+dx[0],i) : (i==k ? 0.0 : ufab(m,i))) - (!xmin ? ufab(m-dx[0],i) : (i==k ? 0.0 : ufab(m,i))))/((xmin || xmax ? 1.0 : 2.0)*DX[0]);,
							     graduU(i,1) = ((!ymax ? ufab(m+dx[1],i) : (i==k ? 0.0 : ufab(m,i))) - (!ymin ? ufab(m-dx[1],i) : (i==k ? 0.0 : ufab(m,i))))/((ymin || ymax ? 1.0 : 2.0)*DX[1]);,
							     graduU(i,2) = ((!zmax ? ufab(m+dx[2],i) : (i==k ? 0.0 : ufab(m,i))) - (!zmin ? ufab(m-dx[2],i) : (i==k ? 0.0 : ufab(m,i))))/((zmin || zmax ? 1.0 : 2.0)*DX[2]););

			
			

						AMREX_D_TERM(gradgraduD[k](0,0) = (i==k ? -2.0 : 0.0)/DX[0]/DX[0];
							     ,// 2D
							     gradgraduD[k](0,1) = 0.0;
							     gradgraduD[k](1,0) = 0.0;
							     gradgraduD[k](1,1) = (i==k ? -2.0 : 0.0)/DX[1]/DX[1];
							     ,// 3D
							     gradgraduD[k](0,2) = 0.0;
							     gradgraduD[k](1,2) = 0.0;
							     gradgraduD[k](2,0) = 0.0;
							     gradgraduD[k](2,1) = 0.0;
							     gradgraduD[k](2,2) = (i==k ? -2.0 : 0.0)/DX[2]/DX[2]);

						AMREX_D_TERM(gradgraduU[i](0,0) = (ufab(m+dx[0],i) - (i==k ? 0.0 : 2.0*ufab(m,i)) + ufab(m-dx[0],i))/DX[0]/DX[0];
							     ,// 2D
							     gradgraduU[i](0,1) = (ufab(m+dx[0]+dx[1],i) + ufab(m-dx[0]-dx[1],i) - ufab(m+dx[0]-dx[1],i) - ufab(m-dx[0]+dx[1],i))/(2.0*DX[0])/(2.0*DX[1]);
							     gradgraduU[i](1,0) = gradgraduU[i](0,1);
							     gradgraduU[i](1,1) = (ufab(m+dx[1],i) - (i==k ? 0.0 : 2.0*ufab(m,i)) + ufab(m-dx[1],i))/DX[1]/DX[1];
							     ,// 3D
							     gradgraduU[i](0,2) = (ufab(m+dx[0]+dx[2],i) + ufab(m-dx[0]-dx[2],i) - ufab(m+dx[0]-dx[2],i) - ufab(m-dx[0]+dx[2],i))/(2.0*DX[0])/(2.0*DX[2]);
							     gradgraduU[i](1,2) = (ufab(m+dx[1]+dx[2],i) + ufab(m-dx[1]-dx[2],i) - ufab(m+dx[1]-dx[2],i) - ufab(m-dx[1]+dx[2],i))/(2.0*DX[1])/(2.0*DX[2]);
							     gradgraduU[i](2,0) = gradgraduU[i](0,2);
							     gradgraduU[i](2,1) = gradgraduU[i](1,2);
							     gradgraduU[i](2,2) = (ufab(m+dx[2],i) - (i==k ? 0.0 : 2.0*ufab(m,i)) + ufab(m-dx[2],i))/DX[2]/DX[2];);

					}

					Set::Matrix epsD = 0.5*(graduD + graduD.transpose());
					Set::Matrix epsU = 0.5*(graduU + graduU.transpose());

					Set::Matrix sigD = C(m)(epsD);
					Set::Matrix sigU = C(m)(epsU);

					if (AMREX_D_TERM(xmax || xmin,
							 || ymax || ymin,
							 || zmax || zmin))
					{
						for (int k = 0; k < AMREX_SPACEDIM; k++) // iterate over DIMENSIONS
						{
							for (int j = 0; j < AMREX_SPACEDIM; j++) // iterate over FACES
							{
								if (m[j] == domain.loVect()[j])
								{
									if (m_bc_lo[j][k] == BC::Displacement)
									{ aa += 1.0; rho += 0.0;}
									else if (m_bc_lo[j][k] == BC::Traction) 
									{ aa -= sigD(k,j); rho -= sigU(k,j); }
									else Util::Abort(INFO, "Invalid BC");
								}
								if (m[j] == domain.hiVect()[j] + 1)
								{
									if (m_bc_hi[j][k] == BC::Displacement)
									{ aa += 1.0; rho += 0.0;}
									else if (m_bc_hi[j][k] == BC::Traction) 
									{ aa += sigD(k,j); rho += sigU(k,j); }
									else Util::Abort(INFO, "Invalid BC");
								}
							}
						}
						if (fabs(aa) < 1E-10) Util::Abort(INFO, "Singular boundary operator in Fsmooth: diagonal = 0");
					}
					else
					{
						Set::Vector fD =
							C(m)(gradgraduD) + 
							AMREX_D_TERM(((C(m+dx[0]) - C(m-dx[0]))/2.0/DX[0])(epsD).col(0),
								     + ((C(m+dx[1]) - C(m-dx[1]))/2.0/DX[1])(epsD).col(1),
								     + ((C(m+dx[2]) - C(m-dx[2]))/2.0/DX[2])(epsD).col(2));
						Set::Vector fU =
							C(m)(gradgraduU) + 
							AMREX_D_TERM(((C(m+dx[0]) - C(m-dx[0]))/2.0/DX[0])(epsU).col(0),
								     + ((C(m+dx[1]) - C(m-dx[1]))/2.0/DX[1])(epsU).col(1),
								     + ((C(m+dx[2]) - C(m-dx[2]))/2.0/DX[2])(epsU).col(2));
						aa += fD(i);
						rho += fU(i);
					}

					if (fabs(aa) < 1E-10) Util::Abort(INFO, "Singular operator in Fsmooth: diagonal = 0");

					ufab(m,i) = (rhsfab(m,i) - rho) / aa;
				}
			}
		}
	}
}

template<class T>
void
Elastic<T>::normalize (int amrlev, int mglev, MultiFab& mf) const
{
	// See FApply for documentation. 

	amrex::Box domain(m_geom[amrlev][mglev].Domain());
	const Real* DX = m_geom[amrlev][mglev].CellSize();
	
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
			bool    AMREX_D_DECL(xmin = (m1 == domain.loVect()[0]),
					     ymin = (m2 == domain.loVect()[1]),
					     zmin = (m3 == domain.loVect()[2])),
				AMREX_D_DECL(xmax = (m1 == domain.hiVect()[0] + 1),
					     ymax = (m2 == domain.hiVect()[1] + 1),
					     zmax = (m3 == domain.hiVect()[2] + 1));

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
						     gradgradu[k](2,2) = (i==k ? -2.0 : 0.0)/DX[2]/DX[2]);
				}

				Set::Matrix eps = 0.5*(gradu + gradu.transpose());
				Set::Matrix sig = C(m)(eps);

				if (AMREX_D_TERM(xmax || xmin,
						 || ymax || ymin,
						 || zmax || zmin)) 
				{
					for (int k = 0; k < AMREX_SPACEDIM; k++) // iterate over DIMENSIONS
					{
						for (int j = 0; j < AMREX_SPACEDIM; j++) // iterate over FACES
						{
							if (m[j] == domain.loVect()[j])
							{
								if (m_bc_lo[j][k] == BC::Displacement)
									aa += 1.0;
								else if (m_bc_lo[j][k] == BC::Traction) 
									aa -= sig(k,j);
								else Util::Abort(INFO, "Invalid BC");
							}
							if (m[j] == domain.hiVect()[j] + 1)
							{
								if (m_bc_hi[j][k] == BC::Displacement)
									aa += 1.0;
								else if (m_bc_hi[j][k] == BC::Traction) 
									aa += sig(k,j);
								else Util::Abort(INFO, "Invalid BC");
							}
						}
					}
					if (fabs(aa) < 1E-10) Util::Abort(INFO, "Singular boundary operator in normalize: diagonal = ", aa);
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

				//Util::Message(INFO, "eps = \n", eps);

				if (fabs(aa) < 1E-10) Util::Abort(INFO, "Singular operator in normalize: diagonal = ", aa);
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
Elastic<T>::Stress (int amrlev,
		    amrex::MultiFab& sigma,
		    const amrex::MultiFab& u,
		    bool voigt) const
{
	amrex::Box domain(m_geom[amrlev][0].Domain());
	if (voigt)
		AMREX_ASSERT(sigma.nComp() == (AMREX_SPACEDIM*(AMREX_SPACEDIM-1)/2));
	else
		AMREX_ASSERT(sigma.nComp() == AMREX_SPACEDIM*AMREX_SPACEDIM);

	const amrex::Real* DX = m_geom[amrlev][0].CellSize();

	for (MFIter mfi(u, true); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();
		amrex::BaseFab<T> &C          = (*(model[amrlev][0]))[mfi];
		amrex::FArrayBox  &sigmafab   = sigma[mfi];
		const amrex::FArrayBox  &ufab = u[mfi];

		AMREX_D_TERM(for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++),
			     for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++),
			     for (int m3 = bx.loVect()[2]; m3<=bx.hiVect()[2]; m3++))
		{
			amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));
			bool    AMREX_D_DECL(xmin = (m1 == domain.loVect()[0]),
					     ymin = (m2 == domain.loVect()[1]),
					     zmin = (m3 == domain.loVect()[2])),
				AMREX_D_DECL(xmax = (m1 == domain.hiVect()[0] + 1),
					     ymax = (m2 == domain.hiVect()[1] + 1),
					     zmax = (m3 == domain.hiVect()[2] + 1));


			Set::Matrix gradu;

			AMREX_D_TERM(gradu(0,0) = ((!xmax ? ufab(m+dx[0],0) : ufab(m,0)) - (!xmin ? ufab(m-dx[0],0) : ufab(m,0)))/((xmin || xmax ? 1.0 : 2.0)*DX[0]);
				     ,
				     gradu(0,1) = ((!ymax ? ufab(m+dx[1],0) : ufab(m,0)) - (!ymin ? ufab(m-dx[1],0) : ufab(m,0)))/((ymin || ymax ? 1.0 : 2.0)*DX[1]);
				     gradu(1,0) = ((!xmax ? ufab(m+dx[0],1) : ufab(m,1)) - (!xmin ? ufab(m-dx[0],1) : ufab(m,1)))/((xmin || xmax ? 1.0 : 2.0)*DX[0]);
				     gradu(1,1) = ((!ymax ? ufab(m+dx[1],1) : ufab(m,1)) - (!ymin ? ufab(m-dx[1],1) : ufab(m,1)))/((ymin || ymax ? 1.0 : 2.0)*DX[1]);
				     ,
				     gradu(0,2) = ((!zmax ? ufab(m+dx[2],0) : ufab(m,0)) - (!zmin ? ufab(m-dx[2],0) : ufab(m,0)))/((zmin || zmax ? 1.0 : 2.0)*DX[2]);
				     gradu(1,2) = ((!zmax ? ufab(m+dx[2],1) : ufab(m,1)) - (!zmin ? ufab(m-dx[2],1) : ufab(m,1)))/((zmin || zmax ? 1.0 : 2.0)*DX[2]);
				     gradu(2,0) = ((!xmax ? ufab(m+dx[0],2) : ufab(m,2)) - (!xmin ? ufab(m-dx[0],2) : ufab(m,2)))/((xmin || xmax ? 1.0 : 2.0)*DX[0]);
				     gradu(2,1) = ((!ymax ? ufab(m+dx[1],2) : ufab(m,2)) - (!ymin ? ufab(m-dx[1],2) : ufab(m,2)))/((ymin || ymax ? 1.0 : 2.0)*DX[1]);
				     gradu(2,2) = ((!zmax ? ufab(m+dx[2],2) : ufab(m,2)) - (!zmin ? ufab(m-dx[2],2) : ufab(m,2)))/((zmin || zmax ? 1.0 : 2.0)*DX[2]););
			
			Set::Matrix eps = 0.5 * (gradu + gradu.transpose());

			Set::Matrix sig = C(m)(eps);

			if (voigt)
			{
#if AMREX_SPACEDIM == 2
				sigmafab(m,0) = sig(0,0); sigmafab(m,1) = sig(1,1); sigmafab(m,2) = sig(0,1); 
#elif AMREX_SPACEDIM == 3
				sigmafab(m,0) = sig(0,0); sigmafab(m,1) = sig(1,1); sigmafab(m,2) = sig(2,2); 
				sigmafab(m,3) = sig(1,2); sigmafab(m,4) = sig(2,0); sigmafab(m,5) = sig(0,1); 
#endif
			}
			else
			{
#if   AMREX_SPACEDIM == 2
				sigmafab(m,0) = sig(0,0); sigmafab(m,1) = sig(0,1); 
				sigmafab(m,2) = sig(1,0); sigmafab(m,3) = sig(1,1); 
#elif AMREX_SPACEDIM == 3
				sigmafab(m,0) = sig(0,0); sigmafab(m,1) = sig(0,1); sigmafab(m,2) = sig(0,2); 
				sigmafab(m,3) = sig(1,0); sigmafab(m,4) = sig(1,1); sigmafab(m,5) = sig(1,2); 
				sigmafab(m,6) = sig(2,0); sigmafab(m,7) = sig(2,1); sigmafab(m,8) = sig(2,2); 
#endif
			}
			}
	}
}


template<class T>
void
Elastic<T>::Energy (int amrlev,
		    amrex::MultiFab& energy,
		    const amrex::MultiFab& u) const
{
	amrex::Box domain(m_geom[amrlev][0].Domain());
	AMREX_ASSERT(energy.nComp() == 1);
	AMREX_ASSERT(u.nComp() == AMREX_SPACEDIM);

	const amrex::Real* DX = m_geom[amrlev][0].CellSize();

	for (MFIter mfi(u, true); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();
		amrex::BaseFab<T> &C          = (*(model[amrlev][0]))[mfi];
		amrex::FArrayBox  &energyfab  = energy[mfi];
		const amrex::FArrayBox  &ufab = u[mfi];

		AMREX_D_TERM(for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++),
			     for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++),
			     for (int m3 = bx.loVect()[2]; m3<=bx.hiVect()[2]; m3++))
		{
			amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));
			bool    AMREX_D_DECL(xmin = (m1 == domain.loVect()[0]),
					     ymin = (m2 == domain.loVect()[1]),
					     zmin = (m3 == domain.loVect()[2])),
				AMREX_D_DECL(xmax = (m1 == domain.hiVect()[0] + 1),
					     ymax = (m2 == domain.hiVect()[1] + 1),
					     zmax = (m3 == domain.hiVect()[2] + 1));


			Set::Matrix gradu;

			AMREX_D_TERM(gradu(0,0) = ((!xmax ? ufab(m+dx[0],0) : ufab(m,0)) - (!xmin ? ufab(m-dx[0],0) : ufab(m,0)))/((xmin || xmax ? 1.0 : 2.0)*DX[0]);
				     ,
				     gradu(0,1) = ((!ymax ? ufab(m+dx[1],0) : ufab(m,0)) - (!ymin ? ufab(m-dx[1],0) : ufab(m,0)))/((ymin || ymax ? 1.0 : 2.0)*DX[1]);
				     gradu(1,0) = ((!xmax ? ufab(m+dx[0],1) : ufab(m,1)) - (!xmin ? ufab(m-dx[0],1) : ufab(m,1)))/((xmin || xmax ? 1.0 : 2.0)*DX[0]);
				     gradu(1,1) = ((!ymax ? ufab(m+dx[1],1) : ufab(m,1)) - (!ymin ? ufab(m-dx[1],1) : ufab(m,1)))/((ymin || ymax ? 1.0 : 2.0)*DX[1]);
				     ,
				     gradu(0,2) = ((!zmax ? ufab(m+dx[2],0) : ufab(m,0)) - (!zmin ? ufab(m-dx[2],0) : ufab(m,0)))/((zmin || zmax ? 1.0 : 2.0)*DX[2]);
				     gradu(1,2) = ((!zmax ? ufab(m+dx[2],1) : ufab(m,1)) - (!zmin ? ufab(m-dx[2],1) : ufab(m,1)))/((zmin || zmax ? 1.0 : 2.0)*DX[2]);
				     gradu(2,0) = ((!xmax ? ufab(m+dx[0],2) : ufab(m,2)) - (!xmin ? ufab(m-dx[0],2) : ufab(m,2)))/((xmin || xmax ? 1.0 : 2.0)*DX[0]);
				     gradu(2,1) = ((!ymax ? ufab(m+dx[1],2) : ufab(m,2)) - (!ymin ? ufab(m-dx[1],2) : ufab(m,2)))/((ymin || ymax ? 1.0 : 2.0)*DX[1]);
				     gradu(2,2) = ((!zmax ? ufab(m+dx[2],2) : ufab(m,2)) - (!zmin ? ufab(m-dx[2],2) : ufab(m,2)))/((zmin || zmax ? 1.0 : 2.0)*DX[2]););
			
			Set::Matrix eps = 0.5 * (gradu + gradu.transpose());
			Set::Matrix sig = C(m)(eps);

			energyfab(m) = (eps.transpose() * sig).trace();

		}
	}
}


template<class T>
void
Elastic<T>::reflux (int crse_amrlev,
		    MultiFab& res, const MultiFab& crse_sol, const MultiFab& crse_rhs,
		    MultiFab& fine_res, MultiFab& fine_sol, const MultiFab& fine_rhs) const
{
	Util::Message(INFO, "reflux implementation in progress");


	int ncomp = AMREX_SPACEDIM;

	const Real* DX = m_geom[crse_amrlev+1][0].CellSize();

	const Geometry& cgeom = m_geom[crse_amrlev  ][0];
 	const Geometry& fgeom = m_geom[crse_amrlev+1][0];
 	const Box& c_cc_domain = cgeom.Domain();
 	const Box& c_nd_domain = amrex::surroundingNodes(c_cc_domain);

	Util::Message(INFO);
 	const BoxArray& fba = fine_sol.boxArray();
 	const DistributionMapping& fdm = fine_sol.DistributionMap();

 	const iMultiFab& fdmsk = *m_dirichlet_mask[crse_amrlev+1][0];

 	MultiFab fine_res_for_coarse(amrex::coarsen(fba, 2), fdm, ncomp, 0);

	Util::Message(INFO);
 	applyBC(crse_amrlev+1, 0, fine_res, BCMode::Inhomogeneous, StateMode::Solution);

	Util::Message(INFO);

	Util::Message(INFO);

	if (AMREX_SPACEDIM != 2) Util::Abort("2D implementation only");

	for (MFIter mfi(fine_res_for_coarse, true); mfi.isValid(); ++mfi)
	{
		Util::Message(INFO);
		const Box& bx = mfi.tilebox();
		const FArrayBox &fine = fine_res[mfi];
		FArrayBox &crse = fine_res_for_coarse[mfi];

		crse.setVal(0.0);
		
		// amrex_mlndlap_restriction
		for (int n = 0 ; n < ncomp; n++)
		{
			AMREX_D_TERM(for (int m1 = bx.loVect()[0] +1; m1<=bx.hiVect()[0] -1; m1++), // Do not include
				     for (int m2 = bx.loVect()[1] +1; m2<=bx.hiVect()[1] -1; m2++), // boundaries, interior
				     for (int m3 = bx.loVect()[2] +1; m3<=bx.hiVect()[2] -1; m3++)) // only.
			{
				amrex::IntVect m_crse(AMREX_D_DECL(m1,m2,m3));
				amrex::IntVect m_fine(AMREX_D_DECL(m1*2,m2*2,m3*2));

				crse(m_crse,n) =
					(fine(m_fine-dx[0]-dx[1],n) + 2.0*fine(m_fine-dx[1],n) + fine(m_fine+dx[0]-dx[1],n)
					 + 2.0*fine(m_fine-dx[0],n) + 4.0*fine(m_fine,n)  + 2.0*fine(m_fine+dx[0],n) 
					 + fine(m_fine-dx[0]+dx[1],n) + 2.0*fine(m_fine+dx[1],n) + fine(m_fine+dx[0]+dx[1],n))/16.0;

			}
		}
	}

	
	res.ParallelCopy(fine_res_for_coarse, cgeom.periodicity());
	MultiFab fine_contrib(amrex::coarsen(fba, 2), fdm, ncomp, 0);
	fine_contrib.setVal(0.0);

	FArrayBox Axfab;
        for (MFIter mfi(fine_contrib, MFItInfo().EnableTiling().SetDynamic(true));
             mfi.isValid(); ++mfi)
	{
		const Box& cbx = mfi.tilebox(); // clo, chi
		const Box& cvbx = mfi.validbox(); // cglo, cghi

		const Box& fbx = amrex::refine(cbx,2); // lo, hi
		const Box& fvbx = amrex::refine(cvbx,2); // glo, ghi


		Util::Message(INFO, "amrlev  = ", crse_amrlev);
		Util::Message(INFO, "cbx  = ", cbx);
		Util::Message(INFO, "cvbx = ", cvbx);
		Util::Message(INFO, "fbx  = ", fbx);
		Util::Message(INFO, "fvbx = ", fvbx);


		const Box& bx_Ax = amrex::grow(fbx,1);
		const Box& b2 = bx_Ax & amrex::grow(fvbx,-1);

		Axfab.resize(bx_Ax,ncomp);
		Axfab.setVal(0.0);
		Axfab.copy(fine_rhs[mfi], b2, 0, b2, 0, 1);
		Axfab.minus(fine_res[mfi], b2, 0, 0, 1);

		amrex::IntVect gtlo(std::min(fbx.loVect()[0]-1, fvbx.loVect()[0]),
				    std::min(fbx.loVect()[1]-1, fvbx.loVect()[1]));

		amrex::IntVect gthi(std::max(fbx.loVect()[0]+1, fvbx.loVect()[0]),
				    std::min(fbx.loVect()[1]+1, fvbx.loVect()[1]));

		FArrayBox &ufab = fine_sol[mfi];
		FArrayBox &ffab = fine_contrib[mfi];
		amrex::BaseFab<T> &C = (*(model[crse_amrlev+1][0]))[mfi];

		for (int m1 = gtlo[0]; m1<=gthi[0]; m1++)
		for (int m2 = gtlo[1]; m2<=gthi[1]; m2++)
		{
			amrex::IntVect m(m1, m2);

			bool    xmin = (m1 == gtlo[0]),
				ymin = (m2 == gtlo[1]),
				xmax = (m1 == gthi[0]),
				ymax = (m2 == gthi[1]);

			// Consider boundary terms only
			if (!(xmin || ymin || xmax || ymin)) continue; 

			
			Set::Matrix gradu; // gradu(i,j) = u_{i,j)

			// The gradient of the displacement gradient tensor
			std::array<Set::Matrix,AMREX_SPACEDIM> gradgradu; // gradgradu[k](l,j) = u_{k,lj}


			// Fill gradu and gradgradu
			for (int i = 0; i < AMREX_SPACEDIM; i++)
			{
				AMREX_D_TERM(gradu(i,0) = ((!xmax ? ufab(m+dx[0],i) : ufab(m,i)) - (!xmin ? ufab(m-dx[0],i) : ufab(m,i)))/((xmin || xmax ? 1.0 : 2.0)*DX[0]);,
					     gradu(i,1) = ((!ymax ? ufab(m+dx[1],i) : ufab(m,i)) - (!ymin ? ufab(m-dx[1],i) : ufab(m,i)))/((ymin || ymax ? 1.0 : 2.0)*DX[1]);,
					     gradu(i,2) = ((!zmax ? ufab(m+dx[2],i) : ufab(m,i)) - (!zmin ? ufab(m-dx[2],i) : ufab(m,i)))/((zmin || zmax ? 1.0 : 2.0)*DX[2]););
			}

			Set::Matrix eps = 0.5*(gradu + gradu.transpose());
			// Stress tensor computed using the model fab
			Set::Matrix sig = C(m)(eps);

			for (int i = 0; i < AMREX_SPACEDIM; i++) // iterate over DIMENSIONS
				for (int j = 0; j < AMREX_SPACEDIM; j++) // iterate over FACES
					if (m[j] == gtlo[j])
						Axfab(m,i) = -sig(i,j);
					else if (m[j] == gthi[j])
						Axfab(m,i) = +sig(i,j);

			for (int i = 0; i < ncomp ; i++)
				for (int m1 = cbx.loVect()[0]; m1<=cbx.loVect()[0]+1; m1++)
					for (int m2 = cbx.loVect()[1]; m2<=cbx.loVect()[1]+1; m2++)
					{
						amrex::IntVect m(m1, m2);
						bool    xmin = (m1 == cbx.loVect()[0]),
							ymin = (m2 == cbx.loVect()[1]),
							xmax = (m1 == cbx.hiVect()[0]),
							ymax = (m2 == cbx.hiVect()[1]);
						if (!(xmin || ymin || xmax || ymax)) continue;
						ffab(m,i) = ffab(m,i) + (0.25)*Axfab(2*m,i)
							+ (0.5)*(Axfab(2*m + dx[0]) +
								 Axfab(2*m - dx[0]) +
								 Axfab(2*m + dx[1]) +
								 Axfab(2*m - dx[1]))
							+ (0.25)*(Axfab(2*m - dx[0] - dx[1]) +
								  Axfab(2*m - dx[0] + dx[1]) +
								  Axfab(2*m + dx[1] - dx[1]) +
								  Axfab(2*m + dx[1] + dx[1]));
				
					
					}
		}
	}
	MultiFab fine_contrib_on_crse(crse_sol.boxArray(), crse_sol.DistributionMap(), 1, 0);
	fine_contrib_on_crse.setVal(0.0);
	fine_contrib_on_crse.ParallelAdd(fine_contrib, cgeom.periodicity());
	     

	//Util::Abort("Exit normally");

// 	const iMultiFab& cdmsk = *m_dirichlet_mask[crse_amrlev][0];
// 	const auto& nd_mask     = m_nd_fine_mask[crse_amrlev];
// 	const auto& cc_mask     = m_cc_fine_mask[crse_amrlev];
// 	const auto& has_fine_bndry = m_has_fine_bndry[crse_amrlev];

// 	const auto& csigma = *m_sigma[crse_amrlev][0][0];

 	for (MFIter mfi(res, MFItInfo().EnableTiling().SetDynamic(true)); mfi.isValid(); ++mfi)
 	{
// 		if ((*has_fine_bndry)[mfi])
// 		{
// 			const Box& bx = mfi.tilebox();
// 			amrex_mlndlap_res_cf_contrib(BL_TO_FORTRAN_BOX(bx),
// 						     BL_TO_FORTRAN_ANYD(res[mfi]),
// 						     BL_TO_FORTRAN_ANYD(crse_sol[mfi]),
// 						     BL_TO_FORTRAN_ANYD(crse_rhs[mfi]),
// 						     BL_TO_FORTRAN_ANYD(csigma[mfi]),
// 						     BL_TO_FORTRAN_ANYD(cdmsk[mfi]),
// 						     BL_TO_FORTRAN_ANYD((*nd_mask)[mfi]),
// 						     BL_TO_FORTRAN_ANYD((*cc_mask)[mfi]),
// 						     BL_TO_FORTRAN_ANYD(fine_contrib_on_crse[mfi]),
// 						     cdxinv, BL_TO_FORTRAN_BOX(c_nd_domain),
// 						     m_lobc.data(), m_hibc.data());
// 		}
 	}

	Util::Message(INFO, "reflux implementation in progress");
}






template class Elastic<Model::Solid::Elastic::Isotropic::Isotropic>;
template class Elastic<Model::Solid::Elastic::Cubic::Cubic>;


}
}
