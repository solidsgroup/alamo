#include "Model/Solid/LinearElastic/Isotropic.H"
#include "Model/Solid/LinearElastic/Cubic.H"
#include "Model/Solid/LinearElastic/Laplacian.H"
#include "Model/Solid/LinearElastic/Degradable/Isotropic.H"
#include "Model/Solid/Viscoelastic/Isotropic.H"
#include "Elastic.H"


namespace Operator
{
template<class T>
Elastic<T>::Elastic (const Vector<Geometry>& a_geom,
		     const Vector<BoxArray>& a_grids,
		     const Vector<DistributionMapping>& a_dmap,
		     const LPInfo& a_info)
{
	BL_PROFILE("Operator::Elastic::Elastic()");

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
	BL_PROFILE("Operator::Elastic::define()");

	Operator::define(a_geom,a_grids,a_dmap,a_info,a_factory);

	int model_nghost = 2;

	model.resize(m_num_amr_levels);
	for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
	{
		model[amrlev].resize(m_num_mg_levels[amrlev]);
		for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
		{
			model[amrlev][mglev].reset(new MultiTab(amrex::convert(m_grids[amrlev][mglev],
									       amrex::IntVect::TheNodeVector()),
								m_dmap[amrlev][mglev], 1, model_nghost));
		}
	}
}

template <class T>
void
Elastic<T>::SetModel (int amrlev, const amrex::FabArray<amrex::BaseFab<T> >& a_model)
{
	BL_PROFILE("Operator::Elastic::SetModel()");

	int nghost = model[amrlev][0]->nGrow();
	for (MFIter mfi(a_model, false); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.validbox();
		amrex::BaseFab<T> &modelfab = (*(model[amrlev][0]))[mfi];
		const amrex::BaseFab<T> &a_modelfab = a_model[mfi];

		AMREX_D_TERM(for (int m1 = bx.loVect()[0]-nghost; m1<=bx.hiVect()[0]+nghost; m1++),
			     for (int m2 = bx.loVect()[1]-nghost; m2<=bx.hiVect()[1]+nghost; m2++),
			     for (int m3 = bx.loVect()[2]-nghost; m3<=bx.hiVect()[2]+nghost; m3++))
		{
			amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));
			modelfab(m) = a_modelfab(m);
		}
	}
}

template<class T>
void
Elastic<T>::Fapply (int amrlev, int mglev, MultiFab& a_f, const MultiFab& a_u) const
{
	BL_PROFILE("Operator::Elastic::Fapply()");

	amrex::Box domain(m_geom[amrlev][mglev].Domain());
	domain.convert(amrex::IntVect::TheNodeVector());

	const Real* DX = m_geom[amrlev][mglev].CellSize();

	for (MFIter mfi(a_f, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		Box bx = mfi.tilebox();
		bx.grow(1);        // Expand to cover first layer of ghost nodes
		bx = bx & domain;  // Take intersection of box and the problem domain
			
		amrex::Array4<T> const& C                 = (*(model[amrlev][mglev])).array(mfi);
		amrex::Array4<const amrex::Real> const& U = a_u.array(mfi);
		amrex::Array4<amrex::Real> const& F       = a_f.array(mfi);

		const Dim3 lo= amrex::lbound(domain), hi = amrex::ubound(domain);
			
		amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {
					
				Set::Vector f = Set::Vector::Zero();

				bool    AMREX_D_DECL(xmin = (i == lo.x), ymin = (j==lo.y), zmin = (k==lo.z)),
					AMREX_D_DECL(xmax = (i == hi.x), ymax = (j==hi.y), zmax = (k==hi.z));

				// The displacement gradient tensor
				Set::Matrix gradu; // gradu(i,j) = u_{i,j)

				// Fill gradu and gradgradu
				for (int p = 0; p < AMREX_SPACEDIM; p++)
				{
					AMREX_D_TERM(gradu(p,0) = ((!xmax ? U(i+1,j,k,p) : U(i,j,k,p)) - (!xmin ? U(i-1,j,k,p) : U(i,j,k,p)))/((xmin || xmax ? 1.0 : 2.0)*DX[0]);,
						     gradu(p,1) = ((!ymax ? U(i,j+1,k,p) : U(i,j,k,p)) - (!ymin ? U(i,j-1,k,p) : U(i,j,k,p)))/((ymin || ymax ? 1.0 : 2.0)*DX[1]);,
						     gradu(p,2) = ((!zmax ? U(i,j,k+1,p) : U(i,j,k,p)) - (!zmin ? U(i,j,k-1,p) : U(i,j,k,p)))/((zmin || zmax ? 1.0 : 2.0)*DX[2]););
				}
					
				// Stress tensor computed using the model fab
				Set::Matrix sig = C(i,j,k)(gradu);

				// Boundary conditions
				amrex::IntVect m(AMREX_D_DECL(i,j,k));
				if (AMREX_D_TERM(xmax || xmin, || ymax || ymin, || zmax || zmin)) 
				{
					for (int p = 0; p < AMREX_SPACEDIM; p++) // iterate over DIMENSIONS
					{
						for (int q = 0; q < AMREX_SPACEDIM; q++) // iterate over FACES
						{
							if (m[q] == domain.loVect()[q])
							{
								if (m_bc_lo[q][p] == BC::Displacement)    f(p) =   U(i,j,k,p);
								else if (m_bc_lo[q][p] == BC::Traction)   f(p) += -sig(p,q);
								else if (m_bc_lo[q][p] == BC::Neumann)    f(p) += -gradu(p,q);
								else Util::Abort(INFO, "Invalid BC");
							}
							if (m[q] == domain.hiVect()[q])
							{
								if (m_bc_hi[q][p] == BC::Displacement)    f(p) = U(i,j,k,p);
								else if (m_bc_hi[q][p] == BC::Traction)   f(p) += +sig(p,q);
								else if (m_bc_hi[j][i] == BC::Neumann)    f(p) += +gradu(p,q);
								else Util::Abort(INFO, "Invalid BC");

							}
						}
					}
				}
				else
				{
					// The gradient of the displacement gradient tensor
					std::array<Set::Matrix,AMREX_SPACEDIM> gradgradu; // gradgradu[k](l,j) = u_{k,lj}

					// Fill gradu and gradgradu
					for (int p = 0; p < AMREX_SPACEDIM; p++)
					{
						// Diagonal terms:
						AMREX_D_TERM(gradgradu[p](0,0) = (U(i+1,j,k,p) - 2.0*U(i,j,k,p) + U(i-1,j,k,p))/DX[0]/DX[0];,
							     gradgradu[p](1,1) = (U(i,j+1,k,p) - 2.0*U(i,j,k,p) + U(i,j-1,k,p))/DX[1]/DX[1];,
							     gradgradu[p](2,2) = (U(i,j,k+1,p) - 2.0*U(i,j,k,p) + U(i,j,k-1,p))/DX[2]/DX[2];);

						// Off-diagonal terms:
						AMREX_D_TERM(,// 2D
							     gradgradu[p](0,1) = (U(i+1,j+1,k,p) + U(i-1,j-1,k,p) - U(i+1,j-1,k,p) - U(i-1,j+1,k,p))/(2.0*DX[0])/(2.0*DX[1]);
							     gradgradu[p](1,0) = gradgradu[i](0,1);
							     ,// 3D
							     gradgradu[p](0,2) = (U(i+1,j,k+1,p) + U(i-1,j,k-1,p) - U(i+1,j,k-1,p) - U(i-1,j,k+1,p))/(2.0*DX[0])/(2.0*DX[2]);
							     gradgradu[p](1,2) = (U(i,j+1,k+1,p) + U(i,j-1,k-1,p) - U(i,j+1,k-1,p) - U(i,j-1,k+1,p))/(2.0*DX[1])/(2.0*DX[2]);
							     gradgradu[p](2,0) = gradgradu[i](0,2);
							     gradgradu[p](2,1) = gradgradu[i](1,2););
					}
	
					//
					// Operator
					//
					// The return value is
					//    f = C(grad grad u) + grad(C)*grad(u)
					// In index notation
					//    f_i = C_{ijkl,j} u_{k,l}  +  C_{ijkl}u_{k,lj}
					//

					f =     C(i,j,k)(gradgradu) + 
						AMREX_D_TERM(( ( C(i+1,j,k) - C(i-1,j,k))/2.0/DX[0])(gradu).col(0),
							     + ((C(i,j+1,k) - C(i,j-1,k))/2.0/DX[1])(gradu).col(1),
							     + ((C(i,j,k+1) - C(i,j,k-1))/2.0/DX[2])(gradu).col(2));

				}

				AMREX_D_TERM(F(i,j,k,0) = f[0];, F(i,j,k,1) = f[1];, F(i,j,k,2) = f[2];);

			});
	}
}



template<class T>
void
Elastic<T>::Diagonal (int amrlev, int mglev, MultiFab& diag)
{
	BL_PROFILE("Operator::Elastic::Diagonal()");

	amrex::Box domain(m_geom[amrlev][mglev].Domain());
	const Real* DX = m_geom[amrlev][mglev].CellSize();
	
	for (MFIter mfi(diag, false); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.validbox();
		amrex::BaseFab<T> &C = (*(model[amrlev][mglev]))[mfi];
		amrex::FArrayBox       &diagfab    = diag[mfi];

		AMREX_D_TERM(for (int m1 = bx.loVect()[0] - 1; m1<=bx.hiVect()[0] + 1; m1++),
		 	     for (int m2 = bx.loVect()[1] - 1; m2<=bx.hiVect()[1] + 1; m2++),
		 	     for (int m3 = bx.loVect()[2] - 1; m3<=bx.hiVect()[2] + 1; m3++))
		{
			amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));


			AMREX_D_TERM(if (m[0] < domain.loVect()[0]) continue;,
				     if (m[1] < domain.loVect()[1]) continue;,
				     if (m[2] < domain.loVect()[2]) continue;);
			AMREX_D_TERM(if (m[0] > domain.hiVect()[0]+1) continue;,
				     if (m[1] > domain.hiVect()[1]+1) continue;,
				     if (m[2] > domain.hiVect()[2]+1) continue;)
					     
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
				diagfab(m,i) = 0.0;
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

				Set::Matrix sig = C(m)(gradu);

				if (AMREX_D_TERM(xmax || xmin, || ymax || ymin, || zmax || zmin)) 
				{
					for (int j = 0; j < AMREX_SPACEDIM; j++) // iterate over FACES
					{
						if (m[j] == domain.loVect()[j])
						{
							if (m_bc_lo[j][i] == BC::Displacement)
								diagfab(m,i) += 1.0;
							else if (m_bc_lo[j][i] == BC::Traction) 
								diagfab(m,i) -= sig(i,j);
							else if (m_bc_lo[j][i] == BC::Neumann) 
								diagfab(m,i) -= gradu(i,j);
							else Util::Abort(INFO, "Invalid BC");
						}
						if (m[j] == domain.hiVect()[j] + 1)
						{
							if (m_bc_hi[j][i] == BC::Displacement)
								diagfab(m,i) += 1.0;
							else if (m_bc_hi[j][i] == BC::Traction) 
								diagfab(m,i) += sig(i,j);
							else if (m_bc_lo[j][i] == BC::Neumann) 
								diagfab(m,i) += gradu(i,j);
							else Util::Abort(INFO, "Invalid BC");
						}
					}
				}
				else
				{
					Set::Vector f =
						C(m)(gradgradu)  + 
						AMREX_D_TERM(((C(m+dx) - C(m-dx))/2.0/DX[0])(gradu).col(0),
						   	     + ((C(m+dy) - C(m-dy))/2.0/DX[1])(gradu).col(1),
						     	     + ((C(m+dz) - C(m-dz))/2.0/DX[2])(gradu).col(2));
					diagfab(m,i) += f(i);
				}
			}
		}
	}
}


template<class T>
void
Elastic<T>::Error0x (int amrlev, int mglev, MultiFab& R0x, const MultiFab& x) const
{
	BL_PROFILE("Operator::Elastic::Error0x()");
	Util::Message(INFO);

	int ncomp = x.nComp();//getNComp();
	int nghost = x.nGrow();

	if (!m_diagonal_computed)
		Util::Abort(INFO,"Operator::Diagonal() must be called before using normalize");

	amrex::MultiFab D0x(x.boxArray(), x.DistributionMap(), ncomp, nghost);
	amrex::MultiFab AD0x(x.boxArray(), x.DistributionMap(), ncomp, nghost);

	amrex::MultiFab::Copy(D0x,x,0,0,ncomp,nghost); // D0x = x
	amrex::MultiFab::Divide(D0x,*m_diag[amrlev][mglev],0,0,ncomp,0); // D0x = x/diag
	amrex::MultiFab::Copy(AD0x,D0x,0,0,ncomp,nghost); // AD0x = D0x

	Fapply(amrlev,mglev,AD0x,D0x);	// AD0x = A * D0 * x
	
	amrex::MultiFab::Copy(R0x,x,0,0,ncomp,nghost); // R0x = x
	amrex::MultiFab::Subtract(R0x,AD0x,0,0,ncomp,nghost); // R0x = x - AD0x
}


template<class T>
void
Elastic<T>::FFlux (int /*amrlev*/, const MFIter& /*mfi*/,
		const std::array<FArrayBox*,AMREX_SPACEDIM>& sigmafab,
		const FArrayBox& /*ufab*/, const int /*face_only*/) const
{
	BL_PROFILE("Operator::Elastic::FFlux()");
	Util::Message(INFO);
	amrex::BaseFab<amrex::Real> AMREX_D_DECL( &fxfab = *sigmafab[0],
	 					  &fyfab = *sigmafab[1],
	 					  &fzfab = *sigmafab[2] ) ;
	AMREX_D_TERM(fxfab.setVal(0.0);,
	 	     fyfab.setVal(0.0);,
	 	     fzfab.setVal(0.0););

}

template<class T>
void
Elastic<T>::Strain  (int amrlev,
		    amrex::MultiFab& eps,
		    const amrex::MultiFab& u,
		    bool voigt) const
{
	BL_PROFILE("Operator::Elastic::Strain()");
	Util::Message(INFO);

	if (voigt)
		AMREX_ASSERT(eps.nComp() == (AMREX_SPACEDIM*(AMREX_SPACEDIM-1)/2));
	else
		AMREX_ASSERT(eps.nComp() == AMREX_SPACEDIM*AMREX_SPACEDIM);
	
	const amrex::Real* DX = m_geom[amrlev][0].CellSize();
	
	for (MFIter mfi(u, true); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();
		amrex::FArrayBox  &epsfab   = eps[mfi];
		const amrex::FArrayBox  &ufab = u[mfi];
		
		AMREX_D_TERM(for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++),
			     for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++),
			     for (int m3 = bx.loVect()[2]; m3<=bx.hiVect()[2]; m3++))
		{
			amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));
			bool	AMREX_D_DECL(xmin = (m1 == bx.loVect()[0]),
					     ymin = (m2 == bx.loVect()[1]),
					     zmin = (m3 == bx.loVect()[2])),
				AMREX_D_DECL(xmax = (m1 == bx.hiVect()[0]),
					     ymax = (m2 == bx.hiVect()[1]),
					     zmax = (m3 == bx.hiVect()[2]));

			Set::Matrix gradu;
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
				     gradu(2,2) = ((!zmax ? ufab(m+dz,2) : ufab(m,2)) - (!zmin ? ufab(m-dz,2) : ufab(m,2)))/((zmin || zmax ? 1.0 : 2.0)*DX[2]););
			
			Set::Matrix strain = 0.5 * (gradu + gradu.transpose());
			
			if (voigt)
			{
#if AMREX_SPACEDIM == 2
				epsfab(m,0) = strain(0,0); epsfab(m,1) = strain(1,1); epsfab(m,2) = strain(0,1); 
#elif AMREX_SPACEDIM == 3
				epsfab(m,0) = strain(0,0); epsfab(m,1) = strain(1,1); epsfab(m,2) = strain(2,2); 
				epsfab(m,3) = strain(1,2); epsfab(m,4) = strain(2,0); epsfab(m,5) = strain(0,1); 
#endif
			}
			else
			{
#if   AMREX_SPACEDIM == 2
				epsfab(m,0) = strain(0,0); epsfab(m,1) = strain(0,1); 
				epsfab(m,2) = strain(1,0); epsfab(m,3) = strain(1,1); 
#elif AMREX_SPACEDIM == 3
				epsfab(m,0) = strain(0,0); epsfab(m,1) = strain(0,1); epsfab(m,2) = strain(0,2); 
				epsfab(m,3) = strain(1,0); epsfab(m,4) = strain(1,1); epsfab(m,5) = strain(1,2); 
				epsfab(m,6) = strain(2,0); epsfab(m,7) = strain(2,1); epsfab(m,8) = strain(2,2); 
#endif
			}
		}
	}
}


template<class T>
void
Elastic<T>::Stress (int amrlev,
		    amrex::MultiFab& a_sigma,
		    const amrex::MultiFab& a_u,
		    bool voigt) const
{
	BL_PROFILE("Operator::Elastic::Stress()");

	const amrex::Real* DX = m_geom[amrlev][0].CellSize();
	amrex::Box domain(m_geom[amrlev][0].Domain());
	domain.convert(amrex::IntVect::TheNodeVector());

	for (MFIter mfi(a_u, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();
		amrex::Array4<T> const& C                 = (*(model[amrlev][0])).array(mfi);
		amrex::Array4<amrex::Real> const& sigma   = a_sigma.array(mfi);
		amrex::Array4<const amrex::Real> const& u = a_u.array(mfi);
		const Dim3 lo= amrex::lbound(domain), hi = amrex::ubound(domain);

		amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k)
				    {
					    Set::Matrix gradu;

					    bool    AMREX_D_DECL(xmin = (i == lo.x), ymin = (j==lo.y), zmin = (k==lo.z)),
						    AMREX_D_DECL(xmax = (i == hi.x), ymax = (j==hi.y), zmax = (k==hi.z));

					    AMREX_D_TERM(gradu(0,0) = ((!xmax ? u(i+1,j,k,0) : u(i,j,k,0)) - (!xmin ? u(i-1,j,k,0) : u(i,j,k,0)))/((xmin || xmax ? 1.0 : 2.0)*DX[0]);
							 ,
							 gradu(0,1) = ((!ymax ? u(i,j+1,k,0) : u(i,j,k,0)) - (!ymin ? u(i,j-1,k,0) : u(i,j,k,0)))/((ymin || ymax ? 1.0 : 2.0)*DX[1]);
							 gradu(1,0) = ((!xmax ? u(i+1,j,k,1) : u(i,j,k,1)) - (!xmin ? u(i-1,j,k,1) : u(i,j,k,1)))/((xmin || xmax ? 1.0 : 2.0)*DX[0]);
							 gradu(1,1) = ((!ymax ? u(i,j+1,k,1) : u(i,j,k,1)) - (!ymin ? u(i,j-1,k,1) : u(i,j,k,1)))/((ymin || ymax ? 1.0 : 2.0)*DX[1]);
							 ,
							 gradu(0,2) = ((!zmax ? u(i,j,k+1,0) : u(i,j,k,0)) - (!zmin ? u(i,j,k-1,0) : u(i,j,k,0)))/((zmin || zmax ? 1.0 : 2.0)*DX[2]);
							 gradu(1,2) = ((!zmax ? u(i,j,k+1,1) : u(i,j,k,1)) - (!zmin ? u(i,j,k-1,1) : u(i,j,k,1)))/((zmin || zmax ? 1.0 : 2.0)*DX[2]);
							 gradu(2,0) = ((!xmax ? u(i+1,j,k,2) : u(i,j,k,2)) - (!xmin ? u(i-1,j,k,2) : u(i,j,k,2)))/((xmin || xmax ? 1.0 : 2.0)*DX[0]);
							 gradu(2,1) = ((!ymax ? u(i,j+1,k,2) : u(i,j,k,2)) - (!ymin ? u(i,j-1,k,2) : u(i,j,k,2)))/((ymin || ymax ? 1.0 : 2.0)*DX[1]);
							 gradu(2,2) = ((!zmax ? u(i,j,k+1,2) : u(i,j,k,2)) - (!zmin ? u(i,j,k-1,2) : u(i,j,k,2)))/((zmin || zmax ? 1.0 : 2.0)*DX[2]););
					 
					    Set::Matrix sig = C(i,j,k)(gradu);

					    if (voigt)
					    {
						    AMREX_D_PICK(sigma(i,j,k,0) = sig(0,0);
								 ,
								 sigma(i,j,k,0) = sig(0,0); sigma(i,j,k,1) = sig(1,1); sigma(i,j,k,2) = sig(0,1); 
								 ,
								 sigma(i,j,k,0) = sig(0,0); sigma(i,j,k,1) = sig(1,1); sigma(i,j,k,2) = sig(2,2); 
								 sigma(i,j,k,3) = sig(1,2); sigma(i,j,k,4) = sig(2,0); sigma(i,j,k,5) = sig(0,1););
					    }
					    else
					    {
						    AMREX_D_PICK(sigma(i,j,k,0) = sig(0,0);
								 ,
								 sigma(i,j,k,0) = sig(0,0); sigma(i,j,k,1) = sig(0,1); 
								 sigma(i,j,k,2) = sig(1,0); sigma(i,j,k,3) = sig(1,1);
								 ,
								 sigma(i,j,k,0) = sig(0,0); sigma(i,j,k,1) = sig(0,1); sigma(i,j,k,2) = sig(0,2); 
								 sigma(i,j,k,3) = sig(1,0); sigma(i,j,k,4) = sig(1,1); sigma(i,j,k,5) = sig(1,2); 
								 sigma(i,j,k,6) = sig(2,0); sigma(i,j,k,7) = sig(2,1); sigma(i,j,k,8) = sig(2,2););
					    }
				    });
	}
}


template<class T>
void
Elastic<T>::Energy (int amrlev,
		    amrex::MultiFab& energy,
		    const amrex::MultiFab& u) const
{
	BL_PROFILE("Operator::Elastic::Energy()");
	Util::Message(INFO);
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
			bool    AMREX_D_DECL(xmin = (m1 == bx.loVect()[0]),
					     ymin = (m2 == bx.loVect()[1]),
					     zmin = (m3 == bx.loVect()[2])),
				AMREX_D_DECL(xmax = (m1 == bx.hiVect()[0]),
					     ymax = (m2 == bx.hiVect()[1]),
					     zmax = (m3 == bx.hiVect()[2]));
			

			Set::Matrix gradu;

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
				     gradu(2,2) = ((!zmax ? ufab(m+dz,2) : ufab(m,2)) - (!zmin ? ufab(m-dz,2) : ufab(m,2)))/((zmin || zmax ? 1.0 : 2.0)*DX[2]););
			
			Set::Matrix eps = 0.5 * (gradu + gradu.transpose());
			Set::Matrix sig = C(m)(eps);

			energyfab(m) = (eps.transpose() * sig).trace();

		}
	}
}



template<class T>
void
Elastic<T>::averageDownCoeffs ()
{
	BL_PROFILE("Elastic::averageDownCoeffs()");
	
	// for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
	// {
	// 	for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
	// 	{
	// 		///\todo replace number of ghost cells with 0
	// 		///\todo I think we can erase this section.
	// 		model[amrlev][mglev].reset(new amrex::FabArray<amrex::BaseFab<T> >(m_grids[amrlev][mglev], m_dmap[amrlev][mglev], 1, 1));
	// 	}
	// }

	for (int amrlev = m_num_amr_levels-1; amrlev > 0; --amrlev)
	{
		averageDownCoeffsSameAmrLevel(amrlev);
	 	averageDownCoeffsToCoarseAmrLevel(amrlev);
	}

	averageDownCoeffsSameAmrLevel(0);
	for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
	{
	 	for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
	 	{
	 		if (model[amrlev][mglev]) {
	 			FillBoundaryCoeff(*model[amrlev][mglev], m_geom[amrlev][mglev]);
	 		}
	 	}
	}
}

template<class T>
void
Elastic<T>::averageDownCoeffsToCoarseAmrLevel (int /*flev*/) 
{
	/*
	BL_PROFILE("Operator::Elastic::averageDownCoeffsToCoarseAmrLevel()");

	//const int mglev = 0;

	// const int idim = 0;  // other dimensions are just aliases

	// // amrex::average_down(*m_sigma[flev][mglev][idim], *m_sigma[flev-1][mglev][idim], 0, 1,
	// // 		    m_amr_ref_ratio[flev-1]);

	int ncomp = AMREX_SPACEDIM;
	MultiTab finemt;
	MultiTab crsemt;
	
	

	amrex::BoxArray fineBACoarsened = finemt.boxArray(); fineBACoarsened.coarsen(m_amr_ref_ratio[flev-1]);


	//MultiFab crse_S_fine(crse_S_fine_BA, S_fine.DistributionMap(), ncomp, 0, MFInfo(), FArrayBoxFactory());
        if ((fineBACoarsened == crsemt.boxArray()) &&
	    (finemt.DistributionMap() == crsemt.DistributionMap()))
	{
		Util::Message(INFO); // this never seems to happen
	}
	else
	{
		//Util::Abort(INFO, "difference in box arrays");

		MultiTab finemtcoarsened(fineBACoarsened, finemt.DistributionMap(), ncomp, 0);

// #ifdef _OPENMP
// #pragma omp parallel
// #endif
		for (MFIter mfi(finemtcoarsened,true); mfi.isValid(); ++mfi)
		{
			//  NOTE: The tilebox is defined at the coarse level.
			const Box& bx = mfi.tilebox();
                
			const TArrayBox &fine = finemtcoarsened[mfi];
			TArrayBox &crse = crsemt[mfi];

			//  NOTE: We copy from component scomp of the fine fab into component 0 of the crse fab
			//        because the crse fab is a temporary which was made starting at comp 0, it is
			//        not part of the actual crse multifab which came in.


			AMREX_D_TERM(for (int m1 = bx.loVect()[0]-1; m1<=bx.hiVect()[0]+1; m1++),
				     for (int m2 = bx.loVect()[1]-1; m2<=bx.hiVect()[1]+1; m2++),
				     for (int m3 = bx.loVect()[2]-1; m3<=bx.hiVect()[2]+1; m3++))
			{
				amrex::IntVect m_crse(AMREX_D_DECL(m1,m2,m3));
				amrex::IntVect m_fine(AMREX_D_DECL(m1*2,m2*2,m3*2));

				crse(m_crse) = fine(m_fine);
			}
		}
            
		crsemt.copy(finemtcoarsened,0,0,ncomp);

	}
	*/
}

template<class T>
void
Elastic<T>::averageDownCoeffsSameAmrLevel (int amrlev)
{
	BL_PROFILE("Elastic::averageDownCoeffsSameAmrLevel()");

// 	if (m_coarsening_strategy != CoarseningStrategy::Sigma) return;

// 	const int nsigma = (m_use_harmonic_average) ? AMREX_SPACEDIM : 1;

 	for (int mglev = 1; mglev < m_num_mg_levels[amrlev]; ++mglev)
 	{
		amrex::Box domain_crse(m_geom[amrlev][mglev].Domain()); domain_crse.convert(amrex::IntVect::TheNodeVector());
		amrex::Box domain_fine(m_geom[amrlev][mglev-1].Domain()); domain_fine.convert(amrex::IntVect::TheNodeVector());

		MultiTab& crse = *model[amrlev][mglev];
		MultiTab& fine = *model[amrlev][mglev-1];
		
		amrex::BoxArray crseba = crse.boxArray();
		amrex::BoxArray fineba = fine.boxArray();

		// Util::Message(INFO,crseba);
		// Util::Message(INFO,fineba);
		//fineba.grow(2);
		//Util::Message(INFO,fineba);
		// bool isMFIterSafe  = (crse.DistributionMap() == fine.DistributionMap()) && BoxArray::SameRefs(crseba,fineba);
		// bool need_parallel_copy = !isMFIterSafe;
		// MultiTab cfine;
		// if (need_parallel_copy) {
		// 	//const BoxArray& ba = amrex::coarsen(fine.boxArray(), 2);
		// 	const BoxArray& ba = amrex::coarsen(fineba, 2);
		// 	cfine.define(ba, fine.DistributionMap(), 1, 2);
		// }


		// if (need_parallel_copy) {
		// 	crse.ParallelCopy(cfine);
		// }
		
		BoxArray newba = crseba;
		newba.refine(2);
		MultiTab fine_on_crseba;
		fine_on_crseba.define(newba,crse.DistributionMap(),1,4);
		fine_on_crseba.ParallelCopy(fine,0,0,1,2,4,m_geom[amrlev][mglev].periodicity());

		for (MFIter mfi(crse, false); mfi.isValid(); ++mfi)
 			{
				const Box& bx = mfi.validbox() & domain_crse;

				// TArrayBox &crsetab = (*pcrse)[mfi];
				// TArrayBox &finetab = fine[mfi];
				TArrayBox &crsetab = crse[mfi];
				TArrayBox &finetab = fine_on_crseba[mfi];

				// Util::Message(INFO,bx);
				// Util::Message(INFO,bx.loVect()[0]," ",bx.hiVect()[0]);
				// Util::Message(INFO,finetab.box());


				// Util::Abort(INFO);
				
				AMREX_D_TERM(for (int m1 = bx.loVect()[0]-2; m1<=bx.hiVect()[0]+2; m1++),
					     for (int m2 = bx.loVect()[1]-2; m2<=bx.hiVect()[1]+2; m2++),
					     for (int m3 = bx.loVect()[2]-2; m3<=bx.hiVect()[2]+2; m3++))
				{

					amrex::IntVect m_crse(AMREX_D_DECL(m1,m2,m3));
					amrex::IntVect m_fine(AMREX_D_DECL(m1*2,m2*2,m3*2));


					bool    AMREX_D_DECL(xmin = (m_crse[0] <= domain_crse.loVect()[0]) || (m_crse[0] == bx.loVect()[0]-2),
							     ymin = (m_crse[1] <= domain_crse.loVect()[1]) || (m_crse[1] == bx.loVect()[1]-2),
							     zmin = (m_crse[2] <= domain_crse.loVect()[2]) || (m_crse[2] == bx.loVect()[2]-2)),
						AMREX_D_DECL(xmax = (m_crse[0] >= domain_crse.hiVect()[0]) || (m_crse[0] == bx.hiVect()[0]+2),
							     ymax = (m_crse[1] >= domain_crse.hiVect()[1]) || (m_crse[1] == bx.hiVect()[1]+2),
							     zmax = (m_crse[2] >= domain_crse.hiVect()[2]) || (m_crse[2] == bx.hiVect()[2]+2));


					
					// AMREX_D_TERM(if (m1 == bx.loVect()[0] - 1) ++m_fine[0];
					// 	     if (m2 == bx.loVect()[1] - 1) ++m_fine[1];
					// 	     if (m3 == bx.loVect()[2] - 1) ++m_fine[1];
					// if (m1 == bx.hiVect()[0] + 1) --m_fine[0];
					// if (m2 == bx.hiVect()[1] + 1) --m_fine[1];
					
#if AMREX_SPACEDIM == 2
					if ((xmin || xmax) && (ymin || ymax)) // corner
					{
						crsetab(m_crse) = finetab(m_fine);
					}
					else if (ymin || ymax) // x edge
					{
						crsetab(m_crse) = finetab(m_fine-dx)*0.25 + finetab(m_fine)*0.5 + finetab(m_fine+dx)*0.25;
					}
					else if (xmin || xmax) // y edge
					{
						crsetab(m_crse) = finetab(m_fine-dy)*0.25 + finetab(m_fine)*0.5 + finetab(m_fine+dy)*0.25;
					}
					else 
					{
						crsetab(m_crse) =
							(finetab(m_fine+dx+dy) + finetab(m_fine+dx-dy) + finetab(m_fine-dx+dy) + finetab(m_fine-dx-dy)) / 16. +
							(finetab(m_fine+dx) + finetab(m_fine-dx) + finetab(m_fine+dy) + finetab(m_fine-dy)) / 8. +
							(finetab(m_fine)) / 4.;
					}

#elif AMREX_SPACEDIM == 3
					//crsetab(m_crse) = finetab(m_fine);
					if ((xmin || xmax) && (ymin || ymax) && (zmin || zmax)) // corner
					{
						crsetab(m_crse) = finetab(m_fine);
					}
					else if ((ymin || ymax) && // x edge
						 (zmin || zmax))
					{
						crsetab(m_crse) = finetab(m_fine-dx)*0.25 + finetab(m_fine)*0.5 + finetab(m_fine+dx)*0.25;
					}
					else if ((xmin || xmax) && // y edge
						 (zmin || zmax))
					{
						crsetab(m_crse) = finetab(m_fine-dy)*0.25 + finetab(m_fine)*0.5 + finetab(m_fine+dy)*0.25;
					}
					else if ((xmin || xmax) && // z edge
						 (ymin || ymax))
					{
						crsetab(m_crse) = finetab(m_fine-dz)*0.25 + finetab(m_fine)*0.5 + finetab(m_fine+dz)*0.25;
					}
					else if ((xmin || xmax)) // x face
					{
						crsetab(m_crse) =
							(finetab(m_fine+dy+dz) + finetab(m_fine+dy-dz) + finetab(m_fine-dy+dz) + finetab(m_fine-dy-dz)) / 16. +
							(finetab(m_fine+dy) + finetab(m_fine-dy) + finetab(m_fine+dz) + finetab(m_fine-dz)) / 8. +
							(finetab(m_fine)) / 4.;							
					}
					else if ((ymin || ymax)) // y face
					{
						crsetab(m_crse) =
							(finetab(m_fine+dz+dx) + finetab(m_fine+dz-dx) + finetab(m_fine-dz+dx) + finetab(m_fine-dz-dx)) / 16. +
							(finetab(m_fine+dz) + finetab(m_fine-dz) + finetab(m_fine+dx) + finetab(m_fine-dx)) / 8. +
							(finetab(m_fine)) / 4.;
					}
					else if ((zmin || zmax)) // z face
					{
						crsetab(m_crse) =
							(finetab(m_fine+dx+dy) + finetab(m_fine+dx-dy) + finetab(m_fine-dx+dy) + finetab(m_fine-dx-dy)) / 16. +
							(finetab(m_fine+dx) + finetab(m_fine-dx) + finetab(m_fine+dy) + finetab(m_fine-dy)) / 8. +
							(finetab(m_fine)) / 4.;
					}
					else
					{
						crsetab(m_crse) =
							(finetab(m_fine-dx-dy-dz) + finetab(m_fine-dx-dy+dz) + finetab(m_fine-dx+dy-dz) + finetab(m_fine-dx+dy+dz) +
							 finetab(m_fine+dx-dy-dz) + finetab(m_fine+dx-dy+dz) + finetab(m_fine+dx+dy-dz) + finetab(m_fine+dx+dy+dz)) / 64.0
							+
							(finetab(m_fine-dy-dz) + finetab(m_fine-dy+dz) + finetab(m_fine+dy-dz) + finetab(m_fine+dy+dz) +
							 finetab(m_fine-dz-dx) + finetab(m_fine-dz+dx) + finetab(m_fine+dz-dx) + finetab(m_fine+dz+dx) +
							 finetab(m_fine-dx-dy) + finetab(m_fine-dx+dy) + finetab(m_fine+dx-dy) + finetab(m_fine+dx+dy)) / 32.0
							+
							(finetab(m_fine-dx) + finetab(m_fine-dy) + finetab(m_fine-dz) +
							 finetab(m_fine+dx) + finetab(m_fine+dy) + finetab(m_fine+dz)) / 16.0
							+
							finetab(m_fine) / 8.0;
					}
#endif
					

				}
 			}

		// if (need_parallel_copy) {
		// 	crse.ParallelCopy(cfine);
		// }
 	}
}

template<class T>
void
Elastic<T>::FillBoundaryCoeff (MultiTab& sigma, const Geometry& geom)
{
	BL_PROFILE("Elastic::FillBoundaryCoeff()");

	//sigma.FillBoundary(geom.periodicity());
	for (int i = 0; i < 4; i++)
	{
		MultiTab & mf = sigma;
		mf.FillBoundary(geom.periodicity());
		const int ncomp = mf.nComp();
		const int ng1 = 1;
		const int ng2 = 2;
		MultiTab tmpmf(mf.boxArray(), mf.DistributionMap(), ncomp, ng1);
		//MultiTab::Copy(tmpmf, mf, 0, 0, ncomp, ng1); 
	  	tmpmf.copy(mf,0,0,ncomp,ng2,ng1,geom.periodicity());

		mf.ParallelCopy   (tmpmf, 0, 0, ncomp, ng1, ng2, geom.periodicity());
	}


	//const Box& domain = geom.Domain();

// #ifdef _OPENMP
// #pragma omp parallel
// #endif
	// for (MFIter mfi(sigma, MFItInfo().SetDynamic(true)); mfi.isValid(); ++mfi)
	// {
	// 	if (!domain.contains(mfi.fabbox()))
	// 	{
			

	// 	}
	// }
	///////Util::Warning(INFO, "FillBoundaryCoeff not fully implemented");
}



template class Elastic<Model::Solid::LinearElastic::Isotropic>;
template class Elastic<Model::Solid::LinearElastic::Cubic>;
template class Elastic<Model::Solid::LinearElastic::Laplacian>;
template class Elastic<Model::Solid::LinearElastic::Degradable::Isotropic>;
template class Elastic<Model::Solid::Viscoelastic::Isotropic>;
}

