// TODO: Remove these 

#include "Elastic.H"
#include "Set/Set.H"

#include "Numeric/Stencil.H"
namespace Operator
{
template<int SYM>
Elastic<SYM>::Elastic (const Vector<Geometry>& a_geom,
		     const Vector<BoxArray>& a_grids,
		     const Vector<DistributionMapping>& a_dmap,
		     const LPInfo& a_info)
{
	BL_PROFILE("Operator::Elastic::Elastic()");

	define(a_geom, a_grids, a_dmap, a_info);
}

template<int SYM>
Elastic<SYM>::~Elastic ()
{}

template<int SYM>
void
Elastic<SYM>::define (const Vector<Geometry>& a_geom,
		    const Vector<BoxArray>& a_grids,
		    const Vector<DistributionMapping>& a_dmap,
		    const LPInfo& a_info,
		    const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
	BL_PROFILE("Operator::Elastic::define()");

	Operator::define(a_geom,a_grids,a_dmap,a_info,a_factory);

	int model_nghost = 2;

	m_ddw_mf.resize(m_num_amr_levels);
	for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
	{
		m_ddw_mf[amrlev].resize(m_num_mg_levels[amrlev]);
		for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
		{
			m_ddw_mf[amrlev][mglev].reset(new MultiTab(amrex::convert(m_grids[amrlev][mglev],
									       amrex::IntVect::TheNodeVector()),
								m_dmap[amrlev][mglev], 1, model_nghost));
		}
	}
}

template <int SYM>
void 
Elastic<SYM>::SetModel (MATRIX4 &a_model)
{
	for (int amrlev = 0; amrlev < m_num_amr_levels; amrlev++)
	{
		amrex::Box domain(m_geom[amrlev][0].Domain());
		domain.convert(amrex::IntVect::TheNodeVector());

		int nghost = m_ddw_mf[amrlev][0]->nGrow();

		for (MFIter mfi(*m_ddw_mf[amrlev][0], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
		{
			Box bx = mfi.tilebox();
			bx.grow(nghost);   // Expand to cover first layer of ghost nodes
			bx = bx & domain;  // Take intersection of box and the problem domain
				
			amrex::Array4<MATRIX4> const& ddw         = (*(m_ddw_mf[amrlev][0])).array(mfi);
	
			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {
					ddw(i,j,k) = a_model;
				});
		}
	}
	m_model_set = true;
}

template <int SYM>
void
Elastic<SYM>::SetModel (int amrlev, const amrex::FabArray<amrex::BaseFab<MATRIX4> >& a_model)
{
	BL_PROFILE("Operator::Elastic::SetModel()");

	amrex::Box domain(m_geom[amrlev][0].Domain());
	domain.convert(amrex::IntVect::TheNodeVector());

	if (a_model.boxArray()        != m_ddw_mf[amrlev][0]->boxArray()) Util::Abort(INFO,"Inconsistent box arrays\n","a_model.boxArray()=\n",a_model.boxArray(),"\n but the current box array is \n",m_ddw_mf[amrlev][0]->boxArray());
	if (a_model.DistributionMap() != m_ddw_mf[amrlev][0]->DistributionMap()) Util::Abort(INFO,"Inconsistent distribution maps");
	if (a_model.nComp()           != m_ddw_mf[amrlev][0]->nComp()) Util::Abort(INFO,"Inconsistent # of components - should be ",m_ddw_mf[amrlev][0]->nComp());
	if (a_model.nGrow()           != m_ddw_mf[amrlev][0]->nGrow()) Util::Abort(INFO,"Inconsistent # of ghost nodes, should be ",m_ddw_mf[amrlev][0]->nGrow());


	int nghost = m_ddw_mf[amrlev][0]->nGrow();

	for (MFIter mfi(a_model, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		Box bx = mfi.tilebox();
		bx.grow(nghost);   // Expand to cover first layer of ghost nodes
		bx = bx & domain;  // Take intersection of box and the problem domain
			
		amrex::Array4<MATRIX4> const& C         = (*(m_ddw_mf[amrlev][0])).array(mfi);
		amrex::Array4<const MATRIX4> const& a_C = a_model.array(mfi);

		amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {
				C(i,j,k) = a_C(i,j,k);
			});
	}
	//FillBoundaryCoeff(*model[amrlev][0], m_geom[amrlev][0]);


	m_model_set = true;
}

template<int SYM>
void
Elastic<SYM>::Fapply (int amrlev, int mglev, MultiFab& a_f, const MultiFab& a_u) const
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
			
		amrex::Array4<MATRIX4> const& DDW                 = (*(m_ddw_mf[amrlev][mglev])).array(mfi);
		amrex::Array4<const amrex::Real> const& U = a_u.array(mfi);
		amrex::Array4<amrex::Real> const& F       = a_f.array(mfi);

		const Dim3 lo= amrex::lbound(domain), hi = amrex::ubound(domain);
			
		amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {
					
				Set::Vector f = Set::Vector::Zero();

				Set::Vector u;
				for (int p = 0; p < AMREX_SPACEDIM; p++) u(p) = U(i,j,k,p);
				

				bool AMREX_D_DECL(xmin = (i == lo.x), ymin = (j==lo.y), zmin = (k==lo.z)),
					 AMREX_D_DECL(xmax = (i == hi.x), ymax = (j==hi.y), zmax = (k==hi.z));

				// Determine if a special stencil will be necessary for first derivatives
				std::array<Numeric::StencilType,AMREX_SPACEDIM>
					sten = Numeric::GetStencil(i,j,k,domain);

				// The displacement gradient tensor
				Set::Matrix gradu; // gradu(i,j) = u_{i,j)

				// Fill gradu
				for (int p = 0; p < AMREX_SPACEDIM; p++)
				{
 					AMREX_D_TERM(gradu(p,0) = (Numeric::Stencil<Set::Scalar,1,0,0>::D(U,i,j,k,p,DX,sten));,
					 	     gradu(p,1) = (Numeric::Stencil<Set::Scalar,0,1,0>::D(U,i,j,k,p,DX,sten));,
					 	     gradu(p,2) = (Numeric::Stencil<Set::Scalar,0,0,1>::D(U,i,j,k,p,DX,sten)););
				}
					
				// Stress tensor computed using the model fab
				Set::Matrix sig = DDW(i,j,k)*gradu;

				// Boundary conditions
				/// \todo Important: we need a way to handle corners and edges.
				amrex::IntVect m(AMREX_D_DECL(i,j,k));
				if (AMREX_D_TERM(xmax || xmin, || ymax || ymin, || zmax || zmin)) 
				{
					f = (*m_bc)(u,gradu,sig,i,j,k,domain);
				}
				else
				{
					

					// The gradient of the displacement gradient tensor
					Set::Matrix3 gradgradu; // gradgradu[k](l,j) = u_{k,lj}

					// Fill gradu and gradgradu
					for (int p = 0; p < AMREX_SPACEDIM; p++)
					{
						// Diagonal terms:
						AMREX_D_TERM(gradgradu(p,0,0) = (Numeric::Stencil<Set::Scalar,2,0,0>::D(U,i,j,k,p,DX));,
							     gradgradu(p,1,1) = (Numeric::Stencil<Set::Scalar,0,2,0>::D(U,i,j,k,p,DX));,
							     gradgradu(p,2,2) = (Numeric::Stencil<Set::Scalar,0,0,2>::D(U,i,j,k,p,DX)););

						// Off-diagonal terms:
						AMREX_D_TERM(,// 2D
							     gradgradu(p,0,1) = (Numeric::Stencil<Set::Scalar,1,1,0>::D(U, i,j,k,p, DX));
							     gradgradu(p,1,0) = gradgradu(p,0,1);
							     ,// 3D
							     gradgradu(p,0,2) = (Numeric::Stencil<Set::Scalar,1,0,1>::D(U, i,j,k,p, DX));
							     gradgradu(p,1,2) = (Numeric::Stencil<Set::Scalar,0,1,1>::D(U, i,j,k,p, DX));
							     gradgradu(p,2,0) = gradgradu(p,0,2);
							     gradgradu(p,2,1) = gradgradu(p,1,2););
					}
	
					//
					// Operator
					//
					// The return value is
					//    f = C(grad grad u) + grad(C)*grad(u)
					// In index notation
					//    f_i = C_{ijkl,j} u_{k,l}  +  C_{ijkl}u_{k,lj}
					//

					f = DDW(i,j,k)*gradgradu;

					if (!m_uniform)
					{
						MATRIX4
						AMREX_D_DECL(Cgrad1 = (Numeric::Stencil<MATRIX4,1,0,0>::D(DDW,i,j,k,0,DX,sten)),
							    	 Cgrad2 = (Numeric::Stencil<MATRIX4,0,1,0>::D(DDW,i,j,k,0,DX,sten)),
							         Cgrad3 = (Numeric::Stencil<MATRIX4,0,0,1>::D(DDW,i,j,k,0,DX,sten)));
						f += AMREX_D_TERM((Cgrad1*gradu).col(0),
										 +(Cgrad2*gradu).col(1),
										 +(Cgrad3*gradu).col(2));
					}
				}
				AMREX_D_TERM(F(i,j,k,0) = f[0];, F(i,j,k,1) = f[1];, F(i,j,k,2) = f[2];);
			});
	}
}



template<int SYM>
void
Elastic<SYM>::Diagonal (int amrlev, int mglev, MultiFab& a_diag)
{
	BL_PROFILE("Operator::Elastic::Diagonal()");

	amrex::Box domain(m_geom[amrlev][mglev].Domain());
	domain.convert(amrex::IntVect::TheNodeVector());
	const Real* DX = m_geom[amrlev][mglev].CellSize();
	
	for (MFIter mfi(a_diag, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		Box bx = mfi.validbox();
		bx.grow(1);        // Expand to cover first layer of ghost nodes
		bx = bx & domain;  // Take intersection of box and the problem domain

		amrex::Array4<MATRIX4> const& DDW         = (*(m_ddw_mf[amrlev][mglev])).array(mfi);
		amrex::Array4<amrex::Real> const& diag    = a_diag.array(mfi);

		const Dim3 lo= amrex::lbound(domain), hi = amrex::ubound(domain);
			
		amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {

				Set::Vector f = Set::Vector::Zero();

				bool    AMREX_D_DECL(xmin = (i == lo.x), ymin = (j==lo.y), zmin = (k==lo.z)),
					    AMREX_D_DECL(xmax = (i == hi.x), ymax = (j==hi.y), zmax = (k==hi.z));

				std::array<Numeric::StencilType,AMREX_SPACEDIM> sten
					= Numeric::GetStencil(i,j,k,domain);


				

				Set::Matrix gradu; // gradu(i,j) = u_{i,j)
				Set::Matrix3 gradgradu; // gradgradu[k](l,j) = u_{k,lj}

				for (int p = 0; p < AMREX_SPACEDIM; p++)
				{

					diag(i,j,k,p) = 0.0;
					for (int q = 0; q < AMREX_SPACEDIM; q++)
					{
						AMREX_D_TERM(gradu(q,0) = ((!xmax ? 0.0 : (p==q ? 1.0 : 0.0)) - (!xmin ? 0.0 : (p==q ? 1.0 : 0.0)))/((xmin || xmax ? 1.0 : 2.0)*DX[0]);,
							     gradu(q,1) = ((!ymax ? 0.0 : (p==q ? 1.0 : 0.0)) - (!ymin ? 0.0 : (p==q ? 1.0 : 0.0)))/((ymin || ymax ? 1.0 : 2.0)*DX[1]);,
							     gradu(q,2) = ((!zmax ? 0.0 : (p==q ? 1.0 : 0.0)) - (!zmin ? 0.0 : (p==q ? 1.0 : 0.0)))/((zmin || zmax ? 1.0 : 2.0)*DX[2]););
			
						AMREX_D_TERM(gradgradu(q,0,0) = (p==q ? -2.0 : 0.0)/DX[0]/DX[0];
							     ,// 2D
							     gradgradu(q,0,1) = 0.0;
							     gradgradu(q,1,0) = 0.0;
							     gradgradu(q,1,1) = (p==q ? -2.0 : 0.0)/DX[1]/DX[1];
							     ,// 3D
							     gradgradu(q,0,2) = 0.0;
							     gradgradu(q,1,2) = 0.0;
							     gradgradu(q,2,0) = 0.0;
							     gradgradu(q,2,1) = 0.0;
							     gradgradu(q,2,2) = (p==q ? -2.0 : 0.0)/DX[2]/DX[2]);
					}

					Set::Matrix sig = DDW(i,j,k)*gradu;

					amrex::IntVect m(AMREX_D_DECL(i,j,k));
					if (AMREX_D_TERM(xmax || xmin, || ymax || ymin, || zmax || zmin)) 
					{
						Set::Vector u = Set::Vector::Zero();
						u(p) = 1.0;
						f = (*m_bc)(u,gradu,sig,i,j,k,domain);
						diag(i,j,k,p) = f(p);
					}
					else
					{
						Set::Matrix4<AMREX_SPACEDIM,SYM>
						AMREX_D_DECL(Cgrad1 = (Numeric::Stencil<Set::Matrix4<AMREX_SPACEDIM,SYM>,1,0,0>::D(DDW,i,j,k,0,DX,sten)),
						           	 Cgrad2 = (Numeric::Stencil<Set::Matrix4<AMREX_SPACEDIM,SYM>,0,1,0>::D(DDW,i,j,k,0,DX,sten)),
							         Cgrad3 = (Numeric::Stencil<Set::Matrix4<AMREX_SPACEDIM,SYM>,0,0,1>::D(DDW,i,j,k,0,DX,sten)));

						Set::Vector f = DDW(i,j,k)*gradgradu + 
							AMREX_D_TERM((Cgrad1*gradu).col(0),
										+(Cgrad2*gradu).col(1),
										+(Cgrad3*gradu).col(2));

						diag(i,j,k,p) += f(p);
					}
					if (std::isnan(diag(i,j,k,p))) Util::Abort(INFO,"diagonal is nan at (", i, ",", j , ",",k,"), amrlev=",amrlev,", mglev=",mglev);

				}
			});
	}
}


template<int SYM>
void
Elastic<SYM>::Error0x (int amrlev, int mglev, MultiFab& R0x, const MultiFab& x) const
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


template<int SYM>
void
Elastic<SYM>::FFlux (int /*amrlev*/, const MFIter& /*mfi*/,
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

template<int SYM>
void
Elastic<SYM>::Strain  (int amrlev,
		    amrex::MultiFab& a_eps,
		    const amrex::MultiFab& a_u,
		    bool voigt) const
{
	BL_PROFILE("Operator::Elastic::Strain()");

	const amrex::Real* DX = m_geom[amrlev][0].CellSize();
	amrex::Box domain(m_geom[amrlev][0].Domain());
	domain.convert(amrex::IntVect::TheNodeVector());


	for (MFIter mfi(a_u, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();
		amrex::Array4<amrex::Real> const& epsilon = a_eps.array(mfi);
		amrex::Array4<const amrex::Real> const& u = a_u.array(mfi);
		amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k)
				    {
					    Set::Matrix gradu;

					    std::array<Numeric::StencilType,AMREX_SPACEDIM> sten
						    = Numeric::GetStencil(i,j,k,domain);

					    // Fill gradu
					    for (int p = 0; p < AMREX_SPACEDIM; p++)
					    {
						    AMREX_D_TERM(gradu(p,0) = (Numeric::Stencil<Set::Scalar,1,0,0>::D(u, i,j,k,p, DX, sten));,
						    		 gradu(p,1) = (Numeric::Stencil<Set::Scalar,0,1,0>::D(u, i,j,k,p, DX, sten));,
						     		 gradu(p,2) = (Numeric::Stencil<Set::Scalar,0,0,1>::D(u, i,j,k,p, DX, sten)););
					    }

					    Set::Matrix eps = 0.5 * (gradu + gradu.transpose());

					    if (voigt)
					    {
						    AMREX_D_PICK(epsilon(i,j,k,0) = eps(0,0);
								 ,
								 epsilon(i,j,k,0) = eps(0,0); epsilon(i,j,k,1) = eps(1,1); epsilon(i,j,k,2) = eps(0,1); 
								 ,
								 epsilon(i,j,k,0) = eps(0,0); epsilon(i,j,k,1) = eps(1,1); epsilon(i,j,k,2) = eps(2,2); 
								 epsilon(i,j,k,3) = eps(1,2); epsilon(i,j,k,4) = eps(2,0); epsilon(i,j,k,5) = eps(0,1););
					    }
					    else
					    {
						    AMREX_D_PICK(epsilon(i,j,k,0) = eps(0,0);
								 ,
								 epsilon(i,j,k,0) = eps(0,0); epsilon(i,j,k,1) = eps(0,1); 
								 epsilon(i,j,k,2) = eps(1,0); epsilon(i,j,k,3) = eps(1,1);
								 ,
								 epsilon(i,j,k,0) = eps(0,0); epsilon(i,j,k,1) = eps(0,1); epsilon(i,j,k,2) = eps(0,2); 
								 epsilon(i,j,k,3) = eps(1,0); epsilon(i,j,k,4) = eps(1,1); epsilon(i,j,k,5) = eps(1,2); 
								 epsilon(i,j,k,6) = eps(2,0); epsilon(i,j,k,7) = eps(2,1); epsilon(i,j,k,8) = eps(2,2););
					    }
				    });
	}
}


template<int SYM>
void
Elastic<SYM>::Stress (int amrlev,
		    amrex::MultiFab& a_sigma,
		    const amrex::MultiFab& a_u,
		    bool voigt, bool a_homogeneous) 
{
	BL_PROFILE("Operator::Elastic::Stress()");
	SetHomogeneous(a_homogeneous);

	const amrex::Real* DX = m_geom[amrlev][0].CellSize();
	amrex::Box domain(m_geom[amrlev][0].Domain());
	domain.convert(amrex::IntVect::TheNodeVector());

	for (MFIter mfi(a_u, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();
		amrex::Array4<Set::Matrix4<AMREX_SPACEDIM,SYM>> const& DDW                 = (*(m_ddw_mf[amrlev][0])).array(mfi);
		amrex::Array4<amrex::Real> const& sigma   = a_sigma.array(mfi);
		amrex::Array4<const amrex::Real> const& u = a_u.array(mfi);
		amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k)
				    {
					    Set::Matrix gradu;

					    std::array<Numeric::StencilType,AMREX_SPACEDIM> sten
						    = Numeric::GetStencil(i,j,k,domain);

					    // Fill gradu
					    for (int p = 0; p < AMREX_SPACEDIM; p++)
					    {
						    AMREX_D_TERM(gradu(p,0) = (Numeric::Stencil<Set::Scalar,1,0,0>::D(u, i,j,k,p, DX, sten));,
						     		 gradu(p,1) = (Numeric::Stencil<Set::Scalar,0,1,0>::D(u, i,j,k,p, DX, sten));,
						      		 gradu(p,2) = (Numeric::Stencil<Set::Scalar,0,0,1>::D(u, i,j,k,p, DX, sten)););
					    }
					 
					    Set::Matrix sig = DDW(i,j,k)*gradu;

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


template<int SYM>
void
Elastic<SYM>::Energy (int amrlev,
		    amrex::MultiFab& a_energy,
		    const amrex::MultiFab& a_u, bool a_homogeneous)
{
	BL_PROFILE("Operator::Elastic::Energy()");
	SetHomogeneous(a_homogeneous);

	amrex::Box domain(m_geom[amrlev][0].Domain());
	domain.convert(amrex::IntVect::TheNodeVector());

	const amrex::Real* DX = m_geom[amrlev][0].CellSize();

	for (MFIter mfi(a_u, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();
		amrex::Array4<Set::Matrix4<AMREX_SPACEDIM,SYM>> const& DDW                  = (*(m_ddw_mf[amrlev][0])).array(mfi);
		amrex::Array4<amrex::Real> const& energy   = a_energy.array(mfi);
		amrex::Array4<const amrex::Real> const& u  = a_u.array(mfi);
		amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k)
				    {
					    Set::Matrix gradu;

					    std::array<Numeric::StencilType,AMREX_SPACEDIM> sten
						    = Numeric::GetStencil(i,j,k,domain);

					    // Fill gradu
					    for (int p = 0; p < AMREX_SPACEDIM; p++)
					    {
						    AMREX_D_TERM(gradu(p,0) = (Numeric::Stencil<Set::Scalar,1,0,0>::D(u, i,j,k,p, DX, sten));,
						     		 gradu(p,1) = (Numeric::Stencil<Set::Scalar,0,1,0>::D(u, i,j,k,p, DX, sten));,
						     		 gradu(p,2) = (Numeric::Stencil<Set::Scalar,0,0,1>::D(u, i,j,k,p, DX, sten)););
					    }

					 	Set::Matrix eps = .5 * (gradu + gradu.transpose());
					    Set::Matrix sig = DDW(i,j,k)*gradu;

					    // energy(i,j,k) = (gradu.transpose() * sig).trace();
						
						//Util::Abort(INFO,"Fix this"); //
						//energy(i,j,k) = C(i,j,k).W(gradu);
						for (int m = 0; m < AMREX_SPACEDIM; m++)
						{
							for(int n = 0; n < AMREX_SPACEDIM; n++)
							{
							energy(i,j,k) += .5 * sig(m,n) * eps(m,n);
							}
						}
				    });
	}
}

/*
template <int SYM>
void 
Elastic<SYM>::Energy (int amrlev, amrex::MultiFab& a_energies, const amrex::MultiFab& a_u, std::vector<T> a_models, bool a_homogeneous)
{
	BL_PROFILE("Operator::Elastic::Energy()");
	SetHomogeneous(a_homogeneous);

	if ((unsigned int)a_energies.nComp() != a_models.size())
	{
		Util::Abort(INFO,"Number of energy components (",a_energies.nComp(), ") does not equal number of models (",a_models.size(),")");
	}

	amrex::Box domain(m_geom[amrlev][0].Domain());
	domain.convert(amrex::IntVect::TheNodeVector());

	const amrex::Real* DX = m_geom[amrlev][0].CellSize();

	for (MFIter mfi(a_u, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();
		amrex::Array4<amrex::Real> const& energies   = a_energies.array(mfi);
		//amrex::Array4<amrex::Real> const& u  = a_u.array(mfi);
		amrex::Array4<const amrex::Real> const& u  = a_u.array(mfi);
		amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k)
	    {
		    Set::Matrix gradu;
		    std::array<Numeric::StencilType,AMREX_SPACEDIM> sten
			    		= Numeric::GetStencil(i,j,k,domain);
		    // Fill gradu
		    for (int p = 0; p < AMREX_SPACEDIM; p++)
		    {
			    AMREX_D_TERM(gradu(p,0) = (Numeric::Stencil<Set::Scalar,1,0,0>::D(u, i,j,k,p, DX, sten));,
			     		 gradu(p,1) = (Numeric::Stencil<Set::Scalar,0,1,0>::D(u, i,j,k,p, DX, sten));,
			     		 gradu(p,2) = (Numeric::Stencil<Set::Scalar,0,0,1>::D(u, i,j,k,p, DX, sten)););
		    }
		 	 
			for (unsigned int p = 0; p < a_models.size(); p++)
			{
				energies(i,j,k,p) = a_models[p].W(gradu);
			}
	    });
	}
}
*/



template<int SYM>
void
Elastic<SYM>::averageDownCoeffs ()
{
	BL_PROFILE("Elastic::averageDownCoeffs()");
	
	if (m_average_down_coeffs)
		for (int amrlev = m_num_amr_levels-1; amrlev > 0; --amrlev)
			averageDownCoeffsDifferentAmrLevels(amrlev);

	averageDownCoeffsSameAmrLevel(0);
	for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
	{
	 	for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
	 	{
	 		if (m_ddw_mf[amrlev][mglev]) {
	 			FillBoundaryCoeff(*m_ddw_mf[amrlev][mglev], m_geom[amrlev][mglev]);
	 		}
	 	}
	}
}

template<int SYM>
void
Elastic<SYM>::averageDownCoeffsDifferentAmrLevels (int fine_amrlev)
{
	BL_PROFILE("Operator::Elastic::averageDownCoeffsDifferentAmrLevels()");
	Util::Assert(INFO,TEST(fine_amrlev > 0));
	
	const int crse_amrlev = fine_amrlev - 1;
	const int ncomp = 1;

	MultiTab & crse_ddw = *m_ddw_mf[crse_amrlev][0];
	MultiTab & fine_ddw = *m_ddw_mf[fine_amrlev][0];

	amrex::Box cdomain(m_geom[crse_amrlev][0].Domain());
	cdomain.convert(amrex::IntVect::TheNodeVector());

	const Geometry& cgeom = m_geom[crse_amrlev  ][0];

 	const BoxArray&            fba = fine_ddw.boxArray();
 	const DistributionMapping& fdm = fine_ddw.DistributionMap();

 	MultiTab fine_ddw_for_coarse(amrex::coarsen(fba, 2), fdm, ncomp, 2);
	fine_ddw_for_coarse.ParallelCopy(crse_ddw,0,0,ncomp,0,0,cgeom.periodicity());

	const int coarse_fine_node = 1;
	const int fine_fine_node = 2;

	amrex::iMultiFab nodemask(amrex::coarsen(fba,2), fdm, 1, 2);
	nodemask.ParallelCopy(*m_nd_fine_mask[crse_amrlev],0,0,1,0,0,cgeom.periodicity());

	amrex::iMultiFab cellmask(amrex::convert(amrex::coarsen(fba,2),amrex::IntVect::TheCellVector()), fdm, 1, 2);
	cellmask.ParallelCopy(*m_cc_fine_mask[crse_amrlev],0,0,1,1,1,cgeom.periodicity());
	
	for (MFIter mfi(fine_ddw_for_coarse, false); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.validbox();

		amrex::Array4<const int> const& nmask = nodemask.array(mfi);
		//amrex::Array4<const int> const& cmask = cellmask.array(mfi);

		amrex::Array4<MATRIX4> const& cdata = fine_ddw_for_coarse.array(mfi);
		amrex::Array4<const MATRIX4> const& fdata       = fine_ddw.array(mfi);

		const Dim3 lo= amrex::lbound(cdomain), hi = amrex::ubound(cdomain);

		for (int n = 0; n < fine_ddw.nComp(); n++)
		{
			// I,J,K == coarse coordinates
			// i,j,k == fine coordinates
			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int I, int J, int K) {
					int i=I*2, j=J*2, k=K*2;
					
					if (nmask(I,J,K) == fine_fine_node || nmask(I,J,K) == coarse_fine_node)
						{
							if ((I == lo.x || I == hi.x) &&
							    (J == lo.y || J == hi.y) &&
							    (K == lo.z || K == hi.z)) // Corner
								cdata(I,J,K,n) = fdata(i,j,k,n);
							else if ((J == lo.y || J == hi.y) &&
								 (K == lo.z || K == hi.z)) // X edge
								cdata(I,J,K,n) = fdata(i-1,j,k,n)*0.25 + fdata(i,j,k,n)*0.5 + fdata(i+1,j,k,n)*0.25;
							else if ((K == lo.z || K == hi.z) &&
								 (I == lo.x || I == hi.x)) // Y edge
								cdata(I,J,K,n) = fdata(i,j-1,k,n)*0.25 + fdata(i,j,k,n)*0.5 + fdata(i,j+1,k,n)*0.25;
							else if ((I == lo.x || I == hi.x) &&
								 (J == lo.y || J == hi.y)) // Z edge
								cdata(I,J,K,n) = fdata(i,j,k-1,n)*0.25 + fdata(i,j,k,n)*0.5 + fdata(i,j,k+1,n)*0.25;
							else if (I == lo.x || I == hi.x) // X face
								cdata(I,J,K,n) =
									(  fdata(i,j-1,k-1,n)     + fdata(i,j,k-1,n)*2.0 + fdata(i,j+1,k-1,n)
									 + fdata(i,j-1,k  ,n)*2.0 + fdata(i,j,k  ,n)*4.0 + fdata(i,j+1,k  ,n)*2.0 
									 + fdata(i,j-1,k+1,n)     + fdata(i,j,k+1,n)*2.0 + fdata(i,j+1,k+1,n)     )/16.0;
							else if (J == lo.y || J == hi.y) // Y face
								cdata(I,J,K,n) =
									(  fdata(i-1,j,k-1,n)     + fdata(i-1,j,k,n)*2.0 + fdata(i-1,j,k+1,n)
									 + fdata(i  ,j,k-1,n)*2.0 + fdata(i  ,j,k,n)*4.0 + fdata(i  ,j,k+1,n)*2.0 
									 + fdata(i+1,j,k-1,n)     + fdata(i+1,j,k,n)*2.0 + fdata(i+1,j,k+1,n)     )/16.0;
							else if (K == lo.z || K == hi.z) // Z face
								cdata(I,J,K,n) =
									(  fdata(i-1,j-1,k,n)     + fdata(i,j-1,k,n)*2.0 + fdata(i+1,j-1,k,n)
									 + fdata(i-1,j  ,k,n)*2.0 + fdata(i,j  ,k,n)*4.0 + fdata(i+1,j  ,k,n)*2.0 
									 + fdata(i-1,j+1,k,n)     + fdata(i,j+1,k,n)*2.0 + fdata(i+1,j+1,k,n)     )/16.0;
							else // Interior
								cdata(I,J,K,n) =
									(fdata(i-1,j-1,k-1,n) + fdata(i-1,j-1,k+1,n) + fdata(i-1,j+1,k-1,n) + fdata(i-1,j+1,k+1,n) +
									 fdata(i+1,j-1,k-1,n) + fdata(i+1,j-1,k+1,n) + fdata(i+1,j+1,k-1,n) + fdata(i+1,j+1,k+1,n)) / 64.0
									+
									(fdata(i,j-1,k-1,n) + fdata(i,j-1,k+1,n) + fdata(i,j+1,k-1,n) + fdata(i,j+1,k+1,n) +
									 fdata(i-1,j,k-1,n) + fdata(i+1,j,k-1,n) + fdata(i-1,j,k+1,n) + fdata(i+1,j,k+1,n) +
									 fdata(i-1,j-1,k,n) + fdata(i-1,j+1,k,n) + fdata(i+1,j-1,k,n) + fdata(i+1,j+1,k,n)) / 32.0
									+
									(fdata(i-1,j,k,n) + fdata(i,j-1,k,n) + fdata(i,j,k-1,n) +
									 fdata(i+1,j,k,n) + fdata(i,j+1,k,n) + fdata(i,j,k+1,n)) / 16.0
									+
									fdata(i,j,k,n) / 8.0;
						}

				});
		}
	}

	// Copy the fine residual restricted onto the coarse grid
	// into the final residual.
	crse_ddw.ParallelCopy(fine_ddw_for_coarse,0,0,ncomp,0,0,cgeom.periodicity());
	const int mglev = 0;
	Util::RealFillBoundary(crse_ddw,m_geom[crse_amrlev][mglev]);
	return;
}



template<int SYM>
void
Elastic<SYM>::averageDownCoeffsSameAmrLevel (int amrlev)
{
	BL_PROFILE("Elastic::averageDownCoeffsSameAmrLevel()");

 	for (int mglev = 1; mglev < m_num_mg_levels[amrlev]; ++mglev)
 	{
		amrex::Box cdomain(m_geom[amrlev][mglev].Domain());
		cdomain.convert(amrex::IntVect::TheNodeVector());
		amrex::Box fdomain(m_geom[amrlev][mglev-1].Domain());
		fdomain.convert(amrex::IntVect::TheNodeVector());

		MultiTab& crse = *m_ddw_mf[amrlev][mglev];
		MultiTab& fine = *m_ddw_mf[amrlev][mglev-1];
		
		amrex::BoxArray crseba = crse.boxArray();
		amrex::BoxArray fineba = fine.boxArray();
		
		BoxArray newba = crseba;
		newba.refine(2);
		MultiTab fine_on_crseba;
		fine_on_crseba.define(newba,crse.DistributionMap(),1,4);
		fine_on_crseba.ParallelCopy(fine,0,0,1,2,4,m_geom[amrlev][mglev].periodicity());

		for (MFIter mfi(crse, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
		{
			Box bx = mfi.tilebox();
			bx = bx & cdomain;

			amrex::Array4<const Set::Matrix4<AMREX_SPACEDIM,SYM>> const& fdata = fine_on_crseba.array(mfi);
			amrex::Array4<Set::Matrix4<AMREX_SPACEDIM,SYM>> const& cdata       = crse.array(mfi);

			const Dim3 lo= amrex::lbound(cdomain), hi = amrex::ubound(cdomain);

			// I,J,K == coarse coordinates
			// i,j,k == fine coordinates
			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int I, int J, int K) {
					int i=2*I, j=2*J, k=2*K;

					if ((I == lo.x || I == hi.x) &&
					    (J == lo.y || J == hi.y) &&
					    (K == lo.z || K == hi.z)) // Corner
						cdata(I,J,K) = fdata(i,j,k);
					else if ((J == lo.y || J == hi.y) &&
						 (K == lo.z || K == hi.z)) // X edge
						cdata(I,J,K) = fdata(i-1,j,k)*0.25 + fdata(i,j,k)*0.5 + fdata(i+1,j,k)*0.25;
					else if ((K == lo.z || K == hi.z) &&
					 	 (I == lo.x || I == hi.x)) // Y edge
					 	cdata(I,J,K) = fdata(i,j-1,k)*0.25 + fdata(i,j,k)*0.5 + fdata(i,j+1,k)*0.25;
					else if ((I == lo.x || I == hi.x) &&
					 	 (J == lo.y || J == hi.y)) // Z edge
					 	cdata(I,J,K) = fdata(i,j,k-1)*0.25 + fdata(i,j,k)*0.5 + fdata(i,j,k+1)*0.25;
					else if (I == lo.x || I == hi.x) // X face
					 	cdata(I,J,K) =
					 		(  fdata(i,j-1,k-1)     + fdata(i,j,k-1)*2.0 + fdata(i,j+1,k-1)
					 		 + fdata(i,j-1,k  )*2.0 + fdata(i,j,k  )*4.0 + fdata(i,j+1,k  )*2.0 
					 		 + fdata(i,j-1,k+1)     + fdata(i,j,k+1)*2.0 + fdata(i,j+1,k+1)    )/16.0;
					else if (J == lo.y || J == hi.y) // Y face
					 	cdata(I,J,K) =
					 		(  fdata(i-1,j,k-1)     + fdata(i-1,j,k)*2.0 + fdata(i-1,j,k+1)
					 		 + fdata(i  ,j,k-1)*2.0 + fdata(i  ,j,k)*4.0 + fdata(i  ,j,k+1)*2.0 
					 		 + fdata(i+1,j,k-1)     + fdata(i+1,j,k)*2.0 + fdata(i+1,j,k+1))/16.0;
					 else if (K == lo.z || K == hi.z) // Z face
					 	cdata(I,J,K) =
					 		(  fdata(i-1,j-1,k)     + fdata(i,j-1,k)*2.0 + fdata(i+1,j-1,k)
					 		 + fdata(i-1,j  ,k)*2.0 + fdata(i,j  ,k)*4.0 + fdata(i+1,j  ,k)*2.0 
					 		 + fdata(i-1,j+1,k)     + fdata(i,j+1,k)*2.0 + fdata(i+1,j+1,k))/16.0;
					 else // Interior
						 cdata(I,J,K) =
							 (fdata(i-1,j-1,k-1) + fdata(i-1,j-1,k+1) + fdata(i-1,j+1,k-1) + fdata(i-1,j+1,k+1) +
							  fdata(i+1,j-1,k-1) + fdata(i+1,j-1,k+1) + fdata(i+1,j+1,k-1) + fdata(i+1,j+1,k+1)) / 64.0
							 +
							 (fdata(i,j-1,k-1) + fdata(i,j-1,k+1) + fdata(i,j+1,k-1) + fdata(i,j+1,k+1) +
							  fdata(i-1,j,k-1) + fdata(i+1,j,k-1) + fdata(i-1,j,k+1) + fdata(i+1,j,k+1) +
							  fdata(i-1,j-1,k) + fdata(i-1,j+1,k) + fdata(i+1,j-1,k) + fdata(i+1,j+1,k)) / 32.0
							 +
							 (fdata(i-1,j,k) + fdata(i,j-1,k) + fdata(i,j,k-1) +
							  fdata(i+1,j,k) + fdata(i,j+1,k) + fdata(i,j,k+1)) / 16.0
							 +
							 fdata(i,j,k) / 8.0;
				});
		}
		FillBoundaryCoeff(crse,m_geom[amrlev][mglev]);
	}
}

template<int SYM>
void
Elastic<SYM>::FillBoundaryCoeff (MultiTab& sigma, const Geometry& geom)
{
	BL_PROFILE("Elastic::FillBoundaryCoeff()");
	for (int i = 0; i < 2; i++)
	{
		MultiTab & mf = sigma;
		mf.FillBoundary(geom.periodicity());
		const int ncomp = mf.nComp();
		const int ng1 = 1;
		const int ng2 = 2;
		MultiTab tmpmf(mf.boxArray(), mf.DistributionMap(), ncomp, ng1);
	  	tmpmf.copy(mf,0,0,ncomp,ng2,ng1,geom.periodicity());
		mf.ParallelCopy   (tmpmf, 0, 0, ncomp, ng1, ng2, geom.periodicity());
	}
}

template class Elastic<Set::Sym::Major>;
template class Elastic<Set::Sym::Isotropic>;
template class Elastic<Set::Sym::MajorMinor>;
template class Elastic<Set::Sym::Diagonal>;

}

