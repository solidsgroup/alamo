#include "Fracture.H"
#include "BC/Constant.H"
#include "Operator/Elastic.H"
#include "Solver/Linear.H"
#include "Model/Solid/LinearElastic/Isotropic.H"
#include "BC/Operator/Elastic.H"
//#include "IC/X.H"
#include "IC/Sphere.H"
//#include "IC/PS.H"
#include "Numeric/Stencil.H"

namespace Integrator
{
Fracture::Fracture() :
	Integrator()
{
	/*
	amrex::ParmParse pp("heat");
	pp.query("alpha", alpha);
	pp.query("refinement_threshold", refinement_threshold);
	pp.query("ic_type", ic_type);

	// // Determine initial condition
	if (ic_type == "cylinder")
		ic = new IC::Cylinder(geom);
	else if (ic_type == "sphere")
		ic = new IC::Sphere(geom);
	else if (ic_type == "constant")
		ic = new IC::Constant(geom);
//	else if (ic_type == "packed_spheres")
//		ic = new IC::PS(geom,20,0.0,1.0);
	else
		ic = new IC::Constant(geom);
    */
	{
		amrex::ParmParse pp("bc");
		amrex::Vector<std::string> bc_hi_str(AMREX_SPACEDIM);
		amrex::Vector<std::string> bc_lo_str(AMREX_SPACEDIM);
		pp.queryarr("lo",bc_lo_str,0,BL_SPACEDIM);
		pp.queryarr("hi",bc_hi_str,0,BL_SPACEDIM);
		amrex::Vector<amrex::Real> bc_lo_1, bc_hi_1;
		if (pp.countval("lo_1")) pp.getarr("lo_1",bc_lo_1);
		if (pp.countval("hi_1")) pp.getarr("hi_1",bc_hi_1);
		amrex::Vector<amrex::Real> bc_lo_2, bc_hi_2;
		if (pp.countval("lo_2")) pp.getarr("lo_2",bc_lo_2);
		if (pp.countval("hi_2")) pp.getarr("hi_2",bc_hi_2);
		amrex::Vector<amrex::Real> bc_lo_3, bc_hi_3;
		if (pp.countval("lo_3")) pp.getarr("lo_3",bc_lo_3);
		if (pp.countval("hi_3")) pp.getarr("hi_3",bc_hi_3);

		mybc = new BC::Constant(bc_hi_str, bc_lo_str,
					AMREX_D_DECL(bc_lo_1, bc_lo_2, bc_lo_3),
					AMREX_D_DECL(bc_hi_1, bc_hi_2, bc_hi_3));
	}
	

	RegisterNewFab(m_c,     mybc, number_of_components, number_of_ghost_cells, "c");
	RegisterNewFab(m_cold, mybc, number_of_components, number_of_ghost_cells, "cold");
	RegisterNodalFab(m_disp, AMREX_SPACEDIM, number_of_ghost_cells, "Disp");
	RegisterNodalFab(m_rhs,  AMREX_SPACEDIM, number_of_ghost_cells, "RHS");
}

Fracture::~Fracture()
{
}

void
Fracture::Initialize (int lev)
{
	//TODO: initiaize

	//ic->Initialize(lev,TempFab);
	//m_c[lev]->setVal(0.0);
	//m_cold[lev]->setVal(0.0);
	//m_disp[lev]->setVal(0.0);
	//m_rhs[lev]->setVal(0.0);
}

void 
Fracture::TimeStepBegin(amrex::Real /*time*/, int iter)
{
	Util::Message(INFO);
	if (iter % plot_int) return;
	Set::Scalar lame = 2.6, shear = 6.0;
	Model::Solid::LinearElastic::Isotropic model(lame,shear);
	Operator::Elastic<Model::Solid::LinearElastic::Isotropic> elastic;
	elastic.SetHomogeneous(false);
	elastic.define(geom,grids,dmap);
	elastic.SetModel(model);
	Util::Message(INFO);

	BC::Operator::Elastic<Model::Solid::LinearElastic::Isotropic> bc;
	elastic.SetBC(&bc);
	Util::Message(INFO);
	/*
	for (int lev = 0; lev < m_c.size(); lev++) 
	{
		m_rhs[lev]->setVal(0.0);
		m_disp[lev]->setVal(0.0);
		const amrex::Real* DX = geom[lev].CellSize();
		Util::Message(INFO);	
		for (MFIter mfi(*m_rhs[lev],amrex::TilingIfNotGPU());mfi.isValid();++mfi)
		{
			Util::Message(INFO);	
			amrex::Box bx = mfi.tilebox();
			bx.grow(2);
			Util::Message(INFO);	
			amrex::Array4<Set::Scalar> const & Rhs = m_rhs[lev]->array(mfi);
			amrex::Array4<const Set::Scalar> const & temp = m_cold[lev]->array(mfi);
			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {
				Util::Message(INFO,"A");	
				std::array<Numeric::StencilType,AMREX_SPACEDIM> stencil 
				= {{AMREX_D_DECL(Numeric::StencilType::CellToNode,Numeric::StencilType::CellToNode,Numeric::StencilType::CellToNode)}};
				Util::Message(INFO,"B");	
				Set::Vector GradT = Numeric::Gradient(temp,i,j,k,0,DX,stencil);

				Util::Message(INFO,"C");	
				Set::Matrix eps = Set::Matrix::Identity();
				Util::Message(INFO,"D");	
				Set::Matrix sig = model(eps);

				Util::Message(INFO,"E");	
				Set::Vector f = sig*GradT;
				AMREX_D_TERM(Rhs(i,j,k,0) = f(0);,
							 Rhs(i,j,k,1) = f(1);,	
							 Rhs(i,j,k,2) = f(2););	
			});
			Util::Message(INFO);	
		}
		Util::Message(INFO);	
	}
	Util::Message(INFO);	*/

	Solver::Linear solver(elastic);
	solver.setVerbose(3);
	Set::Scalar tol_rel = 1E-8, tol_abs = 1E-8;
	//solver.apply(GetVecOfPtrs(rhs),GetVecOfPtrs(eigendef));
	solver.solve(GetVecOfPtrs(m_disp),GetVecOfConstPtrs(m_rhs),tol_rel,tol_abs);
	//for (int lev = 0; lev < sigma.size(); lev++)
	//{
		//elastic.Stress(lev,*sigma[lev],*disp[lev]);
		//disp[lev]->setVal(0.0);
	//}
	Util::Message(INFO);

}


/// \fn    Integrator::Fracture::Advance
///
/// Integrate the heat diffusion equation
/// \f[\nabla^2T = \alpha \frac{\partial T}{\partial t}\f]
/// using an explicit forward Euler method.
/// \f$\alpha\f$ is stored in #alpha
void
Fracture::Advance (int lev, amrex::Real /*time*/, amrex::Real dt)
{
	/*
	static amrex::IntVect AMREX_D_DECL(dx(AMREX_D_DECL(1,0,0)),
									  dy(AMREX_D_DECL(0,1,0)),
									  dz(AMREX_D_DECL(0,0,1)));

	std::swap(*TempFab[lev], *TempOldFab[lev]);

	const amrex::Real* DX = geom[lev].CellSize();

	for ( amrex::MFIter mfi(*TempFab[lev],true); mfi.isValid(); ++mfi )
		{
			amrex::Box bx = mfi.tilebox();
			bx.grow(1);
			bx = bx & geom[lev].Domain();

			amrex::Array4<const Set::Scalar> const & TempOld = (*TempOldFab[lev]).array(mfi);
			amrex::Array4<      Set::Scalar> const & Temp    = (*TempFab[lev]).array(mfi);

			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {
				Temp(i,j,k) = TempOld(i,j,k) 
						      + dt*alpha*(AMREX_D_TERM((Numeric::Stencil<Set::Scalar,2,0,0>::D(TempOld,i,j,k,0,DX)),
							  			  			   + (Numeric::Stencil<Set::Scalar,0,2,0>::D(TempOld,i,j,k,0,DX)) ,
							  			  			   + (Numeric::Stencil<Set::Scalar,0,0,2>::D(TempOld,i,j,k,0,DX)) ));
			});
		}
		*/
}

void
Fracture::TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/)
{
	/*
	const amrex::Real* DX      = geom[lev].CellSize();

	static amrex::IntVect AMREX_D_DECL(dx(AMREX_D_DECL(1,0,0)),
												  dy(AMREX_D_DECL(0,1,0)),
												  dz(AMREX_D_DECL(0,0,1)));

	for (amrex::MFIter mfi(*TempFab[lev],true); mfi.isValid(); ++mfi)
		{
			const amrex::Box&  bx  = mfi.tilebox();
			amrex::TagBox&     tag  = tags[mfi];
 	    
			amrex::FArrayBox &Temp = (*TempFab[lev])[mfi];

			AMREX_D_TERM(for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++),
							 for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++),
							 for (int k = bx.loVect()[2]; k<=bx.hiVect()[2]; k++))
						{
							amrex::IntVect m(AMREX_D_DECL(i,j,k));

							AMREX_D_TERM(amrex::Real grad1 = (Temp(m+dx) - Temp(m-dx));,
											 amrex::Real grad2 = (Temp(m+dy) - Temp(m-dy));,
											 amrex::Real grad3 = (Temp(m+dz) - Temp(m-dz));)

							amrex::Real grad = sqrt(AMREX_D_TERM(grad1*grad1, + grad2*grad2, + grad3*grad3));

							amrex::Real dr = sqrt(AMREX_D_TERM(DX[0]*DX[0], + DX[1]*DX[1], + DX[2]*DX[2]));

							if (grad*dr > refinement_threshold)
								tag(amrex::IntVect(AMREX_D_DECL(i,j,k))) = amrex::TagBox::SET;
						}
		}
		*/
}

}
