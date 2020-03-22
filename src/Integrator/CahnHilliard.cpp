#include <AMReX_MLPoisson.H>

#include "CahnHilliard.H"
#include "BC/Nothing.H"
#include "Numeric/Stencil.H"

namespace Integrator
{
CahnHilliard::CahnHilliard() : Integrator()
{
	bc = new BC::Nothing();
	ic = new IC::Random(geom,2.0);
	RegisterNewFab(etanewmf, bc, ncomp, nghost, "Eta",true);
	RegisterNewFab(etaoldmf, bc, ncomp, nghost, "EtaOld",false);
	RegisterNewFab(intermediate, bc, ncomp, nghost, "int",false);
	LPInfo info;
	op.define(geom,grids,dmap,*bc,info);
}


void
CahnHilliard::TimeStepBegin(amrex::Real /*time*/, int /*iter*/)
{
	// amrex::MLPoisson myop(geom,grids,dmap);
	// amrex::MLMG solver(myop);
	// solver.setMaxIter(elastic.max_iter);
	// solver.setMaxFmgIter(elastic.max_fmg_iter);
	// solver.setVerbose(elastic.verbose);
	// solver.setCGVerbose(elastic.cgverbose);

	// etanewmf[0]->setVal(0.0);
	// etanewmf[0]->setVal(0.0);

	// Set::Scalar tol_rel = 1E-8;
	// Set::Scalar tol_abs = 0.0;
	// solver.solve(GetVecOfPtrs(etanewmf),
	// 	     GetVecOfConstPtrs(etaoldmf),
	// 	     tol_rel,
	// 	     tol_abs);
	
}

void
CahnHilliard::Advance (int lev, Set::Scalar /*time*/, Set::Scalar dt)
{
	std::swap(etaoldmf[lev], etanewmf[lev]);
	const amrex::Real* DX = geom[lev].CellSize();
	for ( amrex::MFIter mfi(*etanewmf[lev],true); mfi.isValid(); ++mfi )
	{
		const amrex::Box& bx = mfi.tilebox();
		amrex::Array4<const amrex::Real> const& eta = etaoldmf[lev]->array(mfi);
		amrex::Array4<amrex::Real> const& inter    = intermediate[lev]->array(mfi);
		amrex::Array4<amrex::Real> const& etanew    = etanewmf[lev]->array(mfi);

		amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k){
				Set::Scalar lap =
				 	Numeric::Stencil<Set::Scalar,2,0,0>::D(eta,i,j,k,0,DX) +
				 	Numeric::Stencil<Set::Scalar,2,0,0>::D(eta,i,j,k,0,DX);

				inter(i,j,k) =
				 	eta(i,j,k)*eta(i,j,k)*eta(i,j,k)
				 	- eta(i,j,k)
				 	- gamma*lap;

				etanew(i,j,k) = eta(i,j,k) - dt*inter(i,j,k); // Allen Cahn
			});

		amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k){
				Set::Scalar lap = 
				 	Numeric::Stencil<Set::Scalar,2,0,0>::D(inter,i,j,k,0,DX) +
				 	Numeric::Stencil<Set::Scalar,2,0,0>::D(inter,i,j,k,0,DX);

				etanew(i,j,k) = eta(i,j,k) + dt*lap;
			});
	}
}

void
CahnHilliard::Initialize (int lev)
{
	etanewmf[lev]->setVal(-1.);
	etaoldmf[lev]->setVal(-1.);
	ic->Add(lev,etanewmf);
	ic->Add(lev,etaoldmf);
}


void
CahnHilliard::TagCellsForRefinement (int /*lev*/, amrex::TagBoxArray& /*tags*/, amrex::Real /*time*/, int /*ngrow*/)
{
}


}
