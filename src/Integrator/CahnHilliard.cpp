#include "CahnHilliard.H"
#include "BC/Nothing.H"
#include "Numeric/Stencil.H"
//#include "Operator/Implicit/Implicit.H"
namespace Integrator
{
CahnHilliard::CahnHilliard() : Integrator()
{
	bc = new BC::Nothing();
	ic = new IC::Random(geom,2.0);
	RegisterNewFab(etanewmf, bc, ncomp, nghost, "Eta");
	RegisterNewFab(etaoldmf, bc, ncomp, nghost, "EtaOld");
	RegisterNewFab(intermediate, bc, ncomp, nghost, "int");
}


void
CahnHilliard::Advance (int lev, Set::Scalar time, Set::Scalar dt)
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
					Numeric::Stencil::D<2,0,0>(eta,i,j,k,0,DX) +
					Numeric::Stencil::D<0,2,0>(eta,i,j,k,0,DX);

				inter(i,j,k) =
					eta(i,j,k)*eta(i,j,k)*eta(i,j,k)
					- eta(i,j,k)
					- gamma*lap;


				//etanew(i,j,k) = eta(i,j,k) - dt*inter(i,j,k); // Allen Cahn
			});

		amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k){

				Set::Scalar lap = 
					Numeric::Stencil::D<2,0,0>(inter,i,j,k,0,DX) +
					Numeric::Stencil::D<0,2,0>(inter,i,j,k,0,DX);

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
CahnHilliard::TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/)
{
}


}
