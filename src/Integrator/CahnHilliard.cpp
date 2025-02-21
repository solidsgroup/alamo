#include <AMReX_MLPoisson.H>

#include "CahnHilliard.H"
#include "BC/Nothing.H"
#include "BC/Constant.H"
#include "IO/ParmParse.H"
#include "Numeric/Stencil.H"
#include "Set/Set.H"

namespace Integrator
{
CahnHilliard::CahnHilliard() : Integrator()
{
}
CahnHilliard::~CahnHilliard()
{
    delete ic;
    delete bc;
}

void CahnHilliard::Parse(CahnHilliard &value, IO::ParmParse &pp)
{
    // Interface energy
    pp.query_default("gamma",value.gamma, 0.0005);
    // Mobility
    pp.query_default("L",    value.L,     1.0);
    // Regridding criterion
    pp.query_default("refinement_threshold",value.refinement_threshold, 1E100);

    // initial condition for :math:`\eta`
    pp.select_default<IC::Random>("eta.ic", value.ic, value.geom);
    // boundary condition for :math:`\eta`
    pp.select_default<BC::Constant>("eta.bc", value.bc, 1);

    value.RegisterNewFab(value.etanew_mf, value.bc, 1, 1, "eta",true);
    value.RegisterNewFab(value.etaold_mf, value.bc, 1, 1, "eta_old",false);
    value.RegisterNewFab(value.intermediate, value.bc, 1, 1, "int",false);
}


void
CahnHilliard::Advance (int lev, Set::Scalar /*time*/, Set::Scalar dt)
{
    std::swap(etaold_mf[lev], etanew_mf[lev]);
    const amrex::Real* DX = geom[lev].CellSize();
    for ( amrex::MFIter mfi(*etanew_mf[lev],true); mfi.isValid(); ++mfi )
    {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<const amrex::Real> const& eta = etaold_mf[lev]->array(mfi);
        amrex::Array4<amrex::Real> const& inter    = intermediate[lev]->array(mfi);
        amrex::Array4<amrex::Real> const& etanew    = etanew_mf[lev]->array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Scalar lap_eta = Numeric::Laplacian(eta,i,j,k,0,DX);
            

            inter(i,j,k) =
                eta(i,j,k)*eta(i,j,k)*eta(i,j,k)
                - eta(i,j,k)
                - gamma*lap_eta;


            etanew(i,j,k) = eta(i,j,k) - dt*inter(i,j,k); // Allen Cahn
        });

        amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k){
            Set::Scalar lap_inter = Numeric::Laplacian(inter,i,j,k,0,DX);

            etanew(i,j,k) = eta(i,j,k) + dt*lap_inter;
        });
    }
}

void
CahnHilliard::Initialize (int lev)
{
    intermediate[lev]->setVal(0.0);
    ic->Initialize(lev,etanew_mf);
    ic->Initialize(lev,etaold_mf);
}


void
CahnHilliard::TagCellsForRefinement (int lev, amrex::TagBoxArray& a_tags, Set::Scalar /*time*/, int /*ngrow*/)
{
        const Set::Scalar* DX = geom[lev].CellSize();
        Set::Scalar dr = sqrt(AMREX_D_TERM(DX[0] * DX[0], +DX[1] * DX[1], +DX[2] * DX[2]));

        for (amrex::MFIter mfi(*etanew_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<char> const&     tags = a_tags.array(mfi);
            Set::Patch<const Set::Scalar>   eta = (*etanew_mf[lev]).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Vector grad = Numeric::Gradient(eta, i, j, k, 0, DX);
                if (grad.lpNorm<2>() * dr > refinement_threshold)
                    tags(i, j, k) = amrex::TagBox::SET;
            });
        }
}


}
