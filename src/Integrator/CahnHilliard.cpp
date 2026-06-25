#include <AMReX_MLPoisson.H>

#include "CahnHilliard.H"
#include "BC/Constant.H"
#include "IO/ParmParse.H"
#include "IC/Random.H"
#include "Numeric/Stencil.H"
#include "Operator/Spectral/FFT.H"
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

    // Which method to use - realspace or spectral method.
    pp.query_validate("method",value.method,{"realspace","spectral"});

    value.RegisterNewFab(value.etanew_mf, value.bc, 1, 1, "eta",true);
    value.RegisterNewFab(value.intermediate, value.bc, 1, 1, "int",true);

    if (value.method == "realspace")
        value.RegisterNewFab(value.etaold_mf, value.bc, 1, 1, "eta_old",false);
}

void
CahnHilliard::Advance(int lev, Set::Scalar time, Set::Scalar dt)
{
    if (method == "realspace")
        AdvanceReal(lev, time, dt);
    else if (method == "spectral")
    {
        if (lev == finest_level) AdvanceSpectral(lev, time, dt);
    }
    else
        Util::Abort(INFO,"Invalid method: ",method);
}


void
CahnHilliard::AdvanceReal (int lev, Set::Scalar /*time*/, Set::Scalar dt)
{
    std::swap(etaold_mf[lev], etanew_mf[lev]);
    const Set::Scalar* DX = geom[lev].CellSize();
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

#ifdef ALAMO_FFT
void
CahnHilliard::AdvanceSpectral (int lev, Set::Scalar time, Set::Scalar dt)
{
    Operator::Spectral::FFT fft(geom, refRatio(), lev);

    //
    // Compute the gradient of the chemical potential in realspace
    //
    for ( amrex::MFIter mfi(*etanew_mf[lev],true); mfi.isValid(); ++mfi )
    {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<const amrex::Real> const& eta = etanew_mf[lev]->array(mfi);
        amrex::Array4<amrex::Real> const& inter    = intermediate[lev]->array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            inter(i,j,k) = eta(i,j,k)*eta(i,j,k)*eta(i,j,k) - eta(i,j,k);
        });
    }

    intermediate[lev]->FillBoundary();

    //
    // FFT of eta
    // 
    amrex::FabArray<amrex::BaseFab<Set::Complex> > eta_hat_mf = fft.MakeSpectralFab();
    fft.Forward(etanew_mf, lev, eta_hat_mf, 0, 0, time);

    //
    // FFT of chemical potential gradient
    //
    amrex::FabArray<amrex::BaseFab<Set::Complex> > chempot_hat_mf = fft.MakeSpectralFab();
    fft.Forward(intermediate, lev, chempot_hat_mf, 0, 0, time);

    //
    // Perform update in spectral coordinatees
    //
    //for (amrex::MFIter mfi(eta_hat_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    for (amrex::MFIter mfi(eta_hat_mf, false); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.tilebox();

        
        amrex::Array4<Set::Complex> const & eta_hat     =  eta_hat_mf.array(mfi);
        amrex::Array4<Set::Complex> const & chempot_hat =  chempot_hat_mf.array(mfi);

        fft.ParallelFor(bx, [=] AMREX_GPU_DEVICE(int m, int n, int p, Set::Scalar omega2) {
            Set::Scalar omega4 = omega2 * omega2;

            eta_hat(m, n, p) =
                (eta_hat(m, n, p) - L * omega2 * chempot_hat(m, n, p) * dt) /
                (1.0 + L * gamma * omega4 * dt);
        });
    }

    //
    // Transform solution back to realspace
    //
    fft.Backward(eta_hat_mf, etanew_mf, lev);
}
#else
void
CahnHilliard::AdvanceSpectral (int, Set::Scalar, Set::Scalar)
{
    Util::Abort(INFO,"Alamo must be compiled with fft");
}
#endif



void
CahnHilliard::Initialize (int lev)
{
    intermediate[lev]->setVal(0.0);
    ic->Initialize(lev,etanew_mf);
    if (method == "realspace")
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
