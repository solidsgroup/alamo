#include <AMReX_MLPoisson.H>

#include <AMReX_FFT.H>

#include "PFC.H"
#include "BC/Constant.H"
#include "IO/ParmParse.H"
#include "IC/Random.H"
#include "Numeric/Stencil.H"
#include "Set/Set.H"

namespace Integrator
{
PFC::PFC() : Integrator()
{
}
PFC::~PFC()
{
    delete ic;
    delete bc;
}

void PFC::Parse(PFC &value, IO::ParmParse &pp)
{
    pp.query_required("q0",value.q0);
    pp.query_required("r",    value.r);
    // Regridding criterion
    pp.query_default("refinement_threshold",value.refinement_threshold, 1E100);

    // initial condition for :math:`\eta`
    pp.select_default<IC::Random>("eta.ic", value.ic, value.geom);
    // boundary condition for :math:`\eta`
    pp.select_default<BC::Constant>("eta.bc", value.bc, 1);

    value.RegisterNewFab(value.eta_mf, value.bc, 1, 1, "eta",true);
    //value.RegisterNewFab(value.etaold_mf, value.bc, 1, 1, "eta_old",false);
    value.RegisterNewFab(value.grad_chempot_mf, value.bc, 1, 1, "grad_chempot",true);
}


void
PFC::Advance (int lev, Set::Scalar /*time*/, Set::Scalar dt)
{

    //
    // FFT Boilerplate
    //
    amrex::FFT::R2C my_fft(this->geom[lev].Domain());
    auto const &[cba, cdm] = my_fft.getSpectralDataLayout();
    const Set::Scalar* DX = geom[lev].CellSize();
    Set::Scalar
        AMREX_D_DECL(
            pi_Lx = 2.0 * Set::Constant::Pi / geom[lev].Domain().length(0) / DX[0],
            pi_Ly = 2.0 * Set::Constant::Pi / geom[lev].Domain().length(1) / DX[1],
            pi_Lz = 2.0 * Set::Constant::Pi / geom[lev].Domain().length(2) / DX[2]);
    Set::Scalar scaling = 1.0 / geom[lev].Domain().d_numPts();
    

    //
    // Compute the gradient of the chemical potential in realspace
    //
    for ( amrex::MFIter mfi(*eta_mf[lev],true); mfi.isValid(); ++mfi )
    {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<const amrex::Real> const& eta = eta_mf[lev]->array(mfi);
        amrex::Array4<amrex::Real> const& grad_chempot    = grad_chempot_mf[lev]->array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            grad_chempot(i,j,k) = -eta(i,j,k)*eta(i,j,k)*eta(i,j,k);
        });
    }

    grad_chempot_mf[lev]->FillBoundary();

    //
    // FFT of eta
    // 
    amrex::FabArray<amrex::BaseFab<amrex::GpuComplex<Set::Scalar> > > eta_hat_mf(cba, cdm, 1, 0);
    my_fft.forward(*eta_mf[lev], eta_hat_mf);

    //
    // FFT of chemical potential gradient
    //
    amrex::FabArray<amrex::BaseFab<amrex::GpuComplex<Set::Scalar> > > chempot_hat_mf(cba, cdm, 1, 0);
    my_fft.forward(*grad_chempot_mf[lev], chempot_hat_mf);

    Set::Scalar small = 1.0E-12;

    //
    // Perform update in spectral coordinatees
    //
    //for (amrex::MFIter mfi(eta_hat_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    for (amrex::MFIter mfi(eta_hat_mf, false); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.tilebox();

        Util::Message(INFO, bx);
        
        amrex::Array4<amrex::GpuComplex<Set::Scalar>> const & eta_hat   =  eta_hat_mf.array(mfi);
        amrex::Array4<amrex::GpuComplex<Set::Scalar>> const & N_hat     =  chempot_hat_mf.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int m, int n, int p) {

            // Get spectral coordinates
            AMREX_D_TERM(
                Set::Scalar k1 = m * pi_Lx;,
                Set::Scalar k2 = (n < bx.length(1)/2 ? n * pi_Ly : (n - bx.length(1)) * pi_Ly);,
                Set::Scalar k3 = (p < bx.length(2)/2 ? p * pi_Lz : (p - bx.length(2)) * pi_Lz););

            Set::Scalar lap = AMREX_D_TERM(k1 * k1, + k2 * k2, + k3*k3);

            Set::Scalar bilap = lap * lap;

            Set::Scalar L = -bilap * (r + (q0*q0 - bilap) * (q0*q0 - bilap));

            Set::Scalar exp_LdT = exp(L * dt);

            eta_hat(m, n, p) = exp_LdT * eta_hat(m, n, p) + ( (exp_LdT - 1.0)/(L+small) ) * N_hat(m,n,p);
                
            eta_hat(m,n,p) *= scaling;
        });
    }

    //
    // Transform solution back to realspace
    //
    my_fft.backward(eta_hat_mf, *eta_mf[lev]);


}

void
PFC::Initialize (int lev)
{
    grad_chempot_mf[lev]->setVal(0.0);
    ic->Initialize(lev,eta_mf);
    //ic->Initialize(lev,etaold_mf);
}


void
PFC::TagCellsForRefinement (int lev, amrex::TagBoxArray& a_tags, Set::Scalar /*time*/, int /*ngrow*/)
{
        const Set::Scalar* DX = geom[lev].CellSize();
        Set::Scalar dr = sqrt(AMREX_D_TERM(DX[0] * DX[0], +DX[1] * DX[1], +DX[2] * DX[2]));

        for (amrex::MFIter mfi(*eta_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<char> const&     tags = a_tags.array(mfi);
            Set::Patch<const Set::Scalar>   eta = (*eta_mf[lev]).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Vector grad = Numeric::Gradient(eta, i, j, k, 0, DX);
                if (grad.lpNorm<2>() * dr > refinement_threshold)
                    tags(i, j, k) = amrex::TagBox::SET;
            });
        }
}


}
