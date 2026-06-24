#include "IC/Expression.H"
#include <AMReX_MLPoisson.H>
#include <algorithm>
#include <cmath>
#include <limits>

#ifdef ALAMO_FFT
#include <AMReX_FFT.H>
#endif

#include "CahnHilliard.H"
#include "BC/Constant.H"
#include "IO/ParmParse.H"
#include "IC/Random.H"
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
    pp.query_validate("mobility", value.mobility, {"constant", "singly_degenerate"});
    pp.query_default("mobility_floor", value.mobility_floor, 0.0);
    pp.query_default("spectral_stabilization", value.spectral_stabilization, std::numeric_limits<Set::Scalar>::quiet_NaN());
    // Regridding criterion
    pp.query_default("refinement_threshold",value.refinement_threshold, 1E100);

    // initial condition for :math:`\eta`
    pp.select_default<IC::Random,IC::Expression>("eta.ic", value.ic, value.geom);
    // boundary condition for :math:`\eta`
    pp.select_default<BC::Constant>("eta.bc", value.bc, 1);

    // Which method to use - realspace or spectral method.
    pp.query_validate("method",value.method,{"realspace","spectral"});

    pp.query_default("tstart",value.tstart,1E100);
    pp.query_default("tstart_timestep",value.tstart_timestep,-1.0);

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
        AdvanceSpectral(lev, time, dt);
    else
        Util::Abort(INFO,"Invalid method: ",method);
}


void
CahnHilliard::AdvanceReal (int lev, Set::Scalar time, Set::Scalar dt)
{
    std::swap(etaold_mf[lev], etanew_mf[lev]);
    etaold_mf[lev]->FillBoundary(geom[lev].periodicity());
    const Set::Scalar* DX = geom[lev].CellSize();

    if (time > tstart && tstart_timestep > 0)
        SetTimestep(tstart_timestep);

    for ( amrex::MFIter mfi(*etanew_mf[lev],true); mfi.isValid(); ++mfi )
    {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<const amrex::Real> const& eta = etaold_mf[lev]->array(mfi);
        amrex::Array4<amrex::Real> const& inter    = intermediate[lev]->array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Scalar lap_eta = Numeric::Laplacian(eta,i,j,k,0,DX);
            inter(i,j,k) = eta(i,j,k)*eta(i,j,k)*eta(i,j,k) - eta(i,j,k) - gamma*lap_eta;
        });
    }

    intermediate[lev]->FillBoundary(geom[lev].periodicity());

    for ( amrex::MFIter mfi(*etanew_mf[lev],true); mfi.isValid(); ++mfi )
    {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<const amrex::Real> const& eta = etaold_mf[lev]->array(mfi);
        amrex::Array4<const amrex::Real> const& inter = intermediate[lev]->array(mfi);
        amrex::Array4<amrex::Real> const& etanew = etanew_mf[lev]->array(mfi);
        const bool singly_degenerate = mobility == "singly_degenerate";
        const Set::Scalar L_local = L;
        const Set::Scalar M_floor = mobility_floor;

        amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k){
            Set::Scalar lap_inter = Numeric::Laplacian(inter,i,j,k,0,DX);
            Set::Scalar rhs = L_local * lap_inter;

            if (singly_degenerate)
            {
                Set::Scalar alpha = 0.5 * (eta(i,j,k) + 1.0);
                Set::Scalar M = L_local * alpha * alpha * (1.0 - alpha) * (1.0 - alpha) + M_floor;
                Set::Scalar dM_deta = 0.5 * L_local * 2.0 * alpha * (1.0 - alpha) * (1.0 - 2.0 * alpha);
                Set::Vector grad_eta = Numeric::Gradient(eta, i, j, k, 0, DX);
                Set::Vector grad_inter = Numeric::Gradient(inter, i, j, k, 0, DX);

                rhs = M * lap_inter + dM_deta * grad_eta.dot(grad_inter);
            }

            etanew(i,j,k) = eta(i,j,k) + dt*rhs;
            etanew(i,j,k) = std::max(-1.0, etanew(i,j,k));
            etanew(i,j,k) = std::min( 1.0, etanew(i,j,k));
        });
    }

    etanew_mf[lev]->FillBoundary(geom[lev].periodicity());
}

#ifdef ALAMO_FFT
void
CahnHilliard::AdvanceSpectral (int lev, Set::Scalar /*time*/, Set::Scalar dt)
{
    Util::Assert(INFO, TEST(lev == 0), "CahnHilliard spectral method currently supports one AMR level");

    amrex::FFT::R2C my_fft(this->geom[lev].Domain());
    auto const &[cba, cdm] = my_fft.getSpectralDataLayout();
    const Set::Scalar* DX = geom[lev].CellSize();
    amrex::Box const & domain = this->geom[lev].Domain();
    Set::Scalar
        AMREX_D_DECL(
            pi_Lx = 2.0 * Set::Constant::Pi / geom[lev].Domain().length(0) / DX[0],
            pi_Ly = 2.0 * Set::Constant::Pi / geom[lev].Domain().length(1) / DX[1],
            pi_Lz = 2.0 * Set::Constant::Pi / geom[lev].Domain().length(2) / DX[2]);
    Set::Scalar scaling = 1.0 / geom[lev].Domain().d_numPts();

    etanew_mf[lev]->FillBoundary(geom[lev].periodicity());

    if (mobility == "constant")
    {
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

        intermediate[lev]->FillBoundary(geom[lev].periodicity());

        amrex::FabArray<amrex::BaseFab<amrex::GpuComplex<Set::Scalar> > > eta_hat_mf(cba, cdm, 1, 0);
        amrex::FabArray<amrex::BaseFab<amrex::GpuComplex<Set::Scalar> > > chempot_hat_mf(cba, cdm, 1, 0);
        my_fft.forward(*etanew_mf[lev], eta_hat_mf);
        my_fft.forward(*intermediate[lev], chempot_hat_mf);

        for (amrex::MFIter mfi(eta_hat_mf, false); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.tilebox();
            amrex::Array4<amrex::GpuComplex<Set::Scalar>> const & eta_hat = eta_hat_mf.array(mfi);
            amrex::Array4<amrex::GpuComplex<Set::Scalar>> const & chempot_hat = chempot_hat_mf.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int m, int n, int p) {
                AMREX_D_TERM(
                    Set::Scalar k1 = m * pi_Lx;,
                    Set::Scalar k2 = (n < domain.length(1)/2 ? n * pi_Ly : (n - domain.length(1)) * pi_Ly);,
                    Set::Scalar k3 = (p < domain.length(2)/2 ? p * pi_Lz : (p - domain.length(2)) * pi_Lz););

                Set::Scalar omega2 = AMREX_D_TERM(k1 * k1, + k2 * k2, + k3*k3);
                Set::Scalar omega4 = omega2 * omega2;

                eta_hat(m, n, p) =
                    (eta_hat(m, n, p) - L * omega2 * chempot_hat(m, n, p) * dt) /
                    (1.0 + L * gamma * omega4 * dt);
                eta_hat(m,n,p) *= scaling;
            });
        }

        my_fft.backward(eta_hat_mf, *etanew_mf[lev]);
    }
    else if (mobility == "singly_degenerate")
    {
        for (amrex::MFIter mfi(*etanew_mf[lev], true); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<const amrex::Real> const& eta = etanew_mf[lev]->array(mfi);
            amrex::Array4<amrex::Real> const& inter = intermediate[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Scalar lap_eta = Numeric::Laplacian(eta, i, j, k, 0, DX);
                inter(i,j,k) = eta(i,j,k)*eta(i,j,k)*eta(i,j,k) - eta(i,j,k) - gamma*lap_eta;
            });
        }

        intermediate[lev]->FillBoundary(geom[lev].periodicity());

        amrex::MultiFab flux_mf(etanew_mf[lev]->boxArray(), etanew_mf[lev]->DistributionMap(), AMREX_SPACEDIM, 0);
        for (amrex::MFIter mfi(*etanew_mf[lev], true); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<const amrex::Real> const& eta = etanew_mf[lev]->array(mfi);
            amrex::Array4<const amrex::Real> const& inter = intermediate[lev]->array(mfi);
            amrex::Array4<amrex::Real> const& flux = flux_mf.array(mfi);
            const Set::Scalar L_local = L;
            const Set::Scalar M_floor = mobility_floor;

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Scalar alpha = 0.5 * (eta(i,j,k) + 1.0);
                Set::Scalar M = L_local * alpha * alpha * (1.0 - alpha) * (1.0 - alpha) + M_floor;
                Set::Vector grad_inter = Numeric::Gradient(inter, i, j, k, 0, DX);
                AMREX_D_TERM(
                    flux(i,j,k,0) = M * grad_inter(0);,
                    flux(i,j,k,1) = M * grad_inter(1);,
                    flux(i,j,k,2) = M * grad_inter(2););
            });
        }

        flux_mf.FillBoundary(geom[lev].periodicity());

        amrex::FabArray<amrex::BaseFab<amrex::GpuComplex<Set::Scalar> > > eta_hat_mf(cba, cdm, 1, 0);
        amrex::FabArray<amrex::BaseFab<amrex::GpuComplex<Set::Scalar> > > flux_hat_mf(cba, cdm, AMREX_SPACEDIM, 0);
        my_fft.forward(*etanew_mf[lev], eta_hat_mf);
        AMREX_D_TERM(
            my_fft.forward(flux_mf, flux_hat_mf, 0, 0);,
            my_fft.forward(flux_mf, flux_hat_mf, 1, 1);,
            my_fft.forward(flux_mf, flux_hat_mf, 2, 2););

        Set::Scalar S = spectral_stabilization;
        if (std::isnan(S))
            S = gamma * (L * 0.0625 + mobility_floor);
        amrex::GpuComplex<Set::Scalar> I(0.0, 1.0);

        for (amrex::MFIter mfi(eta_hat_mf, false); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.tilebox();
            amrex::Array4<amrex::GpuComplex<Set::Scalar>> const & eta_hat = eta_hat_mf.array(mfi);
            amrex::Array4<amrex::GpuComplex<Set::Scalar>> const & flux_hat = flux_hat_mf.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int m, int n, int p) {
                AMREX_D_TERM(
                    Set::Scalar k1 = m * pi_Lx;,
                    Set::Scalar k2 = (n < domain.length(1)/2 ? n * pi_Ly : (n - domain.length(1)) * pi_Ly);,
                    Set::Scalar k3 = (p < domain.length(2)/2 ? p * pi_Lz : (p - domain.length(2)) * pi_Lz););

                Set::Scalar omega2 = AMREX_D_TERM(k1 * k1, + k2 * k2, + k3*k3);
                Set::Scalar omega4 = omega2 * omega2;
                amrex::GpuComplex<Set::Scalar> div_flux_hat = AMREX_D_TERM(
                    I * k1 * flux_hat(m,n,p,0),
                    + I * k2 * flux_hat(m,n,p,1),
                    + I * k3 * flux_hat(m,n,p,2));

                eta_hat(m,n,p) = eta_hat(m,n,p) + dt * div_flux_hat + dt * S * omega4 * eta_hat(m,n,p);
                eta_hat(m,n,p) /= 1.0 + dt * S * omega4;
                eta_hat(m,n,p) *= scaling;
            });
        }

        my_fft.backward(eta_hat_mf, *etanew_mf[lev]);
    }
    else
        Util::Abort(INFO,"Invalid mobility: ",mobility);

    etanew_mf[lev]->FillBoundary(geom[lev].periodicity());
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
    if (method == "spectral") return;

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
