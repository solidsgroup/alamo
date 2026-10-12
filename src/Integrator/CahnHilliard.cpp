#include "IC/Expression.H"
#include <AMReX_MLPoisson.H>
#include <algorithm>
#include <cmath>
#include <limits>

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
    pp.query_validate("mobility", value.mobility, {"constant", "singly_degenerate"});
    pp.query_default("mobility_floor", value.mobility_floor, 0.0);
    pp.query_default("spectral_stabilization", value.spectral_stabilization, std::numeric_limits<Set::Scalar>::quiet_NaN());
    pp.query_default("diagnostics", value.output_diagnostics, false);
    // Regridding criterion
    pp.query_default("refinement_threshold",value.refinement_threshold, 1E100);

    pp.query_default("input_name", value.input_name, value.input_name);

    // initial condition for :math:`\eta`
    pp.select_default<IC::Random,IC::Expression>(value.input_name + ".ic", value.ic, pp.forward_args(value.geom));
    // boundary condition for :math:`\eta`
    pp.select_default<BC::Constant>(value.input_name + ".bc", value.bc, pp.forward_args(1));

    // Which method to use - realspace or spectral method.
    pp.query_validate("method",value.method,{"realspace","spectral"});

    pp.query_default("tstart",value.tstart,1E100);
    pp.query_default("tstart_timestep",value.tstart_timestep,-1.0);
    pp.query_default("field_name", value.field_name, value.field_name);
    pp.query_default("old_field_name", value.old_field_name, value.old_field_name);
    pp.query_default("intermediate_name", value.intermediate_name, value.intermediate_name);

    value.RegisterNewFab(value.etanew_mf, value.bc, 1, 1, value.field_name,true);
    value.RegisterNewFab(value.intermediate, value.bc, 1, 1, value.intermediate_name,true);

    if (value.output_diagnostics)
    {
        std::string prefix;
        const auto separator = value.field_name.find_last_of('.');
        if (separator != std::string::npos)
            prefix = value.field_name.substr(0, separator + 1);

        std::vector<std::string> flux_suffix;
        flux_suffix.push_back("_x");
#if AMREX_SPACEDIM >= 2
        flux_suffix.push_back("_y");
#endif
#if AMREX_SPACEDIM >= 3
        flux_suffix.push_back("_z");
#endif
        value.RegisterNewFab(value.diagnostic_mobility, nullptr, 1, 0,
                             prefix + "mobility", true, false);
        value.RegisterNewFab(value.diagnostic_mobility_mask, nullptr, 1, 0,
                             prefix + "mobility_mask", true, false);
        value.RegisterNewFab(value.diagnostic_flux, nullptr, AMREX_SPACEDIM, 0,
                             prefix + "flux", true, false, flux_suffix);
        value.RegisterNewFab(value.diagnostic_rhs, nullptr, 1, 0,
                             prefix + "rhs", true, false);
        value.RegisterNewFab(value.diagnostic_raw, nullptr, 1, 0,
                             prefix + "alpha_raw", true, false);
        value.RegisterNewFab(value.diagnostic_delta, nullptr, 1, 0,
                             prefix + "alpha_delta", true, false);
    }

    if (value.method == "realspace")
        value.RegisterNewFab(value.etaold_mf, value.bc, 1, 1, value.old_field_name,false);
}


void
CahnHilliard::FillMobilityMask(int /*lev*/, amrex::MultiFab& mask)
{
    mask.setVal(1.0);
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

    std::unique_ptr<amrex::MultiFab> mobility_mask;
    if (mobility == "singly_degenerate")
    {
        mobility_mask = std::make_unique<amrex::MultiFab>(
            etanew_mf[lev]->boxArray(), etanew_mf[lev]->DistributionMap(), 1, 0);
        FillMobilityMask(lev, *mobility_mask);
    }

    for ( amrex::MFIter mfi(*etanew_mf[lev],true); mfi.isValid(); ++mfi )
    {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<const amrex::Real> const& eta = etaold_mf[lev]->array(mfi);
        amrex::Array4<const amrex::Real> const& inter = intermediate[lev]->array(mfi);
        amrex::Array4<amrex::Real> const& etanew = etanew_mf[lev]->array(mfi);
        const bool singly_degenerate = mobility == "singly_degenerate";
        amrex::Array4<const amrex::Real> const mask = singly_degenerate ? mobility_mask->const_array(mfi) : amrex::Array4<const amrex::Real>();
        const Set::Scalar L_local = L;
        const Set::Scalar M_floor = mobility_floor;

        amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k){
            Set::Scalar lap_inter = Numeric::Laplacian(inter,i,j,k,0,DX);
            Set::Scalar rhs = L_local * lap_inter;

            if (singly_degenerate)
            {
                Set::Scalar alpha = std::max(0.0, std::min(1.0, 0.5 * (eta(i,j,k) + 1.0)));
                Set::Scalar active = mask(i,j,k);
                Set::Scalar M = active * (L_local * alpha * alpha * (1.0 - alpha) * (1.0 - alpha) + M_floor);
                Set::Scalar dM_deta = active * 0.5 * L_local * 2.0 * alpha * (1.0 - alpha) * (1.0 - 2.0 * alpha);
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
CahnHilliard::AdvanceSpectral (int lev, Set::Scalar time, Set::Scalar dt)
{
    Operator::Spectral::FFT fft(geom, refRatio(), lev);
    etanew_mf[lev]->FillBoundary(geom[lev].periodicity());

    if (mobility == "constant")
    {
        for (int ilev = 0; ilev <= lev; ++ilev)
        {
            for ( amrex::MFIter mfi(*etanew_mf[ilev],true); mfi.isValid(); ++mfi )
            {
                const amrex::Box& bx = mfi.tilebox();
                amrex::Array4<const amrex::Real> const& eta = etanew_mf[ilev]->array(mfi);
                amrex::Array4<amrex::Real> const& inter = intermediate[ilev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    inter(i,j,k) = eta(i,j,k)*eta(i,j,k)*eta(i,j,k) - eta(i,j,k);
                });
            }
            intermediate[ilev]->FillBoundary(geom[ilev].periodicity());
        }

        amrex::FabArray<amrex::BaseFab<Set::Complex> > eta_hat_mf = fft.MakeSpectralFab();
        amrex::FabArray<amrex::BaseFab<Set::Complex> > chempot_hat_mf = fft.MakeSpectralFab();
        fft.Forward(etanew_mf, lev, eta_hat_mf, 0, 0, time);
        fft.Forward(intermediate, lev, chempot_hat_mf, 0, 0, time);

        for (amrex::MFIter mfi(eta_hat_mf, false); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.tilebox();
            amrex::Array4<Set::Complex> const & eta_hat = eta_hat_mf.array(mfi);
            amrex::Array4<Set::Complex> const & chempot_hat = chempot_hat_mf.array(mfi);

            fft.ParallelFor(bx, [=] AMREX_GPU_DEVICE(int m, int n, int p, Set::Scalar omega2) {
                Set::Scalar omega4 = omega2 * omega2;
                eta_hat(m, n, p) =
                    (eta_hat(m, n, p) - L * omega2 * chempot_hat(m, n, p) * dt) /
                    (1.0 + L * gamma * omega4 * dt);
            });
        }

        fft.Backward(eta_hat_mf, etanew_mf, lev);
    }
    else if (mobility == "singly_degenerate")
    {
        amrex::Vector<std::unique_ptr<amrex::MultiFab> > old_alpha_mf(lev + 1);
        amrex::Vector<std::unique_ptr<amrex::MultiFab> > rhs_mf(lev + 1);
        amrex::Vector<std::unique_ptr<amrex::MultiFab> > mobility_mf(lev + 1);
        amrex::Vector<std::unique_ptr<amrex::MultiFab> > mobility_mask(lev + 1);
        amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab> > > face_flux_mf(lev + 1);

        for (int ilev = 0; ilev <= lev; ++ilev)
        {
            const Set::Scalar* DX = geom[ilev].CellSize();
            etanew_mf[ilev]->FillBoundary(geom[ilev].periodicity());
            old_alpha_mf[ilev] = std::make_unique<amrex::MultiFab>(
                etanew_mf[ilev]->boxArray(), etanew_mf[ilev]->DistributionMap(), 1, 0);
            amrex::MultiFab::Copy(*old_alpha_mf[ilev], *etanew_mf[ilev], 0, 0, 1, 0);

            for (amrex::MFIter mfi(*etanew_mf[ilev], true); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.tilebox();
                amrex::Array4<const amrex::Real> const& eta = etanew_mf[ilev]->array(mfi);
                amrex::Array4<amrex::Real> const& inter = intermediate[ilev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Scalar lap_eta = Numeric::Laplacian(eta, i, j, k, 0, DX);
                    inter(i,j,k) = eta(i,j,k)*eta(i,j,k)*eta(i,j,k) - eta(i,j,k) - gamma*lap_eta;
                });
            }

            intermediate[ilev]->FillBoundary(geom[ilev].periodicity());
            mobility_mask[ilev] = std::make_unique<amrex::MultiFab>(
                etanew_mf[ilev]->boxArray(), etanew_mf[ilev]->DistributionMap(), 1, 1);
            FillMobilityMask(ilev, *mobility_mask[ilev]);
            mobility_mask[ilev]->FillBoundary(geom[ilev].periodicity());
            mobility_mf[ilev] = std::make_unique<amrex::MultiFab>(
                etanew_mf[ilev]->boxArray(), etanew_mf[ilev]->DistributionMap(), 1, 1);

            for (amrex::MFIter mfi(*etanew_mf[ilev], true); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.tilebox();
                const amrex::Array4<const amrex::Real> eta = etanew_mf[ilev]->const_array(mfi);
                const amrex::Array4<const amrex::Real> mask = mobility_mask[ilev]->const_array(mfi);
                const amrex::Array4<amrex::Real> mobility_arr = mobility_mf[ilev]->array(mfi);
                const Set::Scalar L_local = L;
                const Set::Scalar M_floor = mobility_floor;

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Scalar alpha = std::max(0.0, std::min(1.0, 0.5 * (eta(i,j,k) + 1.0)));
                    mobility_arr(i,j,k) = mask(i,j,k) *
                        (L_local * alpha * alpha * (1.0 - alpha) * (1.0 - alpha) + M_floor);
                });
            }
            mobility_mf[ilev]->FillBoundary(geom[ilev].periodicity());

            if (output_diagnostics)
            {
                amrex::MultiFab::Copy(*diagnostic_mobility[ilev], *mobility_mf[ilev], 0, 0, 1, 0);
                amrex::MultiFab::Copy(*diagnostic_mobility_mask[ilev], *mobility_mask[ilev], 0, 0, 1, 0);
            }

            face_flux_mf[ilev].resize(AMREX_SPACEDIM);
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
            {
                amrex::BoxArray face_boxes = etanew_mf[ilev]->boxArray();
                face_boxes.surroundingNodes(dir);
                face_flux_mf[ilev][dir] = std::make_unique<amrex::MultiFab>(
                    face_boxes, etanew_mf[ilev]->DistributionMap(), 1, 0);

                const Set::Scalar inv_dx = 1.0 / DX[dir];
                for (amrex::MFIter mfi(*face_flux_mf[ilev][dir], true); mfi.isValid(); ++mfi)
                {
                    const amrex::Box& bx = mfi.tilebox();
                    const amrex::Array4<const amrex::Real> mu = intermediate[ilev]->const_array(mfi);
                    const amrex::Array4<const amrex::Real> mobility = mobility_mf[ilev]->const_array(mfi);
                    const amrex::Array4<amrex::Real> flux = face_flux_mf[ilev][dir]->array(mfi);
                    const int direction = dir;
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                    {
                        int ilo = i, jlo = j, klo = k;
                        if (direction == 0) --ilo;
                        else if (direction == 1) --jlo;
                        else --klo;
                        const Set::Scalar M_face = 0.5 *
                            (mobility(ilo,jlo,klo) + mobility(i,j,k));
                        flux(i,j,k) = M_face * (mu(i,j,k) - mu(ilo,jlo,klo)) * inv_dx;
                    });
                }
                face_flux_mf[ilev][dir]->FillBoundary(geom[ilev].periodicity());
            }

            rhs_mf[ilev] = std::make_unique<amrex::MultiFab>(
                etanew_mf[ilev]->boxArray(), etanew_mf[ilev]->DistributionMap(), 1, 0);
            for (amrex::MFIter mfi(*rhs_mf[ilev], true); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.tilebox();
                const amrex::Array4<const amrex::Real> fx = face_flux_mf[ilev][0]->const_array(mfi);
#if AMREX_SPACEDIM >= 2
                const amrex::Array4<const amrex::Real> fy = face_flux_mf[ilev][1]->const_array(mfi);
#endif
#if AMREX_SPACEDIM >= 3
                const amrex::Array4<const amrex::Real> fz = face_flux_mf[ilev][2]->const_array(mfi);
#endif
                const amrex::Array4<amrex::Real> rhs = rhs_mf[ilev]->array(mfi);
                amrex::Array4<amrex::Real> out_flux;
                if (output_diagnostics)
                    out_flux = diagnostic_flux[ilev]->array(mfi);
                const bool save_diagnostics = output_diagnostics;
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Scalar divergence = 0.0;
                    AMREX_D_TERM(
                        divergence += (fx(i+1,j,k) - fx(i,j,k)) / DX[0];,
                        divergence += (fy(i,j+1,k) - fy(i,j,k)) / DX[1];,
                        divergence += (fz(i,j,k+1) - fz(i,j,k)) / DX[2];);
                    rhs(i,j,k) = divergence;
                    if (save_diagnostics)
                    {
                        AMREX_D_TERM(
                            out_flux(i,j,k,0) = 0.5 * (fx(i,j,k) + fx(i+1,j,k));,
                            out_flux(i,j,k,1) = 0.5 * (fy(i,j,k) + fy(i,j+1,k));,
                            out_flux(i,j,k,2) = 0.5 * (fz(i,j,k) + fz(i,j,k+1)););
                    }
                });
            }
            if (output_diagnostics)
                amrex::MultiFab::Copy(*diagnostic_rhs[ilev], *rhs_mf[ilev], 0, 0, 1, 0);
        }

        amrex::FabArray<amrex::BaseFab<Set::Complex> > eta_hat_mf = fft.MakeSpectralFab();
        amrex::FabArray<amrex::BaseFab<Set::Complex> > rhs_hat_mf = fft.MakeSpectralFab();
        fft.Forward(etanew_mf, lev, eta_hat_mf, 0, 0, time);
        fft.Forward(rhs_mf, lev, rhs_hat_mf, 0, 0, time);

        Set::Scalar S = spectral_stabilization;
        if (std::isnan(S))
            S = gamma * (L * 0.0625 + mobility_floor);
        const Set::Scalar* DX = geom[lev].CellSize();
        const Set::Scalar dx0 = DX[0];
#if AMREX_SPACEDIM >= 2
        const Set::Scalar dx1 = DX[1];
#endif
#if AMREX_SPACEDIM >= 3
        const Set::Scalar dx2 = DX[2];
#endif

        for (amrex::MFIter mfi(eta_hat_mf, false); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.tilebox();
            amrex::Array4<Set::Complex> const & eta_hat = eta_hat_mf.array(mfi);
            amrex::Array4<Set::Complex> const & rhs_hat = rhs_hat_mf.array(mfi);

            fft.ParallelForModes(bx, [=] AMREX_GPU_DEVICE(
                int m, int n, int p, Set::Scalar k1, Set::Scalar k2, Set::Scalar k3, Set::Scalar omega2) {
                amrex::ignore_unused(omega2);
#if AMREX_SPACEDIM < 3
                amrex::ignore_unused(k3);
#endif
#if AMREX_SPACEDIM < 2
                amrex::ignore_unused(k2);
#endif
                Set::Scalar lambda = 4.0 * std::sin(0.5 * k1 * dx0) * std::sin(0.5 * k1 * dx0) / (dx0 * dx0);
#if AMREX_SPACEDIM >= 2
                lambda += 4.0 * std::sin(0.5 * k2 * dx1) * std::sin(0.5 * k2 * dx1) / (dx1 * dx1);
#endif
#if AMREX_SPACEDIM >= 3
                lambda += 4.0 * std::sin(0.5 * k3 * dx2) * std::sin(0.5 * k3 * dx2) / (dx2 * dx2);
#endif
                const Set::Scalar lambda4 = lambda * lambda;
                eta_hat(m,n,p) = (eta_hat(m,n,p) + dt * rhs_hat(m,n,p) +
                    dt * S * lambda4 * eta_hat(m,n,p)) / (1.0 + dt * S * lambda4);
            });
        }

        fft.Backward(eta_hat_mf, etanew_mf, lev);

        for (int ilev = 0; ilev <= lev; ++ilev)
        {
            if (output_diagnostics)
                amrex::MultiFab::Copy(*diagnostic_raw[ilev], *etanew_mf[ilev], 0, 0, 1, 0);

            for (amrex::MFIter mfi(*etanew_mf[ilev], true); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.tilebox();
                amrex::Array4<amrex::Real> const& eta = etanew_mf[ilev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    eta(i,j,k) = std::max(-1.0, std::min(1.0, eta(i,j,k)));
                });
            }
            if (output_diagnostics)
            {
                amrex::MultiFab::Copy(*diagnostic_delta[ilev], *etanew_mf[ilev], 0, 0, 1, 0);
                amrex::MultiFab::Subtract(*diagnostic_delta[ilev], *old_alpha_mf[ilev], 0, 0, 1, 0);
            }
        }
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
    InitializeDiagnostics(lev);
}


void
CahnHilliard::InitializeDiagnostics(int lev)
{
    if (!output_diagnostics) return;

    amrex::MultiFab mask(etanew_mf[lev]->boxArray(), etanew_mf[lev]->DistributionMap(), 1, 0);
    FillMobilityMask(lev, mask);
    const Set::Scalar L_local = L;
    const Set::Scalar M_floor = mobility_floor;
    const bool singly_degenerate = mobility == "singly_degenerate";

    for (amrex::MFIter mfi(*etanew_mf[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        const amrex::Array4<const amrex::Real> eta = etanew_mf[lev]->const_array(mfi);
        const amrex::Array4<const amrex::Real> mask_arr = mask.const_array(mfi);
        const amrex::Array4<amrex::Real> out_mobility = diagnostic_mobility[lev]->array(mfi);
        const amrex::Array4<amrex::Real> out_mask = diagnostic_mobility_mask[lev]->array(mfi);
        const amrex::Array4<amrex::Real> out_raw = diagnostic_raw[lev]->array(mfi);
        const amrex::Array4<amrex::Real> out_delta = diagnostic_delta[lev]->array(mfi);
        const amrex::Array4<amrex::Real> out_rhs = diagnostic_rhs[lev]->array(mfi);
        const amrex::Array4<amrex::Real> out_flux = diagnostic_flux[lev]->array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Scalar alpha = singly_degenerate ?
                std::max(0.0, std::min(1.0, 0.5 * (eta(i,j,k) + 1.0))) : 0.5;
            out_mask(i,j,k) = mask_arr(i,j,k);
            out_mobility(i,j,k) = singly_degenerate ?
                mask_arr(i,j,k) * (L_local * alpha * alpha * (1.0 - alpha) * (1.0 - alpha) + M_floor) : L_local;
            out_raw(i,j,k) = eta(i,j,k);
            out_delta(i,j,k) = 0.0;
            out_rhs(i,j,k) = 0.0;
            AMREX_D_TERM(out_flux(i,j,k,0) = 0.0;, out_flux(i,j,k,1) = 0.0;, out_flux(i,j,k,2) = 0.0;);
        });
    }
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
