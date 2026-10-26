#include <AMReX_MLPoisson.H>
#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#include "IC/Expression.H"

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
    // Mobility model (constant preserves the original update).
    pp.query_validate("mobility", value.mobility, {"constant", "singly_degenerate"});
    // Additive mobility in the pure phases.
    pp.query_default("mobility_floor", value.mobility_floor, 0.0);
    // Coefficient of the spectral fourth-order stabilizer.
    pp.query_default("spectral_stabilization", value.spectral_stabilization,
                    value.gamma * (value.L / 16.0 + value.mobility_floor));
    // Regridding criterion
    pp.query_default("refinement_threshold",value.refinement_threshold, 1E100);

    pp.query_default("input_name", value.input_name, value.input_name);
    pp.query_default("field_name", value.field_name, value.field_name);
    pp.query_default("old_field_name", value.old_field_name, value.old_field_name);
    pp.query_default("intermediate_name", value.intermediate_name, value.intermediate_name);
    pp.query_default("clipping", value.clipping, false);
    pp.query_default("diagnostics", value.output_diagnostics, false);
    pp.query_default("tstart", value.tstart, 1E100);
    pp.query_default("tstart_timestep", value.tstart_timestep, -1.0);

    // initial condition for :math:`\eta`
    pp.select_default<IC::Random,IC::Expression>(value.input_name + ".ic", value.ic, pp.forward_args(value.geom));
    // boundary condition for :math:`\eta`
    pp.select_default<BC::Constant>(value.input_name + ".bc", value.bc, pp.forward_args(1));

    // Which method to use - realspace or spectral method.
    pp.query_validate("method",value.method,{"realspace","spectral"});

    if (value.mobility == "singly_degenerate" && !pp.InTraversalMode())
    {
        if (!std::isfinite(value.L) || value.L < 0.0 ||
            !std::isfinite(value.gamma) || value.gamma < 0.0 ||
            !std::isfinite(value.mobility_floor) || value.mobility_floor < 0.0 ||
            !std::isfinite(value.spectral_stabilization) || value.spectral_stabilization < 0.0)
            Util::ParmParseException(INFO, "mobility", "SDCH requires finite nonnegative L, gamma, mobility_floor and spectral_stabilization");
        if (!value.geom[0].isAllPeriodic())
            Util::ParmParseException(INFO, "mobility", "SDCH requires periodic boundaries");
    }
    if (value.output_diagnostics && value.mobility != "singly_degenerate" && !pp.InTraversalMode())
        Util::ParmParseException(INFO, "diagnostics", "CahnHilliard diagnostics require mobility=singly_degenerate");

    value.RegisterNewFab(value.etanew_mf, value.bc, 1, 1, value.field_name,true);
    value.RegisterNewFab(value.intermediate, value.bc, 1, 1, value.intermediate_name,true);

    if (value.output_diagnostics)
    {
        std::string prefix;
        const auto separator = value.field_name.find_last_of('.');
        if (separator != std::string::npos) prefix = value.field_name.substr(0, separator + 1);
        std::vector<std::string> flux_suffix{"_x"};
#if AMREX_SPACEDIM >= 2
        flux_suffix.push_back("_y");
#endif
#if AMREX_SPACEDIM >= 3
        flux_suffix.push_back("_z");
#endif
        value.RegisterNewFab(value.diagnostic_mobility, nullptr, 1, 0, prefix + "mobility", true, false);
        value.RegisterNewFab(value.diagnostic_mobility_mask, nullptr, 1, 0, prefix + "mobility_mask", true, false);
        value.RegisterNewFab(value.diagnostic_flux, nullptr, AMREX_SPACEDIM, 0, prefix + "flux", true, false, flux_suffix);
        value.RegisterNewFab(value.diagnostic_rhs, nullptr, 1, 0, prefix + "rhs", true, false);
        value.RegisterNewFab(value.diagnostic_raw, nullptr, 1, 0, prefix + "alpha_raw", true, false);
        value.RegisterNewFab(value.diagnostic_delta, nullptr, 1, 0, prefix + "alpha_delta", true, false);
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
    if (method == "realspace" && time > tstart && tstart_timestep > 0.0)
        SetTimestep(tstart_timestep);
    if (mobility == "singly_degenerate")
    {
        if (method == "realspace")
            AdvanceDegenerateReal(lev, time, dt);
        else if (lev == finest_level)
            AdvanceDegenerateSpectral(lev, time, dt);
        return;
    }
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

namespace
{
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar DegenerateMobility(Set::Scalar eta, Set::Scalar L, Set::Scalar floor)
{
    Set::Scalar alpha = std::max(0.0, std::min(1.0, 0.5 * (eta + 1.0)));
    return L * alpha * alpha * (1.0 - alpha) * (1.0 - alpha) + floor;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void CellFlux(const amrex::Array4<const amrex::Real>& eta,
              const amrex::Array4<const amrex::Real>& mu,
              const amrex::Array4<const amrex::Real>& mask,
              int i, int j, int k, const Set::Scalar* dx,
              Set::Scalar L, Set::Scalar floor,
              const amrex::Array4<amrex::Real>& flux)
{
    const Set::Scalar M = mask(i,j,k) * DegenerateMobility(eta(i,j,k), L, floor);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
    {
        const int di = dir == 0, dj = dir == 1, dk = dir == 2;
        const Set::Scalar Mp = mask(i+di,j+dj,k+dk) * DegenerateMobility(eta(i+di,j+dj,k+dk),L,floor);
        const Set::Scalar Mm = mask(i-di,j-dj,k-dk) * DegenerateMobility(eta(i-di,j-dj,k-dk),L,floor);
        const Set::Scalar fp = 0.5 * (M + Mp) * (mu(i+di,j+dj,k+dk)-mu(i,j,k)) / dx[dir];
        const Set::Scalar fm = 0.5 * (M + Mm) * (mu(i,j,k)-mu(i-di,j-dj,k-dk)) / dx[dir];
        flux(i,j,k,dir) = 0.5 * (fp + fm);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar ConservativeRHS(const amrex::Array4<const amrex::Real>& eta,
                            const amrex::Array4<const amrex::Real>& mu,
                            const amrex::Array4<const amrex::Real>& mask,
                            int i, int j, int k, const Set::Scalar* dx,
                            Set::Scalar L, Set::Scalar floor)
{
    const Set::Scalar M = mask(i,j,k) * DegenerateMobility(eta(i,j,k), L, floor);
    Set::Scalar rhs = 0.0;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
    {
        const int di = dir == 0, dj = dir == 1, dk = dir == 2;
        const Set::Scalar plus = 0.5 * (M + mask(i+di,j+dj,k+dk) * DegenerateMobility(eta(i+di,j+dj,k+dk), L, floor));
        const Set::Scalar minus = 0.5 * (M + mask(i-di,j-dj,k-dk) * DegenerateMobility(eta(i-di,j-dj,k-dk), L, floor));
        rhs += (plus * (mu(i+di,j+dj,k+dk) - mu(i,j,k))
                - minus * (mu(i,j,k) - mu(i-di,j-dj,k-dk))) / (dx[dir] * dx[dir]);
    }
    return rhs;
}

void ChemicalPotential(const amrex::MultiFab& eta_mf, amrex::MultiFab& mu_mf,
                        const amrex::Geometry& geometry, Set::Scalar gamma)
{
    const auto dx = geometry.CellSizeArray();
    for (amrex::MFIter mfi(eta_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const auto eta = eta_mf.const_array(mfi);
        const auto mu = mu_mf.array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            mu(i,j,k) = eta(i,j,k)*eta(i,j,k)*eta(i,j,k) - eta(i,j,k)
                        - gamma * Numeric::Laplacian(eta, i, j, k, 0, dx.data());
        });
    }
}
}

void
CahnHilliard::AdvanceDegenerateReal (int lev, Set::Scalar time, Set::Scalar dt)
{
    std::swap(etaold_mf[lev], etanew_mf[lev]);
    etaold_mf[lev]->FillBoundary(geom[lev].periodicity());
    ChemicalPotential(*etaold_mf[lev], *intermediate[lev], geom[lev], gamma);
    intermediate[lev]->FillBoundary(geom[lev].periodicity());
    amrex::MultiFab mobility_mask(etanew_mf[lev]->boxArray(), etanew_mf[lev]->DistributionMap(), 1, 1);
    mobility_mask.setVal(1.0);
    FillMobilityMask(lev, mobility_mask);
    mobility_mask.FillBoundary(geom[lev].periodicity());
    std::unique_ptr<amrex::MultiFab> coarse_mobility_mask;
    if (lev > 0)
    {
        coarse_mobility_mask = std::make_unique<amrex::MultiFab>(
            etanew_mf[lev-1]->boxArray(), etanew_mf[lev-1]->DistributionMap(), 1, 1);
        coarse_mobility_mask->setVal(1.0);
        FillMobilityMask(lev-1, *coarse_mobility_mask);
        coarse_mobility_mask->FillBoundary(geom[lev-1].periodicity());
    }
    if (lev > 0)
    {
        amrex::Vector<amrex::MultiFab*> coarse{intermediate[lev-1].get()};
        amrex::Vector<amrex::MultiFab*> fine{intermediate[lev].get()};
        amrex::Vector<amrex::Real> times{time};
        bc->define(geom[lev]);
        amrex::Vector<amrex::BCRec> bcs{bc->GetBCRec()};
        amrex::FillPatchTwoLevels(*intermediate[lev], time, coarse, times, fine, times,
            0, 0, 1, geom[lev-1], geom[lev], *bc, 0, *bc, 0, refRatio(lev-1),
            &amrex::cell_cons_interp, bcs, 0);

        amrex::Vector<amrex::MultiFab*> coarse_mask{coarse_mobility_mask.get()};
        amrex::Vector<amrex::MultiFab*> fine_mask{&mobility_mask};
        amrex::FillPatchTwoLevels(mobility_mask, time, coarse_mask, times, fine_mask, times,
            0, 0, 1, geom[lev-1], geom[lev], *bc, 0, *bc, 0, refRatio(lev-1),
            &amrex::cell_cons_interp, bcs, 0);
    }

    const auto dx = geom[lev].CellSizeArray();
    const Set::Scalar L_local = L;
    const Set::Scalar floor = mobility_floor;
    const bool clipping_local = clipping;
    const bool diagnostics_local = output_diagnostics;
    for (amrex::MFIter mfi(*etanew_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const auto eta = etaold_mf[lev]->const_array(mfi);
        const auto mu = intermediate[lev]->const_array(mfi);
        const auto mask = mobility_mask.const_array(mfi);
        const auto next = etanew_mf[lev]->array(mfi);
        amrex::Array4<amrex::Real> out_m, out_mask, out_raw, out_delta, out_rhs, out_flux;
        if (diagnostics_local)
        {
            out_m = diagnostic_mobility[lev]->array(mfi);
            out_mask = diagnostic_mobility_mask[lev]->array(mfi);
            out_raw = diagnostic_raw[lev]->array(mfi);
            out_delta = diagnostic_delta[lev]->array(mfi);
            out_rhs = diagnostic_rhs[lev]->array(mfi);
            out_flux = diagnostic_flux[lev]->array(mfi);
        }
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            const Set::Scalar rhs = ConservativeRHS(eta, mu, mask, i, j, k, dx.data(), L_local, floor);
            const Set::Scalar raw = eta(i,j,k) + dt * rhs;
            next(i,j,k) = clipping_local ? std::max(-1.0, std::min(1.0, raw)) : raw;
            if (diagnostics_local)
            {
                const Set::Scalar alpha = std::max(0.0, std::min(1.0, 0.5*(eta(i,j,k)+1.0)));
                const Set::Scalar center_m = mask(i,j,k)*(L_local*alpha*alpha*(1.0-alpha)*(1.0-alpha)+floor);
                out_m(i,j,k) = center_m; out_mask(i,j,k) = mask(i,j,k);
                out_raw(i,j,k) = raw; out_delta(i,j,k) = next(i,j,k)-eta(i,j,k); out_rhs(i,j,k) = rhs;
                CellFlux(eta, mu, mask, i, j, k, dx.data(), L_local, floor, out_flux);
            }
        });
    }
    etanew_mf[lev]->FillBoundary(geom[lev].periodicity());
}

#ifdef ALAMO_FFT
void
CahnHilliard::AdvanceDegenerateSpectral (int lev, Set::Scalar time, Set::Scalar dt)
{
    using Operator::Spectral::FFT;
    FFT fft(geom[lev]);
    amrex::Vector<amrex::MultiFab*> hierarchy(lev + 1);
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> old_alpha(lev + 1);
    for (int ilev = 0; ilev <= lev; ++ilev)
    {
        hierarchy[ilev] = etanew_mf[ilev].get();
        if (output_diagnostics)
        {
            old_alpha[ilev] = std::make_unique<amrex::MultiFab>(etanew_mf[ilev]->boxArray(), etanew_mf[ilev]->DistributionMap(), 1, 0);
            amrex::MultiFab::Copy(*old_alpha[ilev], *etanew_mf[ilev], 0, 0, 1, 0);
        }
    }
    // Form derivatives on a complete grid so coarse/fine patch edges have valid data.
    auto eta_mf = FFT::CompositeToUniform(hierarchy, geom, refRatio(), lev, time, *bc, 1, 1);
    eta_mf->FillBoundary(geom[lev].periodicity());
    amrex::MultiFab mu_mf(eta_mf->boxArray(), eta_mf->DistributionMap(), 1, 1);
    ChemicalPotential(*eta_mf, mu_mf, geom[lev], gamma);
    mu_mf.FillBoundary(geom[lev].periodicity());
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> mask_hierarchy(lev + 1);
    amrex::Vector<amrex::MultiFab*> mask_pointers(lev + 1);
    for (int ilev = 0; ilev <= lev; ++ilev)
    {
        mask_hierarchy[ilev] = std::make_unique<amrex::MultiFab>(
            etanew_mf[ilev]->boxArray(), etanew_mf[ilev]->DistributionMap(), 1, 1);
        FillMobilityMask(ilev, *mask_hierarchy[ilev]);
        mask_hierarchy[ilev]->FillBoundary(geom[ilev].periodicity());
        mask_pointers[ilev] = mask_hierarchy[ilev].get();
    }
    auto mobility_mask = FFT::CompositeToUniform(mask_pointers, geom, refRatio(), lev, time, *bc, 1, 1);
    mobility_mask->FillBoundary(geom[lev].periodicity());
    std::unique_ptr<amrex::MultiFab> diagnostic_mobility_mf, diagnostic_flux_mf;
    if (output_diagnostics)
    {
        diagnostic_mobility_mf = std::make_unique<amrex::MultiFab>(eta_mf->boxArray(), eta_mf->DistributionMap(), 1, 0);
        diagnostic_flux_mf = std::make_unique<amrex::MultiFab>(eta_mf->boxArray(), eta_mf->DistributionMap(), AMREX_SPACEDIM, 0);
    }
    amrex::Vector<amrex::MultiFab*> potential_hierarchy(lev + 1);
    for (int ilev = 0; ilev <= lev; ++ilev) potential_hierarchy[ilev] = intermediate[ilev].get();
    FFT::UniformToComposite(mu_mf, potential_hierarchy, geom, refRatio(), lev, 1);
    amrex::MultiFab rhs_mf(eta_mf->boxArray(), eta_mf->DistributionMap(), 1, 0);
    const auto dx = geom[lev].CellSizeArray();
    const Set::Scalar L_local = L;
    const Set::Scalar floor = mobility_floor;
    const bool diagnostics_local = output_diagnostics;
    for (amrex::MFIter mfi(*eta_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const auto eta = eta_mf->const_array(mfi);
        const auto mu = mu_mf.const_array(mfi);
        const auto mask = mobility_mask->const_array(mfi);
        const auto rhs = rhs_mf.array(mfi);
        amrex::Array4<amrex::Real> diagnostic_mobility, diagnostic_flux;
        if (diagnostics_local)
        {
            diagnostic_mobility = diagnostic_mobility_mf->array(mfi);
            diagnostic_flux = diagnostic_flux_mf->array(mfi);
        }
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            // Face gradients and their adjoint divergence dissipate the discrete energy.
            rhs(i,j,k) = ConservativeRHS(eta, mu, mask, i, j, k, dx.data(), L_local, floor);
            if (diagnostics_local)
            {
                const Set::Scalar alpha = std::max(0.0, std::min(1.0, 0.5*(eta(i,j,k)+1.0)));
                const Set::Scalar center_m = mask(i,j,k)*(L_local*alpha*alpha*(1.0-alpha)*(1.0-alpha)+floor);
                diagnostic_mobility(i,j,k) = center_m;
                CellFlux(eta, mu, mask, i, j, k, dx.data(), L_local, floor, diagnostic_flux);
            }
        });
    }
    if (output_diagnostics)
    {
        amrex::Vector<amrex::MultiFab*> dmask(lev+1), dmob(lev+1), dflux(lev+1), drhs(lev+1);
        for (int ilev = 0; ilev <= lev; ++ilev)
        {
            dmask[ilev] = diagnostic_mobility_mask[ilev].get();
            dmob[ilev] = diagnostic_mobility[ilev].get();
            dflux[ilev] = diagnostic_flux[ilev].get();
            drhs[ilev] = diagnostic_rhs[ilev].get();
        }
        FFT::UniformToComposite(*mobility_mask, dmask, geom, refRatio(), lev, 1);
        FFT::UniformToComposite(*diagnostic_mobility_mf, dmob, geom, refRatio(), lev, 1);
        FFT::UniformToComposite(*diagnostic_flux_mf, dflux, geom, refRatio(), lev, AMREX_SPACEDIM);
        FFT::UniformToComposite(rhs_mf, drhs, geom, refRatio(), lev, 1);
    }
    auto eta_hat_mf = fft.MakeSpectralFab();
    auto rhs_hat_mf = fft.MakeSpectralFab();
    fft.Forward(*eta_mf, eta_hat_mf);
    fft.Forward(rhs_mf, rhs_hat_mf);
    const Set::Scalar S = spectral_stabilization;
    const auto length = geom[lev].Domain().length3d();
    for (amrex::MFIter mfi(eta_hat_mf, false); mfi.isValid(); ++mfi)
    {
        const auto eta_hat = eta_hat_mf.array(mfi);
        const auto rhs_hat = rhs_hat_mf.const_array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE(int m, int n, int p)
        {
            // The face-flux Laplacian symbol includes both odd-grid and Nyquist modes.
            const Set::Scalar sx = std::sin(Set::Constant::Pi * m / length[0]);
            const Set::Scalar sy = std::sin(Set::Constant::Pi * n / length[1]);
            amrex::ignore_unused(p);
#if AMREX_SPACEDIM == 3
            const Set::Scalar sz = std::sin(Set::Constant::Pi * p / length[2]);
#endif
            const Set::Scalar lap = AMREX_D_TERM(4.0 * sx*sx / (dx[0]*dx[0]),
                                                + 4.0 * sy*sy / (dx[1]*dx[1]),
                                                + 4.0 * sz*sz / (dx[2]*dx[2]));
            // Set the conserved mode exactly; the face flux sum is zero up to roundoff.
            if (m != 0 || n != 0 || p != 0)
                eta_hat(m,n,p) += dt * rhs_hat(m,n,p) / (1.0 + dt * S * lap * lap);
        });
    }
    fft.Backward(eta_hat_mf, *eta_mf);
    if (output_diagnostics)
    {
        amrex::Vector<amrex::MultiFab*> raw_hierarchy(lev+1);
        for (int ilev = 0; ilev <= lev; ++ilev) raw_hierarchy[ilev] = diagnostic_raw[ilev].get();
        FFT::UniformToComposite(*eta_mf, raw_hierarchy, geom, refRatio(), lev, 1);
    }
    FFT::UniformToComposite(*eta_mf, hierarchy, geom, refRatio(), lev, 1);
    if (clipping)
        for (int ilev = 0; ilev <= lev; ++ilev)
            for (amrex::MFIter mfi(*etanew_mf[ilev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const auto eta = etanew_mf[ilev]->array(mfi);
                amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k)
                { eta(i,j,k) = std::max(-1.0, std::min(1.0, eta(i,j,k))); });
            }
    if (output_diagnostics)
        for (int ilev = 0; ilev <= lev; ++ilev)
        {
            amrex::MultiFab::Copy(*diagnostic_delta[ilev], *etanew_mf[ilev], 0, 0, 1, 0);
            amrex::MultiFab::Subtract(*diagnostic_delta[ilev], *old_alpha[ilev], 0, 0, 1, 0);
        }
    for (int ilev = 0; ilev <= lev; ++ilev) etanew_mf[ilev]->FillBoundary(geom[ilev].periodicity());
}
#else
void
CahnHilliard::AdvanceDegenerateSpectral (int, Set::Scalar, Set::Scalar)
{
    Util::Abort(INFO,"Alamo must be compiled with fft");
}
#endif

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
    InitializeDiagnostics(lev);
}

void
CahnHilliard::InitializeDiagnostics(int lev)
{
    if (!output_diagnostics) return;
    amrex::MultiFab mask(etanew_mf[lev]->boxArray(), etanew_mf[lev]->DistributionMap(), 1, 0);
    FillMobilityMask(lev, mask);
    const bool sdch = mobility == "singly_degenerate";
    const Set::Scalar L_local = L, floor = mobility_floor;
    for (amrex::MFIter mfi(*etanew_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const auto eta = etanew_mf[lev]->const_array(mfi);
        const auto m = mask.const_array(mfi);
        const auto out_m = diagnostic_mobility[lev]->array(mfi);
        const auto out_mask = diagnostic_mobility_mask[lev]->array(mfi);
        const auto out_raw = diagnostic_raw[lev]->array(mfi);
        const auto out_delta = diagnostic_delta[lev]->array(mfi);
        const auto out_rhs = diagnostic_rhs[lev]->array(mfi);
        const auto out_flux = diagnostic_flux[lev]->array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            const Set::Scalar alpha = sdch ? std::max(0.0, std::min(1.0, 0.5*(eta(i,j,k)+1.0))) : 0.5;
            out_mask(i,j,k) = m(i,j,k);
            out_m(i,j,k) = sdch ? m(i,j,k)*(L_local*alpha*alpha*(1.0-alpha)*(1.0-alpha)+floor) : L_local;
            out_raw(i,j,k) = eta(i,j,k); out_delta(i,j,k) = 0.0; out_rhs(i,j,k) = 0.0;
            AMREX_D_TERM(out_flux(i,j,k,0)=0.0;, out_flux(i,j,k,1)=0.0;, out_flux(i,j,k,2)=0.0;);
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
