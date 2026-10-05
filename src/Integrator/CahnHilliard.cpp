#include <AMReX_MLPoisson.H>
#include <algorithm>
#include <cmath>

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

    // initial condition for :math:`\eta`
    pp.select_default<IC::Random,IC::Expression>("eta.ic", value.ic, pp.forward_args(value.geom));
    // boundary condition for :math:`\eta`
    pp.select_default<BC::Constant>("eta.bc", value.bc, pp.forward_args(1));

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

    value.RegisterNewFab(value.etanew_mf, value.bc, 1, 1, "eta",true);
    value.RegisterNewFab(value.intermediate, value.bc, 1, 1, "int",true);

    if (value.method == "realspace")
        value.RegisterNewFab(value.etaold_mf, value.bc, 1, 1, "eta_old",false);
}

void
CahnHilliard::Advance(int lev, Set::Scalar time, Set::Scalar dt)
{
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
Set::Scalar ConservativeRHS(const amrex::Array4<const amrex::Real>& eta,
                            const amrex::Array4<const amrex::Real>& mu,
                            int i, int j, int k, const Set::Scalar* dx,
                            Set::Scalar L, Set::Scalar floor)
{
    const Set::Scalar M = DegenerateMobility(eta(i,j,k), L, floor);
    Set::Scalar rhs = 0.0;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
    {
        const int di = dir == 0, dj = dir == 1, dk = dir == 2;
        const Set::Scalar plus = 0.5 * (M + DegenerateMobility(eta(i+di,j+dj,k+dk), L, floor));
        const Set::Scalar minus = 0.5 * (M + DegenerateMobility(eta(i-di,j-dj,k-dk), L, floor));
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
    }

    const auto dx = geom[lev].CellSizeArray();
    const Set::Scalar L_local = L;
    const Set::Scalar floor = mobility_floor;
    for (amrex::MFIter mfi(*etanew_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const auto eta = etaold_mf[lev]->const_array(mfi);
        const auto mu = intermediate[lev]->const_array(mfi);
        const auto next = etanew_mf[lev]->array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            next(i,j,k) = eta(i,j,k) + dt * ConservativeRHS(eta, mu, i, j, k, dx.data(), L_local, floor);
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
    for (int ilev = 0; ilev <= lev; ++ilev) hierarchy[ilev] = etanew_mf[ilev].get();
    // Form derivatives on a complete grid so coarse/fine patch edges have valid data.
    auto eta_mf = FFT::CompositeToUniform(hierarchy, geom, refRatio(), lev, time, *bc, 1, 1);
    eta_mf->FillBoundary(geom[lev].periodicity());
    amrex::MultiFab mu_mf(eta_mf->boxArray(), eta_mf->DistributionMap(), 1, 1);
    ChemicalPotential(*eta_mf, mu_mf, geom[lev], gamma);
    mu_mf.FillBoundary(geom[lev].periodicity());
    amrex::Vector<amrex::MultiFab*> potential_hierarchy(lev + 1);
    for (int ilev = 0; ilev <= lev; ++ilev) potential_hierarchy[ilev] = intermediate[ilev].get();
    FFT::UniformToComposite(mu_mf, potential_hierarchy, geom, refRatio(), lev, 1);
    amrex::MultiFab rhs_mf(eta_mf->boxArray(), eta_mf->DistributionMap(), 1, 0);
    const auto dx = geom[lev].CellSizeArray();
    const Set::Scalar L_local = L;
    const Set::Scalar floor = mobility_floor;
    for (amrex::MFIter mfi(*eta_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const auto eta = eta_mf->const_array(mfi);
        const auto mu = mu_mf.const_array(mfi);
        const auto rhs = rhs_mf.array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            // Face gradients and their adjoint divergence dissipate the discrete energy.
            rhs(i,j,k) = ConservativeRHS(eta, mu, i, j, k, dx.data(), L_local, floor);
        });
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
    FFT::UniformToComposite(*eta_mf, hierarchy, geom, refRatio(), lev, 1);
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
