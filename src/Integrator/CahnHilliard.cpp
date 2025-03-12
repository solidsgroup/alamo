#include <AMReX_MLPoisson.H>

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
        AdvanceSpectral(lev, time, dt);
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

void
CahnHilliard::AdvanceSpectral (int lev, Set::Scalar /*time*/, Set::Scalar dt)
{
#ifdef ALAMO_FFT
    //
    // FFT Boilerplate
    //
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
    amrex::FabArray<amrex::BaseFab<amrex::GpuComplex<Set::Scalar> > > eta_hat_mf(cba, cdm, 1, 0);
    my_fft.forward(*etanew_mf[lev], eta_hat_mf);

    //
    // FFT of chemical potential gradient
    //
    amrex::FabArray<amrex::BaseFab<amrex::GpuComplex<Set::Scalar> > > chempot_hat_mf(cba, cdm, 1, 0);
    my_fft.forward(*intermediate[lev], chempot_hat_mf);

    //
    // Perform update in spectral coordinatees
    //
    //for (amrex::MFIter mfi(eta_hat_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    for (amrex::MFIter mfi(eta_hat_mf, false); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.tilebox();

        Util::Message(INFO, bx);
        
        amrex::Array4<amrex::GpuComplex<Set::Scalar>> const & eta_hat     =  eta_hat_mf.array(mfi);
        amrex::Array4<amrex::GpuComplex<Set::Scalar>> const & chempot_hat =  chempot_hat_mf.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int m, int n, int p) {


            // Get spectral coordinates
            AMREX_D_TERM(
                Set::Scalar k1 = m * pi_Lx;,
                Set::Scalar k2 = (n < domain.length(1)/2 ? n * pi_Ly : (n - domain.length(1)) * pi_Ly);,
                Set::Scalar k3 = (p < domain.length(2)/2 ? p * pi_Lz : (p - domain.length(2)) * pi_Lz););


            Set::Scalar lap = AMREX_D_TERM(k1 * k1, + k2 * k2, + k3*k3);
                        
            Set::Scalar bilap = lap*lap;

            eta_hat(m, n, p) =
                (eta_hat(m, n, p) - L * lap * chempot_hat(m, n, p) * dt) /
                (1.0 + L * gamma * bilap * dt);
                

            eta_hat(m,n,p) *= scaling;
        });
    }

    //
    // Transform solution back to realspace
    //
    my_fft.backward(eta_hat_mf, *etanew_mf[lev]);
#else

    Util::Abort(INFO,"Alamo must be compiled with fft");

#endif
}


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
