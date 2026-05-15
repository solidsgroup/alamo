#include <AMReX_MLPoisson.H>
#include <AMReX_Reduce.H>

#include "BC/Constant.H"
#include "DoublyDegenerateCahnHilliard.H"
#include "IC/Ellipse.H"
#include "IC/Random.H"
#include "IO/ParmParse.H"
#include "Numeric/Stencil.H"
#include "Set/Set.H"

namespace Integrator
{
DoublyDegenerateCahnHilliard::DoublyDegenerateCahnHilliard() : Integrator()
{
}
DoublyDegenerateCahnHilliard::~DoublyDegenerateCahnHilliard()
{
    delete ic;
    delete bc;
}

void DoublyDegenerateCahnHilliard::Parse(DoublyDegenerateCahnHilliard &value, IO::ParmParse &pp)
{
    // Interface energy
    pp.query_default("epsilon", value.epsilon, 0.0005);
    // Chemical potential parameter
    pp.query_default("omega", value.omega, 72.0);
    // Energy restriction parameter
    pp.query_default("gamma", value.gamma, 6.0);
    // Energy restriction exponent
    // pp.query_default("p", value.p, 1.0);
    // Energy restriction normalization parameter
    pp.query_default("alpha", value.alpha, 1e-4);
    // Mobility parameter
    pp.query_default("mu", value.mu, 36.0);

    // initial condition for :math:`\eta`
    pp.select_default<IC::Ellipse, IC::Random>("eta.ic", value.ic, value.geom);
    // boundary condition for :math:`\eta`
    pp.select_default<BC::Constant>("eta.bc", value.bc, 1);

    value.RegisterNewFab(value.etanew_mf, value.bc, 1, 1, "eta", true);
    value.RegisterNewFab(value.etaold_mf, value.bc, 1, 1, "eta_old", false);
    value.RegisterNewFab(value.deta_mf, value.bc, 1, 1, "deta", true);
    value.RegisterNewFab(value.driving_force_mf, value.bc, 1, 1, "driving_force", true, false);
    value.AddField<Set::Vector, Set::Hypercube::Cell>(value.evolution_driving_force_mf, nullptr, 1, 1, "evolution_driving_force", true, false);

    value.RegisterIntegratedVariable(&value.volume, "volume");
    value.integrate_variables_before_advance = false;
    if (value.thermo.plot_int <= 0 && value.thermo.plot_dt <= 0.0)
    {
        if (value.plot_int > 0)
            value.thermo.plot_int = value.plot_int;
        if (value.plot_dt > 0.0)
            value.thermo.plot_dt = value.plot_dt;
    }
}

void
DoublyDegenerateCahnHilliard::Initialize(int lev)
{
    driving_force_mf[lev]->setVal(0.0);
    deta_mf[lev]->setVal(0.0);
    ic->Initialize(lev, etanew_mf);
    ic->Initialize(lev, etaold_mf);
}

void
DoublyDegenerateCahnHilliard::TimeStepBegin(Set::Scalar time, int step)
{
    if (step == 0)
        IntegrateVariables(time, step);
}

void
DoublyDegenerateCahnHilliard::TimeStepComplete(Set::Scalar time, int step)
{
    IntegrateVariables(time + dt[0], step + 1);
}

void
DoublyDegenerateCahnHilliard::Advance(int lev, Set::Scalar /*time*/, Set::Scalar dt)
{
    std::swap(etaold_mf[lev], etanew_mf[lev]);
    const Set::Scalar *DX = geom[lev].CellSize();
    for (amrex::MFIter mfi(*etanew_mf[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.tilebox();
        Set::Patch<const Set::Scalar> &etaold = etaold_mf[lev]->array(mfi);
        Set::Patch<Set::Scalar> &driving_force = driving_force_mf[lev]->array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            Set::Scalar eta = etaold(i, j, k);
            Set::Scalar g_alpha_denominator = sqrt(gamma * gamma * (1 - eta) * (1 - eta) * eta * eta + alpha * alpha * epsilon * epsilon);
            Set::Scalar g_alpha = 1.0 / g_alpha_denominator;
            Set::Scalar g_alphaprime = -gamma * gamma * (eta - 1) * eta * (2 * eta - 1) / (g_alpha_denominator * g_alpha_denominator * g_alpha_denominator);
            Set::Scalar f = omega / 4 * eta * eta * (1 - eta) * (1 - eta);
            Set::Scalar fprime = omega / 2 * eta * (1 - eta) * (1 - 2 * eta);
            Set::Vector grad_eta = Numeric::Gradient(etaold, i, j, k, 0, DX);
            Set::Scalar norm_squared_grad_eta = grad_eta.squaredNorm();
            Set::Scalar lap_eta = Numeric::Laplacian(etaold, i, j, k, 0, DX);

            driving_force(i, j, k) = g_alpha * fprime / epsilon - epsilon * (g_alphaprime * norm_squared_grad_eta + g_alpha * lap_eta) + g_alphaprime * (1 / epsilon * f + epsilon / 2 * norm_squared_grad_eta);
        });
    }

    driving_force_mf[lev]->FillBoundary(geom[lev].periodicity());

    // Uncomment the next two `for` loops to use the non-product rule formulation
    //
    for (amrex::MFIter mfi(*etanew_mf[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.tilebox();
        Set::Patch<const Set::Scalar> &etaold = etaold_mf[lev]->array(mfi);
        Set::Patch<const Set::Scalar> &driving_force = driving_force_mf[lev]->array(mfi);
        Set::Patch<Set::Vector> evolution_driving_force = evolution_driving_force_mf[lev]->array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            Set::Scalar eta = etaold(i, j, k);
            Set::Scalar M_alpha = mu * eta * eta * (1 - eta) * (1 - eta) + alpha * epsilon;
            Set::Vector grad_driving_force = Numeric::Gradient(driving_force, i, j, k, 0, DX);

            evolution_driving_force(i, j, k) = M_alpha * grad_driving_force;
        });
    }

    evolution_driving_force_mf[lev]->FillBoundary(geom[lev].periodicity());

    for (amrex::MFIter mfi(*etanew_mf[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.tilebox();
        Set::Patch<const Set::Vector> evolution_driving_force = evolution_driving_force_mf[lev]->array(mfi);
        Set::Patch<Set::Scalar> &etanew = etanew_mf[lev]->array(mfi);
        Set::Patch<const Set::Scalar> &etaold = etaold_mf[lev]->array(mfi);
        Set::Patch<Set::Scalar> &deta = deta_mf[lev]->array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            etanew(i, j, k) = etaold(i, j, k) + dt / epsilon * Numeric::Divergence(evolution_driving_force, i, j, k, DX);
            deta(i, j, k) = etanew(i, j, k) - etaold(i, j, k);
        });
    }

    // Uncomment the next `for` loop to use the product rule formulation
    //
    // for (amrex::MFIter mfi(*etanew_mf[lev], true); mfi.isValid(); ++mfi)
    // {
    //     const amrex::Box &bx = mfi.tilebox();
    //     Set::Patch<Set::Scalar> &etanew = etanew_mf[lev]->array(mfi);
    //     Set::Patch<const Set::Scalar> &etaold = etaold_mf[lev]->array(mfi);
    //     Set::Patch<const Set::Scalar> &driving_force = driving_force_mf[lev]->array(mfi);
    //     Set::Patch<Set::Scalar> &deta = deta_mf[lev]->array(mfi);

    //     amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
    //         Set::Scalar eta = etaold(i, j, k);
    //         Set::Vector grad_eta = Numeric::Gradient(etaold, i, j, k, 0, DX);
    //         Set::Scalar M_alpha = mu * eta * eta * (1 - eta) * (1 - eta) + alpha * epsilon;
    //         Set::Scalar M_alphaprime = 2.0 * mu * eta * (eta - 1) * (2 * eta - 1);
    //         Set::Vector grad_driving_force = Numeric::Gradient(driving_force, i, j, k, 0, DX);
    //         Set::Scalar lap_driving_force = Numeric::Laplacian(driving_force, i, j, k, 0, DX);

    //         etanew(i, j, k) = etaold(i, j, k) + dt / epsilon * (M_alphaprime * grad_eta.dot(grad_driving_force) + M_alpha * lap_driving_force);
    //         deta(i, j, k) = etanew(i, j, k) - etaold(i, j, k);
    //     });
    // }
}

void
DoublyDegenerateCahnHilliard::Integrate(int amrlev, Set::Scalar /*time*/, int /*step*/, const amrex::MFIter &mfi, const amrex::Box &box)
{
    BL_PROFILE("DoublyDegenerateCahnHilliard::Integrate");

    const Set::Scalar *DX = geom[amrlev].CellSize();
    Set::Scalar dv = AMREX_D_TERM(DX[0], *DX[1], *DX[2]);
    Set::Patch<const Set::Scalar> eta = etanew_mf[amrlev]->array(mfi);

    amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
    amrex::ReduceData<Set::Scalar> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    reduce_op.eval(box, reduce_data, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
        return { eta(i, j, k, 0) * dv };
    });
    ReduceTuple hv = reduce_data.value(reduce_op);
    volume += amrex::get<0>(hv);
}

void
DoublyDegenerateCahnHilliard::TagCellsForRefinement(int /*lev*/, amrex::TagBoxArray & /*tags*/, amrex::Real /*time*/, int /*ngrow*/)
{
}
}
