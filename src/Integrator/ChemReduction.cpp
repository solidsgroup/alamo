#include "ChemReduction.H"
#include "IO/ParmParse.H"
#include "BC/Constant.H"
#include "Numeric/Stencil.H"
#include "Numeric/Function.H"
#include "IC/Laminate.H"
#include "IC/Constant.H"
#include "IC/PSRead.H"
#include "IC/Expression.H"
#include "IC/BMP.H"
#include <cmath>

namespace Integrator
{

ChemReduction::ChemReduction() : Base::Mechanics<model_type>() {}

ChemReduction::ChemReduction(IO::ParmParse& pp) : ChemReduction()
{
    pp.queryclass(*this);
}

// [parser]
void
ChemReduction::Parse(Flame& value, IO::ParmParse& pp)
{
    BL_PROFILE("Integrator::ChemReduction::ChemReduction()");
    {
        pp.query("geometry.y_len", value.y_len); // Domain y length

        value.bc_eta = new BC::Constant(1);
        pp.queryclass("pf.phi.bc", *static_cast<BC::Constant*>(value.bc_phi)); // See :ref:`BC::Constant`
        value.RegisterNewFab(value.phi_o_mf, value.bc_phi, 1, value.ghost_count + 1, "phiO", true);
        value.RegisterNewFab(value.phi_fe_mf, value.bc_phi, 1, value.ghost_count + 1, "phiFe", true);
        value.RegisterNewFab(value.phi_h_mf, value.bc_phi, 1, value.ghost_count + 1, "phiH", true);

        std::string phi_bc_str = "constant";
        pp.query("pf.phi.ic.type", eta_bc_str); // Eta boundary condition [constant, expression]
        if (phi_bc_str == "constant") value.ic_phi = new IC::Constant(value.geom, pp, "pf.phi.ic.constant");
        else if (phi_bc_str == "expression") value.ic_phi = new IC::Expression(value.geom, pp, "pf.phi.ic.expression");

        std::string phi_ic_type = "constant";
        pp.query("phi.ic.type", eta_ic_type); // Eta initial condition [constant, laminate, expression, bmp]
        else if (phi_ic_type == "constant") value.ic_eta = new IC::Constant(value.geom, pp, "phi.ic.constant");
        else if (phi_ic_type == "expression") value.ic_eta = new IC::Expression(value.geom, pp, "phi.ic.expression");
        else if (phi_ic_type == "bmp") value.ic_eta = new IC::BMP(value.geom, pp, "phi.ic.bmp");
        else if (phi_ic_type == "png") value.ic_eta = new IC::PNG(value.geom, pp, "phi.ic.png");
        else Util::Abort(INFO, "Invalid phi IC type", phi_ic_type);
    }

    pp.queryclass("thermal.temp.bc", *static_cast<BC::Constant*>(value.bc_temp));
    value.RegisterNewFab(value.phi_o_mf, value.bc_phi, 1, value.ghost_count + 1, "phiO", true);
    value.RegisterNewFab(value.phi_fe_mf, value.bc_phi, 1, value.ghost_count + 1, "phiFe", true);
    value.RegisterNewFab(value.phi_h_mf, value.bc_phi, 1, value.ghost_count + 1, "phiH", true);

    //value.RegisterNewFab(value.temp_old_mf, value.bc_temp, 1, value.ghost_count + 1, "temp_old", false);
    //value.RegisterIntegratedVariable(&value.volume, "total_area");

}

void ChemReduction::Initialize(int lev)
{
    BL_PROFILE("Integrator::Flame::Initialize");
    ic_eta->Initialize(lev, eta_mf);
    ic_eta->Initialize(lev, eta_old_mf);
    ic_phi->Initialize(lev, phi_mf);
    psi_mf[lev]->setVal(1.0);
}

void ChemReduction::UpdateModel(int /*a_step*/, Set::Scalar /*a_time*/)
{
}

void ChemReduction::TimeStepBegin(Set::Scalar a_time, int a_iter)
{
    BL_PROFILE("Integrator::ChemReduction::TimeStepBegin");
}

void ChemReduction::TimeStepComplete(Set::Scalar /*a_time*/, int /*a_iter*/)
{
    BL_PROFILE("Integrator::ChemReduction::TimeStepComplete");
}

void ChemReduction::Advance(int lev, Set::Scalar time, Set::Scalar dt)
{
    BL_PROFILE("Integrador::ChemReduction::Advance");

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
    {
        Set::Scalar phi_avg = Numeric::Interpolate::NodeToCellAverage(phi, i, j, k, 0);
        Set::Scalar eta_lap = Numeric::Laplacian(eta, i, j, k, 0, DX);
        Set::Scalar K; // Calculate effective thermal conductivity
        Set::Scalar rho; // No special interface mixure rule is needed here.
        Set::Scalar cp;
        Set::Scalar df_deta = ((pf.lambda / pf.eps) * dw(eta(i, j, k)) - pf.eps * pf.kappa * eta_lap);
        if (isnan(etanew(i, j, k)) || isnan(alpha(i, j, k)) || isnan(mdot(i, j, k))) {
            Util::Message(INFO, etanew(i, j, k), "etanew contains nan (i=", i, " j= ", j, ")");
            Util::Message(INFO, mdot(i, j, k), "mdot contains nan (i=", i, " j= ", j, ")");
            Util::Message(INFO, alpha(i, j, k), "alpha contains nan (i=", i, " j= ", j, ")");
            Util::Abort(INFO);
        }
    });

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
    {
        auto sten = Numeric::GetStencil(i, j, k, bx);
        Set::Vector grad_eta = Numeric::Gradient(eta, i, j, k, 0, DX);
        Set::Vector grad_temp = Numeric::Gradient(temp, i, j, k, 0, DX);
        Set::Scalar lap_temp = Numeric::Laplacian(temp, i, j, k, 0, DX);
        Set::Scalar grad_eta_mag = grad_eta.lpNorm<2>();
        Set::Vector grad_alpha = Numeric::Gradient(alpha, i, j, k, 0, DX, sten);
        Set::Scalar dTdt = 0.0;
        dTdt += grad_eta.dot(grad_temp * alpha(i, j, k));
        dTdt += grad_alpha.dot(eta(i, j, k) * grad_temp);
        dTdt += eta(i, j, k) * alpha(i, j, k) * lap_temp;
        dTdt += alpha(i, j, k) * heatflux(i, j, k) * grad_eta_mag;
        Set::Scalar Tsolid;
        Tsolid = dTdt + temps(i, j, k) * (etanew(i, j, k) - eta(i, j, k)) / dt;
        tempsnew(i, j, k) = temps(i, j, k) + dt * Tsolid;
        tempnew(i, j, k) = etanew(i, j, k) * tempsnew(i, j, k) + (1.0 - etanew(i, j, k)) * thermal.T_fluid;
        if (isnan(tempsnew(i, j, k)) || isnan(temps(i, j, k))) {
            Util::Message(INFO, tempsnew(i, j, k), "tempsnew contains nan (i=", i, " j= ", j, ")");
            Util::Message(INFO, temps(i, j, k), "temps contains nan (i=", i, " j= ", j, ")");
            Util::Abort(INFO);
        }

    });
} //Function


void ChemReduction::TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags, Set::Scalar time, int ngrow)
{
    BL_PROFILE("Integrator::ChemReduction::TagCellsForRefinement");
    const Set::Scalar* DX = geom[lev].CellSize();
    Set::Scalar dr = sqrt(AMREX_D_TERM(DX[0] * DX[0], +DX[1] * DX[1], +DX[2] * DX[2]));

    // Eta criterion for refinement
    for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<char> const& tags = a_tags.array(mfi);
        amrex::Array4<const Set::Scalar> const& eta = (*eta_mf[lev]).array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Vector gradeta = Numeric::Gradient(eta, i, j, k, 0, DX);
            if (gradeta.lpNorm<2>() * dr * 2 > m_refinement_criterion && eta(i, j, k) >= t_refinement_restriction)
                tags(i, j, k) = amrex::TagBox::SET;
        });
    }

}

void ChemReduction::Regrid(int lev, Set::Scalar time)
{
    BL_PROFILE("Integrator::ChemReduction::Regrid");
    ic_phi->Initialize(lev, phi_mf, time);
}

void ChemReduction::Integrate(int amrlev, Set::Scalar /*time*/, int /*step*/,
    const amrex::MFIter& mfi, const amrex::Box& box)
{
    BL_PROFILE("ChemReduction::Integrate");
    const Set::Scalar* DX = geom[amrlev].CellSize();
    Set::Scalar dv = AMREX_D_TERM(DX[0], *DX[1], *DX[2]);
    amrex::Array4<amrex::Real> const& eta = (*eta_mf[amrlev]).array(mfi);
    amrex::Array4<amrex::Real> const& mdot = (*mdot_mf[amrlev]).array(mfi);
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
    {
        volume += eta(i, j, k, 0) * dv;
        Set::Vector grad = Numeric::Gradient(eta, i, j, k, 0, DX);
        Set::Scalar normgrad = grad.lpNorm<2>();
        Set::Scalar da = normgrad * dv;
        area += da;

        Set::Vector mgrad = Numeric::Gradient(mdot, i, j, k, 0, DX);
        Set::Scalar mnormgrad = mgrad.lpNorm<2>();
        Set::Scalar dm = mnormgrad * dv;
        massflux += dm;

    });

}

} // namespace Integrator


