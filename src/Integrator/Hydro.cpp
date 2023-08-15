#include "Hydro.H"
#include "IO/ParmParse.H"
#include "BC/Constant.H"
#include "Numeric/Stencil.H"
#include "IC/Laminate.H"
#include "IC/PSRead.H"
#include "IC/Expression.H"
#include "Solver/Local/Riemann/Roe.H"

namespace Integrator
{

Hydro::Hydro(IO::ParmParse& pp) : Hydro()
{
    pp.queryclass(*this);
}

void
Hydro::Parse(Hydro& value, IO::ParmParse& pp)
{
    BL_PROFILE("Integrator::Hydro::Hydro()");
    {
        pp.query("r_refinement_criterion", value.r_refinement_criterion);
        pp.query("e_refinement_criterion", value.e_refinement_criterion);
        pp.query("m_refinement_criterion", value.m_refinement_criterion);
        pp.query("eta_refinement_criterion", value.eta_refinement_criterion);

        pp.query("gamma", value.gamma);
        pp.query("cfl", value.cfl);

        pp.query("rho_solid", value.rho_solid);
        pp.query("rho_fluid", value.rho_fluid);
        pp.query("E_solid", value.E_solid);
        pp.query("E_fluid", value.E_fluid);

        pp.query("Mx_init", value.Mx_init);
        pp.query("My_init", value.My_init);

        pp.query("num_cells_x", value.num_cells_x);

        pp.query("eps", value.eps);

        pp.query("mdot", value.mdot);
        pp.query("Pdot_x", value.Pdot_x);
        pp.query("Pdot_y", value.Pdot_y);
        pp.query("Qdot", value.Qdot);

        value.bc_eta = new BC::Constant(1, pp, "pf.eta.bc");
        value.bc_rho = new BC::Constant(1, pp, "rho.bc");
        value.bc_E = new BC::Constant(1, pp, "E.bc");
        value.bc_M = new BC::Constant(2, pp, "M.bc");
        value.bc_v = new BC::Constant(2, pp, "v.bc");
        value.bc_p = new BC::Constant(1, pp, "p.bc");

    }
    // Register FabFields:
    {
        int nghost = 2;
        value.RegisterNewFab(value.eta_mf, value.bc_eta, 1, nghost, "eta", true);

        value.RegisterNewFab(value.etadot_mf, value.bc_eta, 1, nghost, "etadot", false);

        value.RegisterNewFab(value.Density_mf, value.bc_rho, 1, nghost, "Density", true);
        value.RegisterNewFab(value.Density_old_mf, value.bc_rho, 1, nghost, "rho_old", false);

        value.RegisterNewFab(value.Energy_mf, value.bc_E, 1, nghost, "Energy", true);
        value.RegisterNewFab(value.Energy_old_mf, value.bc_E, 1, nghost, "E_old", false);

        value.RegisterNewFab(value.Momentum_mf, value.bc_M, 2, nghost, "Momentum", true);
        value.RegisterNewFab(value.Momentum_old_mf, value.bc_M, 2, nghost, "M_old", false);

        value.RegisterNewFab(value.Velocity_mf, value.bc_v, 2, nghost, "Velocity", true);

        value.RegisterNewFab(value.Pressure_mf, value.bc_p, 1, nghost, "Pressure", true);
    }
    {
        BL_PROFILE("Integrator::Hydro::Hydro()");
        std::string type = "constant";
        pp.query("eta.ic.type", type);
        if (type == "constant") value.ic_eta = new IC::Constant(value.geom, pp, "eta.ic.constant");
        else if (type == "laminate") value.ic_eta = new IC::Laminate(value.geom, pp, "eta.ic.laminate");
        else Util::Abort(INFO, "Invalid eta.ic: ", type);
    }

}


void Hydro::Initialize(int lev)
{
    BL_PROFILE("Integrator::Hydro::Initialize");

    //ic_eta->Initialize(lev, eta_mf);

    for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<Set::Scalar> const& eta = (*eta_mf[lev]).array(mfi);

        amrex::Array4<Set::Scalar> const& etadot = (*etadot_mf[lev]).array(mfi);

        amrex::Array4<Set::Scalar> const& E_new = (*Energy_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const& E = (*Energy_old_mf[lev]).array(mfi);

        amrex::Array4<Set::Scalar> const& rho_new = (*Density_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const& rho = (*Density_old_mf[lev]).array(mfi);

        amrex::Array4<Set::Scalar> const& M_new = (*Momentum_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const& M = (*Momentum_old_mf[lev]).array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            const Set::Scalar* DX = geom[lev].CellSize();
            Set::Scalar pos_x0 = (i - num_cells_x / 2.) * DX[0];

            eta(i, j, k) = 0.5 - 0.5 * std::erf(pos_x0 / eps);
            etadot(i, j, k) = 0.0;

            rho(i, j, k) = rho_solid * (1 - eta(i, j, k)) + rho_fluid * eta(i, j, k);
            rho_new(i, j, k) = rho(i, j, k);

            M(i, j, k, 0) = 0.0 * (1 - eta(i, j, k)) + Mx_init * eta(i, j, k);
            M_new(i, j, k, 0) = M(i, j, k, 0);
            ///
            M(i, j, k, 1) = 0.0 * (1 - eta(i, j, k)) + My_init * eta(i, j, k);
            M_new(i, j, k, 1) = M(i, j, k, 1);

            E(i, j, k) = E_solid * (1 - eta(i, j, k)) + E_fluid * eta(i, j, k);
            E_new(i, j, k) = E(i, j, k);
        });
    }

    c_max = 0.0;
    vx_max = 0.0;
    vy_max = 0.0;
}

void Hydro::TimeStepBegin(Set::Scalar, int)
{
    BL_PROFILE("Integrator::Hydro::TimeStepBegin");
}


void Hydro::TimeStepComplete(Set::Scalar, int lev)
{
    return;

    const Set::Scalar* DX = geom[lev].CellSize();

    amrex::ParallelDescriptor::ReduceRealMax(c_max);
    amrex::ParallelDescriptor::ReduceRealMax(vx_max);
    amrex::ParallelDescriptor::ReduceRealMax(vy_max);

    Set::Scalar new_timestep = cfl / ((c_max + vx_max) / DX[0] + (c_max + vy_max) / DX[1]);

    Util::Assert(INFO, TEST(AMREX_SPACEDIM == 2));

    SetTimestep(new_timestep);
}



void Hydro::Advance(int lev, Set::Scalar, Set::Scalar dt)
{
    std::swap(Momentum_old_mf[lev], Momentum_mf[lev]);
    ///
    std::swap(Energy_old_mf[lev], Energy_mf[lev]);
    ///
    std::swap(Density_old_mf[lev], Density_mf[lev]);

    const Set::Scalar* DX = geom[lev].CellSize();

    for (amrex::MFIter mfi(*eta_mf[lev], false); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        amrex::Array4<const Set::Scalar> const& E = (*Energy_old_mf[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& rho = (*Density_old_mf[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& M = (*Momentum_old_mf[lev]).array(mfi);

        amrex::Array4<Set::Scalar> const& v = (*Velocity_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const& p = (*Pressure_mf[lev]).array(mfi);

        //Compute primitive variables

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            v(i, j, k, 0) = M(i, j, k, 0) / rho(i, j, k);
            v(i, j, k, 1) = M(i, j, k, 1) / rho(i, j, k);

            Set::Scalar ke = 0.5 * (v(i, j, k, 0) * v(i, j, k, 0) + v(i, j, k, 1) * v(i, j, k, 1));

            p(i, j, k, 0) = (gamma - 1.0) * rho(i, j, k) * (E(i, j, k) / rho(i, j, k) - ke);

            // Set::Scalar c = sqrt(gamma * p(i, j, k) / rho(i, j, k));

            // if (c > c_max) { c_max = c; }
            // if (v(i, j, k, 0) > vx_max) { vx_max = v(i, j, k, 0); }
            // if (v(i, j, k, 1) > vy_max) { vy_max = v(i, j, k, 1); }
        });
    }

    for (amrex::MFIter mfi(*eta_mf[lev], false); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        Set::Scalar step_int = 0;

        amrex::Array4<Set::Scalar> const& eta = (*eta_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const& etadot = (*etadot_mf[lev]).array(mfi);

        amrex::Array4<const Set::Scalar> const& E = (*Energy_old_mf[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& rho = (*Density_old_mf[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& M = (*Momentum_old_mf[lev]).array(mfi);

        amrex::Array4<const Set::Scalar> const& v = (*Velocity_mf[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& p = (*Pressure_mf[lev]).array(mfi);

        amrex::Array4<Set::Scalar> const& E_new = (*Energy_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const& rho_new = (*Density_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const& M_new = (*Momentum_mf[lev]).array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {

            //Advance eta
            Set::Vector grad_eta = Numeric::Gradient(eta, i, j, k, 0, DX);
            Set::Scalar grad_eta_mag = grad_eta.lpNorm<2>();

            Set::Scalar interface_pos_x = (i - num_cells_x / 2.) * DX[0];
            eta(i, j, k) = 0.5 + 0.5 * std::erf((interface_pos_x) / eps);
            etadot(i, j, k) = 0.0;

            ///////////////////////////
            /////GODUNOV ADVECTION/////
            ///////////////////////////

            // slopes in current cell
            Set::Vector drho = Numeric::Gradient(rho, i, j, k, 0, DX);
            Set::Vector dp = Numeric::Gradient(p, i, j, k, 0, DX);
            Set::Vector dvx = Numeric::Gradient(v, i, j, k, 0, DX);
            Set::Vector dvy = Numeric::Gradient(v, i, j, k, 0, DX);

            Set::Vector drho_n = Numeric::Gradient(rho, i - 1, j, k, 0, DX);
            Set::Vector dvx_n = Numeric::Gradient(v, i - 1, j, k, 0, DX);
            Set::Vector dvy_n = Numeric::Gradient(v, i - 1, j, k, 1, DX);
            Set::Vector dp_n = Numeric::Gradient(p, i - 1, j, k, 0, DX);

            Set::Scalar rho_slope_right, vx_slope_right, vy_slope_right, p_slope_right;
            Set::Scalar rho_right, vx_right, vy_right, p_right;
            Set::Scalar rho_slope_left, vx_slope_left, vy_slope_left, p_slope_left;
            Set::Scalar rho_left, vx_left, vy_left, p_left;
            std::array<Set::Scalar, 4> flux_x, flux_y;

            // compute reconstructed states at left interface along x in current cell
            // left interface: right state
            rho_slope_right = (-v(i, j, k, 0) * drho[0] - dvx[0] * rho(i, j, k)) * dt + (-v(i, j, k, 1) * drho[1] - dvy[1] * rho(i, j, k)) * dt;
            vx_slope_right = (-v(i, j, k, 0) * dvx[0] - dp[0] / rho(i, j, k)) * dt + (-v(i, j, k, 1) * dvy[1]) * dt / DX[1];
            vy_slope_right = (-v(i, j, k, 0) * dvy[0]) * dt + (-v(i, j, k, 1) * dvy[1] - dp[1] / rho(i, j, k)) * dt;
            p_slope_right = (-v(i, j, k, 0) * dp[0] - dvx[0] * gamma * p(i, j, k)) * dt + (-v(i, j, k, 1) * dp[1] - dvy[1] * gamma * p(i, j, k)) * dt;

            rho_right = rho(i, j, k) + 0.5 * rho_slope_right - drho[0] * DX[0];
            vx_right = v(i, j, k, 0) + 0.5 * vx_slope_right - dvx[0] * DX[0];
            vy_right = v(i, j, k, 1) + 0.5 * vy_slope_right - dvy[0] * DX[0];
            p_right = p(i, j, k) + 0.5 * p_slope_right - dp[0] * DX[0];

            // left interface: left state
            rho_slope_left = (-v(i - 1, j, k, 0) * drho[0] - dvx[0] * rho(i - 1, j, k)) * dt + (-v(i - 1, j, k, 1) * drho[1] - dvy[1] * rho(i - 1, j, k)) * dt;
            vx_slope_left = (-v(i - 1, j, k, 0) * dvx[0] - dp[0] / rho(i - 1, j, k)) * dt + (-v(i - 1, j, k, 1) * dvy[1]) * dt;
            vy_slope_left = (-v(i - 1, j, k, 0) * dvy[0]) * dt + (-v(i - 1, j, k, 1) * dvy[1] - dp[1] / rho(i - 1, j, k)) * dt;
            p_slope_left = (-v(i - 1, j, k, 0) * dp[0] - dvx[0] * gamma * p(i - 1, j, k)) * dt + (-v(i - 1, j, k, 1) * dp[1] - dvy[1] * gamma * p(i - 1, j, k)) * dt;

            rho_left = rho(i - 1, j, k) + 0.5 * rho_slope_left + drho_n[0] * DX[0];
            vx_left = v(i - 1, j, k, 0) + 0.5 * vx_slope_left + dvx_n[0] * DX[0];
            vy_left = v(i - 1, j, k, 1) + 0.5 * vy_slope_left + dvy_n[0] * DX[0];
            p_left = p(i - 1, j, k) + 0.5 * p_slope_left + dp_n[0] * DX[0];

            double left_state[5] = { rho_left, vx_left, vy_left, p_left, eta(i - 1, j, k) };
            double right_state[5] = { rho_right, vx_right, vy_right, p_right, eta(i, j, k) };

            flux_x = Solver::Local::Riemann::Roe(left_state, right_state, eta(i, j, k), gamma);

            //
            // left interface along y direction
            //

            // compute slopes in left neighbor
            drho_n = Numeric::Gradient(rho, i, j - 1, k, 0, DX);
            dvx_n = Numeric::Gradient(v, i, j - 1, k, 0, DX);
            dvy_n = Numeric::Gradient(v, i, j - 1, k, 1, DX);
            dp_n = Numeric::Gradient(p, i, j - 1, k, 0, DX);

            // compute reconstructed states at left interface along y in current cell
            // left interface: right state
            rho_right = rho(i, j, k) + 0.5 * rho_slope_right - drho[1] * DX[1];
            vx_right = v(i, j, k, 0) + 0.5 * vx_slope_right - dvx[1] * DX[1];
            vy_right = v(i, j, k, 1) + 0.5 * vy_slope_right - dvy[1] * DX[1];
            p_right = p(i, j, k) + 0.5 * p_slope_right - dp[1] * DX[1];

            // left interface: left state
            rho_slope_left = (-v(i, j - 1, k, 0) * drho[0] - dvx[0] * rho(i, j - 1, k)) * dt + (-v(i, j - 1, k, 1) * drho[1] - dvy[1] * rho(i, j - 1, k)) * dt;
            vx_slope_left = (-v(i, j - 1, k, 0) * dvx[0] - dp[0] / rho(i, j - 1, k)) * dt + (-v(i, j - 1, k, 1) * dvy[1]) * dt;
            vy_slope_left = (-v(i, j - 1, k, 0) * dvy[0]) * dt + (-v(i, j - 1, k, 1) * dvy[1] - dp[1] / rho(i, j - 1, k)) * dt;
            p_slope_left = (-v(i, j - 1, k, 0) * dp[0] - dvx[0] * gamma * p(i, j - 1, k)) * dt + (-v(i, j - 1, k, 1) * dp[1] - dvy[1] * gamma * p(i, j - 1, k)) * dt;

            rho_left = rho(i, j - 1, k) + 0.5 * rho_slope_left + drho_n[1] * DX[1];
            vx_left = v(i, j - 1, k, 0) + 0.5 * vx_slope_left + dvx_n[1] * DX[1];
            vy_left = v(i, j - 1, k, 1) + 0.5 * vy_slope_left + dvy_n[1] * DX[1];
            p_left = p(i, j - 1, k) + 0.5 * p_slope_left + dp_n[1] * DX[1];

            // x, y permutations -> y direction flux is computed
            std::swap(vx_left, vy_left);
            std::swap(vx_right, vy_right);

            left_state[0] = rho_left;
            left_state[1] = vx_left;
            left_state[2] = vy_left;
            left_state[3] = p_left;
            left_state[4] = eta(i, j - 1, k);

            right_state[0] = rho_right;
            right_state[1] = vx_right;
            right_state[2] = vy_right;
            right_state[3] = p_right;
            right_state[4] = eta(i, j - 1, k);

            flux_y = Solver::Local::Riemann::Roe(left_state, right_state, eta(i, j, k), gamma);

            // swap flux_y components 
            std::swap(flux_y[2], flux_y[3]);

            //Godunov fluxes
            E_new(i - 1, j, k) = E(i - 1, j, k) - flux_x[1] * dt / DX[0];
            E_new(i, j, k) = E(i, j, k) + flux_x[1] * dt / DX[0];
            E_new(i, j - 1, k) = E(i, j - 1, k) - flux_y[1] * dt / DX[1];
            E_new(i, j, k) = E(i, j, k) + flux_y[1] * dt / DX[1];

            rho_new(i - 1, j, k) = rho(i - 1, j, k) - flux_x[0] * dt / DX[0];
            rho_new(i, j, k) = rho(i, j, k) + flux_x[0] * dt / DX[0];
            rho_new(i, j - 1, k) = rho(i, j - 1, k) - flux_y[0] * dt / DX[1];
            rho_new(i, j, k) = rho(i, j, k) + flux_y[0] * dt / DX[1];

            M_new(i - 1, j, k, 0) = M(i - 1, j, k, 0) - flux_x[2] * dt / DX[0];
            M_new(i, j, k, 0) = M(i, j, k, 0) + flux_x[2] * dt / DX[0];
            M_new(i, j - 1, k, 0) = M(i, j - 1, k, 0) - flux_y[2] * dt / DX[1];
            M_new(i, j, k, 0) = M(i, j, k, 0) + flux_y[2] * dt / DX[1];

            M_new(i - 1, j, k, 1) = M(i - 1, j, k, 1) - flux_x[3] * dt / DX[0];
            M_new(i, j, k, 1) = M(i, j, k, 1) + flux_x[3] * dt / DX[0];
            M_new(i, j - 1, k, 1) = M(i, j - 1, k, 1) - flux_y[3] * dt / DX[1];
            M_new(i, j, k, 1) = M(i, j, k, 1) + flux_y[3] * dt / DX[1];

            ///////////////////////////
            //////DIFFUSE SOURCES//////
            ///////////////////////////

            std::array<Set::Scalar, 3> source;
            source[0] = mdot * grad_eta_mag;
            source[1] = Pdot_x * grad_eta_mag;
            source[2] = Pdot_y * grad_eta_mag;
            source[3] = Qdot * grad_eta_mag;

            E_new(i, j, k) += source[3] + E(i, j, k) * etadot(i, j, k);
            rho_new(i, j, k) += source[0] + rho(i, j, k) * etadot(i, j, k);
            M_new(i, j, k, 0) += source[1] + M(i, j, k, 0) * etadot(i, j, k);
            M_new(i, j, k, 1) += source[2] + M(i, j, k, 1) * etadot(i, j, k);

            //Advance total fields
            rho_new(i, j, k) = rho_solid * (1 - eta(i, j, k)) + rho_new(i, j, k) * eta(i, j, k);
            M_new(i, j, k, 0) = 0.0 * (1 - eta(i, j, k)) + M_new(i, j, k, 0) * eta(i, j, k);
            M_new(i, j, k, 1) = 0.0 * (1 - eta(i, j, k)) + M_new(i, j, k, 1) * eta(i, j, k);
            E_new(i, j, k) = E_solid * (1 - eta(i, j, k)) + E_new(i, j, k) * eta(i, j, k);

        });

        step_int += 1;
    }
}//end Advance

void Hydro::Regrid(int lev, Set::Scalar /* time */)
{
    BL_PROFILE("Integrator::Hydro::Regrid");
    if (lev < finest_level) return;

    Util::Message(INFO, "Regridding on level", lev);
}//end regrid

//void Hydro::TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, Set::Scalar time, int ngrow)
void Hydro::TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags, Set::Scalar, int)
{
    BL_PROFILE("Integrator::Flame::TagCellsForRefinement");

    const Set::Scalar* DX = geom[lev].CellSize();
    Set::Scalar dr = sqrt(AMREX_D_TERM(DX[0] * DX[0], +DX[1] * DX[1], +DX[2] * DX[2]));

    // Eta criterion for refinement
    for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<char> const& tags = a_tags.array(mfi);
        amrex::Array4<const Set::Scalar> const& eta = (*eta_mf[lev]).array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            Set::Vector grad_eta = Numeric::Gradient(eta, i, j, k, 0, DX);
            if (grad_eta.lpNorm<2>() * dr * 2 > eta_refinement_criterion) tags(i, j, k) = amrex::TagBox::SET;
        });
    }



}//end TagCells

// void Hydro::Integrate(int amrlev, Set::Scalar /*time*/, int /*step*/, const amrex::MFIter &mfi, const amrex::Box &box)
// {
//   BL_PROFILE("Hydro::Integrate");
//   const Set::Scalar *DX = geom[amrlev].CellSize();
//   Set::Scalar dv = AMREX_D_TERM(DX[0], *DX[1], *DX[2]);
//   amrex::Array4<amrex::Real> const &eta = (*eta_mf[amrlev]).array(mfi);
//   amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k){
//   volume += eta(i, j, k, 0) * dv;
//   Set::Vector grad = Numeric::Gradient(eta, i, j, k, 0, DX);
//   Set::Scalar normgrad = grad.lpNorm<2>();
//   Set::Scalar da = normgrad * dv;
//   area += da;
//   });


//  }//end Integrate

//  void Hydro::UpdateModel(int /*a_step*/)
//  {
//    for (int lev = 0; lev <= finest_level; ++lev)
//    {
//      eta_mf[lev] -> FillBoundary();
//      Density_mf[lev] -> FillBoundary();
//      Energy_mf[lev] -> FillBoundary();
//      Momentum[lev] -> FillBoundary();
//      
//      for (MFIter mfi(*model_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
//      {
//amrex::Box bx = mfi.nodaltilebox();
//amrex::Array4<model_type> const &model = model_mf[lev]->array(mfi);
//
//amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k){
//// TODO
//
//});
//
//
//     } // end For2

//     Util::RealFillBoundary(*model_mf[lev], geom[lev]);
//     amrex::MultiFab::Copy(*psi_mf[lev], *eta_mf[lev], 0, 0, 1, psi_mf[lev]-> nGrow());

//    } //end For1
//}//end update


}//end code
