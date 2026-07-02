#include "LowMach.H"

#include "AMReX_MLABecLaplacian.H"
#include "AMReX_MLMG.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_TimeIntegrator.H"
#include "Numeric/Stencil.H"

#include "Model/Gas/Thermo/Thermo.H"
#include "Model/Gas/Thermo/CpConstant.H"
#include "Model/Gas/Transport/Transport.H"
#include "Model/Gas/Transport/Mixture_Averaged.H"
#include "Model/Gas/EOS/EOS.H"
#include "Model/Gas/EOS/CPG.H"

namespace Integrator
{
namespace
{
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar lowmach_max(Set::Scalar a, Set::Scalar b)
{
    return a > b ? a : b;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar lowmach_finite_or(Set::Scalar value, Set::Scalar fallback)
{
    return (value == value && value < 1.0e300 && value > -1.0e300) ? value : fallback;
}
}

LowMach::LowMach(IO::ParmParse& pp) : LowMach()
{
    pp_queryclass(*this);
}

void
LowMach::Parse(LowMach& value, IO::ParmParse& pp)
{
    BL_PROFILE("Integrator::LowMach::Parse");

    pp.query_required("cfl", value.cfl);
    pp.query_default("cfl_v", value.cfl_v, 1.0e100);
    pp.query_default("small", value.small, 1.0e-12);
    pp.query_default("density_floor", value.density_floor, value.small);
    pp.query_default("temperature_floor", value.temperature_floor, value.small);
    pp.query_default("pressure_floor", value.pressure_floor, value.small);
    pp.query_default("pressure_scale", value.pressure_scale, 1.0);
    pp.query_default("projection.enabled", value.projection_enabled, true);
    pp.query_default("projection.tol_rel", value.projection_tol_rel, 1.0e-11);
    pp.query_default("projection.tol_abs", value.projection_tol_abs, 1.0e-12);
    pp.query_default("projection.verbose", value.projection_verbose, 0);
    pp.query_default("projection.update_pressure", value.projection_update_pressure, false);
    pp.query_default("include_viscosity", value.include_viscosity, true);
    pp.query_default("include_conduction", value.include_conduction, true);
    pp.query_default("advect_temperature", value.advect_temperature, true);

    pp.query_default("velocity_refinement_criterion", value.velocity_refinement_criterion, 1.0e100);
    pp.query_default("pressure_refinement_criterion", value.pressure_refinement_criterion, 1.0e100);
    pp.query_default("temperature_refinement_criterion", value.temperature_refinement_criterion, 1.0e100);
    pp.queryarr_default("g", value.g, Set::Vector::Zero());
    for (int face = 0; face < 2 * AMREX_SPACEDIM; ++face)
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
            value.velocity_dirichlet_value[face][d] = Numeric::Interpolator::Linear<Set::Scalar>(NAN);
    auto read_velocity_wall_value = [&](const std::string& name, int face)
    {
        std::vector<std::string> vals;
        pp.queryarr_default(name, vals, std::vector<std::string>{"nan"});
        if (vals.size() == 1)
        {
            for (int d = 0; d < AMREX_SPACEDIM; ++d) value.velocity_dirichlet_value[face][d].define(vals[0], Unit::Time(), Unit::Less());
        }
        else if (vals.size() == AMREX_SPACEDIM)
        {
            for (int d = 0; d < AMREX_SPACEDIM; ++d) value.velocity_dirichlet_value[face][d].define(vals[d], Unit::Time(), Unit::Less());
        }
    };
    read_velocity_wall_value("velocity.bc.constant.val.xlo", 0);
    read_velocity_wall_value("velocity.bc.constant.val.ylo", 1);
#if AMREX_SPACEDIM == 3
    read_velocity_wall_value("velocity.bc.constant.val.zlo", 2);
    read_velocity_wall_value("velocity.bc.constant.val.xhi", 3);
    read_velocity_wall_value("velocity.bc.constant.val.yhi", 4);
    read_velocity_wall_value("velocity.bc.constant.val.zhi", 5);
#else
    read_velocity_wall_value("velocity.bc.constant.val.xhi", 2);
    read_velocity_wall_value("velocity.bc.constant.val.yhi", 3);
#endif

    pp.queryclass<Model::Gas::Gas>("gas", value.gas);
    value.nspecies = value.gas.nspecies;

    int nghost = 2;
    pp.select_default<BC::Constant,BC::Expression>("velocity.bc", value.velocity_bc, AMREX_SPACEDIM);
    pp.select_default<BC::Constant,BC::Expression>("temperature.bc", value.temperature_bc, 1);
    pp.select_default<BC::Constant,BC::Expression>("mass_fraction.bc", value.mass_fraction_bc, value.nspecies);
    pp.select_default<BC::Constant,BC::Expression>("pressure.bc", value.pressure_bc, 1);

    pp.select_default<IC::Constant,IC::Expression>("velocity.ic", value.velocity_ic, value.geom);
    pp.select_default<IC::Constant,IC::Expression>("temperature.ic", value.temperature_ic, value.geom);
    pp.select_default<IC::Constant,IC::Expression>("mass_fraction.ic", value.mass_fraction_ic, value.geom);
    pp.select_default<IC::Constant,IC::Expression>("pressure.ic", value.pressure_ic, value.geom);

    value.RegisterNewFab(value.velocity_mf,          value.velocity_bc,      AMREX_SPACEDIM, nghost, "velocity",          true,  true, {"x","y"});
    value.RegisterNewFab(value.velocity_old_mf,      value.velocity_bc,      AMREX_SPACEDIM, nghost, "velocity_old",      false, true, {"x","y"});
    value.RegisterNewFab(value.temperature_mf,       value.temperature_bc,   1,              nghost, "temperature",       true,  true);
    value.RegisterNewFab(value.temperature_old_mf,   value.temperature_bc,   1,              nghost, "temperature_old",   false, true);
    value.RegisterNewFab(value.mass_fraction_mf,     value.mass_fraction_bc, value.nspecies, nghost, "mass_fraction",     true,  true);
    value.RegisterNewFab(value.mass_fraction_old_mf, value.mass_fraction_bc, value.nspecies, nghost, "mass_fraction_old", false, true);

    value.RegisterNewFab(value.density_mf,             &value.bc_nothing, 1,              0,      "density",             true,  false);
    value.RegisterNewFab(value.momentum_mf,            &value.bc_nothing, AMREX_SPACEDIM, 0,      "momentum",            true,  false, {"x","y"});
    value.RegisterNewFab(value.pressure_mf,            value.pressure_bc, 1,              nghost, "pressure",            true,  false);
    value.RegisterNewFab(value.pressure_correction_mf, &value.bc_nothing, 1,              nghost, "pressure_correction", false, false);
    value.RegisterNewFab(value.projection_rhs_mf,      &value.bc_nothing, 1,              0,      "projection_rhs",      false, false);
    value.RegisterNewFab(value.energy_mf,              &value.bc_nothing, 1,              0,      "energy",              true,  false);
    value.RegisterNewFab(value.mole_fraction_mf, &value.bc_nothing, value.nspecies, 0,      "mole_fraction", true, false);
    value.RegisterNewFab(value.vorticity_mf,     &value.bc_nothing, 1,              0,      "vorticity",     true, false);

    bool allow_unused;
    pp.query_default("allow_unused", allow_unused, false);
    if (!allow_unused && pp.AnyUnusedInputs(true, false))
    {
        Util::Warning(INFO, "The following inputs were specified but not used:");
        pp.AllUnusedInputs();
        Util::Exception(INFO, "Aborting. Specify 'allow_unused=True` to ignore this error.");
    }
}

void
LowMach::FillStateBoundaries(int lev, amrex::MultiFab& u_mf, amrex::MultiFab& T_mf, amrex::MultiFab& Y_mf, Set::Scalar time)
{
    auto preserve_valid = [](amrex::MultiFab& mf, auto&& fill)
    {
        amrex::MultiFab valid(mf.boxArray(), mf.DistributionMap(), mf.nComp(), 0);
        amrex::MultiFab::Copy(valid, mf, 0, 0, mf.nComp(), 0);
        fill();
        amrex::MultiFab::Copy(mf, valid, 0, 0, mf.nComp(), 0);
    };

    preserve_valid(u_mf, [&]()
    {
        velocity_bc->FillBoundary(u_mf, 0, AMREX_SPACEDIM, time, 0);
        u_mf.FillBoundary(geom[lev].periodicity());
    });
    ApplyVelocityDirichletCells(lev, u_mf, time);

    preserve_valid(T_mf, [&]()
    {
        temperature_bc->FillBoundary(T_mf, 0, 1, time, 0);
        T_mf.FillBoundary(geom[lev].periodicity());
    });

    preserve_valid(Y_mf, [&]()
    {
        mass_fraction_bc->FillBoundary(Y_mf, 0, nspecies, time, 0);
        Y_mf.FillBoundary(geom[lev].periodicity());
    });

    FillPressureBoundary(lev, time);
}



void
LowMach::ApplyVelocityDirichletCells(int lev, amrex::MultiFab& u_mf, Set::Scalar time)
{
    const amrex::Box domain = geom[lev].Domain();
    const amrex::Dim3 lo = amrex::lbound(domain);
    const amrex::Dim3 hi = amrex::ubound(domain);
    const amrex::BCRec bc = velocity_bc->GetBCRec();
    const bool xlo_dirichlet = !geom[lev].isPeriodic(0) && BC::BCUtil::IsDirichlet(bc.lo(0));
    const bool xhi_dirichlet = !geom[lev].isPeriodic(0) && BC::BCUtil::IsDirichlet(bc.hi(0));
    const bool ylo_dirichlet = !geom[lev].isPeriodic(1) && BC::BCUtil::IsDirichlet(bc.lo(1));
    const bool yhi_dirichlet = !geom[lev].isPeriodic(1) && BC::BCUtil::IsDirichlet(bc.hi(1));
#if AMREX_SPACEDIM == 3
    const bool zlo_dirichlet = !geom[lev].isPeriodic(2) && BC::BCUtil::IsDirichlet(bc.lo(2));
    const bool zhi_dirichlet = !geom[lev].isPeriodic(2) && BC::BCUtil::IsDirichlet(bc.hi(2));
#endif

    Set::Vector xlo_val = Set::Vector::Zero();
    Set::Vector ylo_val = Set::Vector::Zero();
#if AMREX_SPACEDIM == 3
    Set::Vector zlo_val = Set::Vector::Zero();
    Set::Vector xhi_val = Set::Vector::Zero();
    Set::Vector yhi_val = Set::Vector::Zero();
    Set::Vector zhi_val = Set::Vector::Zero();
#else
    Set::Vector xhi_val = Set::Vector::Zero();
    Set::Vector yhi_val = Set::Vector::Zero();
#endif
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        xlo_val(d) = velocity_dirichlet_value[0][d](time);
        ylo_val(d) = velocity_dirichlet_value[1][d](time);
#if AMREX_SPACEDIM == 3
        zlo_val(d) = velocity_dirichlet_value[2][d](time);
        xhi_val(d) = velocity_dirichlet_value[3][d](time);
        yhi_val(d) = velocity_dirichlet_value[4][d](time);
        zhi_val(d) = velocity_dirichlet_value[5][d](time);
#else
        xhi_val(d) = velocity_dirichlet_value[2][d](time);
        yhi_val(d) = velocity_dirichlet_value[3][d](time);
#endif
    }

    for (amrex::MFIter mfi(u_mf, true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        Set::Patch<Set::Scalar> u = u_mf.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
            {
                if (ylo_dirichlet && j == lo.y && ylo_val(d) == ylo_val(d)) u(i,j,k,d) = ylo_val(d);
                if (yhi_dirichlet && j == hi.y && yhi_val(d) == yhi_val(d)) u(i,j,k,d) = yhi_val(d);
                if (xlo_dirichlet && i == lo.x && xlo_val(d) == xlo_val(d)) u(i,j,k,d) = xlo_val(d);
                if (xhi_dirichlet && i == hi.x && xhi_val(d) == xhi_val(d)) u(i,j,k,d) = xhi_val(d);
#if AMREX_SPACEDIM == 3
                if (zlo_dirichlet && k == lo.z && zlo_val(d) == zlo_val(d)) u(i,j,k,d) = zlo_val(d);
                if (zhi_dirichlet && k == hi.z && zhi_val(d) == zhi_val(d)) u(i,j,k,d) = zhi_val(d);
#endif
            }
        });
    }
}


void
LowMach::FillPressureBoundary(int lev, Set::Scalar time)
{
    amrex::MultiFab valid(pressure_mf[lev]->boxArray(), pressure_mf[lev]->DistributionMap(), pressure_mf[lev]->nComp(), 0);
    amrex::MultiFab::Copy(valid, *pressure_mf[lev], 0, 0, pressure_mf[lev]->nComp(), 0);
    pressure_bc->FillBoundary(*pressure_mf[lev], 0, 1, time, 0);
    pressure_mf[lev]->FillBoundary(geom[lev].periodicity());
    amrex::MultiFab::Copy(*pressure_mf[lev], valid, 0, 0, pressure_mf[lev]->nComp(), 0);
}

void
LowMach::EnforceStateBounds(amrex::MultiFab& u_mf, amrex::MultiFab& T_mf, amrex::MultiFab& Y_mf)
{
    for (amrex::MFIter mfi(Y_mf, false); mfi.isValid(); ++mfi)
    {
        const amrex::Box bx = mfi.fabbox();
        Set::Patch<Set::Scalar> u = u_mf.array(mfi);
        Set::Patch<Set::Scalar> T = T_mf.array(mfi);
        Set::Patch<Set::Scalar> Y = Y_mf.array(mfi);
        const int nsp = nspecies;
        const Set::Scalar T_floor = temperature_floor;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Scalar sumY = 0.0;
            for (int n = 0; n < nsp; ++n)
            {
                if (!(Y(i,j,k,n) > 0.0)) Y(i,j,k,n) = 0.0;
                sumY += Y(i,j,k,n);
            }
            if (!(sumY > 0.0))
            {
                Y(i,j,k,0) = 1.0;
                for (int n = 1; n < nsp; ++n) Y(i,j,k,n) = 0.0;
            }
            else
            {
                for (int n = 0; n < nsp; ++n) Y(i,j,k,n) /= sumY;
            }

            if (!(T(i,j,k) > T_floor)) T(i,j,k) = T_floor;
            u(i,j,k,0) = lowmach_finite_or(u(i,j,k,0), 0.0);
            u(i,j,k,1) = lowmach_finite_or(u(i,j,k,1), 0.0);
#if AMREX_SPACEDIM == 3
            u(i,j,k,2) = lowmach_finite_or(u(i,j,k,2), 0.0);
#endif
        });
    }
}

void
LowMach::UpdateDerived(int lev, const amrex::MultiFab& u_mf, const amrex::MultiFab& T_mf, const amrex::MultiFab& Y_mf)
{
    const Set::Scalar* DX = geom[lev].CellSize();
    amrex::Box domain = geom[lev].Domain();

    density_mf[lev]->setVal(0.0);
    momentum_mf[lev]->setVal(0.0);
    energy_mf[lev]->setVal(0.0);
    mole_fraction_mf[lev]->setVal(0.0);
    vorticity_mf[lev]->setVal(0.0);

    for (amrex::MFIter mfi(Y_mf, true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        Set::Patch<const Set::Scalar> u = u_mf.array(mfi);
        Set::Patch<const Set::Scalar> T = T_mf.array(mfi);
        Set::Patch<const Set::Scalar> Y = Y_mf.array(mfi);
        Set::Patch<const Set::Scalar> pressure = pressure_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> rho = density_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> M = momentum_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> E = energy_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> X = mole_fraction_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> omega = vorticity_mf.Patch(lev,mfi);
        const Set::Scalar rho_floor = density_floor;
        const Set::Scalar p_floor = pressure_floor;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Scalar moles = 0.0;
            for (int n = 0; n < gas.nspecies; ++n) moles += Y(i,j,k,n) / gas.MW[n];
            if (!(moles > 0.0)) moles = 1.0 / gas.MW[0];
            for (int n = 0; n < gas.nspecies; ++n) X(i,j,k,n) = (Y(i,j,k,n) / gas.MW[n]) / moles;

            Set::Scalar p = lowmach_max(pressure(i,j,k), p_floor);
            Set::Scalar density = p / (gas.R(X, i, j, k) * T(i,j,k));
            density = lowmach_max(density, rho_floor);
            rho(i,j,k) = density;

            M(i,j,k,0) = density * u(i,j,k,0);
            M(i,j,k,1) = density * u(i,j,k,1);
#if AMREX_SPACEDIM == 3
            M(i,j,k,2) = density * u(i,j,k,2);
#endif
            E(i,j,k) = gas.ComputeE(density, M(i,j,k,0), M(i,j,k,1), T(i,j,k), X, i, j, k);

            auto sten = Numeric::GetStencil(i, j, k, domain);
            Set::Matrix grad_u = Numeric::Gradient(u, i, j, k, DX, sten);
            omega(i,j,k) = lowmach_finite_or(grad_u(1,0) - grad_u(0,1), 0.0);
        });
    }

    density_mf[lev]->FillBoundary(geom[lev].periodicity());
    momentum_mf[lev]->FillBoundary(geom[lev].periodicity());
    energy_mf[lev]->FillBoundary(geom[lev].periodicity());
    mole_fraction_mf[lev]->FillBoundary(geom[lev].periodicity());
    vorticity_mf[lev]->FillBoundary(geom[lev].periodicity());
}


void
LowMach::ProjectVelocity(int lev, Set::Scalar time, Set::Scalar dt)
{
    BL_PROFILE("Integrator::LowMach::ProjectVelocity");
    if (!projection_enabled || !(dt > 0.0)) return;

    amrex::MultiFab& phi_mf = *pressure_correction_mf[lev];
    amrex::MultiFab& rhs_mf = *projection_rhs_mf[lev];
    amrex::MultiFab& u_mf = *velocity_mf[lev];
    amrex::MultiFab& rho_mf = *density_mf[lev];
    amrex::MultiFab& p_mf = *pressure_mf[lev];

    const Set::Scalar* DX = geom[lev].CellSize();
    amrex::Box domain = geom[lev].Domain();
    rhs_mf.setVal(0.0);
    phi_mf.setVal(0.0);

    for (amrex::MFIter mfi(u_mf, true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        Set::Patch<const Set::Scalar> u = u_mf.array(mfi);
        Set::Patch<Set::Scalar> rhs = rhs_mf.array(mfi);
        const Set::Scalar inv_dt = 1.0 / dt;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            auto sten = Numeric::GetStencil(i, j, k, domain);
            Set::Matrix grad_u = Numeric::Gradient(u, i, j, k, DX, sten);
            Set::Scalar div_u = 0.0;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) div_u += grad_u(d,d);
            rhs(i,j,k) = div_u * inv_dt;
        });
    }

    const Set::Scalar rhs_mean = rhs_mf.sum(0, false) / rhs_mf.boxArray().d_numPts();
    for (amrex::MFIter mfi(rhs_mf, true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        Set::Patch<Set::Scalar> rhs = rhs_mf.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            rhs(i,j,k) -= rhs_mean;
        });
    }

    amrex::MultiFab beta_cc(u_mf.boxArray(), u_mf.DistributionMap(), 1, 1);
    beta_cc.setVal(0.0);
    for (amrex::MFIter mfi(beta_cc, true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        Set::Patch<const Set::Scalar> rho = rho_mf.array(mfi);
        Set::Patch<Set::Scalar> beta = beta_cc.array(mfi);
        const Set::Scalar rho_floor = density_floor;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            beta(i,j,k) = 1.0 / lowmach_max(rho(i,j,k), rho_floor);
        });
    }
    beta_cc.FillBoundary(geom[lev].periodicity());

    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> beta_face;
    amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> beta_face_ptr;
    amrex::Array<amrex::MultiFab const*, AMREX_SPACEDIM> beta_face_const_ptr;
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        amrex::BoxArray face_ba = u_mf.boxArray();
        face_ba.surroundingNodes(d);
        beta_face[d].define(face_ba, u_mf.DistributionMap(), 1, 0);
        beta_face_ptr[d] = &beta_face[d];
        beta_face_const_ptr[d] = &beta_face[d];
    }
    amrex::average_cellcenter_to_face(beta_face_ptr, beta_cc, geom[lev], 1, true, 0);

    amrex::LPInfo info;
    amrex::MLABecLaplacian mlabec({geom[lev]}, {u_mf.boxArray()}, {u_mf.DistributionMap()}, info);
    amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> lobc;
    amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> hibc;
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        lobc[d] = geom[lev].isPeriodic(d) ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Neumann;
        hibc[d] = geom[lev].isPeriodic(d) ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Neumann;
    }
    mlabec.setDomainBC(lobc, hibc);
    mlabec.setLevelBC(0, nullptr);
    mlabec.setScalars(0.0, -1.0);
    mlabec.setACoeffs(0, 0.0);
    mlabec.setBCoeffs(0, beta_face_const_ptr);

    amrex::MLMG mlmg(mlabec);
    mlmg.setVerbose(projection_verbose);
    mlmg.solve({&phi_mf}, {&rhs_mf}, projection_tol_rel, projection_tol_abs);
    const Set::Scalar phi_mean = phi_mf.sum(0, false) / phi_mf.boxArray().d_numPts();
    for (amrex::MFIter mfi(phi_mf, true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        Set::Patch<Set::Scalar> phi = phi_mf.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            phi(i,j,k) -= phi_mean;
        });
    }
    phi_mf.FillBoundary(geom[lev].periodicity());

    for (amrex::MFIter mfi(u_mf, true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        Set::Patch<Set::Scalar> u = u_mf.array(mfi);
        Set::Patch<Set::Scalar> p = p_mf.array(mfi);
        Set::Patch<const Set::Scalar> phi = phi_mf.array(mfi);
        Set::Patch<const Set::Scalar> rho = rho_mf.array(mfi);
        const Set::Scalar rho_floor = density_floor;
        const Set::Scalar p_floor = pressure_floor;
        const Set::Scalar p_scale = pressure_scale;
        const Set::Scalar p_scale_inv = (p_scale == p_scale && std::abs(p_scale) > 0.0) ? 1.0 / p_scale : 1.0;
        const bool update_pressure = projection_update_pressure;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            auto sten = Numeric::GetStencil(i, j, k, domain);
            Set::Vector grad_phi = Numeric::Gradient(phi, i, j, k, 0, DX, sten);
            Set::Scalar beta = 1.0 / lowmach_max(rho(i,j,k), rho_floor);
            for (int d = 0; d < AMREX_SPACEDIM; ++d) u(i,j,k,d) -= dt * beta * grad_phi(d);
            if (update_pressure) p(i,j,k) = lowmach_max(p(i,j,k) + p_scale_inv * phi(i,j,k), p_floor);
        });
    }

    velocity_bc->FillBoundary(u_mf, 0, AMREX_SPACEDIM, time, 0);
    u_mf.FillBoundary(geom[lev].periodicity());
    ApplyVelocityDirichletCells(lev, u_mf, time);
    FillPressureBoundary(lev, time);
}

void
LowMach::Initialize(int lev)
{
    BL_PROFILE("Integrator::LowMach::Initialize");

    velocity_ic->Initialize(lev, velocity_mf, 0.0);
    velocity_ic->Initialize(lev, velocity_old_mf, 0.0);
    temperature_ic->Initialize(lev, temperature_mf, 0.0);
    temperature_ic->Initialize(lev, temperature_old_mf, 0.0);
    mass_fraction_ic->Initialize(lev, mass_fraction_mf, 0.0);
    mass_fraction_ic->Initialize(lev, mass_fraction_old_mf, 0.0);
    pressure_ic->Initialize(lev, pressure_mf, 0.0);
    pressure_correction_mf[lev]->setVal(0.0);
    projection_rhs_mf[lev]->setVal(0.0);

    EnforceStateBounds(*velocity_mf[lev], *temperature_mf[lev], *mass_fraction_mf[lev]);
    EnforceStateBounds(*velocity_old_mf[lev], *temperature_old_mf[lev], *mass_fraction_old_mf[lev]);
    FillStateBoundaries(lev, *velocity_mf[lev], *temperature_mf[lev], *mass_fraction_mf[lev], 0.0);
    FillStateBoundaries(lev, *velocity_old_mf[lev], *temperature_old_mf[lev], *mass_fraction_old_mf[lev], 0.0);
    UpdateDerived(lev, *velocity_mf[lev], *temperature_mf[lev], *mass_fraction_mf[lev]);
}

void
LowMach::RHS(int lev, Set::Scalar /*time*/,
             amrex::MultiFab& u_rhs_mf,
             amrex::MultiFab& T_rhs_mf,
             amrex::MultiFab& Y_rhs_mf,
             const amrex::MultiFab& u_mf,
             const amrex::MultiFab& T_mf,
             const amrex::MultiFab& Y_mf)
{
    UpdateDerived(lev, u_mf, T_mf, Y_mf);

    const Set::Scalar* DX = geom[lev].CellSize();
    amrex::Box domain = geom[lev].Domain();

    for (amrex::MFIter mfi(u_mf, false); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        Set::Patch<const Set::Scalar> u = u_mf.array(mfi);
        Set::Patch<const Set::Scalar> T = T_mf.array(mfi);
        Set::Patch<const Set::Scalar> Y = Y_mf.array(mfi);
        Set::Patch<const Set::Scalar> rho = density_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> pressure = pressure_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> X = mole_fraction_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> u_rhs = u_rhs_mf.array(mfi);
        Set::Patch<Set::Scalar> T_rhs = T_rhs_mf.array(mfi);
        Set::Patch<Set::Scalar> Y_rhs = Y_rhs_mf.array(mfi);
        const int nsp = nspecies;
        const bool viscous = include_viscosity;
        const bool conductive = include_conduction;
        const bool advect_T = advect_temperature;
        const Set::Vector gravity = g;
        const Set::Scalar p_scale = pressure_scale;
        const Set::Scalar rho_floor = density_floor;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            auto sten = Numeric::GetStencil(i, j, k, domain);
            Set::Vector vel = Set::Vector::Zero();
            vel(0) = u(i,j,k,0);
            vel(1) = u(i,j,k,1);
#if AMREX_SPACEDIM == 3
            vel(2) = u(i,j,k,2);
#endif
            Set::Scalar density = lowmach_max(rho(i,j,k), rho_floor);
            Set::Matrix grad_u = Numeric::Gradient(u, i, j, k, DX, sten);
            Set::Vector grad_p = Numeric::Gradient(pressure, i, j, k, 0, DX, sten);
            Set::Scalar mu = viscous ? gas.dynamic_viscosity(T(i,j,k), X, i, j, k) : 0.0;

            for (int d = 0; d < AMREX_SPACEDIM; ++d)
            {
                u_rhs(i,j,k,d) = -vel.dot(grad_u.row(d)) - p_scale * grad_p(d) / density + gravity(d);
                if (viscous && mu == mu)
                    u_rhs(i,j,k,d) += (mu / density) * Numeric::Laplacian(u, i, j, k, d, DX);
            }

            T_rhs(i,j,k) = 0.0;
            if (advect_T)
            {
                Set::Vector grad_T = Numeric::Gradient(T, i, j, k, 0, DX, sten);
                T_rhs(i,j,k) -= vel.dot(grad_T);
            }
            if (conductive)
            {
                Set::Scalar kappa = gas.thermal_conductivity(T(i,j,k), X, i, j, k);
                Set::Scalar cp = gas.cp_mass(T(i,j,k), X, i, j, k);
                if (kappa == kappa && cp == cp && cp > 0.0)
                    T_rhs(i,j,k) += kappa / (density * cp) * Numeric::Laplacian(T, i, j, k, 0, DX);
            }

            for (int n = 0; n < nsp; ++n)
            {
                Set::Vector grad_Y = Numeric::Gradient(Y, i, j, k, n, DX, sten);
                Y_rhs(i,j,k,n) = -vel.dot(grad_Y);
            }
        });
    }
}

void
LowMach::Advance(int lev, Set::Scalar time, Set::Scalar dt)
{
    std::swap(velocity_old_mf[lev], velocity_mf[lev]);
    std::swap(temperature_old_mf[lev], temperature_mf[lev]);
    std::swap(mass_fraction_old_mf[lev], mass_fraction_mf[lev]);

    amrex::Vector<amrex::MultiFab> solution_new;
    solution_new.emplace_back(*velocity_mf[lev].get(), amrex::MakeType::make_alias, 0, AMREX_SPACEDIM);
    solution_new.emplace_back(*temperature_mf[lev].get(), amrex::MakeType::make_alias, 0, 1);
    solution_new.emplace_back(*mass_fraction_mf[lev].get(), amrex::MakeType::make_alias, 0, nspecies);

    amrex::Vector<amrex::MultiFab> solution_old;
    solution_old.emplace_back(*velocity_old_mf[lev].get(), amrex::MakeType::make_alias, 0, AMREX_SPACEDIM);
    solution_old.emplace_back(*temperature_old_mf[lev].get(), amrex::MakeType::make_alias, 0, 1);
    solution_old.emplace_back(*mass_fraction_old_mf[lev].get(), amrex::MakeType::make_alias, 0, nspecies);

    amrex::TimeIntegrator timeintegrator(solution_new, time);
    timeintegrator.set_rhs([&](amrex::Vector<amrex::MultiFab>& rhs_mf,
                               amrex::Vector<amrex::MultiFab>& state_mf,
                               const Set::Scalar rhs_time)
    {
        RHS(lev, rhs_time, rhs_mf[0], rhs_mf[1], rhs_mf[2], state_mf[0], state_mf[1], state_mf[2]);
    });

    timeintegrator.set_post_stage_action([&](amrex::Vector<amrex::MultiFab>& stage_mf, Set::Scalar stage_time)
    {
        EnforceStateBounds(stage_mf[0], stage_mf[1], stage_mf[2]);
        FillStateBoundaries(lev, stage_mf[0], stage_mf[1], stage_mf[2], stage_time);
        UpdateDerived(lev, stage_mf[0], stage_mf[1], stage_mf[2]);
    });

    timeintegrator.advance(solution_old, solution_new, time, dt);
    EnforceStateBounds(*velocity_mf[lev], *temperature_mf[lev], *mass_fraction_mf[lev]);
    FillStateBoundaries(lev, *velocity_mf[lev], *temperature_mf[lev], *mass_fraction_mf[lev], time + dt);
    UpdateDerived(lev, *velocity_mf[lev], *temperature_mf[lev], *mass_fraction_mf[lev]);
    ProjectVelocity(lev, time + dt, dt);
    EnforceStateBounds(*velocity_mf[lev], *temperature_mf[lev], *mass_fraction_mf[lev]);
    FillStateBoundaries(lev, *velocity_mf[lev], *temperature_mf[lev], *mass_fraction_mf[lev], time + dt);
    UpdateDerived(lev, *velocity_mf[lev], *temperature_mf[lev], *mass_fraction_mf[lev]);
}

void
LowMach::TimeStepBegin(Set::Scalar /*time*/, int /*iter*/)
{
    if (!dynamictimestep.on) return;

    Set::Scalar vmax = 0.0;
    Set::Scalar viscmax = 0.0;
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        const Set::Scalar* DX = geom[lev].CellSize();
        Set::Scalar dxmin = std::min(DX[0], DX[1]);

        for (amrex::MFIter mfi(*velocity_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            Set::Patch<const Set::Scalar> u = velocity_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> T = temperature_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> rho = density_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> X = mole_fraction_mf.Patch(lev,mfi);
            const Set::Scalar rho_floor = density_floor;
            const bool viscous = include_viscosity;

            amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpMax> reduce_op;
            amrex::ReduceData<Set::Scalar, Set::Scalar> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;
            reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple
            {
                Set::Scalar speed = std::sqrt(u(i,j,k,0)*u(i,j,k,0) + u(i,j,k,1)*u(i,j,k,1));
                Set::Scalar nu = 0.0;
                if (viscous)
                {
                    Set::Scalar mu = gas.dynamic_viscosity(T(i,j,k), X, i, j, k);
                    if (mu == mu) nu = mu / lowmach_max(rho(i,j,k), rho_floor);
                }
                return {speed, nu / (dxmin * dxmin)};
            });
            ReduceTuple hv = reduce_data.value();
            vmax = std::max(vmax, amrex::get<0>(hv));
            viscmax = std::max(viscmax, amrex::get<1>(hv));
        }
    }
    amrex::ParallelDescriptor::ReduceRealMax(vmax);
    amrex::ParallelDescriptor::ReduceRealMax(viscmax);

    Set::Scalar adv_dt = cfl_v;
    if (vmax > 0.0)
    {
        const Set::Scalar* DX = geom[0].CellSize();
        adv_dt = cfl * std::min(DX[0], DX[1]) / vmax;
    }
    Set::Scalar visc_dt = viscmax > 0.0 ? 0.5 * cfl / viscmax : cfl_v;
    DynamicTimestep_SyncTimeStep(0, std::min(adv_dt, visc_dt));
}

void
LowMach::TimeStepComplete(Set::Scalar /*time*/, int /*iter*/)
{
    if (dynamictimestep.on) DynamicTimestep_Update();
}

void
LowMach::TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/)
{
    const Set::Scalar* DX = geom[lev].CellSize();
    Set::Scalar dr = std::sqrt(DX[0] * DX[0] + DX[1] * DX[1]);

    for (amrex::MFIter mfi(*temperature_mf[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<char> const& tag = tags.array(mfi);
        Set::Patch<const Set::Scalar> u = velocity_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> pressure = pressure_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> T = temperature_mf.Patch(lev,mfi);
        const Set::Scalar vcrit = velocity_refinement_criterion;
        const Set::Scalar pcrit = pressure_refinement_criterion;
        const Set::Scalar Tcrit = temperature_refinement_criterion;
        amrex::Box domain = geom[lev].Domain();

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            auto sten = Numeric::GetStencil(i, j, k, domain);
            Set::Matrix grad_u = Numeric::Gradient(u, i, j, k, DX, sten);
            Set::Vector grad_p = Numeric::Gradient(pressure, i, j, k, 0, DX, sten);
            Set::Vector grad_T = Numeric::Gradient(T, i, j, k, 0, DX, sten);
            if (grad_u.norm() * dr > vcrit || grad_p.lpNorm<2>() * dr > pcrit || grad_T.lpNorm<2>() * dr > Tcrit)
                tag(i,j,k) = amrex::TagBox::SET;
        });
    }
}

void
LowMach::Regrid(int lev, Set::Scalar /*time*/)
{
    if (lev < finest_level) return;
    for (int ilev = 0; ilev <= finest_level; ++ilev)
    {
        EnforceStateBounds(*velocity_mf[ilev], *temperature_mf[ilev], *mass_fraction_mf[ilev]);
        FillStateBoundaries(ilev, *velocity_mf[ilev], *temperature_mf[ilev], *mass_fraction_mf[ilev], 0.0);
        UpdateDerived(ilev, *velocity_mf[ilev], *temperature_mf[ilev], *mass_fraction_mf[ilev]);
    }
}
}
