#include "LowMach.H"

#include "AMReX_MLABecLaplacian.H"
#include "AMReX_MLMG.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_TimeIntegrator.H"
#include "Numeric/Stencil.H"
#include "Operator/Dynamic/ElasticLowMach.H"

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
Set::Scalar lowmach_abs(Set::Scalar a)
{
    return a < 0.0 ? -a : a;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar lowmach_min(Set::Scalar a, Set::Scalar b)
{
    return a < b ? a : b;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar lowmach_minmod(Set::Scalar a, Set::Scalar b)
{
    if (a * b <= 0.0) return 0.0;
    Set::Scalar mag = lowmach_min(lowmach_abs(a), lowmach_abs(b));
    return a < 0.0 ? -mag : mag;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar lowmach_mc_slope(Set::Scalar lo, Set::Scalar center, Set::Scalar hi)
{
    Set::Scalar dl = center - lo;
    Set::Scalar dr = hi - center;
    return lowmach_minmod(0.5 * (dl + dr), lowmach_minmod(2.0 * dl, 2.0 * dr));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar lowmach_clamp(Set::Scalar value, Set::Scalar lo, Set::Scalar hi)
{
    return lowmach_min(lowmach_max(value, lowmach_min(lo, hi)), lowmach_max(lo, hi));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool lowmach_is_finite(Set::Scalar value)
{
    return value == value && value < 1.0e300 && value > -1.0e300;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar lowmach_finite_or(Set::Scalar value, Set::Scalar fallback)
{
    return lowmach_is_finite(value) ? value : fallback;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar lowmach_smootherstep(Set::Scalar value)
{
    Set::Scalar x = lowmach_clamp(value, 0.0, 1.0);
    return x * x * x * (x * (x * 6.0 - 15.0) + 10.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Matrix lowmach_regularize_stress(Set::Matrix sigma, Set::Scalar cap)
{
    Set::Scalar norm2 = 0.0;
    for (int a = 0; a < AMREX_SPACEDIM; ++a)
        for (int b = 0; b < AMREX_SPACEDIM; ++b)
        {
            if (!lowmach_is_finite(sigma(a,b))) return Set::Matrix::Zero();
            norm2 += sigma(a,b) * sigma(a,b);
        }

    if (cap > 0.0 && lowmach_is_finite(cap) && norm2 > cap * cap)
    {
        Set::Scalar norm = std::sqrt(norm2);
        if (norm > 0.0 && lowmach_is_finite(norm)) sigma *= cap / norm;
        else return Set::Matrix::Zero();
    }
    return sigma;
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
    pp.query_default("projection.amr_enabled", value.projection_amr_enabled, false);
    pp.query_default("projection.tol_rel", value.projection_tol_rel, 1.0e-11);
    pp.query_default("projection.tol_abs", value.projection_tol_abs, 1.0e-12);
    pp.query_default("projection.verbose", value.projection_verbose, 0);
    pp.query_default("projection.update_pressure", value.projection_update_pressure, false);
    pp.query_default("include_viscosity", value.include_viscosity, true);
    pp.query_default("include_conduction", value.include_conduction, true);
    pp.query_default("advect_temperature", value.advect_temperature, true);
    pp.query_default("eta.initial_value", value.eta_initial_value, 0.0);
    pp.query_default("reference_map.eta_cutoff", value.reference_map_eta_cutoff, 0.5);
    pp.query_default("reference_map.eta_core", value.reference_map_eta_core, value.reference_map_eta_cutoff);
    pp.query_default("reference_map.extrapolation_sweeps", value.reference_map_extrapolation_sweeps, 4);
    pp.query_default("reference_map.smoothing_sweeps", value.reference_map_smoothing_sweeps, 2);

    std::string solid_model_type;
    pp.query_default("solid.model.type", solid_model_type, "none");
    pp.query_default("solid.model.eta_threshold", value.finite_solid_eta_threshold, 0.5);
    pp.query_default("solid.model.J_floor", value.finite_solid_J_floor, 1.0e-6);
    pp.query_default("solid.model.viscosity", value.finite_solid_viscosity, 0.0);
    pp.query_default("solid.model.bulk_viscosity", value.finite_solid_bulk_viscosity, 0.0);
    pp.query_default("solid.model.stress_cap", value.finite_solid_stress_cap, 1.0e100);
    pp.query_default("solid.model.stress_rhs_sign", value.finite_solid_rhs_sign, 0.0);
    pp.query_default("solid.model.implicit", value.finite_solid_implicit, true);
    pp.query_default("solid.model.implicit_coeff_scale", value.finite_solid_implicit_coeff_scale, 1.0e-4);
    pp.query_default("solid.model.implicit_preserve_rigid_modes", value.finite_solid_implicit_preserve_rigid_modes, true);
    pp.query_default("solid.model.implicit_tol_rel", value.finite_solid_implicit_tol_rel, 1.0e-8);
    pp.query_default("solid.model.implicit_tol_abs", value.finite_solid_implicit_tol_abs, 1.0e-10);
    pp.query_default("solid.model.implicit_verbose", value.finite_solid_implicit_verbose, 0);
    if (solid_model_type == "none")
    {
        value.finite_solid_enabled = false;
    }
    else if (solid_model_type == "finite.neohookean")
    {
        value.finite_solid_enabled = true;
        pp.queryclass<Model::Solid::Finite::NeoHookean>("solid.model.finite.neohookean", value.finite_solid_model);
    }
    else
    {
        Util::Exception(INFO, solid_model_type, " is not a valid LowMach solid.model.type");
    }

    pp.query_default("velocity_refinement_criterion", value.velocity_refinement_criterion, 1.0e100);
    pp.query_default("pressure_refinement_criterion", value.pressure_refinement_criterion, 1.0e100);
    pp.query_default("temperature_refinement_criterion", value.temperature_refinement_criterion, 1.0e100);
    pp.query_default("eta_refinement_criterion", value.eta_refinement_criterion, 1.0e100);
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
    if (pp.contains("eta.bc.type")) pp.select_default<BC::Constant,BC::Expression>("eta.bc", value.eta_bc, 1);
    else value.eta_bc = new BC::Constant(BC::Constant::ZeroNeumann(1));
    if (pp.contains("xi.bc.type")) pp.select_default<BC::Constant,BC::Expression>("xi.bc", value.xi_bc, AMREX_SPACEDIM);
    else value.xi_bc = new BC::Constant(BC::Constant::ZeroNeumann(AMREX_SPACEDIM));

    pp.select_default<IC::Constant,IC::Expression>("velocity.ic", value.velocity_ic, value.geom);
    pp.select_default<IC::Constant,IC::Expression>("temperature.ic", value.temperature_ic, value.geom);
    pp.select_default<IC::Constant,IC::Expression>("mass_fraction.ic", value.mass_fraction_ic, value.geom);
    pp.select_default<IC::Constant,IC::Expression>("pressure.ic", value.pressure_ic, value.geom);
    pp.select_default<IC::Constant,IC::Expression>("eta.ic", value.eta_ic, value.geom);

    value.RegisterNewFab(value.velocity_mf,          value.velocity_bc,      AMREX_SPACEDIM, nghost, "velocity",          true,  true, {"x","y"});
    value.RegisterNewFab(value.velocity_old_mf,      value.velocity_bc,      AMREX_SPACEDIM, nghost, "velocity_old",      false, true, {"x","y"});
    value.RegisterNewFab(value.temperature_mf,       value.temperature_bc,   1,              nghost, "temperature",       true,  true);
    value.RegisterNewFab(value.temperature_old_mf,   value.temperature_bc,   1,              nghost, "temperature_old",   false, true);
    value.RegisterNewFab(value.mass_fraction_mf,     value.mass_fraction_bc, value.nspecies, nghost, "mass_fraction",     true,  true);
    value.RegisterNewFab(value.mass_fraction_old_mf, value.mass_fraction_bc, value.nspecies, nghost, "mass_fraction_old", false, true);
    value.RegisterNewFab(value.eta_mf,               value.eta_bc,           1,              nghost, "eta",               true,  true);
    value.RegisterNewFab(value.eta_old_mf,           value.eta_bc,           1,              nghost, "eta_old",           false, true);
    value.RegisterNewFab(value.xi_mf,                value.xi_bc,            AMREX_SPACEDIM, nghost, "xi",                true,  true, {"x","y"});
    value.RegisterNewFab(value.xi_old_mf,            value.xi_bc,            AMREX_SPACEDIM, nghost, "xi_old",            false, true, {"x","y"});

    value.RegisterNewFab(value.density_mf,             &value.bc_nothing, 1,              0,      "density",             true,  false);
    value.RegisterNewFab(value.momentum_mf,            &value.bc_nothing, AMREX_SPACEDIM, 0,      "momentum",            true,  false, {"x","y"});
    value.RegisterNewFab(value.pressure_mf,            value.pressure_bc, 1,              nghost, "pressure",            true,  false);
    value.RegisterNewFab(value.pressure_correction_mf, &value.bc_nothing, 1,              nghost, "pressure_correction", false, false);
    value.RegisterNewFab(value.projection_rhs_mf,      &value.bc_nothing, 1,              0,      "projection_rhs",      false, false);
    value.RegisterNewFab(value.energy_mf,              &value.bc_nothing, 1,              0,      "energy",              true,  false);
    value.RegisterNewFab(value.mole_fraction_mf, &value.bc_nothing, value.nspecies, 0,      "mole_fraction", true, false);
    value.RegisterNewFab(value.vorticity_mf,     &value.bc_nothing, 1,              0,      "vorticity",     true, false);
    value.RegisterNewFab(value.deformation_gradient_mf, &value.bc_nothing, AMREX_SPACEDIM * AMREX_SPACEDIM, 1, "F", true, false, {"_xx","_xy","_yx","_yy"});
    value.RegisterNewFab(value.piola_stress_mf,         &value.bc_nothing, AMREX_SPACEDIM * AMREX_SPACEDIM, 1, "P", true, false, {"_xx","_xy","_yx","_yy"});
    value.RegisterNewFab(value.elastic_stress_mf,       &value.bc_nothing, AMREX_SPACEDIM * AMREX_SPACEDIM, 1, "elastic_stress", true, false, {"_xx","_xy","_yx","_yy"});
    value.RegisterGeneralFab<Set::Matrix>(value.cauchy_stress_mf, 1, 1, true, "stress", false);

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
LowMach::FillStateBoundaries(int lev,
                             amrex::MultiFab& u_mf,
                             amrex::MultiFab& T_mf,
                             amrex::MultiFab& Y_mf,
                             amrex::MultiFab& eta_mf,
                             amrex::MultiFab& xi_mf,
                             Set::Scalar time)
{
    auto preserve_valid = [](amrex::MultiFab& mf, auto&& fill)
    {
        amrex::MultiFab valid(mf.boxArray(), mf.DistributionMap(), mf.nComp(), 0);
        amrex::MultiFab::Copy(valid, mf, 0, 0, mf.nComp(), 0);
        mf.setVal(0.0, 0, mf.nComp(), mf.nGrow());
        amrex::MultiFab::Copy(mf, valid, 0, 0, mf.nComp(), 0);
        fill();
        amrex::MultiFab::Copy(mf, valid, 0, 0, mf.nComp(), 0);
    };

    auto fill_coarse_fine = [&](amrex::MultiFab& mf, Set::Field<Set::Scalar>& source_mf, BC::BC<Set::Scalar>* physbc)
    {
        if (lev == 0) return;

        amrex::Vector<amrex::MultiFab*> cmf;
        amrex::Vector<amrex::MultiFab*> fmf;
        cmf.push_back(source_mf[lev - 1].get());
        fmf.push_back(&mf);

        amrex::Vector<amrex::Real> ctime;
        amrex::Vector<amrex::Real> ftime;
        ctime.push_back(time);
        ftime.push_back(time);

        physbc->define(geom[lev]);
        amrex::Interpolater* mapper = mf.boxArray().ixType() == amrex::IndexType::TheNodeType() ?
            static_cast<amrex::Interpolater*>(&amrex::node_bilinear_interp) :
            static_cast<amrex::Interpolater*>(&amrex::cell_cons_interp);
        amrex::Vector<amrex::BCRec> bcs(mf.nComp(), physbc->GetBCRec());
        amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
            0, 0, mf.nComp(), geom[lev - 1], geom[lev],
            *physbc, 0,
            *physbc, 0,
            refRatio(lev - 1),
            mapper, bcs, 0);
    };

    preserve_valid(u_mf, [&]()
    {
        fill_coarse_fine(u_mf, velocity_mf, velocity_bc);
        velocity_bc->FillBoundary(u_mf, 0, AMREX_SPACEDIM, time, 0);
        u_mf.FillBoundary(geom[lev].periodicity());
    });
    ApplyVelocityDirichletCells(lev, u_mf, time);

    preserve_valid(T_mf, [&]()
    {
        fill_coarse_fine(T_mf, temperature_mf, temperature_bc);
        temperature_bc->FillBoundary(T_mf, 0, 1, time, 0);
        T_mf.FillBoundary(geom[lev].periodicity());
    });

    preserve_valid(Y_mf, [&]()
    {
        fill_coarse_fine(Y_mf, mass_fraction_mf, mass_fraction_bc);
        mass_fraction_bc->FillBoundary(Y_mf, 0, nspecies, time, 0);
        Y_mf.FillBoundary(geom[lev].periodicity());
    });

    preserve_valid(eta_mf, [&]()
    {
        fill_coarse_fine(eta_mf, this->eta_mf, eta_bc);
        eta_bc->FillBoundary(eta_mf, 0, 1, time, 0);
        eta_mf.FillBoundary(geom[lev].periodicity());
    });
    SanitizeEta(eta_mf);

    preserve_valid(xi_mf, [&]()
    {
        fill_coarse_fine(xi_mf, this->xi_mf, xi_bc);
        xi_bc->FillBoundary(xi_mf, 0, AMREX_SPACEDIM, time, 0);
        xi_mf.FillBoundary(geom[lev].periodicity());
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
    pressure_mf[lev]->setVal(0.0, 0, pressure_mf[lev]->nComp(), pressure_mf[lev]->nGrow());
    amrex::MultiFab::Copy(*pressure_mf[lev], valid, 0, 0, pressure_mf[lev]->nComp(), 0);
    if (lev > 0)
    {
        amrex::Vector<amrex::MultiFab*> cmf;
        amrex::Vector<amrex::MultiFab*> fmf;
        cmf.push_back(pressure_mf[lev - 1].get());
        fmf.push_back(pressure_mf[lev].get());
        amrex::Vector<amrex::Real> ctime;
        amrex::Vector<amrex::Real> ftime;
        ctime.push_back(time);
        ftime.push_back(time);
        pressure_bc->define(geom[lev]);
        amrex::Vector<amrex::BCRec> bcs(pressure_mf[lev]->nComp(), pressure_bc->GetBCRec());
        amrex::FillPatchTwoLevels(*pressure_mf[lev], time, cmf, ctime, fmf, ftime,
            0, 0, pressure_mf[lev]->nComp(), geom[lev - 1], geom[lev],
            *pressure_bc, 0,
            *pressure_bc, 0,
            refRatio(lev - 1),
            &amrex::cell_cons_interp, bcs, 0);
        amrex::MultiFab::Copy(*pressure_mf[lev], valid, 0, 0, pressure_mf[lev]->nComp(), 0);
    }
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
LowMach::UpdateSolidStress(int lev,
                           const amrex::MultiFab& u_mf,
                           const amrex::MultiFab& eta_mf,
                           const amrex::MultiFab& xi_mf)
{
    BL_PROFILE("Integrator::LowMach::UpdateSolidStress");

    deformation_gradient_mf[lev]->setVal(0.0, 0, deformation_gradient_mf[lev]->nComp(), deformation_gradient_mf[lev]->nGrow());
    piola_stress_mf[lev]->setVal(0.0, 0, piola_stress_mf[lev]->nComp(), piola_stress_mf[lev]->nGrow());
    elastic_stress_mf[lev]->setVal(0.0, 0, elastic_stress_mf[lev]->nComp(), elastic_stress_mf[lev]->nGrow());

    for (amrex::MFIter mfi(*cauchy_stress_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox();
        amrex::Array4<Set::Matrix> const& stress = cauchy_stress_mf[lev]->array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            stress(i,j,k) = Set::Matrix::Zero();
        });
    }

    if (!finite_solid_enabled)
    {
        cauchy_stress_mf[lev]->FillBoundary(geom[lev].periodicity());
        return;
    }

    const Set::Scalar* DX = geom[lev].CellSize();
    amrex::Box domain = geom[lev].Domain();
    const Model::Solid::Finite::NeoHookean solid_model = finite_solid_model;
    const Set::Scalar eta_threshold = lowmach_clamp(finite_solid_eta_threshold, 0.0, 1.0);
    const Set::Scalar det_floor = lowmach_max(small, finite_solid_J_floor);
    const Set::Scalar solid_viscosity = finite_solid_viscosity;
    const Set::Scalar solid_bulk_viscosity = finite_solid_bulk_viscosity;
    const Set::Scalar stress_limit = finite_solid_stress_cap;

    for (amrex::MFIter mfi(xi_mf, true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        Set::Patch<const Set::Scalar> u = u_mf.array(mfi);
        Set::Patch<const Set::Scalar> eta = eta_mf.array(mfi);
        Set::Patch<const Set::Scalar> xi = xi_mf.array(mfi);
        Set::Patch<Set::Scalar> F_field = deformation_gradient_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> P_field = piola_stress_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> elastic_stress = elastic_stress_mf.Patch(lev,mfi);
        amrex::Array4<Set::Matrix> const& cauchy_stress = cauchy_stress_mf[lev]->array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            auto sten = Numeric::GetStencil(i, j, k, domain);
            Set::Matrix F = Set::Matrix::Identity();
            Set::Matrix P = Set::Matrix::Zero();
            Set::Matrix sigma = Set::Matrix::Zero();

            Set::Scalar eta_val = lowmach_clamp(eta(i,j,k), 0.0, 1.0);
            Set::Scalar solid_weight = 0.0;
            if (eta_val > eta_threshold)
            {
                Set::Scalar denom = lowmach_max(1.0 - eta_threshold, 1.0e-12);
                solid_weight = lowmach_smootherstep((eta_val - eta_threshold) / denom);
            }

            if (solid_weight > 0.0)
            {
                Set::Matrix grad_xi = Set::Matrix::Zero();
                bool grad_valid = true;
                auto eta_is_solid = [=] AMREX_GPU_DEVICE (int ii, int jj, int kk) -> bool
                {
                    return lowmach_clamp(eta(ii,jj,kk), 0.0, 1.0) >= eta_threshold;
                };
                auto xi_derivative = [=, &grad_valid] AMREX_GPU_DEVICE (int n, int d) -> Set::Scalar
                {
                    int lo_i = i;
                    int lo_j = j;
                    int lo_k = k;
                    int hi_i = i;
                    int hi_j = j;
                    int hi_k = k;
                    if (d == 0)
                    {
                        lo_i -= 1;
                        hi_i += 1;
                    }
#if AMREX_SPACEDIM >= 2
                    else if (d == 1)
                    {
                        lo_j -= 1;
                        hi_j += 1;
                    }
#endif
#if AMREX_SPACEDIM == 3
                    else
                    {
                        lo_k -= 1;
                        hi_k += 1;
                    }
#endif

                    bool lo_valid = eta_is_solid(lo_i, lo_j, lo_k);
                    bool hi_valid = eta_is_solid(hi_i, hi_j, hi_k);
                    if (lo_valid && hi_valid)
                        return (xi(hi_i,hi_j,hi_k,n) - xi(lo_i,lo_j,lo_k,n)) / (2.0 * DX[d]);
                    if (hi_valid)
                        return (xi(hi_i,hi_j,hi_k,n) - xi(i,j,k,n)) / DX[d];
                    if (lo_valid)
                        return (xi(i,j,k,n) - xi(lo_i,lo_j,lo_k,n)) / DX[d];

                    grad_valid = false;
                    return n == d ? 1.0 : 0.0;
                };

                for (int n = 0; n < AMREX_SPACEDIM; ++n)
                    for (int d = 0; d < AMREX_SPACEDIM; ++d)
                        grad_xi(n,d) = xi_derivative(n, d);

                Set::Scalar det_grad_xi = grad_xi.determinant();
                bool valid = grad_valid && lowmach_is_finite(det_grad_xi) && lowmach_abs(det_grad_xi) > det_floor;
                if (valid)
                {
#if AMREX_SPACEDIM == 2
                    Set::Scalar inv_det = 1.0 / det_grad_xi;
                    F(0,0) =  grad_xi(1,1) * inv_det;
                    F(0,1) = -grad_xi(0,1) * inv_det;
                    F(1,0) = -grad_xi(1,0) * inv_det;
                    F(1,1) =  grad_xi(0,0) * inv_det;
#else
                    F = grad_xi.inverse();
#endif
                    for (int a = 0; a < AMREX_SPACEDIM; ++a)
                        for (int b = 0; b < AMREX_SPACEDIM; ++b)
                            valid = valid && lowmach_is_finite(F(a,b));
                }

                Set::Scalar J = F.determinant();
                if (valid && lowmach_is_finite(J) && lowmach_abs(J) > det_floor)
                {
                    P = solid_model.DW(F);
                    sigma = (P * F.transpose()) / J;
                }

                if (solid_viscosity != 0.0 || solid_bulk_viscosity != 0.0)
                {
                    Set::Matrix grad_u = Numeric::Gradient(u, i, j, k, DX, sten);
                    Set::Scalar div_u = grad_u.trace();
                    sigma += solid_viscosity * (grad_u + grad_u.transpose()) +
                             solid_bulk_viscosity * div_u * Set::Matrix::Identity();
                }
            }

            sigma = lowmach_regularize_stress(solid_weight * sigma, stress_limit);

            F_field(i,j,k,0) = lowmach_finite_or(F(0,0), 1.0);
            F_field(i,j,k,1) = lowmach_finite_or(F(0,1), 0.0);
            F_field(i,j,k,2) = lowmach_finite_or(F(1,0), 0.0);
            F_field(i,j,k,3) = lowmach_finite_or(F(1,1), 1.0);
            P_field(i,j,k,0) = lowmach_finite_or(P(0,0), 0.0);
            P_field(i,j,k,1) = lowmach_finite_or(P(0,1), 0.0);
            P_field(i,j,k,2) = lowmach_finite_or(P(1,0), 0.0);
            P_field(i,j,k,3) = lowmach_finite_or(P(1,1), 0.0);
            elastic_stress(i,j,k,0) = lowmach_finite_or(sigma(0,0), 0.0);
            elastic_stress(i,j,k,1) = lowmach_finite_or(sigma(0,1), 0.0);
            elastic_stress(i,j,k,2) = lowmach_finite_or(sigma(1,0), 0.0);
            elastic_stress(i,j,k,3) = lowmach_finite_or(sigma(1,1), 0.0);
            cauchy_stress(i,j,k) = sigma;
        });
    }

    deformation_gradient_mf[lev]->FillBoundary(geom[lev].periodicity());
    piola_stress_mf[lev]->FillBoundary(geom[lev].periodicity());
    elastic_stress_mf[lev]->FillBoundary(geom[lev].periodicity());
    cauchy_stress_mf[lev]->FillBoundary(geom[lev].periodicity());
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
LowMach::ImplicitElasticVelocitySolve(int lev, Set::Scalar time, Set::Scalar dt)
{
    BL_PROFILE("Integrator::LowMach::ImplicitElasticVelocitySolve");
    if (!finite_solid_enabled || !finite_solid_implicit || !(dt > 0.0)) return;
    if (!(finite_solid_rhs_sign != 0.0)) return;

    amrex::MultiFab& u_mf = *velocity_mf[lev];
    amrex::MultiFab& eta = *eta_mf[lev];
    amrex::MultiFab& rho_mf = *density_mf[lev];

    const Set::Scalar eta_threshold = lowmach_clamp(finite_solid_eta_threshold, 0.0, 1.0);
    const Set::Scalar rho_floor = density_floor;
    const Set::Scalar stress_sign = finite_solid_rhs_sign;
    if (!(finite_solid_model.mu == finite_solid_model.mu) || !(finite_solid_model.kappa == finite_solid_model.kappa)) return;
    if (!(finite_solid_implicit_coeff_scale == finite_solid_implicit_coeff_scale) || !(finite_solid_implicit_coeff_scale > 0.0)) return;

    Set::Scalar rigid_mass = 0.0;
    Set::Scalar rigid_xcm = 0.0;
    Set::Scalar rigid_ycm = 0.0;
    Set::Scalar rigid_ux = 0.0;
    Set::Scalar rigid_uy = 0.0;
    Set::Scalar rigid_omega = 0.0;
    if (finite_solid_implicit_preserve_rigid_modes)
    {
        Set::Scalar sum_m = 0.0;
        Set::Scalar sum_mx = 0.0;
        Set::Scalar sum_my = 0.0;
        Set::Scalar sum_mux = 0.0;
        Set::Scalar sum_muy = 0.0;
        const Set::Scalar* DX = geom[lev].CellSize();
        const Set::Scalar cell_vol = AMREX_D_TERM(DX[0], * DX[1], * DX[2]);
        for (amrex::MFIter mfi(u_mf, false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            const amrex::Dim3 lo = amrex::lbound(bx);
            const amrex::Dim3 hi = amrex::ubound(bx);
            Set::Patch<const Set::Scalar> eta_patch = eta.array(mfi);
            Set::Patch<const Set::Scalar> rho = rho_mf.array(mfi);
            Set::Patch<const Set::Scalar> u = u_mf.array(mfi);
            for (int k = lo.z; k <= hi.z; ++k)
                for (int j = lo.y; j <= hi.y; ++j)
                    for (int i = lo.x; i <= hi.x; ++i)
                    {
                        Set::Scalar eta_val = lowmach_clamp(eta_patch(i,j,k), 0.0, 1.0);
                        if (eta_val <= eta_threshold) continue;
                        Set::Scalar denom = lowmach_max(1.0 - eta_threshold, 1.0e-12);
                        Set::Scalar solid_weight = lowmach_smootherstep((eta_val - eta_threshold) / denom);
                        Set::Scalar m = solid_weight * lowmach_max(rho(i,j,k), rho_floor) * cell_vol;
                        Set::Vector x = Set::Position(i, j, k, geom[lev], amrex::IndexType::TheCellType());
                        sum_m += m;
                        sum_mx += m * x(0);
                        sum_my += m * x(1);
                        sum_mux += m * u(i,j,k,0);
                        sum_muy += m * u(i,j,k,1);
                    }
        }
        amrex::ParallelDescriptor::ReduceRealSum(sum_m);
        amrex::ParallelDescriptor::ReduceRealSum(sum_mx);
        amrex::ParallelDescriptor::ReduceRealSum(sum_my);
        amrex::ParallelDescriptor::ReduceRealSum(sum_mux);
        amrex::ParallelDescriptor::ReduceRealSum(sum_muy);

        rigid_mass = sum_m;
        if (rigid_mass > 0.0)
        {
            rigid_xcm = sum_mx / rigid_mass;
            rigid_ycm = sum_my / rigid_mass;
            rigid_ux = sum_mux / rigid_mass;
            rigid_uy = sum_muy / rigid_mass;

#if AMREX_SPACEDIM == 2
            Set::Scalar sum_I = 0.0;
            Set::Scalar sum_L = 0.0;
            for (amrex::MFIter mfi(u_mf, false); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                const amrex::Dim3 lo = amrex::lbound(bx);
                const amrex::Dim3 hi = amrex::ubound(bx);
                Set::Patch<const Set::Scalar> eta_patch = eta.array(mfi);
                Set::Patch<const Set::Scalar> rho = rho_mf.array(mfi);
                Set::Patch<const Set::Scalar> u = u_mf.array(mfi);
                for (int k = lo.z; k <= hi.z; ++k)
                    for (int j = lo.y; j <= hi.y; ++j)
                        for (int i = lo.x; i <= hi.x; ++i)
                        {
                            Set::Scalar eta_val = lowmach_clamp(eta_patch(i,j,k), 0.0, 1.0);
                            if (eta_val <= eta_threshold) continue;
                            Set::Scalar denom = lowmach_max(1.0 - eta_threshold, 1.0e-12);
                            Set::Scalar solid_weight = lowmach_smootherstep((eta_val - eta_threshold) / denom);
                            Set::Scalar m = solid_weight * lowmach_max(rho(i,j,k), rho_floor) * cell_vol;
                            Set::Vector x = Set::Position(i, j, k, geom[lev], amrex::IndexType::TheCellType());
                            Set::Scalar rx = x(0) - rigid_xcm;
                            Set::Scalar ry = x(1) - rigid_ycm;
                            Set::Scalar vx = u(i,j,k,0) - rigid_ux;
                            Set::Scalar vy = u(i,j,k,1) - rigid_uy;
                            sum_I += m * (rx * rx + ry * ry);
                            sum_L += m * (rx * vy - ry * vx);
                        }
            }
            amrex::ParallelDescriptor::ReduceRealSum(sum_I);
            amrex::ParallelDescriptor::ReduceRealSum(sum_L);
            if (sum_I > 0.0) rigid_omega = sum_L / sum_I;
#endif
        }
    }

    auto rigid_velocity = [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k, int comp) -> Set::Scalar
    {
        if (!(rigid_mass > 0.0)) return 0.0;
        Set::Vector x = Set::Position(i, j, k, geom[lev], amrex::IndexType::TheCellType());
        if (comp == 0) return rigid_ux - rigid_omega * (x(1) - rigid_ycm);
        if (comp == 1) return rigid_uy + rigid_omega * (x(0) - rigid_xcm);
#if AMREX_SPACEDIM == 3
        if (comp == 2) return 0.0;
#endif
        return 0.0;
    };

    amrex::MultiFab rhs(u_mf.boxArray(), u_mf.DistributionMap(), AMREX_SPACEDIM, 0);
    amrex::MultiFab sol(u_mf.boxArray(), u_mf.DistributionMap(), AMREX_SPACEDIM, u_mf.nGrow());
    amrex::MultiFab::Copy(rhs, u_mf, 0, 0, AMREX_SPACEDIM, 0);
    amrex::MultiFab::Copy(sol, u_mf, 0, 0, AMREX_SPACEDIM, u_mf.nGrow());

    if (rigid_mass > 0.0)
    {
        for (amrex::MFIter mfi(rhs, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.growntilebox();
            amrex::Array4<Set::Scalar> const& rhs_arr = rhs.array(mfi);
            amrex::Array4<Set::Scalar> const& sol_arr = sol.array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                for (int comp = 0; comp < AMREX_SPACEDIM; ++comp)
                {
                    Set::Scalar ur = rigid_velocity(i, j, k, comp);
                    rhs_arr(i,j,k,comp) -= ur;
                    sol_arr(i,j,k,comp) -= ur;
                }
            });
        }
        amrex::Gpu::streamSynchronize();
    }

    amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> lobc;
    amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> hibc;
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        lobc[d] = geom[lev].isPeriodic(d) ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Neumann;
        hibc[d] = geom[lev].isPeriodic(d) ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Neumann;
    }

    amrex::LPInfo info;
    Operator::Dynamic::ElasticLowMach elastic_op({geom[lev]}, {u_mf.boxArray()}, {u_mf.DistributionMap()}, info);
    elastic_op.setDomainBC(lobc, hibc);
    elastic_op.setMaxOrder(2);
    elastic_op.setLevelBC(0, nullptr);
    elastic_op.setScalars(1.0, dt * stress_sign);
    elastic_op.setACoeffs(0, 1.0);
    elastic_op.SetCoefficients(0, eta, rho_mf, finite_solid_model, eta_threshold, rho_floor, finite_solid_implicit_coeff_scale);

    amrex::MLMG mlmg(elastic_op);
    mlmg.setVerbose(finite_solid_implicit_verbose);
    mlmg.solve({&sol}, {&rhs}, finite_solid_implicit_tol_rel, finite_solid_implicit_tol_abs);

    if (rigid_mass > 0.0)
    {
        for (amrex::MFIter mfi(sol, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<Set::Scalar> const& sol_arr = sol.array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                for (int comp = 0; comp < AMREX_SPACEDIM; ++comp)
                    sol_arr(i,j,k,comp) += rigid_velocity(i, j, k, comp);
            });
        }
        amrex::Gpu::streamSynchronize();
    }

    amrex::MultiFab::Copy(u_mf, sol, 0, 0, AMREX_SPACEDIM, 0);
    FillStateBoundaries(lev, *velocity_mf[lev], *temperature_mf[lev], *mass_fraction_mf[lev], *eta_mf[lev], *xi_mf[lev], time);
}

void
LowMach::ProjectVelocity(int lev, Set::Scalar time, Set::Scalar dt)
{
    BL_PROFILE("Integrator::LowMach::ProjectVelocity");
    if (!projection_enabled || !(dt > 0.0)) return;
    if (finest_level > 0 && !projection_amr_enabled) return;

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
    if (lev > 0)
    {
        amrex::MultiFab valid_phi(phi_mf.boxArray(), phi_mf.DistributionMap(), phi_mf.nComp(), 0);
        amrex::MultiFab::Copy(valid_phi, phi_mf, 0, 0, phi_mf.nComp(), 0);
        amrex::Vector<amrex::MultiFab*> cmf;
        amrex::Vector<amrex::MultiFab*> fmf;
        cmf.push_back(pressure_correction_mf[lev - 1].get());
        fmf.push_back(&phi_mf);
        amrex::Vector<amrex::Real> ctime;
        amrex::Vector<amrex::Real> ftime;
        ctime.push_back(time);
        ftime.push_back(time);
        bc_nothing.define(geom[lev]);
        amrex::Vector<amrex::BCRec> bcs(phi_mf.nComp(), bc_nothing.GetBCRec());
        amrex::FillPatchTwoLevels(phi_mf, time, cmf, ctime, fmf, ftime,
            0, 0, phi_mf.nComp(), geom[lev - 1], geom[lev],
            bc_nothing, 0,
            bc_nothing, 0,
            refRatio(lev - 1),
            &amrex::cell_cons_interp, bcs, 0);
        amrex::MultiFab::Copy(phi_mf, valid_phi, 0, 0, phi_mf.nComp(), 0);
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
LowMach::SanitizeEta(amrex::MultiFab& eta_mf)
{
    const amrex::IndexType type = eta_mf.ixType();
    for (amrex::MFIter mfi(eta_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        amrex::Box bx;
        if (type == amrex::IndexType::TheCellType()) bx = mfi.growntilebox();
        else if (type == amrex::IndexType::TheNodeType()) bx = mfi.grownnodaltilebox();
        else Util::Abort(INFO, "Unknown eta index type");

        amrex::Array4<Set::Scalar> const& eta = eta_mf.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Scalar val = eta(i,j,k);
            if (!(val == val && val < 1.0e300 && val > -1.0e300)) val = 0.0;
            eta(i,j,k) = lowmach_clamp(val, 0.0, 1.0);
        });
    }
    amrex::Gpu::streamSynchronize();
}

void
LowMach::FillReferenceMapBoundary(int lev, amrex::MultiFab& xi_mf)
{
    amrex::MultiFab valid(xi_mf.boxArray(), xi_mf.DistributionMap(), xi_mf.nComp(), 0);
    amrex::MultiFab::Copy(valid, xi_mf, 0, 0, xi_mf.nComp(), 0);

    const amrex::IndexType type = xi_mf.ixType();
    const amrex::Geometry geom_lev = geom[lev];
    for (amrex::MFIter mfi(xi_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        amrex::Box bx;
        if (type == amrex::IndexType::TheCellType()) bx = mfi.growntilebox();
        else if (type == amrex::IndexType::TheNodeType()) bx = mfi.grownnodaltilebox();
        else Util::Abort(INFO, "Unknown xi index type");

        Set::Patch<Set::Scalar> xi = xi_mf.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Vector x = Set::Position(i, j, k, geom_lev, type);
            for (int d = 0; d < AMREX_SPACEDIM; ++d) xi(i,j,k,d) = x(d);
        });
    }
    amrex::Gpu::streamSynchronize();

    amrex::MultiFab::Copy(xi_mf, valid, 0, 0, xi_mf.nComp(), 0);
    xi_mf.FillBoundary(geom[lev].periodicity());
}

void
LowMach::InitializeReferenceMap(int lev, amrex::MultiFab& xi_mf)
{
    const amrex::IndexType type = xi_mf.ixType();
    for (amrex::MFIter mfi(xi_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        amrex::Box bx;
        if (type == amrex::IndexType::TheCellType()) bx = mfi.growntilebox();
        else if (type == amrex::IndexType::TheNodeType()) bx = mfi.grownnodaltilebox();
        else Util::Abort(INFO, "Unknown xi index type");

        Set::Patch<Set::Scalar> xi = xi_mf.array(mfi);
        const amrex::Geometry geom_lev = geom[lev];
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Vector x = Set::Position(i, j, k, geom_lev, type);
            for (int d = 0; d < AMREX_SPACEDIM; ++d) xi(i,j,k,d) = x(d);
        });
    }
    FillReferenceMapBoundary(lev, xi_mf);
}

void
LowMach::RebuildReferenceMapOutsideEta(int lev, const amrex::MultiFab& eta_stage_mf, amrex::MultiFab& xi_stage_mf, Set::Scalar time)
{
    const Set::Scalar extension_eta = lowmach_clamp(reference_map_eta_cutoff, 0.0, 1.0);
    const Set::Scalar core_eta = lowmach_clamp(lowmach_max(reference_map_eta_core, extension_eta), 0.0, 1.0);
    const int extrap_sweeps = reference_map_extrapolation_sweeps < 0 ? 0 : reference_map_extrapolation_sweeps;
    const int smooth_sweeps = reference_map_smoothing_sweeps < 0 ? 0 : reference_map_smoothing_sweeps;
    const int xi_ngrow = xi_stage_mf.nGrow();
    amrex::Geometry const geom_lev = geom[lev];

    amrex::MultiFab eta_work(eta_stage_mf.boxArray(), eta_stage_mf.DistributionMap(), 1, eta_stage_mf.nGrow());
    amrex::MultiFab::Copy(eta_work, eta_stage_mf, 0, 0, 1, eta_stage_mf.nGrow());
    {
        amrex::MultiFab eta_valid(eta_work.boxArray(), eta_work.DistributionMap(), eta_work.nComp(), 0);
        amrex::MultiFab::Copy(eta_valid, eta_work, 0, 0, eta_work.nComp(), 0);
        eta_bc->FillBoundary(eta_work, 0, 1, time, 0);
        eta_work.FillBoundary(geom[lev].periodicity());
        amrex::MultiFab::Copy(eta_work, eta_valid, 0, 0, eta_work.nComp(), 0);
    }
    SanitizeEta(eta_work);

    auto fill_xi_boundary_preserve_valid = [&]()
    {
        FillReferenceMapBoundary(lev, xi_stage_mf);
    };

    fill_xi_boundary_preserve_valid();

    amrex::MultiFab xi_protected(xi_stage_mf.boxArray(), xi_stage_mf.DistributionMap(), AMREX_SPACEDIM, xi_ngrow);
    amrex::MultiFab::Copy(xi_protected, xi_stage_mf, 0, 0, AMREX_SPACEDIM, xi_ngrow);

    auto restore_protected_xi = [&]()
    {
        for (amrex::MFIter mfi(xi_stage_mf, true); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& eta = eta_work.const_array(mfi);
            amrex::Array4<const Set::Scalar> const& xi_in = xi_protected.const_array(mfi);
            amrex::Array4<Set::Scalar> const& xi = xi_stage_mf.array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                if (eta(i,j,k) < core_eta) return;
                for (int n = 0; n < AMREX_SPACEDIM; ++n) xi(i,j,k,n) = xi_in(i,j,k,n);
            });
        }
        amrex::Gpu::streamSynchronize();
    };

    amrex::MultiFab xi_known(xi_stage_mf.boxArray(), xi_stage_mf.DistributionMap(), 1, xi_ngrow);
    xi_known.setVal(0.0, 0, 1, xi_ngrow);

    for (amrex::MFIter mfi(xi_stage_mf, true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        amrex::Array4<const Set::Scalar> const& eta = eta_work.const_array(mfi);
        amrex::Array4<Set::Scalar> const& known = xi_known.array(mfi);
        amrex::Array4<Set::Scalar> const& xi = xi_stage_mf.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            if (eta(i,j,k) >= core_eta)
            {
                known(i,j,k) = 1.0;
                return;
            }

            Set::Vector pos = Set::Position(i, j, k, geom_lev, amrex::IndexType::TheCellType());
            for (int n = 0; n < AMREX_SPACEDIM; ++n) xi(i,j,k,n) = pos(n);
        });
    }
    amrex::Gpu::streamSynchronize();
    xi_known.FillBoundary(geom[lev].periodicity());
    restore_protected_xi();
    fill_xi_boundary_preserve_valid();

    for (int sweep = 0; sweep < extrap_sweeps; ++sweep)
    {
        amrex::MultiFab xi_old(xi_stage_mf.boxArray(), xi_stage_mf.DistributionMap(), AMREX_SPACEDIM, xi_ngrow);
        amrex::MultiFab known_old(xi_stage_mf.boxArray(), xi_stage_mf.DistributionMap(), 1, xi_ngrow);
        amrex::MultiFab::Copy(xi_old, xi_stage_mf, 0, 0, AMREX_SPACEDIM, xi_ngrow);
        amrex::MultiFab::Copy(known_old, xi_known, 0, 0, 1, xi_ngrow);
        known_old.FillBoundary(geom[lev].periodicity());

        for (amrex::MFIter mfi(xi_stage_mf, true); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& eta = eta_work.const_array(mfi);
            amrex::Array4<const Set::Scalar> const& xi_in = xi_old.const_array(mfi);
            amrex::Array4<const Set::Scalar> const& known_in = known_old.const_array(mfi);
            amrex::Array4<Set::Scalar> const& xi = xi_stage_mf.array(mfi);
            amrex::Array4<Set::Scalar> const& known = xi_known.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                if (known_in(i,j,k) > 0.5) return;
                if (eta(i,j,k) < extension_eta) return;

                Set::Vector pos = Set::Position(i, j, k, geom_lev, amrex::IndexType::TheCellType());
                for (int n = 0; n < AMREX_SPACEDIM; ++n)
                {
                    Set::Scalar sum = 0.0;
                    Set::Scalar count = 0.0;
                    if (known_in(i-1,j,k) > 0.5)
                    {
                        Set::Vector npos = Set::Position(i-1, j, k, geom_lev, amrex::IndexType::TheCellType());
                        sum += xi_in(i-1,j,k,n) - npos(n);
                        count += 1.0;
                    }
                    if (known_in(i+1,j,k) > 0.5)
                    {
                        Set::Vector npos = Set::Position(i+1, j, k, geom_lev, amrex::IndexType::TheCellType());
                        sum += xi_in(i+1,j,k,n) - npos(n);
                        count += 1.0;
                    }
#if AMREX_SPACEDIM >= 2
                    if (known_in(i,j-1,k) > 0.5)
                    {
                        Set::Vector npos = Set::Position(i, j-1, k, geom_lev, amrex::IndexType::TheCellType());
                        sum += xi_in(i,j-1,k,n) - npos(n);
                        count += 1.0;
                    }
                    if (known_in(i,j+1,k) > 0.5)
                    {
                        Set::Vector npos = Set::Position(i, j+1, k, geom_lev, amrex::IndexType::TheCellType());
                        sum += xi_in(i,j+1,k,n) - npos(n);
                        count += 1.0;
                    }
#endif
#if AMREX_SPACEDIM == 3
                    if (known_in(i,j,k-1) > 0.5)
                    {
                        Set::Vector npos = Set::Position(i, j, k-1, geom_lev, amrex::IndexType::TheCellType());
                        sum += xi_in(i,j,k-1,n) - npos(n);
                        count += 1.0;
                    }
                    if (known_in(i,j,k+1) > 0.5)
                    {
                        Set::Vector npos = Set::Position(i, j, k+1, geom_lev, amrex::IndexType::TheCellType());
                        sum += xi_in(i,j,k+1,n) - npos(n);
                        count += 1.0;
                    }
#endif
                    if (count > 0.0) xi(i,j,k,n) = pos(n) + sum / count;
                }

                if (known_in(i-1,j,k) > 0.5 || known_in(i+1,j,k) > 0.5
#if AMREX_SPACEDIM >= 2
                    || known_in(i,j-1,k) > 0.5 || known_in(i,j+1,k) > 0.5
#endif
#if AMREX_SPACEDIM == 3
                    || known_in(i,j,k-1) > 0.5 || known_in(i,j,k+1) > 0.5
#endif
                    )
                {
                    known(i,j,k) = 1.0;
                }
            });
        }
        amrex::Gpu::streamSynchronize();
        xi_known.FillBoundary(geom[lev].periodicity());
        restore_protected_xi();
        fill_xi_boundary_preserve_valid();
    }

    for (int sweep = 0; sweep < smooth_sweeps; ++sweep)
    {
        amrex::MultiFab xi_old(xi_stage_mf.boxArray(), xi_stage_mf.DistributionMap(), AMREX_SPACEDIM, xi_ngrow);
        amrex::MultiFab known_old(xi_stage_mf.boxArray(), xi_stage_mf.DistributionMap(), 1, xi_ngrow);
        amrex::MultiFab::Copy(xi_old, xi_stage_mf, 0, 0, AMREX_SPACEDIM, xi_ngrow);
        amrex::MultiFab::Copy(known_old, xi_known, 0, 0, 1, xi_ngrow);
        known_old.FillBoundary(geom[lev].periodicity());

        for (amrex::MFIter mfi(xi_stage_mf, true); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& eta = eta_work.const_array(mfi);
            amrex::Array4<const Set::Scalar> const& xi_in = xi_old.const_array(mfi);
            amrex::Array4<const Set::Scalar> const& known = known_old.const_array(mfi);
            amrex::Array4<Set::Scalar> const& xi = xi_stage_mf.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                if (eta(i,j,k) >= core_eta) return;
                if (eta(i,j,k) < extension_eta) return;
                if (known(i,j,k) <= 0.5) return;

                Set::Vector pos = Set::Position(i, j, k, geom_lev, amrex::IndexType::TheCellType());
                for (int n = 0; n < AMREX_SPACEDIM; ++n)
                {
                    Set::Scalar sum = xi_in(i,j,k,n) - pos(n);
                    Set::Scalar count = 1.0;
                    if (known(i-1,j,k) > 0.5)
                    {
                        Set::Vector npos = Set::Position(i-1, j, k, geom_lev, amrex::IndexType::TheCellType());
                        sum += xi_in(i-1,j,k,n) - npos(n);
                        count += 1.0;
                    }
                    if (known(i+1,j,k) > 0.5)
                    {
                        Set::Vector npos = Set::Position(i+1, j, k, geom_lev, amrex::IndexType::TheCellType());
                        sum += xi_in(i+1,j,k,n) - npos(n);
                        count += 1.0;
                    }
#if AMREX_SPACEDIM >= 2
                    if (known(i,j-1,k) > 0.5)
                    {
                        Set::Vector npos = Set::Position(i, j-1, k, geom_lev, amrex::IndexType::TheCellType());
                        sum += xi_in(i,j-1,k,n) - npos(n);
                        count += 1.0;
                    }
                    if (known(i,j+1,k) > 0.5)
                    {
                        Set::Vector npos = Set::Position(i, j+1, k, geom_lev, amrex::IndexType::TheCellType());
                        sum += xi_in(i,j+1,k,n) - npos(n);
                        count += 1.0;
                    }
#endif
#if AMREX_SPACEDIM == 3
                    if (known(i,j,k-1) > 0.5)
                    {
                        Set::Vector npos = Set::Position(i, j, k-1, geom_lev, amrex::IndexType::TheCellType());
                        sum += xi_in(i,j,k-1,n) - npos(n);
                        count += 1.0;
                    }
                    if (known(i,j,k+1) > 0.5)
                    {
                        Set::Vector npos = Set::Position(i, j, k+1, geom_lev, amrex::IndexType::TheCellType());
                        sum += xi_in(i,j,k+1,n) - npos(n);
                        count += 1.0;
                    }
#endif
                    xi(i,j,k,n) = pos(n) + sum / count;
                }
            });
        }
        amrex::Gpu::streamSynchronize();
        restore_protected_xi();
        fill_xi_boundary_preserve_valid();
    }

    restore_protected_xi();

    for (amrex::MFIter mfi(xi_stage_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox();
        amrex::Array4<Set::Scalar> const& xi = xi_stage_mf.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Vector pos = Set::Position(i, j, k, geom_lev, amrex::IndexType::TheCellType());
            for (int n = 0; n < AMREX_SPACEDIM; ++n)
            {
                Set::Scalar val = xi(i,j,k,n);
                if (!(val == val && val < 1.0e300 && val > -1.0e300)) xi(i,j,k,n) = pos(n);
            }
        });
    }
    amrex::Gpu::streamSynchronize();
    FillReferenceMapBoundary(lev, xi_stage_mf);
}

void
LowMach::Initialize(int lev)
{
    BL_PROFILE("Integrator::LowMach::Initialize");

    velocity_mf[lev]->setVal(0.0, 0, velocity_mf[lev]->nComp(), velocity_mf[lev]->nGrow());
    velocity_old_mf[lev]->setVal(0.0, 0, velocity_old_mf[lev]->nComp(), velocity_old_mf[lev]->nGrow());
    temperature_mf[lev]->setVal(0.0, 0, temperature_mf[lev]->nComp(), temperature_mf[lev]->nGrow());
    temperature_old_mf[lev]->setVal(0.0, 0, temperature_old_mf[lev]->nComp(), temperature_old_mf[lev]->nGrow());
    mass_fraction_mf[lev]->setVal(0.0, 0, mass_fraction_mf[lev]->nComp(), mass_fraction_mf[lev]->nGrow());
    mass_fraction_old_mf[lev]->setVal(0.0, 0, mass_fraction_old_mf[lev]->nComp(), mass_fraction_old_mf[lev]->nGrow());
    eta_mf[lev]->setVal(0.0, 0, eta_mf[lev]->nComp(), eta_mf[lev]->nGrow());
    eta_old_mf[lev]->setVal(0.0, 0, eta_old_mf[lev]->nComp(), eta_old_mf[lev]->nGrow());
    xi_mf[lev]->setVal(0.0, 0, xi_mf[lev]->nComp(), xi_mf[lev]->nGrow());
    xi_old_mf[lev]->setVal(0.0, 0, xi_old_mf[lev]->nComp(), xi_old_mf[lev]->nGrow());
    pressure_mf[lev]->setVal(0.0, 0, pressure_mf[lev]->nComp(), pressure_mf[lev]->nGrow());

    velocity_ic->Initialize(lev, velocity_mf, 0.0);
    velocity_ic->Initialize(lev, velocity_old_mf, 0.0);
    temperature_ic->Initialize(lev, temperature_mf, 0.0);
    temperature_ic->Initialize(lev, temperature_old_mf, 0.0);
    mass_fraction_ic->Initialize(lev, mass_fraction_mf, 0.0);
    mass_fraction_ic->Initialize(lev, mass_fraction_old_mf, 0.0);
    if (eta_ic)
    {
        eta_ic->Initialize(lev, eta_mf, 0.0);
        eta_ic->Initialize(lev, eta_old_mf, 0.0);
    }
    else
    {
        eta_mf[lev]->setVal(eta_initial_value);
        eta_old_mf[lev]->setVal(eta_initial_value);
    }
    SanitizeEta(*eta_mf[lev]);
    SanitizeEta(*eta_old_mf[lev]);
    InitializeReferenceMap(lev, *xi_mf[lev]);
    InitializeReferenceMap(lev, *xi_old_mf[lev]);
    pressure_ic->Initialize(lev, pressure_mf, 0.0);
    pressure_correction_mf[lev]->setVal(0.0);
    projection_rhs_mf[lev]->setVal(0.0);

    EnforceStateBounds(*velocity_mf[lev], *temperature_mf[lev], *mass_fraction_mf[lev]);
    EnforceStateBounds(*velocity_old_mf[lev], *temperature_old_mf[lev], *mass_fraction_old_mf[lev]);
    FillStateBoundaries(lev, *velocity_mf[lev], *temperature_mf[lev], *mass_fraction_mf[lev], *eta_mf[lev], *xi_mf[lev], 0.0);
    FillStateBoundaries(lev, *velocity_old_mf[lev], *temperature_old_mf[lev], *mass_fraction_old_mf[lev], *eta_old_mf[lev], *xi_old_mf[lev], 0.0);
    RebuildReferenceMapOutsideEta(lev, *eta_mf[lev], *xi_mf[lev], 0.0);
    RebuildReferenceMapOutsideEta(lev, *eta_old_mf[lev], *xi_old_mf[lev], 0.0);
    EnforceStateBounds(*velocity_mf[lev], *temperature_mf[lev], *mass_fraction_mf[lev]);
    EnforceStateBounds(*velocity_old_mf[lev], *temperature_old_mf[lev], *mass_fraction_old_mf[lev]);
    UpdateDerived(lev, *velocity_mf[lev], *temperature_mf[lev], *mass_fraction_mf[lev]);
    UpdateSolidStress(lev, *velocity_mf[lev], *eta_mf[lev], *xi_mf[lev]);
}

void
LowMach::RHS(int lev, Set::Scalar /*time*/,
             amrex::MultiFab& u_rhs_mf,
             amrex::MultiFab& T_rhs_mf,
             amrex::MultiFab& Y_rhs_mf,
             amrex::MultiFab& eta_rhs_mf,
             amrex::MultiFab& xi_rhs_mf,
             const amrex::MultiFab& u_mf,
             const amrex::MultiFab& T_mf,
             const amrex::MultiFab& Y_mf,
             const amrex::MultiFab& eta_mf,
             const amrex::MultiFab& xi_mf)
{
    UpdateDerived(lev, u_mf, T_mf, Y_mf);
    UpdateSolidStress(lev, u_mf, eta_mf, xi_mf);

    const Set::Scalar* DX = geom[lev].CellSize();
    amrex::Box domain = geom[lev].Domain();

    for (amrex::MFIter mfi(u_mf, false); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        Set::Patch<const Set::Scalar> u = u_mf.array(mfi);
        Set::Patch<const Set::Scalar> T = T_mf.array(mfi);
        Set::Patch<const Set::Scalar> Y = Y_mf.array(mfi);
        Set::Patch<const Set::Scalar> eta = eta_mf.array(mfi);
        Set::Patch<const Set::Scalar> xi = xi_mf.array(mfi);
        Set::Patch<const Set::Scalar> rho = density_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> pressure = pressure_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> X = mole_fraction_mf.Patch(lev,mfi);
        amrex::Array4<const Set::Matrix> const& stress = cauchy_stress_mf[lev]->const_array(mfi);
        Set::Patch<Set::Scalar> u_rhs = u_rhs_mf.array(mfi);
        Set::Patch<Set::Scalar> T_rhs = T_rhs_mf.array(mfi);
        Set::Patch<Set::Scalar> Y_rhs = Y_rhs_mf.array(mfi);
        Set::Patch<Set::Scalar> eta_rhs = eta_rhs_mf.array(mfi);
        Set::Patch<Set::Scalar> xi_rhs = xi_rhs_mf.array(mfi);
        const int nsp = nspecies;
        const bool viscous = include_viscosity;
        const bool conductive = include_conduction;
        const bool advect_T = advect_temperature;
        const Set::Scalar xi_eta_cutoff = lowmach_clamp(lowmach_max(reference_map_eta_core, reference_map_eta_cutoff), 0.0, 1.0);
        const Set::Vector gravity = g;
        const Set::Scalar p_scale = pressure_scale;
        const Set::Scalar rho_floor = density_floor;
        const Set::Scalar stress_sign = finite_solid_implicit ? 0.0 : finite_solid_rhs_sign;
        amrex::Box const cell_domain = domain;

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

            auto clamp_i = [=] AMREX_GPU_DEVICE (int ii) -> int
            {
                return ii < cell_domain.smallEnd(0) ? cell_domain.smallEnd(0) :
                       (ii > cell_domain.bigEnd(0) ? cell_domain.bigEnd(0) : ii);
            };
            auto clamp_j = [=] AMREX_GPU_DEVICE (int jj) -> int
            {
                return jj < cell_domain.smallEnd(1) ? cell_domain.smallEnd(1) :
                       (jj > cell_domain.bigEnd(1) ? cell_domain.bigEnd(1) : jj);
            };
            int im = clamp_i(i - 1);
            int ip = clamp_i(i + 1);
            int jm = clamp_j(j - 1);
            int jp = clamp_j(j + 1);
            Set::Scalar dx_den = lowmach_max(static_cast<Set::Scalar>(ip - im) * DX[0], DX[0]);
            Set::Scalar dy_den = lowmach_max(static_cast<Set::Scalar>(jp - jm) * DX[1], DX[1]);
            Set::Vector div_sigma = Set::Vector::Zero();
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
            {
                div_sigma(d) += (stress(ip,j,k)(d,0) - stress(im,j,k)(d,0)) / dx_den;
#if AMREX_SPACEDIM >= 2
                div_sigma(d) += (stress(i,jp,k)(d,1) - stress(i,jm,k)(d,1)) / dy_den;
#endif
#if AMREX_SPACEDIM == 3
                int km = k - 1 < cell_domain.smallEnd(2) ? cell_domain.smallEnd(2) : k - 1;
                int kp = k + 1 > cell_domain.bigEnd(2) ? cell_domain.bigEnd(2) : k + 1;
                Set::Scalar dz_den = lowmach_max(static_cast<Set::Scalar>(kp - km) * DX[2], DX[2]);
                div_sigma(d) += (stress(i,j,kp)(d,2) - stress(i,j,km)(d,2)) / dz_den;
#endif
                div_sigma(d) = lowmach_finite_or(div_sigma(d), 0.0);
            }

            for (int d = 0; d < AMREX_SPACEDIM; ++d)
            {
                u_rhs(i,j,k,d) = -vel.dot(grad_u.row(d)) - p_scale * grad_p(d) / density + gravity(d) + stress_sign * div_sigma(d) / density;
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

            auto advect_scalar = [=] AMREX_GPU_DEVICE (auto const& phi, int n) -> Set::Scalar
            {
                Set::Scalar phi_c = phi(i,j,k,n);

                Set::Scalar phi_xm = phi(i-1,j,k,n);
                Set::Scalar phi_xp = phi(i+1,j,k,n);
                Set::Scalar slope_xm = lowmach_mc_slope(phi(i-2,j,k,n), phi_xm, phi_c);
                Set::Scalar slope_x  = lowmach_mc_slope(phi_xm, phi_c, phi_xp);
                Set::Scalar slope_xp = lowmach_mc_slope(phi_c, phi_xp, phi(i+2,j,k,n));

                Set::Scalar ux_xlo = 0.5 * (u(i-1,j,k,0) + u(i,j,k,0));
                Set::Scalar ux_xhi = 0.5 * (u(i,j,k,0) + u(i+1,j,k,0));

                Set::Scalar phi_xlo_l = lowmach_clamp(phi_xm + 0.5 * slope_xm, phi_xm, phi_c);
                Set::Scalar phi_xlo_r = lowmach_clamp(phi_c  - 0.5 * slope_x,  phi_xm, phi_c);
                Set::Scalar phi_xhi_l = lowmach_clamp(phi_c  + 0.5 * slope_x,  phi_c, phi_xp);
                Set::Scalar phi_xhi_r = lowmach_clamp(phi_xp - 0.5 * slope_xp, phi_c, phi_xp);

                Set::Scalar flux_xlo = ux_xlo * (ux_xlo >= 0.0 ? phi_xlo_l : phi_xlo_r);
                Set::Scalar flux_xhi = ux_xhi * (ux_xhi >= 0.0 ? phi_xhi_l : phi_xhi_r);
                Set::Scalar div_phi_u = (flux_xhi - flux_xlo) / DX[0];
                Set::Scalar div_u = (ux_xhi - ux_xlo) / DX[0];

#if AMREX_SPACEDIM >= 2
                Set::Scalar phi_ym = phi(i,j-1,k,n);
                Set::Scalar phi_yp = phi(i,j+1,k,n);
                Set::Scalar slope_ym = lowmach_mc_slope(phi(i,j-2,k,n), phi_ym, phi_c);
                Set::Scalar slope_y  = lowmach_mc_slope(phi_ym, phi_c, phi_yp);
                Set::Scalar slope_yp = lowmach_mc_slope(phi_c, phi_yp, phi(i,j+2,k,n));

                Set::Scalar uy_ylo = 0.5 * (u(i,j-1,k,1) + u(i,j,k,1));
                Set::Scalar uy_yhi = 0.5 * (u(i,j,k,1) + u(i,j+1,k,1));

                Set::Scalar phi_ylo_l = lowmach_clamp(phi_ym + 0.5 * slope_ym, phi_ym, phi_c);
                Set::Scalar phi_ylo_r = lowmach_clamp(phi_c  - 0.5 * slope_y,  phi_ym, phi_c);
                Set::Scalar phi_yhi_l = lowmach_clamp(phi_c  + 0.5 * slope_y,  phi_c, phi_yp);
                Set::Scalar phi_yhi_r = lowmach_clamp(phi_yp - 0.5 * slope_yp, phi_c, phi_yp);

                Set::Scalar flux_ylo = uy_ylo * (uy_ylo >= 0.0 ? phi_ylo_l : phi_ylo_r);
                Set::Scalar flux_yhi = uy_yhi * (uy_yhi >= 0.0 ? phi_yhi_l : phi_yhi_r);
                div_phi_u += (flux_yhi - flux_ylo) / DX[1];
                div_u += (uy_yhi - uy_ylo) / DX[1];
#endif

#if AMREX_SPACEDIM == 3
                Set::Scalar phi_zm = phi(i,j,k-1,n);
                Set::Scalar phi_zp = phi(i,j,k+1,n);
                Set::Scalar slope_zm = lowmach_mc_slope(phi(i,j,k-2,n), phi_zm, phi_c);
                Set::Scalar slope_z  = lowmach_mc_slope(phi_zm, phi_c, phi_zp);
                Set::Scalar slope_zp = lowmach_mc_slope(phi_c, phi_zp, phi(i,j,k+2,n));

                Set::Scalar uz_zlo = 0.5 * (u(i,j,k-1,2) + u(i,j,k,2));
                Set::Scalar uz_zhi = 0.5 * (u(i,j,k,2) + u(i,j,k+1,2));

                Set::Scalar phi_zlo_l = lowmach_clamp(phi_zm + 0.5 * slope_zm, phi_zm, phi_c);
                Set::Scalar phi_zlo_r = lowmach_clamp(phi_c  - 0.5 * slope_z,  phi_zm, phi_c);
                Set::Scalar phi_zhi_l = lowmach_clamp(phi_c  + 0.5 * slope_z,  phi_c, phi_zp);
                Set::Scalar phi_zhi_r = lowmach_clamp(phi_zp - 0.5 * slope_zp, phi_c, phi_zp);

                Set::Scalar flux_zlo = uz_zlo * (uz_zlo >= 0.0 ? phi_zlo_l : phi_zlo_r);
                Set::Scalar flux_zhi = uz_zhi * (uz_zhi >= 0.0 ? phi_zhi_l : phi_zhi_r);
                div_phi_u += (flux_zhi - flux_zlo) / DX[2];
                div_u += (uz_zhi - uz_zlo) / DX[2];
#endif

                return -div_phi_u + phi_c * div_u;
            };

            eta_rhs(i,j,k) = advect_scalar(eta, 0);
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
                xi_rhs(i,j,k,d) = eta(i,j,k) >= xi_eta_cutoff ? advect_scalar(xi, d) : 0.0;
        });
    }
}

void
LowMach::Advance(int lev, Set::Scalar time, Set::Scalar dt)
{
    std::swap(velocity_old_mf[lev], velocity_mf[lev]);
    std::swap(temperature_old_mf[lev], temperature_mf[lev]);
    std::swap(mass_fraction_old_mf[lev], mass_fraction_mf[lev]);
    std::swap(eta_old_mf[lev], eta_mf[lev]);
    std::swap(xi_old_mf[lev], xi_mf[lev]);

    amrex::Vector<amrex::MultiFab> solution_new;
    solution_new.emplace_back(*velocity_mf[lev].get(), amrex::MakeType::make_alias, 0, AMREX_SPACEDIM);
    solution_new.emplace_back(*temperature_mf[lev].get(), amrex::MakeType::make_alias, 0, 1);
    solution_new.emplace_back(*mass_fraction_mf[lev].get(), amrex::MakeType::make_alias, 0, nspecies);
    solution_new.emplace_back(*eta_mf[lev].get(), amrex::MakeType::make_alias, 0, 1);
    solution_new.emplace_back(*xi_mf[lev].get(), amrex::MakeType::make_alias, 0, AMREX_SPACEDIM);

    amrex::Vector<amrex::MultiFab> solution_old;
    solution_old.emplace_back(*velocity_old_mf[lev].get(), amrex::MakeType::make_alias, 0, AMREX_SPACEDIM);
    solution_old.emplace_back(*temperature_old_mf[lev].get(), amrex::MakeType::make_alias, 0, 1);
    solution_old.emplace_back(*mass_fraction_old_mf[lev].get(), amrex::MakeType::make_alias, 0, nspecies);
    solution_old.emplace_back(*eta_old_mf[lev].get(), amrex::MakeType::make_alias, 0, 1);
    solution_old.emplace_back(*xi_old_mf[lev].get(), amrex::MakeType::make_alias, 0, AMREX_SPACEDIM);

    amrex::TimeIntegrator timeintegrator(solution_new, time);
    timeintegrator.set_rhs([&](amrex::Vector<amrex::MultiFab>& rhs_mf,
                               amrex::Vector<amrex::MultiFab>& state_mf,
                               const Set::Scalar rhs_time)
    {
        RHS(lev, rhs_time, rhs_mf[0], rhs_mf[1], rhs_mf[2], rhs_mf[3], rhs_mf[4],
            state_mf[0], state_mf[1], state_mf[2], state_mf[3], state_mf[4]);
    });

    timeintegrator.set_post_stage_action([&](amrex::Vector<amrex::MultiFab>& stage_mf, Set::Scalar stage_time)
    {
        EnforceStateBounds(stage_mf[0], stage_mf[1], stage_mf[2]);
        FillStateBoundaries(lev, stage_mf[0], stage_mf[1], stage_mf[2], stage_mf[3], stage_mf[4], stage_time);
        RebuildReferenceMapOutsideEta(lev, stage_mf[3], stage_mf[4], stage_time);
        UpdateDerived(lev, stage_mf[0], stage_mf[1], stage_mf[2]);
        UpdateSolidStress(lev, stage_mf[0], stage_mf[3], stage_mf[4]);
    });

    timeintegrator.advance(solution_old, solution_new, time, dt);
    EnforceStateBounds(*velocity_mf[lev], *temperature_mf[lev], *mass_fraction_mf[lev]);
    FillStateBoundaries(lev, *velocity_mf[lev], *temperature_mf[lev], *mass_fraction_mf[lev], *eta_mf[lev], *xi_mf[lev], time + dt);
    RebuildReferenceMapOutsideEta(lev, *eta_mf[lev], *xi_mf[lev], time + dt);
    UpdateDerived(lev, *velocity_mf[lev], *temperature_mf[lev], *mass_fraction_mf[lev]);
    UpdateSolidStress(lev, *velocity_mf[lev], *eta_mf[lev], *xi_mf[lev]);
    ImplicitElasticVelocitySolve(lev, time + dt, dt);
    UpdateDerived(lev, *velocity_mf[lev], *temperature_mf[lev], *mass_fraction_mf[lev]);
    UpdateSolidStress(lev, *velocity_mf[lev], *eta_mf[lev], *xi_mf[lev]);
    ProjectVelocity(lev, time + dt, dt);
    EnforceStateBounds(*velocity_mf[lev], *temperature_mf[lev], *mass_fraction_mf[lev]);
    FillStateBoundaries(lev, *velocity_mf[lev], *temperature_mf[lev], *mass_fraction_mf[lev], *eta_mf[lev], *xi_mf[lev], time + dt);
    RebuildReferenceMapOutsideEta(lev, *eta_mf[lev], *xi_mf[lev], time + dt);
    UpdateDerived(lev, *velocity_mf[lev], *temperature_mf[lev], *mass_fraction_mf[lev]);
    UpdateSolidStress(lev, *velocity_mf[lev], *eta_mf[lev], *xi_mf[lev]);
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
        Set::Patch<const Set::Scalar> eta = eta_mf.Patch(lev,mfi);
        const Set::Scalar vcrit = velocity_refinement_criterion;
        const Set::Scalar pcrit = pressure_refinement_criterion;
        const Set::Scalar Tcrit = temperature_refinement_criterion;
        const Set::Scalar etacrit = eta_refinement_criterion;
        amrex::Box domain = geom[lev].Domain();

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            auto sten = Numeric::GetStencil(i, j, k, domain);
            Set::Matrix grad_u = Numeric::Gradient(u, i, j, k, DX, sten);
            Set::Vector grad_p = Numeric::Gradient(pressure, i, j, k, 0, DX, sten);
            Set::Vector grad_T = Numeric::Gradient(T, i, j, k, 0, DX, sten);
            Set::Vector grad_eta = Numeric::Gradient(eta, i, j, k, 0, DX, sten);
            if (grad_u.norm() * dr > vcrit ||
                grad_p.lpNorm<2>() * dr > pcrit ||
                grad_T.lpNorm<2>() * dr > Tcrit ||
                grad_eta.lpNorm<2>() * dr * 2.0 > etacrit)
                tag(i,j,k) = amrex::TagBox::SET;
        });
    }
}

void
LowMach::Regrid(int lev, Set::Scalar time)
{
    if (lev < finest_level) return;
    for (int ilev = 0; ilev <= finest_level; ++ilev)
    {
        EnforceStateBounds(*velocity_mf[ilev], *temperature_mf[ilev], *mass_fraction_mf[ilev]);
        FillStateBoundaries(ilev, *velocity_mf[ilev], *temperature_mf[ilev], *mass_fraction_mf[ilev], *eta_mf[ilev], *xi_mf[ilev], 0.0);
        RebuildReferenceMapOutsideEta(ilev, *eta_mf[ilev], *xi_mf[ilev], time);
        UpdateDerived(ilev, *velocity_mf[ilev], *temperature_mf[ilev], *mass_fraction_mf[ilev]);
        UpdateSolidStress(ilev, *velocity_mf[ilev], *eta_mf[ilev], *xi_mf[ilev]);
    }
}
}
