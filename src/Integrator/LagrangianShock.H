#ifndef INTEGRATOR_LAGRANGIANSHOCK_H
#define INTEGRATOR_LAGRANGIANSHOCK_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>

#include "AMReX.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"
#include "AMReX_MultiFabUtil.H"

#include "IO/ParmParse.H"
#include "Integrator/Base/Mechanics.H"

#include "IC/IC.H"
#include "IC/Constant.H"
#include "IC/Expression.H"
#include "BC/BC.H"
#include "BC/Nothing.H"
#include "BC/Constant.H"
#include "BC/Operator/Elastic/Elastic.H"

#include "Set/Set.H"

namespace Integrator
{

struct MieGruneisenEOS
{
    Set::Scalar C0     = 3430.0;
    Set::Scalar s1     = 1.19;
    Set::Scalar Gamma0 = 1.67;
    Set::Scalar rho0   = 16580.0;
    Set::Scalar Cv     = 140.0;
    Set::Scalar J_min  = 0.8;   

    AMREX_GPU_HOST_DEVICE
    Set::Scalar HugoniotPressure(Set::Scalar rho) const
    {
        Set::Scalar eta = 1.0 - rho0 / rho;
        const Set::Scalar eta_max = 1.0 - J_min;
        if (eta > eta_max) eta = eta_max;          
        if (eta >= 0.0)
        {
            Set::Scalar denom = 1.0 - s1 * eta;
            if (denom < 0.01) denom = 0.01;
            return rho0 * C0 * C0 * eta / (denom * denom);
        }
        else
        {

            return rho0 * C0 * C0 * eta;
        }
    }

    AMREX_GPU_HOST_DEVICE
    Set::Scalar HugoniotEnergy(Set::Scalar rho) const
    {
        Set::Scalar eta = 1.0 - rho0 / rho;
        Set::Scalar pH  = HugoniotPressure(rho);
        if (std::abs(eta) < 1e-12) return 0.0;
        return 0.5 * pH * (1.0 / rho0 - 1.0 / rho);
    }

    AMREX_GPU_HOST_DEVICE
    Set::Scalar Pressure(Set::Scalar rho, Set::Scalar e) const
    {
        Set::Scalar pH = HugoniotPressure(rho);
        Set::Scalar eH = HugoniotEnergy(rho);
        return pH + Gamma0 * rho0 * (e - eH);
    }

    AMREX_GPU_HOST_DEVICE
    Set::Scalar SoundSpeed(Set::Scalar rho, Set::Scalar p) const
    {
        Set::Scalar rho_safe  = std::max(rho, 0.1 * rho0);
        Set::Scalar eta       = 1.0 - rho0 / rho_safe;
        const Set::Scalar eta_max = 1.0 - J_min;
        if (eta > eta_max) eta = eta_max;          
        Set::Scalar deta_drho = rho0 / (rho_safe * rho_safe);

        // dp_H/drho : frozen-energy Hugoniot pressure derivative
        Set::Scalar dpH_drho;
        if (eta >= 0.0)
        {
            Set::Scalar denom    = 1.0 - s1 * eta;
            if (denom < 0.01) denom = 0.01;
            Set::Scalar dpH_deta = rho0 * C0 * C0 * (1.0 + s1 * eta)
                                 / (denom * denom * denom);
            dpH_drho = dpH_deta * deta_drho;
        }
        else
        {

            dpH_drho = rho0 * C0 * C0 * deta_drho;
        }

        // de_H/drho with  e_H = 1/2 p_H (1/rho0 - 1/rho)
        //   de_H/drho = 1/2 [ p_H'(rho) (1/rho0 - 1/rho) + p_H/rho^2 ]
        Set::Scalar pH       = HugoniotPressure(rho_safe);
        Set::Scalar deH_drho = 0.5 * ( dpH_drho * (1.0 / rho0 - 1.0 / rho_safe)
                                     + pH / (rho_safe * rho_safe) );

        // (dp/drho)_e = dp_H/drho - Gamma0 rho0 de_H/drho
        Set::Scalar dp_drho_e = dpH_drho - Gamma0 * rho0 * deH_drho;

        // c^2 = (dp/drho)_e + (p/rho^2)(dp/de)_rho,  with (dp/de)_rho = Gamma0 rho0.
        
        Set::Scalar c2 = dp_drho_e + Gamma0 * rho0 * p / (rho_safe * rho_safe);
        c2 = std::max(c2, 0.25 * C0 * C0);
        return std::sqrt(c2);
    }

    static void Parse(MieGruneisenEOS& value, IO::ParmParse& pp)
    {
        pp_query_default("C0",     value.C0,     3430.0);
        pp_query_default("s1",     value.s1,     1.19);
        pp_query_default("Gamma0", value.Gamma0, 1.67);
        pp_query_default("rho0",   value.rho0,   16580.0);
        pp_query_default("Cv",     value.Cv,     140.0);
        pp_query_default("J_min",  value.J_min,  0.8);   
    }
};

template<class MODEL>
class LagrangianShock : virtual public Base::Mechanics<MODEL>
{
public:
    static constexpr const char* name = "lagrangian_shock";

    LagrangianShock() : Base::Mechanics<MODEL>() {}
    LagrangianShock(IO::ParmParse& pp) : Base::Mechanics<MODEL>()
    { Parse(*this, pp); }

    ~LagrangianShock() {}

    static void Parse(LagrangianShock& value, IO::ParmParse& pp)
    {
        BL_PROFILE("Integrator::LagrangianShock::Parse");

        value.m_type = Base::Mechanics<MODEL>::Type::Dynamic;

        pp.queryclass<Base::Mechanics<MODEL>>(value);

        pp_queryclass("eos", value.eos);
        pp_query_default("hydro.rho0", value.hydro_rho0, value.eos.rho0);

        pp_query_default("av.CQ", value.av_CQ,     2.0);
        pp_query_default("av.CL", value.av_CL,     0.5);
        pp_query_default("cfl",   value.hydro_cfl, 0.4);

        pp_query_default("e_init", value.e_init, 0.0);

        pp_query_default("diag_file", value.diag_filename, std::string("lagshock_diag.dat"));
        pp_query_default("diag_dt",   value.diag_dt,       1e-10);

        Util::Message(INFO, "EOS: C0=", value.eos.C0, " s1=", value.eos.s1,
                      " Gamma0=", value.eos.Gamma0, " rho0=", value.eos.rho0);

        int nghost = 2;
        value.RegisterNewFab(value.density_mf,    &value.bc_nothing, 1, nghost, "density",    true,  true);
        value.RegisterNewFab(value.pressure_mf,   &value.bc_nothing, 1, nghost, "pressure",   true,  false);
        value.RegisterNewFab(value.energy_mf,     &value.bc_nothing, 1, nghost, "energy",     true,  true);
        value.RegisterNewFab(value.artvisc_mf,    &value.bc_nothing, 1, nghost, "artvisc",    true,  false);
        value.RegisterNewFab(value.soundspeed_mf, &value.bc_nothing, 1, nghost, "soundspeed", true,  false);

        value.RegisterNewFab(value.density_half_mf,    &value.bc_nothing, 1, nghost, "density_half",    false, false);
        value.RegisterNewFab(value.pressure_half_mf,   &value.bc_nothing, 1, nghost, "pressure_half",   false, false);
        value.RegisterNewFab(value.energy_half_mf,     &value.bc_nothing, 1, nghost, "energy_half",     false, false);
        value.RegisterNewFab(value.artvisc_half_mf,    &value.bc_nothing, 1, nghost, "artvisc_half",    false, false);
        value.RegisterNewFab(value.soundspeed_half_mf, &value.bc_nothing, 1, nghost, "soundspeed_half", false, false);

        value.RegisterIntegratedVariable(&value.diag_fs_vel,       "fs_velocity",  false);
        value.RegisterIntegratedVariable(&value.diag_max_pressure, "max_pressure", false);
        value.RegisterIntegratedVariable(&value.diag_min_pressure, "min_pressure", false);

        value.rho = value.hydro_rho0;
    }

    void Initialize(int lev) override
    {
        BL_PROFILE("Integrator::LagrangianShock::Initialize");

        Base::Mechanics<MODEL>::Initialize(lev);

        density_mf[lev]->setVal(hydro_rho0);
        pressure_mf[lev]->setVal(0.0);
        energy_mf[lev]->setVal(e_init);
        artvisc_mf[lev]->setVal(0.0);
        soundspeed_mf[lev]->setVal(eos.C0);

        density_half_mf[lev]->setVal(hydro_rho0);
        pressure_half_mf[lev]->setVal(0.0);
        energy_half_mf[lev]->setVal(e_init);
        artvisc_half_mf[lev]->setVal(0.0);
        soundspeed_half_mf[lev]->setVal(eos.C0);

        density_mf[lev]->FillBoundary(geom[lev].periodicity());
        energy_mf[lev]->FillBoundary(geom[lev].periodicity());
        vel_mf[lev]->FillBoundary(geom[lev].periodicity());
    }

    void UpdateModel(int , Set::Scalar ) override {}

    void Advance(int lev, Set::Scalar time, Set::Scalar a_dt) override
    {
        BL_PROFILE("Integrator::LagrangianShock::Advance");
        if (m_type == Base::Mechanics<MODEL>::Type::Disable) return;

        const amrex::Real* DX = geom[lev].CellSize();
        Set::Scalar DX0 = DX[0];

        amrex::Box cdomain = geom[lev].Domain();
        amrex::Box ndomain = amrex::convert(cdomain, amrex::IntVect::TheNodeVector());
        const amrex::Dim3 nlo = amrex::lbound(ndomain);
        const amrex::Dim3 nhi = amrex::ubound(ndomain);

        MieGruneisenEOS local_eos        = eos;
        Set::Scalar local_rho0           = hydro_rho0;
        Set::Scalar local_CQ             = av_CQ;
        Set::Scalar local_CL             = av_CL;

        Set::Scalar half_dt  = 0.5 * a_dt;
        Set::Scalar t_full   = time + a_dt;

        std::swap(*disp_mf[lev], *disp_old_mf[lev]);
        std::swap(*vel_mf[lev],  *vel_old_mf[lev]);

        vel_old_mf[lev]->FillBoundary(geom[lev].periodicity());
        disp_old_mf[lev]->FillBoundary(geom[lev].periodicity());

        ComputePressureAndViscosity(lev, DX0, local_eos, local_rho0,
                                    local_CQ, local_CL,
                                    disp_old_mf, vel_old_mf);

        pressure_mf[lev]->FillBoundary(geom[lev].periodicity());
        artvisc_mf[lev]->FillBoundary(geom[lev].periodicity());
        soundspeed_mf[lev]->FillBoundary(geom[lev].periodicity());
        density_mf[lev]->FillBoundary(geom[lev].periodicity());

        UpdateVelocityAndDisplacement(lev, DX0, half_dt, nlo, nhi, ndomain,
                                       local_rho0,
                                       vel_old_mf, disp_old_mf,
                                       vel_old_mf,
                                       pressure_mf, artvisc_mf,
                                       vel_mf, disp_mf);

        vel_mf[lev]->FillBoundary(geom[lev].periodicity());
        disp_mf[lev]->FillBoundary(geom[lev].periodicity());

        ComputeHalfStepEnergy(lev, DX0, local_eos, local_rho0,
                               local_CQ, local_CL);

        pressure_half_mf[lev]->FillBoundary(geom[lev].periodicity());
        artvisc_half_mf[lev]->FillBoundary(geom[lev].periodicity());
        soundspeed_half_mf[lev]->FillBoundary(geom[lev].periodicity());
        density_half_mf[lev]->FillBoundary(geom[lev].periodicity());

        UpdateVelocityAndDisplacement(lev, DX0, a_dt, nlo, nhi, ndomain,
                                       local_rho0,
                                       vel_old_mf, disp_old_mf,
                                       vel_mf,
                                       pressure_half_mf, artvisc_half_mf,
                                       vel_mf, disp_mf);

        vel_mf[lev]->FillBoundary(geom[lev].periodicity());
        disp_mf[lev]->FillBoundary(geom[lev].periodicity());

        for (amrex::MFIter mfi(*density_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box cbx = mfi.tilebox();
            amrex::Array4<Set::Scalar>       const& rho_arr = density_mf[lev]->array(mfi);
            amrex::Array4<const Set::Vector> const& u_arr   = disp_mf[lev]->array(mfi);
            amrex::ParallelFor(cbx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Scalar dx_new = DX0 + (u_arr(i+1,j,k)(0) - u_arr(i,j,k)(0));
                if (dx_new < 0.01 * DX0) dx_new = 0.01 * DX0;
                rho_arr(i,j,k) = local_rho0 * DX0 / dx_new;
            });
        }

        for (amrex::MFIter mfi(*energy_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box cbx = mfi.tilebox();
            amrex::Array4<Set::Scalar>       const& e_arr     = energy_mf[lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const& ph_arr    = pressure_half_mf[lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const& qh_arr    = artvisc_half_mf[lev]->array(mfi);
            amrex::Array4<const Set::Vector> const& u_old_arr = disp_old_mf[lev]->array(mfi);
            amrex::Array4<const Set::Vector> const& u_new_arr = disp_mf[lev]->array(mfi);

            amrex::ParallelFor(cbx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Scalar dx_old = DX0 + (u_old_arr(i+1,j,k)(0) - u_old_arr(i,j,k)(0));
                Set::Scalar dx_new = DX0 + (u_new_arr(i+1,j,k)(0) - u_new_arr(i,j,k)(0));
                if (dx_old < 0.01 * DX0) dx_old = 0.01 * DX0;
                if (dx_new < 0.01 * DX0) dx_new = 0.01 * DX0;

                Set::Scalar dV        = dx_new - dx_old;
                Set::Scalar cell_mass = local_rho0 * DX0;
                Set::Scalar ptot_half = ph_arr(i,j,k) + qh_arr(i,j,k);
                Set::Scalar de        = -ptot_half * dV / cell_mass;
                if (std::isfinite(de)) e_arr(i,j,k) += de;
            });
        }
        energy_mf[lev]->FillBoundary(geom[lev].periodicity());

        for (amrex::MFIter mfi(*disp_mf[lev], false); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.grownnodaltilebox() & ndomain;
            amrex::Array4<Set::Matrix>       const& eps     = (*strain_mf[lev]).array(mfi);
            amrex::Array4<Set::Matrix>       const& sig     = (*stress_mf[lev]).array(mfi);
            amrex::Array4<MODEL>             const& model   = (*model_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& p_arr   = pressure_mf[lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const& q_arr   = artvisc_mf[lev]->array(mfi);
            amrex::Array4<const Set::Vector> const& u_arr   = disp_mf[lev]->array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                eps(i,j,k) = Set::Matrix::Zero();
                sig(i,j,k) = Set::Matrix::Zero();

                int i_left  = (i > nlo.x) ? i - 1 : i;
                int i_right = (i < nhi.x) ? i     : i - 1;

                Set::Scalar p_node = 0.5 * (p_arr(i_left, j, k) + p_arr(i_right, j, k));
                Set::Scalar q_node = 0.5 * (q_arr(i_left, j, k) + q_arr(i_right, j, k));

                Set::Scalar du_left  = (i > nlo.x) ?
                    (u_arr(i,j,k)(0) - u_arr(i-1,j,k)(0)) : 0.0;
                Set::Scalar du_right = (i < nhi.x) ?
                    (u_arr(i+1,j,k)(0) - u_arr(i,j,k)(0)) : 0.0;
                Set::Scalar eps_xx;
                if (i > nlo.x && i < nhi.x)
                    eps_xx = 0.5 * (du_left + du_right) / DX0;
                else
                    eps_xx = (du_left + du_right) / DX0;

                eps(i,j,k)(0,0) = eps_xx;
                sig(i,j,k)(0,0) = -(p_node + q_node);

                model(i, j, k).Advance(a_dt, eps(i, j, k), sig(i, j, k), time);
            });
        }
        model_mf[lev]->setMultiGhost(true);
        model_mf[lev]->FillBoundaryAndSync(geom[lev].periodicity());

        Set::Scalar dt_local = 1e100;
        Set::Scalar max_p    = -1e100;
        Set::Scalar min_p    =  1e100;
        Set::Scalar fs_vel   = 0.0;

        for (amrex::MFIter mfi(*vel_mf[lev], false); mfi.isValid(); ++mfi)
        {
            amrex::Box nbx = mfi.nodaltilebox() & ndomain;
            amrex::Array4<const Set::Vector> const& v_arr = vel_mf[lev]->array(mfi);
            if (nbx.contains(amrex::IntVect(AMREX_D_DECL(nhi.x, nlo.y, 0))))
                fs_vel = v_arr(nhi.x, nlo.y, 0)(0);
        }

        ComputePressureAndViscosity(lev, DX0, local_eos, local_rho0,
                                    local_CQ, local_CL,
                                    disp_mf, vel_mf);

        for (amrex::MFIter mfi(*density_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box cbx = mfi.tilebox();
            amrex::Array4<const Set::Scalar> const& p_arr = pressure_mf[lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const& c_arr = soundspeed_mf[lev]->array(mfi);
            amrex::Array4<const Set::Vector> const& u_arr = disp_mf[lev]->array(mfi);
            amrex::Array4<const Set::Vector> const& v_arr = vel_mf[lev]->array(mfi);

            amrex::Loop(cbx, [&](int i, int j, int k)
            {
                Set::Scalar c      = c_arr(i,j,k);
                Set::Scalar dx_cur = DX0 + (u_arr(i+1,j,k)(0) - u_arr(i,j,k)(0));
                if (dx_cur < 0.01 * DX0) dx_cur = 0.01 * DX0;
                Set::Scalar dv     = std::abs(v_arr(i+1,j,k)(0) - v_arr(i,j,k)(0));
                Set::Scalar dt_cfl = hydro_cfl * dx_cur / (c + local_CQ * dv + 1e-20);
                dt_local = std::min(dt_local, dt_cfl);
                if (p_arr(i,j,k) > max_p) max_p = p_arr(i,j,k);
                if (p_arr(i,j,k) < min_p) min_p = p_arr(i,j,k);
            });
        }

        
        amrex::ParallelDescriptor::ReduceRealMax(max_p);
        amrex::ParallelDescriptor::ReduceRealMin(min_p);
        amrex::ParallelDescriptor::ReduceRealSum(fs_vel);

        diag_fs_vel       = fs_vel;
        diag_max_pressure = max_p;
        diag_min_pressure = min_p;

 
        if (dt_local > 0.0 && dt_local < 1e90)
            this->DynamicTimestep_SyncTimeStep(lev, dt_local);

        if (lev == 0 && amrex::ParallelDescriptor::IOProcessor()
            && t_full >= next_diag_time)
        {
            WriteDiagnostics(t_full, a_dt);
            next_diag_time = t_full + diag_dt;
        }
    }

    void TagCellsForRefinement(int, amrex::TagBoxArray&, Set::Scalar, int) override {}
    void Integrate(int, Set::Scalar, int, const amrex::MFIter&, const amrex::Box&) override {}

protected:

    void ComputePressureAndViscosity(
        int lev, Set::Scalar DX0,
        MieGruneisenEOS local_eos, Set::Scalar local_rho0,
        Set::Scalar local_CQ, Set::Scalar local_CL,
        Set::Field<Set::Vector>& u_field, Set::Field<Set::Vector>& v_field)
    {
        for (amrex::MFIter mfi(*density_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box cbx = mfi.tilebox();
            amrex::Array4<Set::Scalar>       const& rho_arr = density_mf[lev]->array(mfi);
            amrex::Array4<Set::Scalar>       const& p_arr   = pressure_mf[lev]->array(mfi);
            amrex::Array4<Set::Scalar>       const& e_arr   = energy_mf[lev]->array(mfi);
            amrex::Array4<Set::Scalar>       const& q_arr   = artvisc_mf[lev]->array(mfi);
            amrex::Array4<Set::Scalar>       const& c_arr   = soundspeed_mf[lev]->array(mfi);
            amrex::Array4<const Set::Vector> const& u_arr   = u_field[lev]->array(mfi);
            amrex::Array4<const Set::Vector> const& v_arr   = v_field[lev]->array(mfi);

            amrex::ParallelFor(cbx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Scalar dx_cur = DX0 + (u_arr(i+1,j,k)(0) - u_arr(i,j,k)(0));
                if (dx_cur < 0.01 * DX0) dx_cur = 0.01 * DX0;
                Set::Scalar rho = local_rho0 * DX0 / dx_cur;
                rho_arr(i,j,k) = rho;

                Set::Scalar p = local_eos.Pressure(rho, e_arr(i,j,k));
                if (!std::isfinite(p)) p = 0.0;

                p_arr(i,j,k) = p;
                c_arr(i,j,k) = local_eos.SoundSpeed(rho, p);

                Set::Scalar dv = v_arr(i+1,j,k)(0) - v_arr(i,j,k)(0);
                Set::Scalar q  = 0.0;
                if (dv < 0.0)
                    q = local_CQ * rho * dv * dv
                      + local_CL * rho * c_arr(i,j,k) * std::abs(dv);
                q_arr(i,j,k) = q;
            });
        }
    }

    void UpdateVelocityAndDisplacement(
        int lev, Set::Scalar DX0, Set::Scalar dt,
        amrex::Dim3 nlo, amrex::Dim3 nhi, amrex::Box ndomain,
        Set::Scalar local_rho0,
        Set::Field<Set::Vector>& v_in,  Set::Field<Set::Vector>& u_in,
        Set::Field<Set::Vector>& v_drift,
        Set::Field<Set::Scalar>& p_field, Set::Field<Set::Scalar>& q_field,
        Set::Field<Set::Vector>& v_out, Set::Field<Set::Vector>& u_out)
    {

        auto bc_ptr = bc;

        for (amrex::MFIter mfi(*v_out[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box nbx = mfi.nodaltilebox() & ndomain;

            amrex::Array4<const Set::Vector> const& v_old   = v_in[lev]->array(mfi);
            amrex::Array4<const Set::Vector> const& u_old   = u_in[lev]->array(mfi);
            amrex::Array4<const Set::Vector> const& v_drf   = v_drift[lev]->array(mfi);
            amrex::Array4<Set::Vector>       const& v_new   = v_out[lev]->array(mfi);
            amrex::Array4<Set::Vector>       const& u_new   = u_out[lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const& p_arr   = p_field[lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const& q_arr   = q_field[lev]->array(mfi);
            amrex::Array4<const Set::Vector> const& b_arr   = rhs_mf[lev]->array(mfi);

            amrex::ParallelFor(nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Scalar v_i = v_old(i,j,k)(0);
                Set::Scalar u_i = u_old(i,j,k)(0);
                Set::Scalar vd_i = v_drf(i,j,k)(0);
                Set::Vector vnew_vec = v_old(i,j,k);
                Set::Vector unew_vec = u_old(i,j,k);

                bool AMREX_D_DECL(xmin = (i == nlo.x), ymin = (j == nlo.y), zmin = (k == nlo.z));
                bool AMREX_D_DECL(xmax = (i == nhi.x), ymax = (j == nhi.y), zmax = (k == nhi.z));

                if (AMREX_D_TERM(xmax || xmin, || ymax || ymin, || zmax || zmin))
                {

                    auto bctype = bc_ptr->getType(i, j, k, ndomain);

                    if (bctype[0] == BC::Operator::Elastic::Elastic::Type::Displacement)
                    {

                        unew_vec(0) = b_arr(i,j,k)(0);
                        vnew_vec(0) = (unew_vec(0) - u_i) / dt;
                    }
                    else
                    {

                        Set::Scalar half_mass = 0.5 * local_rho0 * DX0;

                        if (xmin)
                        {
                            Set::Scalar ptot = p_arr(i, j, k) + q_arr(i, j, k);
                            Set::Scalar accel = -ptot / half_mass;
                            vnew_vec(0) = v_i + dt * accel;
                            unew_vec(0) = u_i + dt * vd_i;
                        }
                        else if (xmax)
                        {
                            Set::Scalar ptot = p_arr(i-1, j, k) + q_arr(i-1, j, k);
                            Set::Scalar accel = ptot / half_mass;
                            vnew_vec(0) = v_i + dt * accel;
                            unew_vec(0) = u_i + dt * vd_i;
                        }
                        else
                        {

                            Set::Scalar ptot_right = p_arr(i,   j, k) + q_arr(i,   j, k);
                            Set::Scalar ptot_left  = p_arr(i-1, j, k) + q_arr(i-1, j, k);
                            Set::Scalar accel = -(ptot_right - ptot_left)
                                               / (local_rho0 * DX0);
                            vnew_vec(0) = v_i + dt * accel;
                            unew_vec(0) = u_i + dt * vd_i;
                        }
                    }

                    v_new(i,j,k) = vnew_vec;
                    u_new(i,j,k) = unew_vec;
                    return;
                }

                {
                    Set::Scalar ptot_right = p_arr(i,   j, k) + q_arr(i,   j, k);
                    Set::Scalar ptot_left  = p_arr(i-1, j, k) + q_arr(i-1, j, k);
                    Set::Scalar accel      = -(ptot_right - ptot_left)
                                           / (local_rho0 * DX0);
                    vnew_vec(0) = v_i + dt * accel;
                    unew_vec(0) = u_i + dt * vd_i;
                }
                v_new(i,j,k) = vnew_vec;
                u_new(i,j,k) = unew_vec;
            });
        }
    }

    void ComputeHalfStepEnergy(
        int lev, Set::Scalar DX0,
        MieGruneisenEOS local_eos, Set::Scalar local_rho0,
        Set::Scalar local_CQ, Set::Scalar local_CL)
    {
        for (amrex::MFIter mfi(*density_half_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box cbx = mfi.tilebox();
            amrex::Array4<Set::Scalar>       const& rho_h  = density_half_mf[lev]->array(mfi);
            amrex::Array4<Set::Scalar>       const& p_h    = pressure_half_mf[lev]->array(mfi);
            amrex::Array4<Set::Scalar>       const& e_h    = energy_half_mf[lev]->array(mfi);
            amrex::Array4<Set::Scalar>       const& q_h    = artvisc_half_mf[lev]->array(mfi);
            amrex::Array4<Set::Scalar>       const& c_h    = soundspeed_half_mf[lev]->array(mfi);
            amrex::Array4<const Set::Vector> const& u_old  = disp_old_mf[lev]->array(mfi);
            amrex::Array4<const Set::Vector> const& u_half = disp_mf[lev]->array(mfi);
            amrex::Array4<const Set::Vector> const& v_half = vel_mf[lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const& p_old  = pressure_mf[lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const& q_old  = artvisc_mf[lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const& e_old  = energy_mf[lev]->array(mfi);

            amrex::ParallelFor(cbx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Scalar dx_half = DX0 + (u_half(i+1,j,k)(0) - u_half(i,j,k)(0));
                Set::Scalar dx_old  = DX0 + (u_old(i+1,j,k)(0)  - u_old(i,j,k)(0));
                if (dx_half < 0.01 * DX0) dx_half = 0.01 * DX0;
                if (dx_old  < 0.01 * DX0) dx_old  = 0.01 * DX0;

                Set::Scalar rho = local_rho0 * DX0 / dx_half;
                rho_h(i,j,k) = rho;

                Set::Scalar dV        = dx_half - dx_old;
                Set::Scalar cell_mass = local_rho0 * DX0;
                Set::Scalar ptot_old  = p_old(i,j,k) + q_old(i,j,k);

                Set::Scalar e_new = e_old(i,j,k) - ptot_old * dV / cell_mass;
                if (!std::isfinite(e_new)) e_new = e_old(i,j,k);
                e_h(i,j,k) = e_new;

                Set::Scalar p = local_eos.Pressure(rho, e_new);
                if (!std::isfinite(p)) p = 0.0;
                p_h(i,j,k) = p;
                c_h(i,j,k) = local_eos.SoundSpeed(rho, p);

                Set::Scalar dv = v_half(i+1,j,k)(0) - v_half(i,j,k)(0);
                Set::Scalar q  = 0.0;
                if (dv < 0.0)
                    q = local_CQ * rho * dv * dv
                      + local_CL * rho * c_h(i,j,k) * std::abs(dv);
                q_h(i,j,k) = q;
            });
        }
    }

    void WriteDiagnostics(Set::Scalar time, Set::Scalar dt)
    {
        if (!amrex::ParallelDescriptor::IOProcessor()) return;
        std::ofstream out;
        if (!diag_header_written)
        {
            out.open(diag_filename, std::ios::out);
            out << "# Lagrangian Shock Diagnostics\n";
            out << "# Jones et al. (2018) Ta validation\n";
            out << "# Col 1: time (s)\n";
            out << "# Col 2: free-surface velocity (m/s)\n";
            out << "# Col 3: max pressure (Pa)\n";
            out << "# Col 4: min pressure (Pa)\n";
            out << "# Col 5: timestep (s)\n";
            out << std::scientific << std::setprecision(8);
            diag_header_written = true;
        }
        else
        {
            out.open(diag_filename, std::ios::app);
            out << std::scientific << std::setprecision(8);
        }
        out << time              << "  "
            << diag_fs_vel       << "  "
            << diag_max_pressure << "  "
            << diag_min_pressure << "  "
            << dt                << "\n";
        out.close();
    }

    MieGruneisenEOS eos;
    Set::Scalar hydro_rho0 = 16580.0;
    Set::Scalar e_init     = 0.0;

    Set::Scalar av_CQ     = 2.0;
    Set::Scalar av_CL     = 0.5;
    Set::Scalar hydro_cfl = 0.4;

    Set::Field<Set::Scalar> density_mf;
    Set::Field<Set::Scalar> pressure_mf;
    Set::Field<Set::Scalar> energy_mf;
    Set::Field<Set::Scalar> artvisc_mf;
    Set::Field<Set::Scalar> soundspeed_mf;

    Set::Field<Set::Scalar> density_half_mf;
    Set::Field<Set::Scalar> pressure_half_mf;
    Set::Field<Set::Scalar> energy_half_mf;
    Set::Field<Set::Scalar> artvisc_half_mf;
    Set::Field<Set::Scalar> soundspeed_half_mf;

    BC::Nothing bc_nothing;

    Set::Scalar diag_fs_vel        = 0.0;
    Set::Scalar diag_max_pressure  = 0.0;
    Set::Scalar diag_min_pressure  = 0.0;
    std::string diag_filename      = "lagshock_diag.dat";
    Set::Scalar diag_dt            = 1e-10;
    Set::Scalar next_diag_time     = 0.0;
    bool        diag_header_written = false;

    using Base::Mechanics<MODEL>::m_type;
    using Base::Mechanics<MODEL>::finest_level;
    using Base::Mechanics<MODEL>::geom;
    using Base::Mechanics<MODEL>::disp_mf;
    using Base::Mechanics<MODEL>::vel_mf;
    using Base::Mechanics<MODEL>::disp_old_mf;
    using Base::Mechanics<MODEL>::vel_old_mf;
    using Base::Mechanics<MODEL>::stress_mf;
    using Base::Mechanics<MODEL>::strain_mf;
    using Base::Mechanics<MODEL>::model_mf;
    using Base::Mechanics<MODEL>::rhs_mf;
    using Base::Mechanics<MODEL>::rho;
    using Base::Mechanics<MODEL>::timestep;
    using Base::Mechanics<MODEL>::bc;
};

}

#endif