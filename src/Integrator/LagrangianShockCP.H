#ifndef INTEGRATOR_LAGRANGIANSHOCKCP_H
#define INTEGRATOR_LAGRANGIANSHOCKCP_H

#include "Integrator/LagrangianShock.H"
#include "IC/PNG3.H"
#include "Numeric/Stencil.H"
#include <type_traits>
#include <utility>
namespace Integrator
{

template<class, class = void>
struct has_wp_inc : std::false_type {};
template<class M>
struct has_wp_inc<M, std::void_t<decltype(std::declval<M>().wp_inc)>>
    : std::true_type {};

template<class MODEL>
class LagrangianShockCP : public LagrangianShock<MODEL>
{
public:
    static constexpr const char* name = "lagrangian_shock_cp";
    LagrangianShockCP() : LagrangianShock<MODEL>() {}
    LagrangianShockCP(IO::ParmParse& pp) : LagrangianShock<MODEL>()
    { Parse(*this, pp); }
    ~LagrangianShockCP() {}

    static void Parse(LagrangianShockCP& value, IO::ParmParse& pp)
    {
        BL_PROFILE("Integrator::LagrangianShockCP::Parse");
        LagrangianShock<MODEL>::Parse(value, pp);
        pp_query_default("taylor_quinney", value.taylor_quinney, 0.9);
        pp_query_default("theta_init", value.theta_init, 294.0);
        int nghost = 2;
        value.RegisterNewFab(value.temperature_mf, &value.bc_nothing, 1, nghost, "temperature", true, true);
        value.RegisterNewFab(value.sigma_cp_mf, &value.bc_nothing, 1, nghost, "sigma_cp", true, false);

        value.RegisterNewFab(value.artvisc2d_mf,      &value.bc_nothing, 3, nghost, "artvisc2d",      false, false);
        value.RegisterNewFab(value.artvisc2d_half_mf, &value.bc_nothing, 3, nghost, "artvisc2d_half", false, false);

        value.RegisterNewFab(value.cell_stress_mf,     &value.bc_nothing, 4, nghost, "cell_stress",    false, false);

        value.RegisterNewFab(value.cell_stress_raw_mf, &value.bc_nothing, 4, nghost, "cell_stress_raw", false, false);

        pp_query_default("stress_scale", value.stress_scale, 1e9);

        pp_query_default("cp_substeps", value.cp_substeps_max, 8);

        pp_query_default("hg_coeff", value.hg_coeff, 0.1);

        Util::Message(INFO, "LagrangianShockCP: TQ=", value.taylor_quinney,
                      " theta0=", value.theta_init, " stress_scale=", value.stress_scale,
                      " cp_substeps=", value.cp_substeps_max,
                      " hg_coeff=", value.hg_coeff);

        if (pp.contains("euler.ic.type"))
        {
            value.RegisterNodalFab(value.euler_mf, 3, 3, "euler", true);
            std::string type;
            pp_query("euler.ic.type", type);
            if (type == "eulerangles")
            {
                value.ic_euler = new IC::PNG3(value.geom, pp, "euler.ic.eulerangles");
                value.use_euler_image = true;
            }
            else Util::Abort(INFO, "Invalid euler.ic.type: ", type);
            Util::Message(INFO, "LagrangianShockCP: reading crystal orientations from PNG");
        }
    }

    void Initialize(int lev) override
    {
        LagrangianShock<MODEL>::Initialize(lev);
        temperature_mf[lev]->setVal(theta_init);
        sigma_cp_mf[lev]->setVal(0.0);
        artvisc2d_mf[lev]->setVal(0.0);
        artvisc2d_half_mf[lev]->setVal(0.0);
        cell_stress_mf[lev]->setVal(0.0);
        cell_stress_raw_mf[lev]->setVal(0.0);

        Util::RealFillBoundary(*model_mf[lev], geom[lev]);

        bool needs_manual_init = false;
        for (amrex::MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi) {
            auto const& m = (*model_mf[lev]).array(mfi);
            amrex::Box bx = mfi.validbox();
            auto lo = bx.smallEnd();
            if constexpr (MODEL::kinvar == Model::Solid::KinematicVariable::F) {
                if (m(lo[0],lo[1],0).Fp(0,0) == 0.0) needs_manual_init = true;
            }
            break;
        }

        if (needs_manual_init) {
            MODEL model_template;
            IO::ParmParse pp("model1");
            MODEL::Parse(model_template, pp);
            for (amrex::MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi) {
                auto const& marr = (*model_mf[lev]).array(mfi);
                amrex::Box bx = mfi.fabbox();
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    marr(i,j,k) = model_template;
                });
            }
            Util::Message(INFO, "LagrangianShockCP: manually initialized model_mf from model1.*");
            Util::RealFillBoundary(*model_mf[lev], geom[lev]);
        }

        if (use_euler_image && ic_euler)
        {
            euler_mf[lev]->setVal(0.0);
            ic_euler->Initialize(lev, euler_mf);
            amrex::Box domain = geom[lev].Domain();
            domain.convert(amrex::IntVect::TheNodeVector());
            for (amrex::MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Array4<MODEL> const& model = (*model_mf[lev]).array(mfi);
                amrex::Array4<const Set::Scalar> const& euler = euler_mf[lev]->array(mfi);
                amrex::Box bx = mfi.grownnodaltilebox() & domain;
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    if constexpr (MODEL::kinvar == Model::Solid::KinematicVariable::F)
                    {
                        Set::Scalar qx = euler(i, j, k, 0);
                        Set::Scalar qy = euler(i, j, k, 1);
                        Set::Scalar qz = euler(i, j, k, 2);
                        Set::Scalar sum_sq = qx*qx + qy*qy + qz*qz;
                        Set::Scalar qw = (sum_sq < 1.0) ? std::sqrt(1.0 - sum_sq) : 0.0;
                        Set::Scalar qnorm2 = qw*qw + sum_sq;
                        if (qnorm2 < 1e-12 || !std::isfinite(qnorm2))
                        { qw = 1.0; qx = 0.0; qy = 0.0; qz = 0.0; }
                        Eigen::Quaterniond q(qw, qx, qy, qz);
                        q.normalize();
                        model(i, j, k).q = q;
                    }
                });
            }
            Util::RealFillBoundary(*model_mf[lev], geom[lev]);
            Util::Message(INFO, "LagrangianShockCP: stamped euler orientations on lev ", lev);
        }
    }

    void UpdateModel(int , Set::Scalar ) override
    {
        if (m_type == Base::Mechanics<MODEL>::Type::Disable) return;
        if (!use_euler_image) return;

        for (int lev = 0; lev <= finest_level; ++lev)
        {
            euler_mf[lev]->FillBoundary();
            amrex::Box domain = geom[lev].Domain();
            domain.convert(amrex::IntVect::TheNodeVector());

            for (amrex::MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Array4<MODEL> const& model = model_mf[lev]->array(mfi);
                amrex::Box bx = mfi.grownnodaltilebox() & domain;
                amrex::Array4<const Set::Scalar> const& euler = euler_mf[lev]->array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    if constexpr (MODEL::kinvar == Model::Solid::KinematicVariable::F)
                    {
                        Set::Scalar qx = euler(i, j, k, 0);
                        Set::Scalar qy = euler(i, j, k, 1);
                        Set::Scalar qz = euler(i, j, k, 2);
                        Set::Scalar sum_sq = qx*qx + qy*qy + qz*qz;
                        Set::Scalar qw = (sum_sq < 1.0) ? std::sqrt(1.0 - sum_sq) : 0.0;
                        Set::Scalar qnorm2 = qw*qw + sum_sq;
                        if (qnorm2 < 1e-12 || !std::isfinite(qnorm2))
                        { qw = 1.0; qx = 0.0; qy = 0.0; qz = 0.0; }
                        Eigen::Quaterniond q(qw, qx, qy, qz);
                        q.normalize();
                        model(i, j, k).q = q;
                    }
                });
            }
            Util::RealFillBoundary(*model_mf[lev], geom[lev]);
        }
    }

    void Advance(int lev, Set::Scalar time, Set::Scalar a_dt) override
    {
        BL_PROFILE("Integrator::LagrangianShockCP::Advance");
        if (m_type == Base::Mechanics<MODEL>::Type::Disable) return;

        const amrex::Real* DX = geom[lev].CellSize();
        Set::Scalar DX0 = DX[0];
        Set::Scalar DX1 = DX[1];
        amrex::Box cdomain = geom[lev].Domain();
        amrex::Box ndomain = amrex::convert(cdomain, amrex::IntVect::TheNodeVector());
        const amrex::Dim3 nlo = amrex::lbound(ndomain);
        const amrex::Dim3 nhi = amrex::ubound(ndomain);

        MieGruneisenEOS local_eos = eos;
        Set::Scalar local_rho0 = hydro_rho0;
        Set::Scalar local_CQ = av_CQ, local_CL = av_CL;
        Set::Scalar half_dt = 0.5 * a_dt;
        Set::Scalar t_full = time + a_dt;

        std::swap(*disp_mf[lev], *disp_old_mf[lev]);
        std::swap(*vel_mf[lev],  *vel_old_mf[lev]);
        vel_old_mf[lev]->FillBoundary(geom[lev].periodicity());
        disp_old_mf[lev]->FillBoundary(geom[lev].periodicity());

        ComputePressureAndViscosity2D(lev, DX0, DX1, local_eos, local_rho0,
                                      local_CQ, local_CL, disp_old_mf, vel_old_mf,
                                      pressure_mf, artvisc2d_mf, density_mf, soundspeed_mf);
        pressure_mf[lev]->FillBoundary(geom[lev].periodicity());
        artvisc2d_mf[lev]->FillBoundary(geom[lev].periodicity());
        soundspeed_mf[lev]->FillBoundary(geom[lev].periodicity());
        density_mf[lev]->FillBoundary(geom[lev].periodicity());

        CheckEOSReferenceState(lev);

        ComputeCellStressField(lev, DX0, DX1, nlo, nhi,
                               disp_old_mf, pressure_mf, artvisc2d_mf);
        cell_stress_mf[lev]->FillBoundary(geom[lev].periodicity());

        UpdateVelocityAndDisplacement2D(lev, DX0, DX1, half_dt, nlo, nhi, ndomain,
            local_rho0, vel_old_mf, disp_old_mf, vel_mf, disp_mf);
        vel_mf[lev]->FillBoundary(geom[lev].periodicity());
        disp_mf[lev]->FillBoundary(geom[lev].periodicity());

        ComputeHalfStepEnergy2D(lev, DX0, DX1, local_eos, local_rho0, local_CQ, local_CL);
        pressure_half_mf[lev]->FillBoundary(geom[lev].periodicity());
        artvisc2d_half_mf[lev]->FillBoundary(geom[lev].periodicity());
        soundspeed_half_mf[lev]->FillBoundary(geom[lev].periodicity());
        density_half_mf[lev]->FillBoundary(geom[lev].periodicity());

        ComputeCellStressField(lev, DX0, DX1, nlo, nhi,
                               disp_mf, pressure_half_mf, artvisc2d_half_mf);
        cell_stress_mf[lev]->FillBoundary(geom[lev].periodicity());

        UpdateVelocityAndDisplacement2D(lev, DX0, DX1, a_dt, nlo, nhi, ndomain,
            local_rho0, vel_old_mf, disp_old_mf, vel_mf, disp_mf);
        vel_mf[lev]->FillBoundary(geom[lev].periodicity());
        disp_mf[lev]->FillBoundary(geom[lev].periodicity());

        ApplyHourglassControl(lev, DX0, DX1, a_dt, nlo, nhi, ndomain, local_rho0);

        for (amrex::MFIter mfi(*density_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            amrex::Box cbx = mfi.tilebox();
            auto const& rho_arr = density_mf[lev]->array(mfi);
            auto const& u_arr = disp_mf[lev]->array(mfi);
            Set::Scalar lDX0 = DX0, lDX1 = DX1, lrho0 = local_rho0;
            amrex::ParallelFor(cbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Set::Scalar J = CellJ(u_arr, i, j, k, lDX0, lDX1);
                rho_arr(i,j,k) = lrho0 / J;
            });
        }

        for (amrex::MFIter mfi(*energy_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            amrex::Box cbx = mfi.tilebox();
            auto const& e_arr = energy_mf[lev]->array(mfi);
            auto const& ph_arr = pressure_half_mf[lev]->array(mfi);
            auto const& q2dh_arr = artvisc2d_half_mf[lev]->array(mfi);
            auto const& u_old_arr = disp_old_mf[lev]->array(mfi);
            auto const& u_new_arr = disp_mf[lev]->array(mfi);
            Set::Scalar lDX0 = DX0, lDX1 = DX1, lrho0 = local_rho0;
            amrex::ParallelFor(cbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Set::Scalar J_old = CellJ(u_old_arr, i, j, k, lDX0, lDX1);
                Set::Scalar J_new = CellJ(u_new_arr, i, j, k, lDX0, lDX1);
                Set::Scalar dJ = J_new - J_old;

                Set::Scalar dExx_old = 0.5*((u_old_arr(i+1,j,k)(0)-u_old_arr(i,j,k)(0))
                                           +(u_old_arr(i+1,j+1,k)(0)-u_old_arr(i,j+1,k)(0))) / lDX0;
                Set::Scalar dExx_new = 0.5*((u_new_arr(i+1,j,k)(0)-u_new_arr(i,j,k)(0))
                                           +(u_new_arr(i+1,j+1,k)(0)-u_new_arr(i,j+1,k)(0))) / lDX0;
                Set::Scalar dEyy_old = 0.5*((u_old_arr(i,j+1,k)(1)-u_old_arr(i,j,k)(1))
                                           +(u_old_arr(i+1,j+1,k)(1)-u_old_arr(i+1,j,k)(1))) / lDX1;
                Set::Scalar dEyy_new = 0.5*((u_new_arr(i,j+1,k)(1)-u_new_arr(i,j,k)(1))
                                           +(u_new_arr(i+1,j+1,k)(1)-u_new_arr(i+1,j,k)(1))) / lDX1;
                Set::Scalar d_exx = dExx_new - dExx_old;
                Set::Scalar d_eyy = dEyy_new - dEyy_old;

                Set::Scalar de = -ph_arr(i,j,k) * dJ / lrho0;
                de -= (q2dh_arr(i,j,k,0) * d_exx + q2dh_arr(i,j,k,1) * d_eyy) / lrho0;
                if (std::isfinite(de)) e_arr(i,j,k) += de;
            });
        }
        energy_mf[lev]->FillBoundary(geom[lev].periodicity());

        AdvanceCPModelAllNodes(lev, DX0, DX1, a_dt, t_full, nlo, nhi, ndomain,
                               disp_mf, disp_old_mf);
        Util::RealFillBoundary(*model_mf[lev], geom[lev]);

        UpdateTemperature(lev, DX0, DX1, a_dt, nlo, nhi, disp_old_mf, disp_mf);

        PopulateVisualizationFields(lev, DX, nlo, nhi, ndomain, disp_mf, pressure_mf);

        Set::Scalar dt_local = 1e100, max_p = -1e100, min_p = 1e100, fs_vel = 0.0;
        for (amrex::MFIter mfi(*vel_mf[lev], false); mfi.isValid(); ++mfi) {
            amrex::Box nbx = mfi.nodaltilebox() & ndomain;
            auto const& v_arr = vel_mf[lev]->array(mfi);
            if (nbx.contains(amrex::IntVect(AMREX_D_DECL(nhi.x, nlo.y, 0))))
                fs_vel = v_arr(nhi.x, nlo.y, 0)(0);
        }
        ComputePressureAndViscosity2D(lev, DX0, DX1, local_eos, local_rho0,
                                      local_CQ, local_CL, disp_mf, vel_mf,
                                      pressure_mf, artvisc2d_mf, density_mf, soundspeed_mf);
        Set::Scalar min_dx = std::min(DX0, DX1);
        for (amrex::MFIter mfi(*density_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            amrex::Box cbx = mfi.tilebox();
            auto const& p_arr = pressure_mf[lev]->array(mfi);
            auto const& c_arr = soundspeed_mf[lev]->array(mfi);
            amrex::Loop(cbx, [&](int i, int j, int k) {
                Set::Scalar c = c_arr(i,j,k);
                Set::Scalar dt_cfl = hydro_cfl * min_dx / (c + 1e-20);
                dt_local = std::min(dt_local, dt_cfl);
                if (p_arr(i,j,k) > max_p) max_p = p_arr(i,j,k);
                if (p_arr(i,j,k) < min_p) min_p = p_arr(i,j,k);
            });
        }
        diag_fs_vel = fs_vel; diag_max_pressure = max_p; diag_min_pressure = min_p;
        if (dt_local > 0.0 && dt_local < 1e90) timestep = dt_local;
        if (lev == 0 && amrex::ParallelDescriptor::IOProcessor() && t_full >= next_diag_time) {
            this->WriteDiagnostics(t_full, a_dt);
            next_diag_time = t_full + diag_dt;
        }
    }

protected:

    void CheckEOSReferenceState(int lev)
    {
#ifndef NDEBUG
        static bool checked = false;
        if (checked || lev != 0) return;
        checked = true;
        Set::Scalar e_min = energy_mf[lev]->min(0);
        Set::Scalar p_ref = eos.Pressure(hydro_rho0, e_min);

        if (std::fabs(p_ref) > 1.0e6) {
            Util::Warning(INFO,
                "LagrangianShockCP: EOS(rho0, e_min) = ", p_ref,
                " Pa is not zero.  Eq. (10) stress decomposition will "
                "inject a constant offset into the deviator on step 1.  "
                "Check energy_init in your input deck.");
        }
#else
        (void)lev;
#endif
    }

    AMREX_GPU_HOST_DEVICE static
    Set::Scalar CellJ(amrex::Array4<const Set::Vector> const& u,
                      int i, int j, int k,
                      Set::Scalar DX0, Set::Scalar DX1)
    {
        Set::Scalar dux_dx = 0.5*((u(i+1,j,k)(0)-u(i,j,k)(0))
                                 +(u(i+1,j+1,k)(0)-u(i,j+1,k)(0))) / DX0;
        Set::Scalar duy_dy = 0.5*((u(i,j+1,k)(1)-u(i,j,k)(1))
                                 +(u(i+1,j+1,k)(1)-u(i+1,j,k)(1))) / DX1;
        Set::Scalar dux_dy = 0.5*((u(i,j+1,k)(0)-u(i,j,k)(0))
                                 +(u(i+1,j+1,k)(0)-u(i+1,j,k)(0))) / DX1;
        Set::Scalar duy_dx = 0.5*((u(i+1,j,k)(1)-u(i,j,k)(1))
                                 +(u(i+1,j+1,k)(1)-u(i,j+1,k)(1))) / DX0;
        Set::Scalar J = (1.0+dux_dx)*(1.0+duy_dy) - dux_dy*duy_dx;

        if (!(J > 0.01)) J = 0.01;
        return J;
    }

    AMREX_GPU_HOST_DEVICE static
    void CellF2D(amrex::Array4<const Set::Vector> const& u,
                 int ci, int cj, int k,
                 Set::Scalar DX0, Set::Scalar DX1,
                 Set::Matrix& Fc)
    {
        Set::Scalar dux_dx = 0.5*((u(ci+1,cj,k)(0)-u(ci,cj,k)(0))
                                 +(u(ci+1,cj+1,k)(0)-u(ci,cj+1,k)(0))) / DX0;
        Set::Scalar duy_dy = 0.5*((u(ci,cj+1,k)(1)-u(ci,cj,k)(1))
                                 +(u(ci+1,cj+1,k)(1)-u(ci+1,cj,k)(1))) / DX1;
        Set::Scalar dux_dy = 0.5*((u(ci,cj+1,k)(0)-u(ci,cj,k)(0))
                                 +(u(ci+1,cj+1,k)(0)-u(ci+1,cj,k)(0))) / DX1;
        Set::Scalar duy_dx = 0.5*((u(ci+1,cj,k)(1)-u(ci,cj,k)(1))
                                 +(u(ci+1,cj+1,k)(1)-u(ci,cj+1,k)(1))) / DX0;
        Fc(0,0) = 1.0 + dux_dx;  Fc(0,1) = dux_dy;
        Fc(1,0) = duy_dx;        Fc(1,1) = 1.0 + duy_dy;
    }

    AMREX_GPU_HOST_DEVICE static
    Set::Scalar Compute3DCauchyTrace(const MODEL& mdl, Set::Matrix const& Fc)
    {
        if constexpr (MODEL::kinvar == Model::Solid::KinematicVariable::F) {
            Eigen::Matrix3d F3d = Eigen::Matrix3d::Identity();
            F3d(0,0) = Fc(0,0);  F3d(0,1) = Fc(0,1);
            F3d(1,0) = Fc(1,0);  F3d(1,1) = Fc(1,1);
            Eigen::Matrix3d Fpinv = mdl.Fp.inverse();
            Eigen::Matrix3d Fe3d = F3d * Fpinv;
            Eigen::Matrix3d Ce3d = Fe3d.transpose() * Fe3d;
            Eigen::Matrix3d Ee3d = 0.5 * (Ce3d - Eigen::Matrix3d::Identity());
            mdl.q.normalize();
            Eigen::Matrix3d R = mdl.q.toRotationMatrix();
            Set::Scalar c11 = mdl.C11, c12 = mdl.C12, c44 = mdl.C44;
            Eigen::Matrix3d Ee_cryst = R.transpose() * Ee3d * R;
            Eigen::Matrix3d Te_cryst = Eigen::Matrix3d::Zero();
            Te_cryst(0,0) = c11*Ee_cryst(0,0) + c12*Ee_cryst(1,1) + c12*Ee_cryst(2,2);
            Te_cryst(1,1) = c12*Ee_cryst(0,0) + c11*Ee_cryst(1,1) + c12*Ee_cryst(2,2);
            Te_cryst(2,2) = c12*Ee_cryst(0,0) + c12*Ee_cryst(1,1) + c11*Ee_cryst(2,2);
            Te_cryst(0,1) = 2.0*c44*Ee_cryst(0,1);  Te_cryst(1,0) = Te_cryst(0,1);
            Te_cryst(0,2) = 2.0*c44*Ee_cryst(0,2);  Te_cryst(2,0) = Te_cryst(0,2);
            Te_cryst(1,2) = 2.0*c44*Ee_cryst(1,2);  Te_cryst(2,1) = Te_cryst(1,2);
            Eigen::Matrix3d Te3d = R * Te_cryst * R.transpose();
            Set::Scalar J3d = F3d.determinant();
            if (std::abs(J3d) < 1e-10) J3d = (J3d >= 0.0) ? 1e-10 : -1e-10;
            Eigen::Matrix3d T3d = (1.0/J3d) * Fe3d * Te3d * Fe3d.transpose();
            return T3d(0,0) + T3d(1,1) + T3d(2,2);
        }
        else {
            Set::Matrix Pc = mdl.DW(Fc);
            Set::Scalar Jc = Fc(0,0)*Fc(1,1) - Fc(0,1)*Fc(1,0);
            if (std::abs(Jc) < 1e-10) Jc = (Jc >= 0.0) ? 1e-10 : -1e-10;
            Set::Matrix Tc = (1.0/Jc) * Pc * Fc.transpose();

            return Tc(0,0) + Tc(1,1);
        }
    }

    AMREX_GPU_HOST_DEVICE static
    Set::Matrix CellStressCorrected(
        const MODEL& mdl, Set::Matrix const& Fc,
        Set::Scalar P_eos, Set::Scalar qxx, Set::Scalar qyy,
        Set::Scalar stress_scale = 1.0)
    {
        if constexpr (MODEL::kinvar == Model::Solid::KinematicVariable::F) {
            if (mdl.Fp(0,0) == 0.0 || mdl.C11 == 0.0) return Set::Matrix::Zero();
        }
        Set::Matrix Pc = mdl.DW(Fc);
        Pc *= stress_scale;
        Set::Scalar trT_3d = Compute3DCauchyTrace(mdl, Fc) * stress_scale;
        Set::Scalar pres_cpl = -(1.0/3.0) * trT_3d;
        Set::Scalar dp = P_eos - pres_cpl;

        Set::Scalar Jc = Fc(0,0)*Fc(1,1) - Fc(0,1)*Fc(1,0);
        if (std::abs(Jc) < 1e-10) Jc = (Jc >= 0.0) ? 1e-10 : -1e-10;
        Set::Matrix FinvT;
        FinvT(0,0) =  Fc(1,1) / Jc; FinvT(0,1) = -Fc(1,0) / Jc;
        FinvT(1,0) = -Fc(0,1) / Jc; FinvT(1,1) =  Fc(0,0) / Jc;
        Set::Matrix Pcorr = Pc - dp * Jc * FinvT;
        Pcorr(0,0) -= qxx;
        Pcorr(1,1) -= qyy;
        return Pcorr;
    }

    void ComputePressureAndViscosity2D(
        int lev, Set::Scalar DX0, Set::Scalar DX1,
        MieGruneisenEOS local_eos, Set::Scalar local_rho0,
        Set::Scalar local_CQ, Set::Scalar local_CL,
        Set::Field<Set::Vector>& u_field, Set::Field<Set::Vector>& v_field,
        Set::Field<Set::Scalar>& p_out, Set::Field<Set::Scalar>& q2d_out,
        Set::Field<Set::Scalar>& rho_out, Set::Field<Set::Scalar>& cs_out)
    {
        for (amrex::MFIter mfi(*rho_out[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box cbx = mfi.tilebox();
            auto const& rho_arr = rho_out[lev]->array(mfi);
            auto const& p_arr   = p_out[lev]->array(mfi);
            auto const& e_arr   = energy_mf[lev]->array(mfi);
            auto const& q2d_arr = q2d_out[lev]->array(mfi);
            auto const& c_arr   = cs_out[lev]->array(mfi);
            auto const& q1d_arr = artvisc_mf[lev]->array(mfi);
            auto const& u_arr   = u_field[lev]->array(mfi);
            auto const& v_arr   = v_field[lev]->array(mfi);
            amrex::ParallelFor(cbx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Scalar J = CellJ(u_arr, i, j, k, DX0, DX1);
                Set::Scalar rho = local_rho0 / J;
                rho_arr(i,j,k) = rho;
                Set::Scalar p = local_eos.Pressure(rho, e_arr(i,j,k));
                if (!std::isfinite(p)) p = 0.0;
                p_arr(i,j,k) = p;
                c_arr(i,j,k) = local_eos.SoundSpeed(rho, p);
                Set::Scalar dvx_dx = 0.5*((v_arr(i+1,j,k)(0)-v_arr(i,j,k)(0))
                                         +(v_arr(i+1,j+1,k)(0)-v_arr(i,j+1,k)(0))) / DX0;
                Set::Scalar dvy_dy = 0.5*((v_arr(i,j+1,k)(1)-v_arr(i,j,k)(1))
                                         +(v_arr(i+1,j+1,k)(1)-v_arr(i+1,j,k)(1))) / DX1;
                Set::Scalar qxx = 0.0;
                if (dvx_dx < 0.0) {
                    Set::Scalar dv_abs = std::abs(dvx_dx)*DX0;
                    qxx = local_CQ*rho*dv_abs*dv_abs + local_CL*rho*c_arr(i,j,k)*dv_abs;
                }
                Set::Scalar qyy = 0.0;
                if (dvy_dy < 0.0) {
                    Set::Scalar dv_abs = std::abs(dvy_dy)*DX1;
                    qyy = local_CQ*rho*dv_abs*dv_abs + local_CL*rho*c_arr(i,j,k)*dv_abs;
                }
                q2d_arr(i,j,k,0) = qxx; q2d_arr(i,j,k,1) = qyy; q2d_arr(i,j,k,2) = 0.0;
                q1d_arr(i,j,k) = qxx;
            });
        }
    }

void ComputeCellStressField(
        int lev, Set::Scalar DX0, Set::Scalar DX1,
        amrex::Dim3 nlo, amrex::Dim3 nhi,
        Set::Field<Set::Vector>& u_field,
        Set::Field<Set::Scalar>& p_eos_field,
        Set::Field<Set::Scalar>& q2d_field)
    {
        amrex::Box cell_range(
            amrex::IntVect(AMREX_D_DECL(nlo.x, nlo.y, 0)),
            amrex::IntVect(AMREX_D_DECL(nhi.x-1, nhi.y-1, 0)));
        for (amrex::MFIter mfi(*cell_stress_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box cbx = mfi.tilebox() & cell_range;
            auto const& cs_arr = cell_stress_mf[lev]->array(mfi);
            auto const& mdl    = (*model_mf[lev]).array(mfi);
            auto const& u_arr  = u_field[lev]->array(mfi);
            auto const& p_arr  = p_eos_field[lev]->array(mfi);
            auto const& q2d    = q2d_field[lev]->array(mfi);
            Set::Scalar lDX0 = DX0, lDX1 = DX1, lscale = stress_scale;
            amrex::ParallelFor(cbx, [=] AMREX_GPU_DEVICE(int ci, int cj, int ck)
            {
                Set::Matrix Fc;
                CellF2D(u_arr, ci, cj, ck, lDX0, lDX1, Fc);

                Set::Scalar P_eos_c = p_arr(ci, cj, ck);
                Set::Scalar qxx_c   = q2d(ci, cj, ck, 0);
                Set::Scalar qyy_c   = q2d(ci, cj, ck, 1);

                Set::Matrix Pcorr_avg = Set::Matrix::Zero();
                int n_valid = 0;

                if constexpr (MODEL::kinvar == Model::Solid::KinematicVariable::F)
                {

                    Eigen::Matrix3d Fp_avg = Eigen::Matrix3d::Zero();
                    MODEL mdl_cell;
                    bool have_base = false;
                    for (int dj = 0; dj <= 1; ++dj)
                    {
                        for (int di = 0; di <= 1; ++di)
                        {
                            int ni = ci + di;
                            int nj = cj + dj;
                            if (mdl(ni,nj,ck).Fp(0,0) == 0.0 ||
                                mdl(ni,nj,ck).C11    == 0.0)
                                continue;
                            if (!have_base) { mdl_cell = mdl(ni,nj,ck); have_base = true; }
                            Fp_avg += mdl(ni,nj,ck).Fp;
                            ++n_valid;
                        }
                    }
                    if (n_valid > 0)
                    {
                        Fp_avg *= (1.0 / (Set::Scalar)n_valid);
                        mdl_cell.Fp = Fp_avg;
                        Pcorr_avg = CellStressCorrected(
                            mdl_cell, Fc, P_eos_c, qxx_c, qyy_c, lscale);
                    }
                }
                else
                {

                    for (int dj = 0; dj <= 1; ++dj)
                        for (int di = 0; di <= 1; ++di)
                        {
                            Pcorr_avg += CellStressCorrected(
                                mdl(ci+di, cj+dj, ck), Fc, P_eos_c, qxx_c, qyy_c, lscale);
                            ++n_valid;
                        }
                    if (n_valid > 0)
                        Pcorr_avg *= (1.0 / (Set::Scalar)n_valid);
                }

                cs_arr(ci,cj,ck,0) = Pcorr_avg(0,0);
                cs_arr(ci,cj,ck,1) = Pcorr_avg(1,1);
                cs_arr(ci,cj,ck,2) = Pcorr_avg(0,1);
                cs_arr(ci,cj,ck,3) = Pcorr_avg(1,0);
            });
        }
    }

    void UpdateVelocityAndDisplacement2D(int lev,
        Set::Scalar DX0, Set::Scalar DX1, Set::Scalar dt,
        amrex::Dim3 nlo, amrex::Dim3 nhi, amrex::Box ndomain, Set::Scalar local_rho0,
        Set::Field<Set::Vector>& v_in, Set::Field<Set::Vector>& u_in,
        Set::Field<Set::Vector>& v_out, Set::Field<Set::Vector>& u_out)
    {
        auto bc_ptr = bc;
        bool y_periodic = geom[lev].isPeriodic(1);
        for (amrex::MFIter mfi(*v_out[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            amrex::Box nbx = mfi.nodaltilebox() & ndomain;
            auto const& v_old = v_in[lev]->array(mfi);
            auto const& u_old = u_in[lev]->array(mfi);
            auto const& v_new = v_out[lev]->array(mfi);
            auto const& u_new = u_out[lev]->array(mfi);
            auto const& b_arr = rhs_mf[lev]->array(mfi);
            auto const& sig_vis = (*stress_mf[lev]).array(mfi);
            auto const& cs = cell_stress_mf[lev]->array(mfi);
            Set::Scalar lDX0 = DX0, lDX1 = DX1;
            amrex::ParallelFor(nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Set::Vector vnew_vec = v_old(i,j,k);
                Set::Vector unew_vec = u_old(i,j,k);
                bool AMREX_D_DECL(xmin=(i==nlo.x), ymin=(j==nlo.y), zmin=(k==nlo.z));
                bool AMREX_D_DECL(xmax=(i==nhi.x), ymax=(j==nhi.y), zmax=(k==nhi.z));
                if (y_periodic) { ymin = false; ymax = false; }
                bool on_boundary = AMREX_D_TERM(xmax||xmin, ||ymax||ymin, ||zmax||zmin);
                bool has_left  = (i > nlo.x);
                bool has_right = (i < nhi.x);
                bool has_bot   = (j > nlo.y) || y_periodic;
                bool has_top   = (j < nhi.y) || y_periodic;
                auto ReadCell = [&](int ci, int cj) -> Set::Matrix {
                    Set::Matrix S;
                    S(0,0) = cs(ci,cj,k,0); S(1,1) = cs(ci,cj,k,1);
                    S(0,1) = cs(ci,cj,k,2); S(1,0) = cs(ci,cj,k,3);
                    return S;
                };
                bool have_BL = has_left && has_bot;
                bool have_BR = has_right && has_bot;
                bool have_TL = has_left && has_top;
                bool have_TR = has_right && has_top;
                Set::Matrix S_BL=Set::Matrix::Zero(), S_BR=Set::Matrix::Zero();
                Set::Matrix S_TL=Set::Matrix::Zero(), S_TR=Set::Matrix::Zero();
                if (have_BL) S_BL = ReadCell(i-1,j-1);
                if (have_BR) S_BR = ReadCell(i,  j-1);
                if (have_TL) S_TL = ReadCell(i-1,j);
                if (have_TR) S_TR = ReadCell(i,  j);
                Set::Matrix S_xR=Set::Matrix::Zero(), S_xL=Set::Matrix::Zero();
                Set::Matrix S_yT=Set::Matrix::Zero(), S_yB=Set::Matrix::Zero();
                int n_xR=0, n_xL=0, n_yT=0, n_yB=0;
                if (have_BR){S_xR+=S_BR; n_xR++;} if (have_TR){S_xR+=S_TR; n_xR++;}
                if (have_BL){S_xL+=S_BL; n_xL++;} if (have_TL){S_xL+=S_TL; n_xL++;}
                if (have_TL){S_yT+=S_TL; n_yT++;} if (have_TR){S_yT+=S_TR; n_yT++;}
                if (have_BL){S_yB+=S_BL; n_yB++;} if (have_BR){S_yB+=S_BR; n_yB++;}
                if (n_xR>0) S_xR*=(1.0/n_xR);
                if (n_xL>0) S_xL*=(1.0/n_xL);
                if (n_yT>0) S_yT*=(1.0/n_yT);
                if (n_yB>0) S_yB*=(1.0/n_yB);

                Set::Scalar fy_len = (has_bot && has_top)    ? lDX1 : 0.5*lDX1;
                Set::Scalar fx_len = (has_left && has_right) ? lDX0 : 0.5*lDX0;
                Set::Vector force;
                for (int d = 0; d < AMREX_SPACEDIM; d++)
                    force(d) = (S_xR(d,0) - S_xL(d,0)) * fy_len
                             + (S_yT(d,1) - S_yB(d,1)) * fx_len;
                Set::Scalar mx = (has_left && has_right) ? lDX0 : 0.5*lDX0;
                Set::Scalar my = (has_bot  && has_top)   ? lDX1 : 0.5*lDX1;
                Set::Scalar node_mass = local_rho0 * mx * my;

                Set::Matrix Savg=Set::Matrix::Zero(); int na=0;
                if (have_BL){Savg+=S_BL;na++;} if (have_BR){Savg+=S_BR;na++;}
                if (have_TL){Savg+=S_TL;na++;} if (have_TR){Savg+=S_TR;na++;}
                if (na>0) Savg*=(1.0/na);
                sig_vis(i,j,k) = Savg;
                if (on_boundary) {
                    auto bctype = bc_ptr->getType(i,j,k,ndomain);
                    for (int d=0;d<AMREX_SPACEDIM;d++) {
                        if (bctype[d]==BC::Operator::Elastic::Elastic::Type::Displacement) {
                            unew_vec(d) = b_arr(i,j,k)(d);
                            vnew_vec(d) = (unew_vec(d)-u_old(i,j,k)(d))/dt;
                        } else {
                            Set::Scalar accel = force(d)/node_mass;
                            vnew_vec(d) = v_old(i,j,k)(d)+dt*accel;
                            unew_vec(d) = u_old(i,j,k)(d)+dt*vnew_vec(d);
                        }
                    }
                } else {
                    for (int d=0;d<AMREX_SPACEDIM;d++) {
                        Set::Scalar accel = force(d)/node_mass;
                        vnew_vec(d) = v_old(i,j,k)(d)+dt*accel;
                        unew_vec(d) = u_old(i,j,k)(d)+dt*vnew_vec(d);
                    }
                }
                v_new(i,j,k) = vnew_vec;
                u_new(i,j,k) = unew_vec;
            });
        }
    }

    void AdvanceCPModelAllNodes(int lev,
        Set::Scalar DX0, Set::Scalar DX1,
        Set::Scalar dt, Set::Scalar time,
        amrex::Dim3 nlo, amrex::Dim3 nhi, amrex::Box ndomain,
        Set::Field<Set::Vector>& u_new_field,
        Set::Field<Set::Vector>& u_old_field)
    {
        int nsub_max = cp_substeps_max;
        bool y_per = geom[lev].isPeriodic(1);

        for (amrex::MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi) {
            amrex::Box bx = mfi.nodaltilebox() & ndomain;
            auto const& model     = (*model_mf[lev]).array(mfi);
            auto const& u_new_arr = u_new_field[lev]->array(mfi);
            auto const& u_old_arr = u_old_field[lev]->array(mfi);
            Set::Scalar lDX0 = DX0, lDX1 = DX1;
            int lnsub_max = nsub_max;
            int lny = nhi.y - nlo.y;

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if constexpr (MODEL::kinvar == Model::Solid::KinematicVariable::F) {
                    if (model(i,j,k).Fp(0,0) == 0.0 || model(i,j,k).C11 == 0.0)
                        return;
                }

                bool has_left  = (i > nlo.x);
                bool has_right = (i < nhi.x);
                bool has_bot   = (j > nlo.y) || y_per;
                bool has_top   = (j < nhi.y) || y_per;

                bool have_BL = has_left  && has_bot;
                bool have_BR = has_right && has_bot;
                bool have_TL = has_left  && has_top;
                bool have_TR = has_right && has_top;

                int nvalid = (int)have_BL + (int)have_BR + (int)have_TL + (int)have_TR;
                if (nvalid == 0) return;

                auto wrap_j = [=] AMREX_GPU_HOST_DEVICE (int cj) -> int {
                    if (!y_per) return cj;
                    if (cj < nlo.y)        cj += lny;
                    else if (cj > nhi.y-1) cj -= lny;
                    return cj;
                };

                Set::Matrix Fc_new = Set::Matrix::Zero();
                Set::Matrix Fc_old = Set::Matrix::Zero();
                Set::Matrix Ftmp;

                auto accum = [&](int ci, int cj) {
                    cj = wrap_j(cj);
                    if (ci < nlo.x)   ci = nlo.x;
                    if (ci > nhi.x-1) ci = nhi.x - 1;
                    CellF2D(u_new_arr, ci, cj, k, lDX0, lDX1, Ftmp);
                    Fc_new += Ftmp;
                    CellF2D(u_old_arr, ci, cj, k, lDX0, lDX1, Ftmp);
                    Fc_old += Ftmp;
                };
                if (have_BL) accum(i-1, j-1);
                if (have_BR) accum(i,   j-1);
                if (have_TL) accum(i-1, j);
                if (have_TR) accum(i,   j);

                Set::Scalar inv_n = 1.0 / (Set::Scalar)nvalid;
                Fc_new *= inv_n;
                Fc_old *= inv_n;

                Set::Scalar dF_norm = 0.0;
                for (int ii = 0; ii < AMREX_SPACEDIM; ii++)
                    for (int jj = 0; jj < AMREX_SPACEDIM; jj++) {
                        Set::Scalar d = Fc_new(ii,jj) - Fc_old(ii,jj);
                        dF_norm += d * d;
                    }
                dF_norm = std::sqrt(dF_norm);

                int nsub = 1;
                if (dF_norm > 0.005)
                    nsub = std::min(lnsub_max, (int)(dF_norm / 0.005) + 1);
                if (nsub < 2 && dF_norm > 1e-8) nsub = 2;

                Set::Scalar dt_sub = dt / (Set::Scalar)nsub;

                for (int isub = 0; isub < nsub; isub++) {
                    Set::Scalar s = (Set::Scalar)(isub + 1) / (Set::Scalar)nsub;
                    Set::Matrix Fc_sub;
                    for (int ii = 0; ii < AMREX_SPACEDIM; ii++)
                        for (int jj = 0; jj < AMREX_SPACEDIM; jj++)
                            Fc_sub(ii,jj) = Fc_old(ii,jj)
                                          + s * (Fc_new(ii,jj) - Fc_old(ii,jj));
                    Set::Scalar t_sub = time - dt + s * dt;
                    Set::Matrix sig_dummy = Set::Matrix::Zero();
                    model(i,j,k).Advance(dt_sub, Fc_sub, sig_dummy, t_sub);
                }
            });
        }
    }

    void ApplyHourglassControl(int lev,
        Set::Scalar DX0, Set::Scalar DX1, Set::Scalar dt,
        amrex::Dim3 nlo, amrex::Dim3 nhi, amrex::Box ndomain,
        Set::Scalar local_rho0)
    {
        if (hg_coeff <= 0.0) return;

        amrex::Box cell_range(
            amrex::IntVect(AMREX_D_DECL(nlo.x, nlo.y, 0)),
            amrex::IntVect(AMREX_D_DECL(nhi.x-1, nhi.y-1, 0)));

        Set::Scalar lhg     = hg_coeff;
        Set::Scalar A_cell  = DX0 * DX1;
        bool y_periodic     = geom[lev].isPeriodic(1);

        amrex::BoxArray nodal_ba = amrex::convert(this->boxArray(lev),
                                                  amrex::IntVect::TheNodeVector());
        amrex::MultiFab hg_force(nodal_ba, this->DistributionMap(lev),
                                 AMREX_SPACEDIM, 1);
        hg_force.setVal(0.0);

        for (amrex::MFIter mfi(*cell_stress_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box cbx = mfi.tilebox() & cell_range;
            auto const& v_arr   = vel_mf[lev]->array(mfi);
            auto const& rho_arr = density_mf[lev]->array(mfi);
            auto const& cs_arr  = soundspeed_mf[lev]->array(mfi);
            auto const& f_arr   = hg_force.array(mfi);

            amrex::ParallelFor(cbx, [=] AMREX_GPU_DEVICE(int ci, int cj, int k)
            {
                Set::Scalar rho_c = rho_arr(ci,cj,k);
                Set::Scalar cs_c  = cs_arr(ci,cj,k);
                if (!std::isfinite(rho_c) || !std::isfinite(cs_c)) return;
                if (rho_c <= 0.0 || cs_c <= 0.0) return;

                int iBL = ci,   jBL = cj;
                int iBR = ci+1, jBR = cj;
                int iTR = ci+1, jTR = cj+1;
                int iTL = ci,   jTL = cj+1;

                for (int d = 0; d < AMREX_SPACEDIM; ++d)
                {
                    Set::Scalar h_d = ( v_arr(iBL,jBL,k)(d)
                                       -v_arr(iBR,jBR,k)(d)
                                       +v_arr(iTR,jTR,k)(d)
                                       -v_arr(iTL,jTL,k)(d));
                    Set::Scalar coef = -0.25 * lhg * rho_c * cs_c * A_cell * h_d;

                    amrex::Gpu::Atomic::AddNoRet(&f_arr(iBL,jBL,k,d),  coef * (+1.0));
                    amrex::Gpu::Atomic::AddNoRet(&f_arr(iBR,jBR,k,d),  coef * (-1.0));
                    amrex::Gpu::Atomic::AddNoRet(&f_arr(iTR,jTR,k,d),  coef * (+1.0));
                    amrex::Gpu::Atomic::AddNoRet(&f_arr(iTL,jTL,k,d),  coef * (-1.0));
                }
            });
        }

        hg_force.SumBoundary(geom[lev].periodicity());

        auto bc_ptr = bc;
        for (amrex::MFIter mfi(*vel_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box nbx = mfi.nodaltilebox() & ndomain;
            auto const& v_arr = vel_mf[lev]->array(mfi);
            auto const& f_arr = hg_force.const_array(mfi);
            Set::Scalar lDX0 = DX0, lDX1 = DX1;
            bool y_per = y_periodic;

            amrex::ParallelFor(nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                bool xmin = (i == nlo.x);
                bool xmax = (i == nhi.x);
                bool ymin = (j == nlo.y) && !y_per;
                bool ymax = (j == nhi.y) && !y_per;
                bool on_boundary = xmin || xmax || ymin || ymax;

                Set::Scalar mx = (!xmin && !xmax) ? lDX0 : 0.5*lDX0;
                Set::Scalar my = (!ymin && !ymax) ? lDX1 : 0.5*lDX1;
                Set::Scalar node_mass = local_rho0 * mx * my;
                if (node_mass <= 0.0) return;

                Set::Vector vnew = v_arr(i,j,k);
                if (on_boundary)
                {

                    auto bctype = bc_ptr->getType(i,j,k,ndomain);
                    for (int d = 0; d < AMREX_SPACEDIM; ++d)
                    {
                        if (bctype[d] == BC::Operator::Elastic::Elastic::Type::Displacement)
                            continue;
                        Set::Scalar dv = (f_arr(i,j,k,d) / node_mass) * dt;
                        if (std::isfinite(dv)) vnew(d) += dv;
                    }
                }
                else
                {
                    for (int d = 0; d < AMREX_SPACEDIM; ++d)
                    {
                        Set::Scalar dv = (f_arr(i,j,k,d) / node_mass) * dt;
                        if (std::isfinite(dv)) vnew(d) += dv;
                    }
                }
                v_arr(i,j,k) = vnew;
            });
        }
        vel_mf[lev]->FillBoundary(geom[lev].periodicity());
    }

    void UpdateTemperature(int lev,
        Set::Scalar , Set::Scalar , Set::Scalar ,
        amrex::Dim3 , amrex::Dim3 ,
        Set::Field<Set::Vector>& ,
        Set::Field<Set::Vector>& )
    {
        if (taylor_quinney <= 0.0) return;

        if constexpr (!has_wp_inc<MODEL>::value) {
            return;
        } else {
            Set::Scalar beta     = taylor_quinney;
            Set::Scalar Cv       = eos.Cv;
            Set::Scalar rho0_loc = hydro_rho0;
            Set::Scalar scale    = stress_scale;

            for (amrex::MFIter mfi(*temperature_mf[lev], amrex::TilingIfNotGPU());
                 mfi.isValid(); ++mfi)
            {
                amrex::Box cbx = mfi.tilebox();
                auto const& theta_arr = temperature_mf[lev]->array(mfi);
                auto const& mdl       = (*model_mf[lev]).array(mfi);

                amrex::ParallelFor(cbx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {

                    Set::Scalar wp = 0.0;
                    int cnt = 0;
                    auto accum = [&](int ni, int nj) {
                        if constexpr (MODEL::kinvar == Model::Solid::KinematicVariable::F) {
                            if (mdl(ni,nj,k).Fp(0,0) == 0.0) return;
                        }
                        wp += mdl(ni,nj,k).wp_inc;
                        cnt++;
                    };
                    accum(i,   j);
                    accum(i+1, j);
                    accum(i,   j+1);
                    accum(i+1, j+1);
                    if (cnt == 0) return;
                    wp /= (Set::Scalar)cnt;

                    Set::Scalar W_p_Pa = wp * scale;
                    Set::Scalar dtheta = beta * W_p_Pa / (rho0_loc * Cv);
                    if (std::isfinite(dtheta) && dtheta > 0.0)
                        theta_arr(i,j,k) += dtheta;
                });
            }
            temperature_mf[lev]->FillBoundary(geom[lev].periodicity());
        }
    }

    void PopulateVisualizationFields(int lev, const amrex::Real* DX,
        amrex::Dim3 nlo, amrex::Dim3 nhi, amrex::Box ndomain,
        Set::Field<Set::Vector>& u_field, Set::Field<Set::Scalar>& )
    {

        bool yper_vis = geom[lev].isPeriodic(1);
        for (amrex::MFIter mfi(*stress_mf[lev], false); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.nodaltilebox() & ndomain;
            auto const& sig    = (*stress_mf[lev]).array(mfi);
            auto const& eps    = (*strain_mf[lev]).array(mfi);
            auto const& u_arr  = u_field[lev]->array(mfi);
            auto const& cs_arr = cell_stress_mf[lev]->const_array(mfi);
            bool yper = yper_vis;

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {

                bool has_left  = (i > nlo.x);
                bool has_right = (i < nhi.x);
                bool has_bot   = (j > nlo.y) || yper;
                bool has_top   = (j < nhi.y) || yper;

                Set::Matrix S = Set::Matrix::Zero();
                int n_valid = 0;
                auto add_cell = [&](int ci, int cj) {
                    Set::Matrix Sc;
                    Sc(0,0) = cs_arr(ci,cj,k,0);
                    Sc(1,1) = cs_arr(ci,cj,k,1);
                    Sc(0,1) = cs_arr(ci,cj,k,2);
                    Sc(1,0) = cs_arr(ci,cj,k,3);
                    S += Sc;
                    ++n_valid;
                };
                if (has_left  && has_bot) add_cell(i-1, j-1);
                if (has_right && has_bot) add_cell(i,   j-1);
                if (has_left  && has_top) add_cell(i-1, j);
                if (has_right && has_top) add_cell(i,   j);
                if (n_valid > 0) S *= (1.0 / (Set::Scalar)n_valid);
                sig(i,j,k) = S;

                auto sten = Numeric::GetStencil(i, j, k, bx);
                Set::Matrix gradu = Numeric::Gradient(u_arr, i, j, k, DX, sten);
                eps(i,j,k) = Set::Matrix::Identity() + gradu;
            });
        }
        stress_mf[lev]->setMultiGhost(true);
        stress_mf[lev]->FillBoundaryAndSync(geom[lev].periodicity());

        for (amrex::MFIter mfi(*sigma_cp_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box cbx = mfi.tilebox();
            auto const& sc     = sigma_cp_mf[lev]->array(mfi);
            auto const& cs_arr = cell_stress_mf[lev]->const_array(mfi);
            amrex::ParallelFor(cbx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                sc(i,j,k) = -0.5 * (cs_arr(i,j,k,0) + cs_arr(i,j,k,1));
            });
        }
    }

    void ComputeHalfStepEnergy2D(
        int lev, Set::Scalar DX0, Set::Scalar DX1,
        MieGruneisenEOS local_eos, Set::Scalar local_rho0,
        Set::Scalar local_CQ, Set::Scalar local_CL)
    {
        for (amrex::MFIter mfi(*density_half_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box cbx = mfi.tilebox();
            auto const& rho_h  = density_half_mf[lev]->array(mfi);
            auto const& p_h    = pressure_half_mf[lev]->array(mfi);
            auto const& e_h    = energy_half_mf[lev]->array(mfi);
            auto const& q2d_h  = artvisc2d_half_mf[lev]->array(mfi);
            auto const& c_h    = soundspeed_half_mf[lev]->array(mfi);
            auto const& u_old  = disp_old_mf[lev]->array(mfi);
            auto const& u_half = disp_mf[lev]->array(mfi);
            auto const& v_half = vel_mf[lev]->array(mfi);
            auto const& p_old  = pressure_mf[lev]->array(mfi);
            auto const& q2d_old = artvisc2d_mf[lev]->array(mfi);
            auto const& e_old  = energy_mf[lev]->array(mfi);
            auto const& qh_1d  = artvisc_half_mf[lev]->array(mfi);
            amrex::ParallelFor(cbx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Scalar J_half = CellJ(u_half, i, j, k, DX0, DX1);
                Set::Scalar J_old  = CellJ(u_old,  i, j, k, DX0, DX1);
                Set::Scalar rho = local_rho0 / J_half;
                rho_h(i,j,k) = rho;
                Set::Scalar dJ = J_half - J_old;
                Set::Scalar dExx_old = 0.5*((u_old(i+1,j,k)(0)-u_old(i,j,k)(0))
                                           +(u_old(i+1,j+1,k)(0)-u_old(i,j+1,k)(0))) / DX0;
                Set::Scalar dExx_half = 0.5*((u_half(i+1,j,k)(0)-u_half(i,j,k)(0))
                                            +(u_half(i+1,j+1,k)(0)-u_half(i,j+1,k)(0))) / DX0;
                Set::Scalar dEyy_old = 0.5*((u_old(i,j+1,k)(1)-u_old(i,j,k)(1))
                                           +(u_old(i+1,j+1,k)(1)-u_old(i+1,j,k)(1))) / DX1;
                Set::Scalar dEyy_half = 0.5*((u_half(i,j+1,k)(1)-u_half(i,j,k)(1))
                                            +(u_half(i+1,j+1,k)(1)-u_half(i+1,j,k)(1))) / DX1;
                Set::Scalar d_exx = dExx_half - dExx_old;
                Set::Scalar d_eyy = dEyy_half - dEyy_old;
                Set::Scalar e_new = e_old(i,j,k) - p_old(i,j,k) * dJ / local_rho0
                                  - (q2d_old(i,j,k,0)*d_exx + q2d_old(i,j,k,1)*d_eyy) / local_rho0;
                if (!std::isfinite(e_new)) e_new = e_old(i,j,k);
                e_h(i,j,k) = e_new;
                Set::Scalar p = local_eos.Pressure(rho, e_new);
                if (!std::isfinite(p)) p = 0.0;
                p_h(i,j,k) = p;
                c_h(i,j,k) = local_eos.SoundSpeed(rho, p);
                Set::Scalar dvx_dx = 0.5*((v_half(i+1,j,k)(0)-v_half(i,j,k)(0))
                                         +(v_half(i+1,j+1,k)(0)-v_half(i,j+1,k)(0))) / DX0;
                Set::Scalar dvy_dy = 0.5*((v_half(i,j+1,k)(1)-v_half(i,j,k)(1))
                                         +(v_half(i+1,j+1,k)(1)-v_half(i+1,j,k)(1))) / DX1;
                Set::Scalar qxx = 0.0;
                if (dvx_dx < 0.0) {
                    Set::Scalar dv_abs = std::abs(dvx_dx)*DX0;
                    qxx = local_CQ*rho*dv_abs*dv_abs + local_CL*rho*c_h(i,j,k)*dv_abs;
                }
                Set::Scalar qyy = 0.0;
                if (dvy_dy < 0.0) {
                    Set::Scalar dv_abs = std::abs(dvy_dy)*DX1;
                    qyy = local_CQ*rho*dv_abs*dv_abs + local_CL*rho*c_h(i,j,k)*dv_abs;
                }
                q2d_h(i,j,k,0) = qxx; q2d_h(i,j,k,1) = qyy; q2d_h(i,j,k,2) = 0.0;
                qh_1d(i,j,k) = qxx;
            });
        }
    }

    Set::Field<Set::Scalar> temperature_mf;
    Set::Field<Set::Scalar> sigma_cp_mf;
    Set::Field<Set::Scalar> artvisc2d_mf;
    Set::Field<Set::Scalar> artvisc2d_half_mf;
    Set::Field<Set::Scalar> cell_stress_mf;
    Set::Field<Set::Scalar> cell_stress_raw_mf;
    Set::Scalar taylor_quinney = 0.9;
    Set::Scalar theta_init = 294.0;
    Set::Scalar stress_scale = 1e9;
    int cp_substeps_max = 8;
    Set::Scalar hg_coeff = 0.1;

    Set::Field<Set::Scalar> euler_mf;
    IC::IC<Set::Scalar>* ic_euler = nullptr;
    bool use_euler_image = false;

    void Regrid(int lev, Set::Scalar time) override
    {
        if (use_euler_image && ic_euler)
            ic_euler->Initialize(lev, euler_mf, time);
    }

    using LagrangianShock<MODEL>::eos;
    using LagrangianShock<MODEL>::hydro_rho0;
    using LagrangianShock<MODEL>::av_CQ;
    using LagrangianShock<MODEL>::av_CL;
    using LagrangianShock<MODEL>::hydro_cfl;
    using LagrangianShock<MODEL>::density_mf;
    using LagrangianShock<MODEL>::pressure_mf;
    using LagrangianShock<MODEL>::energy_mf;
    using LagrangianShock<MODEL>::artvisc_mf;
    using LagrangianShock<MODEL>::soundspeed_mf;
    using LagrangianShock<MODEL>::density_half_mf;
    using LagrangianShock<MODEL>::pressure_half_mf;
    using LagrangianShock<MODEL>::energy_half_mf;
    using LagrangianShock<MODEL>::artvisc_half_mf;
    using LagrangianShock<MODEL>::soundspeed_half_mf;
    using LagrangianShock<MODEL>::bc_nothing;
    using LagrangianShock<MODEL>::diag_fs_vel;
    using LagrangianShock<MODEL>::diag_max_pressure;
    using LagrangianShock<MODEL>::diag_min_pressure;
    using LagrangianShock<MODEL>::diag_filename;
    using LagrangianShock<MODEL>::diag_dt;
    using LagrangianShock<MODEL>::next_diag_time;
    using LagrangianShock<MODEL>::diag_header_written;
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