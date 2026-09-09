#ifndef INTEGRATOR_LAGRANGIANSHOCKCPDAMAGE_H
#define INTEGRATOR_LAGRANGIANSHOCKCPDAMAGE_H

#include "Integrator/LagrangianShockCP.H"

#include "IC/IC.H"
#include "IC/PNG3.H"
#include "BC/BC.H"
#include "BC/Constant.H"
#include "BC/Nothing.H"
#include "Numeric/Stencil.H"
#include "Numeric/Function.H"
#include "Util/Util.H"

#include <random>
#include <cstdint>   
#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <string>
#include <limits>
#include <cctype>

namespace Integrator
{

enum class DamageModel { DoubleWell = 0, AT1 = 1 };

enum class XTReduce { Row = 0, Mean = 1, Max = 2 };

template <class MODEL>
class LagrangianShockCPDamage : public LagrangianShockCP<MODEL>
{
public:
    using Base = LagrangianShockCP<MODEL>;

    static constexpr const char* name = "lagrangian_shock_cp_damage";

    LagrangianShockCPDamage() : Base() {}
    LagrangianShockCPDamage(IO::ParmParse& pp) : Base()
    { Parse(*this, pp); }

    ~LagrangianShockCPDamage() override
    {
        delete ic_heatmap;
        delete bc_d;
        delete bc_nothing_scalar;
    }

    static void Parse(LagrangianShockCPDamage& v, IO::ParmParse& pp)
    {
        BL_PROFILE("Integrator::LagrangianShockCPDamage::Parse");

        Base::Parse(v, pp);

        pp_query_default("xt.file", v.xt_filename, std::string("lagshock_xt.dat"));
        pp_query_default("xt.dt",   v.xt_dt,       1.0e-8);
        {
            int xt_on = 0;
            pp_query_default("xt.enable", xt_on, 0);
            v.xt_enable = xt_on;
        }
        {
            int jprobe = -1;
            pp_query_default("xt.jprobe", jprobe, -1);   // -1 => centerline
            v.xt_jprobe = jprobe;
        }
        {
            std::string red = "row";
            pp_query_default("xt.reduce", red, std::string("row"));

            std::vector<std::string> toks;
            {
                std::string cur;
                for (char c : red)
                {
                    if (c == ',') { if (!cur.empty()) toks.push_back(cur); cur.clear(); }
                    else if (!std::isspace((unsigned char)c)) cur.push_back(c);
                }
                if (!cur.empty()) toks.push_back(cur);
            }
            if (toks.empty()) toks.push_back("row");

            for (const std::string& t : toks)
            {
                XTReduce m;
                if      (t == "row")                                  m = XTReduce::Row;
                else if (t == "mean" || t == "mean_y" || t == "avg")  m = XTReduce::Mean;
                else if (t == "max"  || t == "max_y")                 m = XTReduce::Max;
                else
                {
                    Util::Abort(INFO,
                        "xt.reduce token '", t, "' not recognized; use a "
                        "comma-separated subset of 'row', 'mean', 'max' "
                        "(e.g. xt.reduce = row,mean).");
                    continue;
                }
                bool dup = false;
                for (XTReduce e : v.xt_modes) if (e == m) dup = true;
                if (dup)
                {
                    Util::Warning(INFO,
                        "xt.reduce lists '", t, "' more than once; ignoring "
                        "the repeat.");
                    continue;
                }
                v.xt_modes.push_back(m);
            }

            const bool multi = (v.xt_modes.size() > 1);
            for (XTReduce m : v.xt_modes)
            {
                if (!multi) { v.xt_files.push_back(v.xt_filename); continue; }

                const char* sfx = (m == XTReduce::Mean) ? "_mean"
                                : (m == XTReduce::Max)  ? "_max" : "_row";
                const std::size_t dot = v.xt_filename.find_last_of('.');
                if (dot == std::string::npos)
                    v.xt_files.push_back(v.xt_filename + sfx);
                else
                    v.xt_files.push_back(v.xt_filename.substr(0, dot) + sfx
                                       + v.xt_filename.substr(dot));
            }
            v.xt_header_written.assign(v.xt_modes.size(), 0);

            bool has_nonrow = false;
            for (XTReduce m : v.xt_modes) if (m != XTReduce::Row) has_nonrow = true;
            if (has_nonrow && v.xt_jprobe >= 0)
                Util::Warning(INFO,
                    "xt.jprobe = ", v.xt_jprobe, " applies only to the 'row' "
                    "mode; the mean/max reductions span all j.");

            if (amrex::ParallelDescriptor::IOProcessor() && v.xt_enable > 0)
                for (std::size_t q = 0; q < v.xt_modes.size(); ++q)
                    Util::Message(INFO, "[xt] mode ",
                        (v.xt_modes[q] == XTReduce::Mean ? "mean"
                       : v.xt_modes[q] == XTReduce::Max  ? "max" : "row"),
                        " -> ", v.xt_files[q]);
        }

        pp_query_required("damage.M",            v.dmg.M);


        {
            std::string model_str = "double_well";
            pp_query("damage.model", model_str);
            if (model_str == "double_well" || model_str == "doublewell"
                || model_str == "dw")
                v.dmg.model = DamageModel::DoubleWell;
            else if (model_str == "at1" || model_str == "AT1")
                v.dmg.model = DamageModel::AT1;
            else
                Util::Abort(INFO,
                    "damage.model = '", model_str,
                    "' not recognized; use 'double_well' (default) or 'at1'.");
        }

        if (v.dmg.model == DamageModel::DoubleWell)
        {

            bool has_bh = pp.contains("damage.barrier_height");
            bool has_cp = pp.contains("damage.cp");
            if (has_bh)
            {
                pp_query_required("damage.barrier_height", v.dmg.barrier_height);
                if (v.dmg.barrier_height <= 0.0)
                    Util::Abort(INFO, "damage.barrier_height must be > 0");
                v.dmg.cp = 16.0 * v.dmg.barrier_height;
                if (has_cp)
                {
                    Set::Scalar cp_override = 0.0;
                    pp_query_default("damage.cp", cp_override, v.dmg.cp);
                    Util::Warning(INFO,
                        "damage.barrier_height and damage.cp both set; "
                        "barrier_height wins (c_p = 16*barrier_height = ",
                        v.dmg.cp, "), ignoring damage.cp = ", cp_override);
                }
            }
            else if (has_cp)
            {
                pp_query_required("damage.cp", v.dmg.cp);
                if (v.dmg.cp <= 0.0)
                    Util::Abort(INFO, "damage.cp must be > 0");
                v.dmg.barrier_height = v.dmg.cp / 16.0;
                Util::Warning(INFO,
                    "damage.cp set directly without damage.barrier_height; "
                    "double-well barrier has no pinned physical identity. "
                    "Prefer setting damage.barrier_height (the void-nucleation "
                    "activation energy density); implied barrier_height = ",
                    v.dmg.barrier_height);
            }
            else
            {
                Util::Abort(INFO,
                    "damage.model=double_well requires a double-well energy "
                    "scale: set damage.barrier_height (preferred) or "
                    "damage.cp.");
            }
        }
        else // DamageModel::AT1
        {

            v.dmg.cw = 8.0 / 3.0;

            const bool has_Gc  = pp.contains("damage.Gc");
            const bool has_ell = pp.contains("damage.ell");
            const bool has_bh  = pp.contains("damage.barrier_height");
            const bool has_cp  = pp.contains("damage.cp");
            const bool has_kap = pp.contains("damage.kappa");

            if (has_Gc && has_ell)
            {
                pp_query_required("damage.Gc",  v.dmg.Gc);
                pp_query_required("damage.ell", v.dmg.ell);
                if (v.dmg.Gc <= 0.0)
                    Util::Abort(INFO, "damage.Gc must be > 0");
                if (v.dmg.ell <= 0.0)
                    Util::Abort(INFO, "damage.ell must be > 0");
                v.dmg.barrier_height = v.dmg.Gc / (v.dmg.cw * v.dmg.ell);
                v.dmg.cp    = 16.0 * v.dmg.barrier_height;
                v.dmg.kappa = v.dmg.cp * v.dmg.ell * v.dmg.ell;
            }
            else if ((has_bh || has_cp) && has_kap)
            {

                pp_query_required("damage.kappa", v.dmg.kappa);
                if (has_bh)
                {
                    pp_query_required("damage.barrier_height",
                                      v.dmg.barrier_height);
                    if (v.dmg.barrier_height <= 0.0)
                        Util::Abort(INFO,
                            "damage.barrier_height must be > 0");
                    v.dmg.cp = 16.0 * v.dmg.barrier_height;
                }
                else
                {
                    pp_query_required("damage.cp", v.dmg.cp);
                    if (v.dmg.cp <= 0.0)
                        Util::Abort(INFO, "damage.cp must be > 0");
                    v.dmg.barrier_height = v.dmg.cp / 16.0;
                }
                if (v.dmg.kappa <= 0.0)
                    Util::Abort(INFO, "damage.kappa must be > 0");
                v.dmg.ell = std::sqrt(v.dmg.kappa / v.dmg.cp);
                v.dmg.Gc  = v.dmg.barrier_height * v.dmg.cw * v.dmg.ell;
                Util::Warning(INFO,
                    "damage.model=at1 with no damage.Gc/damage.ell: "
                    "back-deriving from the legacy double-well deck. "
                    "ell = sqrt(kappa/cp) = ", v.dmg.ell,
                    ", G_c = barrier_height*c_w*ell = ", v.dmg.Gc,
                    " (AT1 threshold G_c/(c_w*ell) = ",
                    v.dmg.Gc / (v.dmg.cw * v.dmg.ell),
                    " Pa, matched to the old barrier-height energy "
                    "density). Set damage.Gc and damage.ell explicitly "
                    "to override.");
            }
            else
            {
                Util::Abort(INFO,
                    "damage.model=at1 requires either (damage.Gc and "
                    "damage.ell), or a legacy (damage.barrier_height|cp "
                    "and damage.kappa) pair to back-derive them.");
            }
        }

        pp_query_default ("damage.small",        v.dmg.small, 1.0e-3);

        pp_query_default ("damage.d_void",       v.dmg.d_void, 0.95);

        pp_query_default ("damage.void_res",     v.dmg.void_res, 0.02);
        if (v.dmg.void_res <= 0.0 || v.dmg.void_res > 1.0)
            Util::Abort(INFO,
                "damage.void_res must lie in (0, 1]: residual stiffness "
                "fraction of a locked spall void (sound speed, tensile "
                "pressure, cell stress). Use a small positive value (e.g. "
                "0.02); void_res -> 0 recovers the hard-vacuum stamp and its "
                "mesh-refinement pressure runaway.");

        
        pp_query_default ("damage.d_seed_floor", v.dmg.d_seed_floor, 0.55);
        if (v.dmg.model == DamageModel::DoubleWell)
        {

            if (v.dmg.d_seed_floor <= 0.5 || v.dmg.d_seed_floor > 1.0)
                Util::Abort(INFO,
                    "damage.d_seed_floor must lie in (0.5, 1.0]; got ",
                    v.dmg.d_seed_floor, ". Values at or below the double-well "
                    "barrier d = 1/2 leave stamped cells to be driven back to "
                    "zero by the barrier.");
        }
        else // AT1
        {

            if (v.dmg.d_seed_floor < 0.0 || v.dmg.d_seed_floor > 1.0)
                Util::Abort(INFO,
                    "damage.d_seed_floor must lie in [0, 1]; got ",
                    v.dmg.d_seed_floor);
            if (v.dmg.d_seed_floor > 0.0 && v.dmg.d_seed_floor <= 0.5)
                Util::Message(INFO,
                    "damage.model=at1: damage.d_seed_floor = ",
                    v.dmg.d_seed_floor, " (<= 0.5 is allowed under AT1; "
                    "there is no double-well barrier to cross -- it is "
                    "just the stamped void-core initialization amplitude).");
        }

        pp_query_default ("damage.hist_seed_factor", v.dmg.hist_seed_factor, 0.0);

        pp_query_default ("damage.p_tens_floor", v.dmg.p_tens_floor, 0.0);
        if (v.dmg.p_tens_floor < 0.0)
            Util::Abort(INFO,
                "damage.p_tens_floor must be >= 0 (it is a tensile-pressure "
                "floor, in the same units as the EOS pressure).");


        pp_query_default ("growth.M",          v.grow.M, 0.0);
        if (pp.contains("growth.M_plastic"))
        {
            Set::Scalar mp_legacy = 0.0;
            pp_query_default("growth.M_plastic", mp_legacy, 0.0);
            if (!pp.contains("growth.M") && mp_legacy > 0.0)
            {
                v.grow.M = mp_legacy;
                Util::Warning(INFO,
                    "growth.M_plastic is deprecated (it scaled the wrong "
                    "driver -- plastic work rather than plastic strain "
                    "rate). Using its value ", mp_legacy, " for growth.M; "
                    "rename to growth.M in the input deck.");
            }
            else if (mp_legacy > 0.0)
            {
                Util::Warning(INFO,
                    "growth.M_plastic is deprecated and ignored because "
                    "growth.M is also set; remove growth.M_plastic.");
            }
        }

        {
            int prop_int = 0;
            pp_query_default("growth.propagate", prop_int, 0);
            v.grow.propagate = (prop_int != 0);
        }

        pp_query_default ("growth.T_cap",      v.grow.T_cap, 3.0);

        pp_query_default("damage.bypass", v.dbg_bypass, 0);   // ### BISECT SWITCH
        if (v.dbg_bypass)
            Util::Warning(INFO,
                "LagrangianShockCPDamage: damage.bypass = 1 -- ALL damage-specific "
                "work inside Advance is disabled. The time-stepping path is now "
                "identical to LagrangianShockCP. Registrations/Initialize are "
                "unchanged. DIAGNOSTIC ONLY.");

        pp_query_default("eta_ref_threshold", v.m_eta_ref_threshold, 1e-2);

        pp_query_default("nucleation.enable",            v.nuc.enable, 0);
        pp_query_default("nucleation.use_stress",        v.nuc.use_stress, 1);
        pp_query_default("nucleation.base_rate",         v.nuc.base_rate, 0.0);
        pp_query_default("nucleation.p_thresh",          v.nuc.p_thresh, 0.0);
        pp_query_default("nucleation.p_ref",             v.nuc.p_ref, 1.0e9);

        pp_query_default("nucleation.hm_alpha",          v.nuc.hm_alpha, 0.0);
        pp_query_default("nucleation.hm_beta",           v.nuc.hm_beta, 0.0);


        pp_query_default("nucleation.sigma_spall",       v.nuc.sigma_spall, 0.0);
        pp_query_default("nucleation.g1_ref",            v.nuc.g1_ref, 0.0);
        pp_query_default("nucleation.m_exp",             v.nuc.m_exp, 1.0);

        pp_query_default("nucleation.heatmap_threshold", v.nuc.hm_thresh, 0.0);
        pp_query_default("nucleation.void_radius",       v.nuc.void_radius, 2.5e-6);
        pp_query_default("nucleation.eps",               v.nuc.eps, 5.0e-7);
        pp_query_default("nucleation.exclusion_radius",  v.nuc.exclusion_radius, 1.0e-5);
        pp_query_default("nucleation.check_frequency",   v.nuc.check_freq, 10);
        pp_query_default("nucleation.boundary_exclude_cells", v.nuc.bnd_excl_cells, 4);
        pp_query_default("nucleation.max_voids",         v.nuc.max_voids, 10);
        {
            int seed_val = 12345;
            pp_query_default("nucleation.seed", seed_val, 12345);
            v.nuc.seed = (unsigned)seed_val;
        }

        if (v.nuc.use_stress == 3)
        {
            if (v.nuc.p_ref <= 0.0)
                Util::Abort(INFO, "Boltzmann nucleation requires nucleation.p_ref > 0");
            if (v.nuc.hm_alpha < 0.0 || v.nuc.hm_alpha > 1.0)
                Util::Warning(INFO,
                    "nucleation.hm_alpha = ", v.nuc.hm_alpha,
                    " is outside [0,1]; barrier may become negative or grow");
            if (v.nuc.hm_beta < 0.0)
                Util::Warning(INFO,
                    "nucleation.hm_beta = ", v.nuc.hm_beta,
                    " is negative; this inverts heatmap weighting");
        }

        if (v.nuc.use_stress == 4)
        {
            if (v.nuc.sigma_spall <= 0.0)
                Util::Abort(INFO,
                    "use_stress=4 (dual-barrier Boltzmann) requires "
                    "nucleation.sigma_spall > 0 -- it is the spall-strength "
                    "energy scale of the dG2 barrier.");
            if (v.nuc.g1_ref < 0.0)
                Util::Abort(INFO,
                    "use_stress=4 requires nucleation.g1_ref >= 0 "
                    "(it is a barrier height, in kB*theta units).");
            if (v.nuc.m_exp <= 0.0)
                Util::Abort(INFO,
                    "use_stress=4 requires nucleation.m_exp > 0.");

            if (!pp.contains("heatmap.ic.type"))
                Util::Abort(INFO,
                    "use_stress=4 needs a heatmap: the dG1 barrier "
                    "g1_ref*(1-h) reads the ML heatmap field. Set "
                    "heatmap.ic.type, or use a mode that does not "
                    "require the heatmap.");
            if (v.nuc.g1_ref == 0.0)
                Util::Warning(INFO,
                    "use_stress=4 with nucleation.g1_ref = 0: the heatmap "
                    "barrier dG1 is identically zero, so nucleation is "
                    "gated by the spall-strength barrier alone (the "
                    "heatmap still acts as a hard mask via "
                    "heatmap_threshold).");
        }

        if (v.nuc.use_stress == 5)
        {
            if (v.nuc.g1_ref < 0.0)
                Util::Abort(INFO,
                    "use_stress=5 (heatmap-only) requires nucleation.g1_ref "
                    ">= 0 (it is the heatmap barrier height in kB*theta units; "
                    "the per-cell rate is base_rate*exp(-g1_ref*(1-h))).");

            if (v.nuc.p_thresh <= 0.0)
                Util::Abort(INFO,
                    "use_stress=5 (heatmap-only) requires nucleation.p_thresh "
                    "> 0. Mode 5 has no dG2 spall-strength barrier, so the "
                    "binary tension gate p_tens > p_thresh is the only stress "
                    "control; with p_thresh = 0 the model nucleates in any "
                    "weakly-tensile cell, including shock-transit and release "
                    "transients, yielding spurious voids during compression. "
                    "Set nucleation.p_thresh to a physical spall-onset tension "
                    "(a few GPa, e.g. 3.0e9 Pa, below the ~5 GPa spall "
                    "strength but above sub-GPa shock-transit transients).");

            if (!pp.contains("heatmap.ic.type"))
                Util::Abort(INFO,
                    "use_stress=5 (heatmap-only) needs a heatmap: the rate is "
                    "base_rate*exp(-g1_ref*(1-h)) and reads the ML heatmap "
                    "field h. Without a heatmap h defaults to 1, making the "
                    "rate a spatially uniform base_rate -- which defeats the "
                    "purpose of this mode. Set heatmap.ic.type.");

            if (v.nuc.g1_ref == 0.0)
                Util::Warning(INFO,
                    "use_stress=5 with nucleation.g1_ref = 0: the heatmap "
                    "barrier is identically zero, so the rate is a spatially "
                    "uniform base_rate (the heatmap still acts as a hard mask "
                    "via heatmap_threshold, but provides no rate contrast). "
                    "Raise g1_ref to obtain heatmap discrimination: the "
                    "contrast between the highest- and lowest-h eligible cells "
                    "is exp(g1_ref*(h_max-h_min)) on the threshold-shifted "
                    "heatmap.");
        }

        pp_query_default("seed_void.enabled", v.seed.enabled, 0);
        pp_query_default("seed_void.x",       v.seed.x, 0.0);
        pp_query_default("seed_void.y",       v.seed.y, 0.0);
        pp_query_default("seed_void.radius",  v.seed.radius, 0.0);

        v.bc_d             = new BC::Constant(1, pp, "d.bc.constant");
        v.bc_nothing_scalar = new BC::Nothing();

        int nghost = 3;

        v.RegisterNewFab(v.d_mf,         v.bc_d,             1, nghost, "d",         true,  true);
        v.RegisterNewFab(v.d_old_mf,     v.bc_d,             1, nghost, "d_old",     false, false);
        v.RegisterNewFab(v.p_drv_mf,     v.bc_nothing_scalar, 1, nghost, "p_drv",     true,  false);
        v.RegisterNewFab(v.pr_mf,        v.bc_nothing_scalar, 1, nghost, "pr",        true,  false);
        v.RegisterNewFab(v.nucleated_mf, v.bc_nothing_scalar, 1, nghost, "nucleated", true,  false);
        v.RegisterNewFab(v.psi_plus_mf,  v.bc_nothing_scalar, 1, nghost, "psi_plus",  true,  false);

        v.RegisterNewFab(v.hist_mf,      v.bc_nothing_scalar, 1, nghost, "hist",      true,  false);

        v.RegisterNewFab(v.edot_p_eq_mf, v.bc_nothing_scalar, 1, nghost, "edot_p_eq", true,  false);

        v.RegisterNewFab(v.triaxiality_mf, v.bc_nothing_scalar, 1, nghost, "triaxiality", true, false);

        
        v.RegisterNodalFab(v.ddw_diag_mf, 3, nghost, "ddw_diag", true);

        v.RegisterNodalFab(v.d_node_mf, 1, nghost, "d_node", true);

        if (pp.contains("heatmap.ic.type"))
        {
            v.RegisterNodalFab(v.heatmap_mf, 3, nghost, "heatmap", true);
            std::string type;
            pp_query("heatmap.ic.type", type);
            if (type == "hm")
            {
                v.ic_heatmap = new IC::PNG3(v.geom, pp, "heatmap.ic.hm");
                v.use_heatmap = true;
                Util::Message(INFO, "LagrangianShockCPDamage: heatmap loaded from PNG");
            }
            else
            {
                Util::Abort(INFO, "Invalid heatmap.ic.type: ", type);
            }
        }

        {
            IO::ParmParse pp_model("model1");
            pp_model.query("C11", v.C11_pristine);
            pp_model.query("C12", v.C12_pristine);
            pp_model.query("C44", v.C44_pristine);

        }
        if (v.C11_pristine <= 0.0 || v.C44_pristine <= 0.0)
        {
            Util::Warning(INFO,
                "LagrangianShockCPDamage: model1.C11/C44 not found or zero; "
                "damage stiffness degradation will drive moduli to zero. "
                "Set model1.C11/C12/C44 in the input deck.");
        }

        v.RegisterIntegratedVariable(&v.diag_max_d,      "max_damage");
        v.RegisterIntegratedVariable(&v.diag_max_p_tens, "max_tensile_p");
        v.RegisterIntegratedVariable(&v.diag_max_psi_plus, "max_psi_plus");

        if (amrex::ParallelDescriptor::IOProcessor())
        {
            const Set::Scalar dx0 = v.geom[0].CellSize(0);
            const Set::Scalar dy0 = v.geom[0].CellSize(1);
            const Set::Scalar dx_min = std::min(dx0, dy0);


            const bool is_at1 = (v.dmg.model == DamageModel::AT1);
            const Set::Scalar ell = is_at1
                ? v.dmg.ell
                : ((v.dmg.cp > 0.0)
                       ? std::sqrt(v.dmg.kappa / v.dmg.cp) : 0.0);

            Util::Message(INFO,
                "LagrangianShockCPDamage [mesh sanity]: ",
                (is_at1 ? "model=AT1 " : "model=double_well "),
                "dx0=", dx0, " dy0=", dy0,
                " void_radius=", v.nuc.void_radius,
                " (void_radius/dx_min=", v.nuc.void_radius / dx_min, ")",
                (is_at1
                    ? " Gc=" : " barrier_height="),
                (is_at1 ? v.dmg.Gc : v.dmg.barrier_height),
                (is_at1
                    ? " AT1_threshold=Gc/(cw*ell)="
                    : " cp=16*barrier_height="),
                (is_at1
                    ? (v.dmg.Gc / (v.dmg.cw * v.dmg.ell))
                    : v.dmg.cp),
                (is_at1
                    ? " regularization ell="
                    : " AC interface ell=sqrt(kappa/cp)="),
                ell,
                " (ell/dx_min=", ell / dx_min, ")",
                " d_void=", v.dmg.d_void,
                " void_res=", v.dmg.void_res);

            if (v.nuc.void_radius < dx_min)
            {
                Util::Warning(INFO,
                    "nucleation.void_radius (", v.nuc.void_radius,
                    ") is smaller than the coarsest cell size (", dx_min,
                    "). The stamped void core is now guaranteed to lock at "
                    "least the centre-containing cell to d=1 (resolution-"
                    "independent), but the void will still be a single cell. "
                    "For a resolved multi-cell void set void_radius >= 2*dx, "
                    "or refine the mesh (amr.max_level).");
            }

            if (ell > 3.0 * dx_min)
            {
                if (is_at1)
                    Util::Warning(INFO,
                        "AT1 regularization ell=", ell, " spans ",
                        ell / dx_min, " cells. This will visibly enlarge "
                        "stamped voids. Target ell ~ (1-2)*dx by lowering "
                        "damage.ell toward ", 2.0 * dx_min,
                        " (and rescaling damage.Gc to hold the threshold "
                        "Gc/(cw*ell) fixed if a fixed spall strength is "
                        "intended).");
                else
                    Util::Warning(INFO,
                        "Allen-Cahn equilibrium interface ell=", ell,
                        " spans ", ell / dx_min, " cells. This will visibly "
                        "enlarge stamped voids. Target ell ~ (1-2)*dx by "
                        "lowering damage.kappa such that "
                        "kappa = cp*(2*dx)^2 = ",
                        v.dmg.cp * 4.0 * dx_min * dx_min);
            }
            else if (ell > 0.0 && ell < 0.5 * dx_min)
            {
                if (is_at1)
                    Util::Warning(INFO,
                        "AT1 regularization ell=", ell,
                        " is under half a cell (ell/dx_min=", ell / dx_min,
                        "). The diffuse front is unresolved; the "
                        "propagating front and coalescence behaviour will "
                        "be mesh-dominated. Raise damage.ell toward ",
                        dx_min, ".");
                else
                    Util::Warning(INFO,
                        "Allen-Cahn equilibrium interface ell=", ell,
                        " is under half a cell (ell/dx_min=", ell / dx_min,
                        "). The diffuse front is unresolved; the propagating "
                        "front and coalescence behaviour will be mesh-dominated. "
                        "Raise damage.kappa toward kappa = cp*(dx)^2 = ",
                        v.dmg.cp * dx_min * dx_min);
            }
        }

        const char* mode_name = "unknown";
        switch (v.nuc.use_stress)
        {
            case 0: mode_name = "ungated";          break;
            case 1: mode_name = "binary";           break;
            case 2: mode_name = "linear";           break;
            case 3: mode_name = "Boltzmann";        break;
            case 4: mode_name = "dual-barrier";     break;
            case 5: mode_name = "heatmap-only";     break;
        }

        Util::Message(INFO,
            "LagrangianShockCPDamage: model=",
            (v.dmg.model == DamageModel::AT1 ? "at1" : "double_well"),
            " M=",  v.dmg.M,
            " kappa=",                      v.dmg.kappa,
            " barrier_height=",             v.dmg.barrier_height,
            " cp=",                         v.dmg.cp,
            " Gc=",                         v.dmg.Gc,
            " ell=",                        v.dmg.ell,
            " cw=",                         v.dmg.cw,
            " d_void=",                     v.dmg.d_void,
            " void_res=",                   v.dmg.void_res,
            " hist_seed_factor=",           v.dmg.hist_seed_factor,
            " growth.M=",                   v.grow.M,
            " growth.propagate=",           (v.grow.propagate ? 1 : 0),
            " growth.T_cap=",               v.grow.T_cap,
            " nuc.enable=",                 v.nuc.enable,
            " nuc.mode=",                   mode_name,
            "(", v.nuc.use_stress, ")",
            " nuc.base_rate=",              v.nuc.base_rate,
            " nuc.p_thresh=",               v.nuc.p_thresh,
            " nuc.p_ref=",                  v.nuc.p_ref,
            " nuc.hm_alpha=",               v.nuc.hm_alpha,
            " nuc.hm_beta=",                v.nuc.hm_beta,
            " nuc.sigma_spall=",            v.nuc.sigma_spall,
            " nuc.g1_ref=",                 v.nuc.g1_ref,
            " nuc.m_exp=",                  v.nuc.m_exp,
            " nuc.max_voids=",              v.nuc.max_voids);
    }

    void Initialize(int lev) override
    {
        Base::Initialize(lev);

        d_mf[lev]->setVal(0.0);
        d_old_mf[lev]->setVal(0.0);
        p_drv_mf[lev]->setVal(0.0);
        pr_mf[lev]->setVal(0.0);
        nucleated_mf[lev]->setVal(0.0);
        ddw_diag_mf[lev]->setVal(0.0);   // ### SHADOWING FIX: was ddw_mf
        psi_plus_mf[lev]->setVal(0.0);
        hist_mf[lev]->setVal(0.0);
        edot_p_eq_mf[lev]->setVal(0.0);
        triaxiality_mf[lev]->setVal(0.0);

        if (use_heatmap && ic_heatmap)
        {
            heatmap_mf[lev]->setVal(0.0);
            ic_heatmap->Initialize(lev, heatmap_mf);
        }

        if (seed.enabled > 0 && seed.radius > 0.0)
        {
            const Set::Scalar* DX      = geom[lev].CellSize();
            const Set::Scalar* prob_lo = geom[lev].ProbLo();
            const Set::Scalar r0       = seed.radius;
            const Set::Scalar cx       = seed.x;
            const Set::Scalar cy       = seed.y;
            const Set::Scalar eps      = nuc.eps;

            for (amrex::MFIter mfi(*d_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.tilebox();
                auto const& d = d_mf[lev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Scalar x = prob_lo[0] + (i + 0.5) * DX[0];
                    Set::Scalar y = prob_lo[1] + (j + 0.5) * DX[1];
                    Set::Scalar r = std::sqrt((x - cx)*(x - cx) + (y - cy)*(y - cy));
                    if (r <= r0 + 3.0 * eps)
                    {
                        Set::Scalar m = 0.5 * (1.0 + std::erf((-r + r0) / eps));
                        d(i,j,k) = d(i,j,k) + m * (1.0 - d(i,j,k));
                    }
                });
            }
            
        }

        FillDHalo(lev, 0.0);   
        amrex::MultiFab::Copy(*d_old_mf[lev], *d_mf[lev], 0, 0, 1,
                              d_mf[lev]->nGrow());
        d_old_mf[lev]->FillBoundary(geom[lev].periodicity());

        d_node_mf[lev]->setVal(0.0);
        FillDNode(lev, 0.0);
    }

    void TimeStepBegin(Set::Scalar a_time, int a_step) override
    {
        BL_PROFILE("LagrangianShockCPDamage::TimeStepBegin");
        Base::TimeStepBegin(a_time, a_step);

        if (nuc.enable <= 0) return;
        if (nuc.check_freq <= 0) return;
        if (a_step % nuc.check_freq != 0) return;

        if (amrex::ParallelDescriptor::IOProcessor())
            Util::Message(INFO,
                "[NUC] step=", a_step, " t=", a_time,
                " n_voids_total=", n_voids_total, "/", nuc.max_voids);

        NucleationCheck(a_time, a_step);
    }

    void Advance(int lev, Set::Scalar time, Set::Scalar a_dt) override
    {
        BL_PROFILE("LagrangianShockCPDamage::Advance");
        if (m_type == ::Integrator::Base::Mechanics<MODEL>::Type::Disable) return;

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

        this->SyncNodal(vel_old_mf,  lev);
        this->SyncNodal(disp_old_mf, lev);
        this->AssertSynced(disp_old_mf, lev, "disp_old @ entry");

        Base::ComputePressureAndViscosity2D(lev, DX0, DX1, local_eos, local_rho0,
                                            local_CQ, local_CL, disp_old_mf, vel_old_mf,
                                            pressure_mf, artvisc2d_mf, density_mf, soundspeed_mf);

        DegradeEOSPressure(lev, pressure_mf);
        DegradeSoundSpeed(lev, soundspeed_mf);

        pressure_mf[lev]->FillBoundary(geom[lev].periodicity());
        artvisc2d_mf[lev]->FillBoundary(geom[lev].periodicity());
        soundspeed_mf[lev]->FillBoundary(geom[lev].periodicity());
        density_mf[lev]->FillBoundary(geom[lev].periodicity());

        Base::CheckEOSReferenceState(lev);


        Base::ComputeCellStressField(lev, DX0, DX1, nlo, nhi,
                                     disp_old_mf, pressure_mf, artvisc2d_mf);
        DegradeDeviatoricStress(lev, disp_old_mf);
        SuppressStressInVoids(lev);
        cell_stress_mf[lev]->FillBoundary(geom[lev].periodicity());
        AssertGhostConsistent(cell_stress_mf, lev, "cell_stress @ predictor");

        Base::UpdateVelocityAndDisplacement2D(lev, DX0, DX1, half_dt, a_dt, nlo, nhi, ndomain,
            local_rho0, vel_old_mf, disp_old_mf, vel_old_mf, vel_mf, disp_mf);
        this->AssertSynced(vel_mf, lev, "vel after half-step update (pre-sync)");
        this->SyncNodal(vel_mf,  lev);
        this->SyncNodal(disp_mf, lev);

        Base::ComputeHalfStepEnergy2D(lev, DX0, DX1, local_eos, local_rho0, local_CQ, local_CL);

        DegradeEOSPressure(lev, pressure_half_mf);
        DegradeSoundSpeed(lev, soundspeed_half_mf);

        pressure_half_mf[lev]->FillBoundary(geom[lev].periodicity());
        artvisc2d_half_mf[lev]->FillBoundary(geom[lev].periodicity());
        soundspeed_half_mf[lev]->FillBoundary(geom[lev].periodicity());
        density_half_mf[lev]->FillBoundary(geom[lev].periodicity());


        Base::ComputeCellStressField(lev, DX0, DX1, nlo, nhi,
                                     disp_mf, pressure_half_mf, artvisc2d_half_mf);

        if (!dbg_bypass)   // ### BISECT SWITCH
        {
            amrex::MultiFab::Copy(*cell_stress_raw_mf[lev], *cell_stress_mf[lev],
                                  0, 0, cell_stress_mf[lev]->nComp(),
                                  cell_stress_mf[lev]->nGrow());
            cell_stress_raw_mf[lev]->FillBoundary(geom[lev].periodicity());
        }
        DegradeDeviatoricStress(lev, disp_mf);
        SuppressStressInVoids(lev);
        cell_stress_mf[lev]->FillBoundary(geom[lev].periodicity());
        AssertGhostConsistent(cell_stress_mf, lev, "cell_stress @ corrector");

        Base::UpdateVelocityAndDisplacement2D(lev, DX0, DX1, a_dt, a_dt, nlo, nhi, ndomain,
            local_rho0, vel_old_mf, disp_old_mf, vel_mf, vel_mf, disp_mf);
        this->AssertSynced(vel_mf, lev, "vel after full-step update (pre-sync)");
        this->SyncNodal(vel_mf,  lev);
        this->SyncNodal(disp_mf, lev);

        Base::ApplyHourglassControl(lev, DX0, DX1, a_dt, nlo, nhi, ndomain, local_rho0);
        this->AssertSynced(vel_mf, lev, "vel after hourglass control (pre-sync)");
        this->SyncNodal(vel_mf, lev);

        for (amrex::MFIter mfi(*density_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box cbx = mfi.tilebox();
            auto const& rho_arr = density_mf[lev]->array(mfi);
            auto const& u_arr = disp_mf[lev]->array(mfi);
            Set::Scalar lDX0 = DX0, lDX1 = DX1, lrho0 = local_rho0;
            amrex::ParallelFor(cbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Set::Scalar J = Base::CellJ(u_arr, i, j, k, lDX0, lDX1);
                rho_arr(i,j,k) = lrho0 / J;
            });
        }

        for (amrex::MFIter mfi(*energy_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box cbx = mfi.tilebox();
            auto const& e_arr     = energy_mf[lev]->array(mfi);
            auto const& ph_arr    = pressure_half_mf[lev]->array(mfi);
            auto const& q2dh_arr  = artvisc2d_half_mf[lev]->array(mfi);
            auto const& u_old_arr = disp_old_mf[lev]->array(mfi);
            auto const& u_new_arr = disp_mf[lev]->array(mfi);
            Set::Scalar lDX0 = DX0, lDX1 = DX1, lrho0 = local_rho0;
            amrex::ParallelFor(cbx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Scalar J_old = Base::CellJ(u_old_arr, i, j, k, lDX0, lDX1);
                Set::Scalar J_new = Base::CellJ(u_new_arr, i, j, k, lDX0, lDX1);
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

                
                Set::Scalar gam_old = 0.5*((u_old_arr(i,j+1,k)(0)-u_old_arr(i,j,k)(0))
                                          +(u_old_arr(i+1,j+1,k)(0)-u_old_arr(i+1,j,k)(0))) / lDX1
                                    + 0.5*((u_old_arr(i+1,j,k)(1)-u_old_arr(i,j,k)(1))
                                          +(u_old_arr(i+1,j+1,k)(1)-u_old_arr(i,j+1,k)(1))) / lDX0;
                Set::Scalar gam_new = 0.5*((u_new_arr(i,j+1,k)(0)-u_new_arr(i,j,k)(0))
                                          +(u_new_arr(i+1,j+1,k)(0)-u_new_arr(i+1,j,k)(0))) / lDX1
                                    + 0.5*((u_new_arr(i+1,j,k)(1)-u_new_arr(i,j,k)(1))
                                          +(u_new_arr(i+1,j+1,k)(1)-u_new_arr(i,j+1,k)(1))) / lDX0;
                Set::Scalar d_gam = gam_new - gam_old;

                Set::Scalar de = -ph_arr(i,j,k) * dJ / lrho0;
                de -= (q2dh_arr(i,j,k,0) * d_exx + q2dh_arr(i,j,k,1) * d_eyy) / lrho0;
                de += q2dh_arr(i,j,k,2) * d_gam / lrho0;   // ### STALE-FORK FIX
                if (std::isfinite(de)) e_arr(i,j,k) += de;
            });
        }
        energy_mf[lev]->FillBoundary(geom[lev].periodicity());
        AssertGhostConsistent(energy_mf,   lev, "energy @ full step");
        AssertGhostConsistent(pressure_mf, lev, "pressure @ full step");

        SuppressEnergyInVoids(lev);

        Base::AdvanceCPModelAllNodes(lev, DX0, DX1, a_dt, t_full, nlo, nhi, ndomain,
                                     disp_mf, disp_old_mf);
        this->SyncModel(lev);

        Base::UpdateTemperature(lev, DX0, DX1, a_dt, nlo, nhi, disp_old_mf, disp_mf);

        Base::PopulateVisualizationFields(lev, DX, nlo, nhi, ndomain, disp_mf, pressure_mf);

        Set::Scalar dt_local = 1e100, max_p = -1e100, min_p = 1e100;
        Set::Scalar fs_vel = 0.0, fs_cnt = 0.0;

        Base::ComputePressureAndViscosity2D(lev, DX0, DX1, local_eos, local_rho0,
                                            local_CQ, local_CL, disp_mf, vel_mf,
                                            pressure_mf, artvisc2d_mf, density_mf, soundspeed_mf);

        DegradeEOSPressure(lev, pressure_mf);
        DegradeSoundSpeed(lev, soundspeed_mf);

        pressure_mf[lev]->FillBoundary(geom[lev].periodicity());
        soundspeed_mf[lev]->FillBoundary(geom[lev].periodicity());
        this->SyncNodal(vel_mf, lev);

        for (amrex::MFIter mfi(*vel_mf[lev], false); mfi.isValid(); ++mfi)
        {
            amrex::Box nbx = mfi.nodaltilebox() & ndomain;
            auto const& v_arr = vel_mf[lev]->array(mfi);
            if (nbx.contains(amrex::IntVect(AMREX_D_DECL(nhi.x, nlo.y, 0)))) {
                fs_vel += v_arr(nhi.x, nlo.y, 0)(0);
                fs_cnt += 1.0;
            }
        }

        amrex::ParallelDescriptor::ReduceRealSum(fs_vel);
        amrex::ParallelDescriptor::ReduceRealSum(fs_cnt);
        if (fs_cnt > 0.0) fs_vel /= fs_cnt;

        Set::Scalar min_dx = std::min(DX0, DX1);
        for (amrex::MFIter mfi(*density_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box cbx = mfi.tilebox();
            auto const& p_arr = pressure_mf[lev]->array(mfi);
            auto const& c_arr = soundspeed_mf[lev]->array(mfi);
            auto const& v_arr = vel_mf[lev]->array(mfi);
            amrex::Loop(cbx, [&](int i, int j, int k) {
                Set::Scalar c = c_arr(i,j,k);
                Set::Scalar dvx = std::abs(0.5*((v_arr(i+1,j,k)(0)-v_arr(i,j,k)(0))
                                               +(v_arr(i+1,j+1,k)(0)-v_arr(i,j+1,k)(0))));
                Set::Scalar dvy = std::abs(0.5*((v_arr(i,j+1,k)(1)-v_arr(i,j,k)(1))
                                               +(v_arr(i+1,j+1,k)(1)-v_arr(i+1,j,k)(1))));
                Set::Scalar dv  = std::max(dvx, dvy);
                Set::Scalar dt_cfl = hydro_cfl * min_dx / (c + local_CQ * dv + 1e-20);
                dt_local = std::min(dt_local, dt_cfl);
                if (p_arr(i,j,k) > max_p) max_p = p_arr(i,j,k);
                if (p_arr(i,j,k) < min_p) min_p = p_arr(i,j,k);
            });
        }

        amrex::ParallelDescriptor::ReduceRealMax(max_p);
        amrex::ParallelDescriptor::ReduceRealMin(min_p);
        amrex::ParallelDescriptor::ReduceRealMin(dt_local);

        diag_fs_vel = fs_vel; diag_max_pressure = max_p; diag_min_pressure = min_p;

        if (dt_local > 0.0 && dt_local < 1e90)
            this->DynamicTimestep_SyncTimeStep(lev, dt_local);

        if (lev == 0 && t_full >= next_diag_time) {
            if (amrex::ParallelDescriptor::IOProcessor())
                this->WriteDiagnostics(t_full, a_dt);
            next_diag_time = t_full + diag_dt;
        }


        PopulatePlasticStrainRate(lev);
        ComputePsiPlus(lev, disp_mf);
        UpdateHistory(lev);
        ComputeTriaxiality(lev);
        AdvanceDamage(lev, time, a_dt);
        AssertGhostConsistent(d_mf,        lev, "d after AdvanceDamage");
        AssertGhostConsistent(psi_plus_mf, lev, "psi_plus after AdvanceDamage");
        AssertGhostConsistent(hist_mf,     lev, "hist after AdvanceDamage");
        DegradeModuli(lev);
        PopulateDDWDiagnostic(lev);
        UpdateDamageDiagnostics(lev);

        if (lev == 0 && t_full >= next_xt_time)
        {
            for (std::size_t q = 0; q < xt_modes.size(); ++q)
                WriteXTLineout(lev, t_full, xt_modes[q], xt_files[q],
                               xt_header_written[q]);
            next_xt_time = t_full + xt_dt;
        }
    }


    void WriteXTLineout(int lev, Set::Scalar time, XTReduce mode,
                       const std::string& fname, char& header_written)
    {
        BL_PROFILE("LagrangianShockCPDamage::WriteXTLineout");
        if (dbg_bypass) return;   // 

        if (xt_enable <= 0) return;

        const amrex::Box& cdom = geom[lev].Domain();
        const int ilo = cdom.smallEnd(0);
        const int ihi = cdom.bigEnd(0);
        const int Nx  = ihi - ilo + 1;
        const int Ny  = cdom.length(1);
        const int jc  = (xt_jprobe >= 0)
                      ? xt_jprobe
                      : (cdom.smallEnd(1) + cdom.bigEnd(1)) / 2;

        const Set::Scalar DX0 = geom[lev].CellSize(0);
        const Set::Scalar x0  = geom[lev].ProbLo(0);

        const bool is_max = (mode == XTReduce::Max);

        const Set::Scalar ident = is_max
            ? -std::numeric_limits<Set::Scalar>::max()
            : 0.0;

        std::vector<Set::Scalar> vx(Nx, ident), ux(Nx, ident),
                                 pr(Nx, ident), rh(Nx, ident), dm(Nx, ident);

        for (amrex::MFIter mfi(*d_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& cbx = mfi.tilebox();
            auto const& p_arr = pressure_mf[lev]->const_array(mfi);
            auto const& r_arr = density_mf [lev]->const_array(mfi);
            auto const& d_arr = d_mf       [lev]->const_array(mfi);
            auto const& v_arr = vel_mf     [lev]->const_array(mfi);   // nodal
            auto const& u_arr = disp_mf    [lev]->const_array(mfi);   // nodal

            amrex::Loop(cbx, [&](int i, int j, int k)
            {
                if (k != 0) return;
                if (mode == XTReduce::Row && j != jc) return;
                const int q = i - ilo;
                if (q < 0 || q >= Nx) return;

                const Set::Scalar vxc =
                    0.25 * (v_arr(i  , j  , k)(0) + v_arr(i+1, j  , k)(0)
                          + v_arr(i  , j+1, k)(0) + v_arr(i+1, j+1, k)(0));
                const Set::Scalar uxc =
                    0.25 * (u_arr(i  , j  , k)(0) + u_arr(i+1, j  , k)(0)
                          + u_arr(i  , j+1, k)(0) + u_arr(i+1, j+1, k)(0));
                const Set::Scalar pc = p_arr(i, j, k);
                const Set::Scalar rc = r_arr(i, j, k);
                const Set::Scalar dc = d_arr(i, j, k);

                if (is_max)
                {
                    vx[q] = std::max(vx[q], vxc);
                    ux[q] = std::max(ux[q], uxc);
                    pr[q] = std::max(pr[q], pc);
                    rh[q] = std::max(rh[q], rc);
                    dm[q] = std::max(dm[q], dc);
                }
                else
                {
                    vx[q] += vxc;
                    ux[q] += uxc;
                    pr[q] += pc;
                    rh[q] += rc;
                    dm[q] += dc;
                }
            });
        }

        if (is_max)
        {
            amrex::ParallelDescriptor::ReduceRealMax(vx.data(), Nx);
            amrex::ParallelDescriptor::ReduceRealMax(ux.data(), Nx);
            amrex::ParallelDescriptor::ReduceRealMax(pr.data(), Nx);
            amrex::ParallelDescriptor::ReduceRealMax(rh.data(), Nx);
            amrex::ParallelDescriptor::ReduceRealMax(dm.data(), Nx);
        }
        else
        {
            amrex::ParallelDescriptor::ReduceRealSum(vx.data(), Nx);
            amrex::ParallelDescriptor::ReduceRealSum(ux.data(), Nx);
            amrex::ParallelDescriptor::ReduceRealSum(pr.data(), Nx);
            amrex::ParallelDescriptor::ReduceRealSum(rh.data(), Nx);
            amrex::ParallelDescriptor::ReduceRealSum(dm.data(), Nx);
        }

        if (!amrex::ParallelDescriptor::IOProcessor()) return;

        if (mode == XTReduce::Mean)
        {
            const Set::Scalar inv_Ny = 1.0 / (Set::Scalar)Ny;
            for (int q = 0; q < Nx; ++q)
            {
                vx[q] *= inv_Ny;  ux[q] *= inv_Ny;  pr[q] *= inv_Ny;
                rh[q] *= inv_Ny;  dm[q] *= inv_Ny;
            }
        }

        std::ofstream out;
        if (!header_written)
        {
            const char* red_name =
                (mode == XTReduce::Mean) ? "mean_y"
              : (mode == XTReduce::Max)  ? "max_y" : "row";

            out.open(fname, std::ios::out);
            out << "# Lagrangian shock x-t lineout\n";
            out << "# reduce = " << red_name;
            if (mode == XTReduce::Row) out << " (probe row j = " << jc << ")";
            else                            out << " (over Ny = " << Ny << " rows)";
            out << ", level = " << lev << ", Nx = " << Nx << "\n";
            if (mode == XTReduce::Mean)
                out << "# mean_y: reference-area (= mass-weighted) transverse "
                       "average. Col 8 (d) is the VOID AREA FRACTION at depth x.\n";
            if (mode == XTReduce::Max)
                out << "# max_y: cols 5-7 (vx,p,rho) are transverse maxima and "
                       "are NOT physical profiles; col 8 (d) is the void indicator.\n";
            out << "# X_ref = material (reference) coordinate; "
                   "x_spatial = X_ref + u_x = deformed position\n";
            out << "# Col 1: time (s)\n";
            out << "# Col 2: cell index i\n";
            out << "# Col 3: X_ref (m)\n";
            out << "# Col 4: x_spatial (m)\n";
            out << "# Col 5: vx (m/s)\n";
            out << "# Col 6: pressure (Pa)\n";
            out << "# Col 7: density (kg/m^3)\n";
            out << "# Col 8: damage d (-)\n";
            out << std::scientific << std::setprecision(8);
            header_written = 1;
        }
        else
        {
            out.open(fname, std::ios::app);
            out << std::scientific << std::setprecision(8);
        }

        for (int q = 0; q < Nx; ++q)
        {
            const int i = ilo + q;
            const Set::Scalar Xref = x0 + ((Set::Scalar)i + 0.5) * DX0;
            out << time           << "  "
                << i              << "  "
                << Xref           << "  "
                << (Xref + ux[q]) << "  "
                << vx[q]          << "  "
                << pr[q]          << "  "
                << rh[q]          << "  "
                << dm[q]          << "\n";
        }
        out.close();
    }


    void PopulatePlasticStrainRate(int lev)
    {
        BL_PROFILE("LagrangianShockCPDamage::PopulatePlasticStrainRate");
        if (dbg_bypass) return;   // ### BISECT SWITCH

        if constexpr (MODEL::kinvar != Model::Solid::KinematicVariable::F)
        {
            return;
        }
        else
        {
            for (amrex::MFIter mfi(*edot_p_eq_mf[lev], amrex::TilingIfNotGPU());
                 mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.tilebox();
                auto const& ep_arr    = edot_p_eq_mf[lev]->array(mfi);
                auto const& model_arr = model_mf[lev]->const_array(mfi);

                const auto lo = amrex::lbound(bx);
                const auto hi = amrex::ubound(bx);
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    int k = 0;
                    Set::Scalar e = 0.25 * (
                        model_arr(i  , j  , k).edot_p_eq
                      + model_arr(i+1, j  , k).edot_p_eq
                      + model_arr(i  , j+1, k).edot_p_eq
                      + model_arr(i+1, j+1, k).edot_p_eq);
                    if (!std::isfinite(e) || e < 0.0) e = 0.0;
                    ep_arr(i,j,k) = e;
                }}
            }
            edot_p_eq_mf[lev]->FillBoundary(geom[lev].periodicity());

            Set::Scalar max_ep = edot_p_eq_mf[lev]->max(0);
            static int print_counter = 0;
            if (amrex::ParallelDescriptor::IOProcessor()
                && (print_counter++ % 50 == 0))
            {
                Util::Message(INFO,
                    "[edot_p_eq] lev=", lev,
                    " max=", max_ep, " [1/s]");
            }
        }
    }


    void UpdateHistory(int lev)
    {
        BL_PROFILE("LagrangianShockCPDamage::UpdateHistory");
        if (dbg_bypass) return;   // ### BISECT SWITCH

        const Set::Scalar lhsf = dmg.hist_seed_factor;

        for (amrex::MFIter mfi(*hist_mf[lev], amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            auto const& H_arr   = hist_mf[lev]->array(mfi);
            auto const& psi_arr = psi_plus_mf[lev]->const_array(mfi);
            auto const& nuc_arr = nucleated_mf[lev]->const_array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Scalar psi = psi_arr(i,j,k);
                if (!std::isfinite(psi) || psi < 0.0) psi = 0.0;

                Set::Scalar Hnew = H_arr(i,j,k);
                if (psi > Hnew) Hnew = psi;

                if (lhsf > 0.0 && nuc_arr(i,j,k) > 0.4)
                {
                    Set::Scalar Hfloor = lhsf * psi;
                    if (Hfloor > Hnew) Hnew = Hfloor;
                }

                if (!std::isfinite(Hnew) || Hnew < 0.0) Hnew = 0.0;
                H_arr(i,j,k) = Hnew;
            });
        }
        hist_mf[lev]->FillBoundary(geom[lev].periodicity());
    }



    void TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags,
                               Set::Scalar a_time, int a_ngrow) override
    {
        Base::TagCellsForRefinement(lev, a_tags, a_time, a_ngrow);

        if (dbg_bypass) return;   // ### BISECT SWITCH
        if (m_type == ::Integrator::Base::Mechanics<MODEL>::Type::Disable) return;

        Set::Vector DX(geom[lev].CellSize());
        Set::Scalar DXnorm = DX.lpNorm<2>();
        Set::Scalar threshold = m_eta_ref_threshold;

        //
        
        //
        const amrex::Box& cdomain = geom[lev].Domain();

        for (amrex::MFIter mfi(*d_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            auto const& tags = a_tags.array(mfi);
            auto const& d    = d_mf[lev]->const_array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                auto sten = Numeric::GetStencil(i, j, k, cdomain);
                Set::Vector g = Numeric::Gradient(d, i, j, k, 0, DX.data(), sten);
                if (g.lpNorm<2>() * DXnorm > threshold)
                    tags(i, j, k) = amrex::TagBox::SET;
            });
        }
    }

    void Regrid(int lev, Set::Scalar time) override
    {
        Base::Regrid(lev, time);
        if (use_heatmap && ic_heatmap)
            ic_heatmap->Initialize(lev, heatmap_mf, time);
    }

private:


    void DegradeDeviatoricStress(int lev,
                                 Set::Field<Set::Vector>& u_field)
    {
        BL_PROFILE("LagrangianShockCPDamage::DegradeDeviatoricStress");
        if (dbg_bypass) return;   // 

        const Set::Scalar lsmall = dmg.small;
        const Set::Scalar* DX    = geom[lev].CellSize();
        const Set::Scalar lDX0   = DX[0];
        const Set::Scalar lDX1   = DX[1];

        for (amrex::MFIter mfi(*cell_stress_mf[lev], amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            amrex::Box cbx   = mfi.tilebox();
            auto const& cs   = cell_stress_mf[lev]->array(mfi);
            auto const& u    = u_field[lev]->array(mfi);
            auto const& d_arr = d_mf[lev]->const_array(mfi);

            amrex::ParallelFor(cbx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Scalar dv = d_arr(i, j, k);
                if (dv <= 0.0) return;
                if (dv > 1.0) dv = 1.0;

                Set::Scalar g = (1.0 - dv) * (1.0 - dv);
                if (g < lsmall) g = lsmall;


                Set::Matrix Fc;
                Base::CellF2D(u, i, j, k, lDX0, lDX1, Fc);
                Set::Scalar Jc = Fc(0,0)*Fc(1,1) - Fc(0,1)*Fc(1,0);
                if (std::abs(Jc) < 1e-10) Jc = (Jc >= 0.0) ? 1e-10 : -1e-10;


                Set::Matrix Pcorr;
                Pcorr(0,0) = cs(i,j,k,0);
                Pcorr(1,1) = cs(i,j,k,1);
                Pcorr(0,1) = cs(i,j,k,2);
                Pcorr(1,0) = cs(i,j,k,3);


                Set::Matrix sigma = (1.0 / Jc) * Pcorr * Fc.transpose();


                Set::Scalar p_mean_neg = 0.5 * (sigma(0,0) + sigma(1,1));
                Set::Matrix s = sigma;
                s(0,0) -= p_mean_neg;
                s(1,1) -= p_mean_neg;


                s *= g;


                sigma = s;
                sigma(0,0) += p_mean_neg;
                sigma(1,1) += p_mean_neg;


                Set::Matrix FinvT;
                FinvT(0,0) =  Fc(1,1) / Jc; FinvT(0,1) = -Fc(1,0) / Jc;
                FinvT(1,0) = -Fc(0,1) / Jc; FinvT(1,1) =  Fc(0,0) / Jc;
                Set::Matrix Pnew = Jc * sigma * FinvT;

                cs(i,j,k,0) = Pnew(0,0);
                cs(i,j,k,1) = Pnew(1,1);
                cs(i,j,k,2) = Pnew(0,1);
                cs(i,j,k,3) = Pnew(1,0);
            });
        }
    }


    void DegradeEOSPressure(int lev, Set::Field<Set::Scalar>& p_field)
    {
        BL_PROFILE("LagrangianShockCPDamage::DegradeEOSPressure");
        if (dbg_bypass) return;   // 

        const Set::Scalar lsmall    = dmg.small;
        const Set::Scalar lvoid_res = dmg.void_res;   
        const Set::Scalar lnuc_lock = 0.9;

        for (amrex::MFIter mfi(*p_field[lev], amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            auto const& p_arr   = p_field[lev]->array(mfi);
            auto const& d_arr   = d_mf[lev]->const_array(mfi);
            auto const& nuc_arr = nucleated_mf[lev]->const_array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                
                if (nuc_arr(i, j, k) > lnuc_lock)
                {
                    if (p_arr(i, j, k) < 0.0) p_arr(i, j, k) *= lvoid_res;
                    return;
                }

                Set::Scalar dv = d_arr(i, j, k);
                if (dv <= 0.0) return;
                if (dv > 1.0) dv = 1.0;

                // Compression carried unchanged
                if (p_arr(i, j, k) >= 0.0) return;

                // Continuous tensile softening 
                Set::Scalar g = (1.0 - dv) * (1.0 - dv);
                if (g < lsmall) g = lsmall;

                p_arr(i, j, k) *= g;
            });
        }
    }


    void DegradeSoundSpeed(int lev, Set::Field<Set::Scalar>& c_field)
    {
        BL_PROFILE("LagrangianShockCPDamage::DegradeSoundSpeed");
        if (dbg_bypass) return;   // 

        const Set::Scalar lsmall    = dmg.small;
        const Set::Scalar lvoid_res = dmg.void_res;   // 
        const Set::Scalar lnuc_lock = 0.9;

        for (amrex::MFIter mfi(*c_field[lev], amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            auto const& c_arr   = c_field[lev]->array(mfi);
            auto const& d_arr   = d_mf[lev]->const_array(mfi);
            auto const& nuc_arr = nucleated_mf[lev]->const_array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                
                if (nuc_arr(i, j, k) > lnuc_lock)
                {
                    c_arr(i, j, k) *= lvoid_res;
                    return;
                }

                Set::Scalar dv = d_arr(i, j, k);
                if (dv <= 0.0) return;
                if (dv > 1.0) dv = 1.0;
                Set::Scalar g = 1.0 - dv;    // linear, not (1-d)^2
                if (g < lsmall) g = lsmall;
                c_arr(i, j, k) *= g;
            });
        }
    }


    void SuppressStressInVoids(int lev)
    {
        BL_PROFILE("LagrangianShockCPDamage::SuppressStressInVoids");
        if (dbg_bypass) return;   

  
        const Set::Scalar lnuc_lock = 0.9;
        const Set::Scalar lvoid_res = dmg.void_res;   

        for (amrex::MFIter mfi(*cell_stress_mf[lev], amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            amrex::Box cbx = mfi.tilebox();
            auto const& cs      = cell_stress_mf[lev]->array(mfi);
            auto const& nuc_arr = nucleated_mf[lev]->const_array(mfi);

            amrex::ParallelFor(cbx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
               
                if (nuc_arr(i, j, k) > lnuc_lock)
                {
                    cs(i, j, k, 0) *= lvoid_res;
                    cs(i, j, k, 1) *= lvoid_res;
                    cs(i, j, k, 2) *= lvoid_res;
                    cs(i, j, k, 3) *= lvoid_res;
                }
            });
        }
    }

    void SuppressEnergyInVoids(int lev)
    {
        BL_PROFILE("LagrangianShockCPDamage::SuppressEnergyInVoids");
        if (dbg_bypass) return;   // ### BISECT SWITCH

       
        const Set::Scalar lnuc_lock = 0.9;

        for (amrex::MFIter mfi(*energy_mf[lev], amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            auto const& e_arr   = energy_mf[lev]->array(mfi);
            auto const& nuc_arr = nucleated_mf[lev]->const_array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
               
                if (nuc_arr(i, j, k) > lnuc_lock)
                    e_arr(i, j, k) = 0.0;
            });
        }
    }



    void ComputeTriaxiality(int lev)
    {
        BL_PROFILE("LagrangianShockCPDamage::ComputeTriaxiality");
        if (dbg_bypass) return;   // 

        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar lDX0 = DX[0];
        const Set::Scalar lDX1 = DX[1];


        const Set::Scalar svm_floor = 1.0e-6 * stress_scale;

        for (amrex::MFIter mfi(*triaxiality_mf[lev], amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            auto const& T_arr  = triaxiality_mf[lev]->array(mfi);
            auto const& cs_raw = cell_stress_raw_mf[lev]->const_array(mfi);
            auto const& u      = disp_mf[lev]->const_array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {

                Set::Matrix Pcorr;
                Pcorr(0,0) = cs_raw(i,j,k,0);
                Pcorr(1,1) = cs_raw(i,j,k,1);
                Pcorr(0,1) = cs_raw(i,j,k,2);
                Pcorr(1,0) = cs_raw(i,j,k,3);


                Set::Matrix Fc;
                Base::CellF2D(u, i, j, k, lDX0, lDX1, Fc);
                Set::Scalar Jc = Fc(0,0)*Fc(1,1) - Fc(0,1)*Fc(1,0);
                if (std::abs(Jc) < 1e-10) Jc = (Jc >= 0.0) ? 1e-10 : -1e-10;
                Set::Matrix sigma = (1.0 / Jc) * Pcorr * Fc.transpose();

                Set::Scalar p_mean_neg = 0.5 * (sigma(0,0) + sigma(1,1));
                Set::Scalar p_phys = -p_mean_neg;   // p > 0 compressive

                Set::Matrix s = sigma;
                s(0,0) -= p_mean_neg;
                s(1,1) -= p_mean_neg;

                Set::Scalar s_dd = s(0,0)*s(0,0) + s(1,1)*s(1,1)
                                 + s(0,1)*s(0,1) + s(1,0)*s(1,0);
                Set::Scalar sigma_vM = std::sqrt(std::max(1.5 * s_dd, 0.0));
                if (sigma_vM < svm_floor) sigma_vM = svm_floor;


                Set::Scalar T = -p_phys / sigma_vM;
                if (!std::isfinite(T)) T = 0.0;

                T_arr(i,j,k) = T;
            });
        }
        triaxiality_mf[lev]->FillBoundary(geom[lev].periodicity());

        Set::Scalar min_T =  1.0e100;
        Set::Scalar max_T = -1.0e100;
        for (amrex::MFIter mfi(*triaxiality_mf[lev], false);
             mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            auto const& T_a = triaxiality_mf[lev]->const_array(mfi);
            amrex::Loop(bx, [&](int i, int j, int k) {
                Set::Scalar v = T_a(i,j,k);
                if (std::isfinite(v)) {
                    if (v < min_T) min_T = v;
                    if (v > max_T) max_T = v;
                }
            });
        }
        amrex::ParallelDescriptor::ReduceRealMin(min_T);
        amrex::ParallelDescriptor::ReduceRealMax(max_T);
        static int print_counter = 0;
        if (amrex::ParallelDescriptor::IOProcessor()
            && (print_counter++ % 50 == 0))
            Util::Message(INFO,
                "[triaxiality] lev=", lev,
                " T range: [", min_T, ", ", max_T, "]");
    }


    void ComputePsiPlus(int lev,
                        Set::Field<Set::Vector>& u_field)
    {
        BL_PROFILE("LagrangianShockCPDamage::ComputePsiPlus");
        if (dbg_bypass) return;   //

        const Set::Scalar* DX    = geom[lev].CellSize();
        const Set::Scalar lDX0   = DX[0];
        const Set::Scalar lDX1   = DX[1];
        const Set::Scalar lscale = stress_scale;

        pressure_mf[lev]->FillBoundary(geom[lev].periodicity());

        amrex::Box cdomain = geom[lev].Domain();
        amrex::Box ndomain = amrex::convert(cdomain, amrex::IntVect::TheNodeVector());
        const amrex::Dim3 nlo = amrex::lbound(ndomain);
        const amrex::Dim3 nhi = amrex::ubound(ndomain);
        amrex::Box cell_range(
            amrex::IntVect(AMREX_D_DECL(nlo.x, nlo.y, 0)),
            amrex::IntVect(AMREX_D_DECL(nhi.x-1, nhi.y-1, 0)));

        const auto cdlo   = amrex::lbound(cdomain);
        const auto cdhi   = amrex::ubound(cdomain);
        const int  cdlo_x = cdlo.x;
        const int  cdhi_x = cdhi.x;
        const Set::Scalar lp_floor = dmg.p_tens_floor;

        for (amrex::MFIter mfi(*psi_plus_mf[lev], amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            amrex::Box cbx    = mfi.tilebox() & cell_range;
            auto const& psi   = psi_plus_mf[lev]->array(mfi);
            auto const& u     = u_field[lev]->array(mfi);
            auto const& mdl   = model_mf[lev]->const_array(mfi);
            auto const& p_eos = pressure_mf[lev]->const_array(mfi);

            amrex::ParallelFor(cbx, [=] AMREX_GPU_DEVICE(int ci, int cj, int ck)
            {
                Set::Matrix Fc;
                Base::CellF2D(u, ci, cj, ck, lDX0, lDX1, Fc);


                Eigen::Matrix3d F3d = Eigen::Matrix3d::Identity();
                F3d(0,0) = Fc(0,0);  F3d(0,1) = Fc(0,1);
                F3d(1,0) = Fc(1,0);  F3d(1,1) = Fc(1,1);

                Set::Scalar psi_accum = 0.0;
                int n_valid = 0;

                for (int dj = 0; dj <= 1; ++dj)
                {
                    for (int di = 0; di <= 1; ++di)
                    {
                        int ni = ci + di;
                        int nj = cj + dj;

                        if constexpr (MODEL::kinvar ==
                                      Model::Solid::KinematicVariable::F)
                        {
                            const auto& m = mdl(ni, nj, ck);
                            if (m.Fp(0,0) == 0.0 || m.C11 == 0.0) continue;

                            Eigen::Matrix3d Fpinv = m.Fp.inverse();
                            Eigen::Matrix3d Fe3d  = F3d * Fpinv;

                            Eigen::Matrix3d Ce3d = Fe3d.transpose() * Fe3d;
                            Eigen::Matrix3d Ee3d = 0.5 *
                                (Ce3d - Eigen::Matrix3d::Identity());

                            
                            Eigen::Matrix3d R = m.q.normalized().toRotationMatrix();

                            Eigen::Matrix3d Ee_cryst = R.transpose() * Ee3d * R;

                            Set::Scalar c11 = m.C11;
                            Set::Scalar c12 = m.C12;
                            Set::Scalar c44 = m.C44;


                            Set::Scalar trE =
                                Ee_cryst(0,0) + Ee_cryst(1,1) + Ee_cryst(2,2);
                            Set::Scalar K_bulk = (c11 + 2.0 * c12) / 3.0;
                            Set::Scalar trE_plus =
                                (trE > 0.0) ? trE : 0.0;
                            Set::Scalar psi_vol_plus =
                                0.5 * K_bulk * trE_plus * trE_plus;


                            Set::Scalar third_trE = trE / 3.0;
                            Eigen::Matrix3d Ed = Ee_cryst;
                            Ed(0,0) -= third_trE;
                            Ed(1,1) -= third_trE;
                            Ed(2,2) -= third_trE;

                            Set::Scalar psi_dev =
                                0.5 * (c11 - c12) *
                                (Ed(0,0)*Ed(0,0) + Ed(1,1)*Ed(1,1)
                                 + Ed(2,2)*Ed(2,2))
                              + 2.0 * c44 *
                                (Ed(0,1)*Ed(0,1) + Ed(0,2)*Ed(0,2)
                                 + Ed(1,2)*Ed(1,2));

                            Set::Scalar psi_dev_gated =
                                (trE > 0.0) ? psi_dev : 0.0;

                            psi_accum += psi_vol_plus + psi_dev_gated;
                            ++n_valid;
                        }
                        else
                        {
                            (void)F3d; (void)di; (void)dj;
                            (void)ni; (void)nj;
                        }
                    }
                }

                Set::Scalar psi_here = 0.0;
                if (n_valid > 0)
                    psi_here = (psi_accum / (Set::Scalar)n_valid) * lscale;

                if (!std::isfinite(psi_here) || psi_here < 0.0)
                    psi_here = 0.0;

                const int im1 = (ci - 1 < cdlo_x) ? cdlo_x : ci - 1;
                const int ip1 = (ci + 1 > cdhi_x) ? cdhi_x : ci + 1;
                const Set::Scalar p_sm =
                    0.25 * (p_eos(im1, cj, ck)
                          + 2.0 * p_eos(ci, cj, ck)
                          + p_eos(ip1, cj, ck));
                const Set::Scalar p_tens_sm = -p_sm;
                if (p_tens_sm <= lp_floor) psi_here = 0.0;

                psi(ci, cj, ck) = psi_here;
            });
        }
        psi_plus_mf[lev]->FillBoundary(geom[lev].periodicity());
    }

    
    void AssertGhostConsistent(Set::Field<Set::Scalar>& f, int lev, const char* tag)
    {
        if (!check_sync) return;
        if (!f.size() || !f[lev]) return;

        const int            nc = f[lev]->nComp();
        const amrex::IntVect ng = f[lev]->nGrowVect();

        amrex::MultiFab cpy(f[lev]->boxArray(), f[lev]->DistributionMap(), nc, ng);
        amrex::MultiFab::Copy(cpy, *f[lev], 0, 0, nc, ng);
        cpy.FillBoundary(geom[lev].periodicity());
        amrex::MultiFab::Subtract(cpy, *f[lev], 0, 0, nc, ng);

        Set::Scalar diff = 0.0;
        for (int c = 0; c < nc; ++c)
            diff = std::max(diff, cpy.norm0(c, f[lev]->nGrow()));

        if (amrex::ParallelDescriptor::IOProcessor())
            Util::Message(INFO, "ghost residual [", tag, "] lev=", lev, " : ", diff);
    }

    
    void FillDHalo(int lev, Set::Scalar a_time)
    {
        d_mf[lev]->FillBoundary(geom[lev].periodicity());
        if (bc_d)
        {
            bc_d->define(geom[lev]);
            bc_d->FillBoundary(*d_mf[lev], 0, 0, a_time, 0);
        }
    }

    void AdvanceDamage(int lev, Set::Scalar a_time, Set::Scalar dt)
    {
        BL_PROFILE("LagrangianShockCPDamage::AdvanceDamage");
        if (dbg_bypass) return;   //

        pressure_mf[lev]->FillBoundary(geom[lev].periodicity());
        FillDHalo(lev, a_time);   // 

        amrex::MultiFab::Copy(*d_old_mf[lev], *d_mf[lev], 0, 0, 1,
                               d_mf[lev]->nGrow());
        // d_old is a full copy of d INCLUDING the halo, which FillDHalo has just
        // made correct, so d_old needs no separate boundary treatment.
        d_old_mf[lev]->FillBoundary(geom[lev].periodicity());

        const Set::Scalar* DX = geom[lev].CellSize();

        amrex::Box cdomain = geom[lev].Domain();
        const auto cdlo   = amrex::lbound(cdomain);
        const auto cdhi   = amrex::ubound(cdomain);
        const int  cdlo_x = cdlo.x;
        const int  cdhi_x = cdhi.x;
        const Set::Scalar lp_floor = dmg.p_tens_floor;

        const Set::Scalar  lM        = dmg.M;
        const Set::Scalar  lkappa    = dmg.kappa;
        const Set::Scalar  lcp       = dmg.cp;
        const DamageModel  lmodel    = dmg.model;
        const Set::Scalar  lGc       = dmg.Gc;
        const Set::Scalar  lell      = dmg.ell;
        const Set::Scalar  lcw       = dmg.cw;
        const Set::Scalar  lMg       = grow.M;
        const Set::Scalar  lT_cap    = grow.T_cap;
        const bool         lpropagate = grow.propagate;

        for (amrex::MFIter mfi(*d_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            auto const& d         = d_mf[lev]->array(mfi);
            auto const& dold      = d_old_mf[lev]->const_array(mfi);
            auto const& p_eos     = pressure_mf[lev]->const_array(mfi);
            auto const& p_drv     = p_drv_mf[lev]->array(mfi);
            auto const& nuc_arr   = nucleated_mf[lev]->const_array(mfi);
            auto const& ep_arr    = edot_p_eq_mf[lev]->const_array(mfi);
            auto const& T_arr     = triaxiality_mf[lev]->const_array(mfi);
            auto const& hist_arr  = hist_mf[lev]->const_array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Scalar p_tens = -p_eos(i,j,k);
                p_drv(i,j,k) = p_tens;

                if (nuc_arr(i,j,k) > 0.9)
                {
                    d(i,j,k) = 1.0;
                    return;
                }

                Set::Scalar dd    = dold(i,j,k);
                Set::Scalar lap_d = Numeric::Laplacian(dold, i, j, k, 0, DX);

                const int im1 = (i - 1 < cdlo_x) ? cdlo_x : i - 1;
                const int ip1 = (i + 1 > cdhi_x) ? cdhi_x : i + 1;
                const Set::Scalar p_tens_sm =
                    -0.25 * (p_eos(im1, j, k)
                           + 2.0 * p_eos(i, j, k)
                           + p_eos(ip1, j, k));
                const bool in_tension = (p_tens_sm > lp_floor);

                Set::Scalar H_hp    = hist_arr(i, j, k);
                Set::Scalar Y_elast = 2.0 * (1.0 - dd) * H_hp;

                Set::Scalar Y_AC;
                Set::Scalar d_dot = 0.0;

                if (lmodel == DamageModel::DoubleWell)
                {

                    Set::Scalar Wprime =
                        2.0 * dd * (1.0 - 2.0 * dd) * (1.0 - dd);
                    Y_AC = -lcp * Wprime + lkappa * lap_d;
                    Set::Scalar Y_tot = Y_elast + Y_AC;

                    if (nuc_arr(i,j,k) > 0.4 && in_tension)
                    {
                        d_dot = lM * std::max(Y_tot, 0.0);
                    }
                }
                else // DamageModel-AT1
                {

                    Y_AC = -(lGc / (lcw * lell))
                         + (2.0 * lGc * lell / lcw) * lap_d;
                    Set::Scalar Y_tot = Y_elast + Y_AC;

                    if (in_tension)
                        d_dot = lM * std::max(Y_tot, 0.0);
                }


                if (lMg > 0.0)
                {
                    bool growth_active = in_tension
                                       && (lpropagate
                                           ? true
                                           : (nuc_arr(i,j,k) > 0.4));
                    if (growth_active)
                    {
                        Set::Scalar edp = ep_arr(i,j,k);
                        Set::Scalar T   = T_arr(i,j,k);
                        if (std::isfinite(edp) && edp > 0.0
                            && std::isfinite(T))
                        {

                            Set::Scalar T_gate = (T > 0.0) ? T : 0.0;
                            if (T_gate > 0.0)
                            {
                                Set::Scalar T_amp =
                                    (T < lT_cap) ? T : lT_cap;
                                Set::Scalar amp =
                                    std::exp(1.5 * T_amp);

                                Set::Scalar d_dot_grow =
                                    lMg * (1.0 - dd) * edp
                                        * amp * T_gate;


                                if (std::isfinite(d_dot_grow)
                                    && d_dot_grow > 0.0)
                                    d_dot += d_dot_grow;
                            }
                        }
                    }
                }

                Set::Scalar dnew = dd + dt * d_dot;
                if (dnew < 0.0) dnew = 0.0;
                if (dnew > 1.0) dnew = 1.0;
                if (!std::isfinite(dnew)) dnew = dd;
                d(i,j,k) = dnew;
            });
        }

        FillDHalo(lev, a_time + dt);
    }


    void FillDNode(int lev, Set::Scalar a_time)
    {
        if (dbg_bypass) return;  
        FillDHalo(lev, a_time);

        amrex::Box domain = geom[lev].Domain();
        domain.convert(amrex::IntVect::TheNodeVector());


        const int ng = std::max(0, d_mf[lev]->nGrow() - 1);

        for (amrex::MFIter mfi(*d_node_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.grownnodaltilebox(-1, ng) & domain;
            auto const& dn = d_node_mf[lev]->array(mfi);
            auto const& dc = d_mf[lev]->const_array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Scalar m = 0.0;
#if AMREX_SPACEDIM == 3
                for (int dk = -1; dk <= 0; ++dk)
#else
                const int dk = 0;
#endif
                for (int dj = -1; dj <= 0; ++dj)
                for (int di = -1; di <= 0; ++di)
                    m = std::max(m, dc(i + di, j + dj, k + dk));

                dn(i, j, k) = m;
            });
        }
        
        this->SyncNodal(d_node_mf, lev);
    }

    void TimeStepComplete(Set::Scalar a_time, int a_step) override
    {
        BL_PROFILE("LagrangianShockCPDamage::TimeStepComplete");
        Base::TimeStepComplete(a_time, a_step);
        for (int lev = 0; lev <= finest_level; ++lev)
            FillDNode(lev, a_time);
    }


    void DegradeModuli(int /*lev*/)
    {
        BL_PROFILE("LagrangianShockCPDamage::DegradeModuli");

    }



    void PopulateDDWDiagnostic(int lev)
    {
        BL_PROFILE("LagrangianShockCPDamage::PopulateDDWDiagnostic");
        if (dbg_bypass) return;   

        amrex::Box domain = geom[lev].Domain();
        domain.convert(amrex::IntVect::TheNodeVector());

       
        for (amrex::MFIter mfi(*ddw_diag_mf[lev], false); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.grownnodaltilebox() & domain;
            auto const& ddw_arr = ddw_diag_mf[lev]->array(mfi);
            auto const& model   = model_mf[lev]->const_array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                auto D = model(i,j,k).DDW(Set::Matrix::Identity());
                ddw_arr(i,j,k,0) = D(0,0,0,0);
                ddw_arr(i,j,k,1) = D(0,0,1,1);
                ddw_arr(i,j,k,2) = D(0,1,0,1);
            });
        }

        this->SyncNodal(ddw_diag_mf, lev);
    }

    
    static AMREX_FORCE_INLINE
    Set::Scalar CellRand(std::uint64_t a_seed, int i, int j, int k, int a_step)
    {
        auto mix = [](std::uint64_t x) -> std::uint64_t
        {
            x += 0x9e3779b97f4a7c15ULL;
            x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
            x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
            return x ^ (x >> 31);
        };
        std::uint64_t h = mix(a_seed);
        h = mix(h ^ (std::uint64_t)(std::uint32_t)i);
        h = mix(h ^ (std::uint64_t)(std::uint32_t)j);
        h = mix(h ^ (std::uint64_t)(std::uint32_t)k);
        h = mix(h ^ (std::uint64_t)(std::uint32_t)a_step);
        // 53-bit mantissa -> [0,1)
        return (Set::Scalar)(h >> 11) * (1.0 / 9007199254740992.0);
    }

    void NucleationCheck(Set::Scalar a_time, int a_step)
    {
        
#ifndef NDEBUG
        {
            int n_lo = n_voids_total, n_hi = n_voids_total;
            amrex::ParallelDescriptor::ReduceIntMin(n_lo);
            amrex::ParallelDescriptor::ReduceIntMax(n_hi);
            Set::Scalar dt_lo = this->timestep, dt_hi = this->timestep;
            amrex::ParallelDescriptor::ReduceRealMin(dt_lo);
            amrex::ParallelDescriptor::ReduceRealMax(dt_hi);
            Util::Assert(INFO, TEST(n_lo == n_hi),
                "n_voids_total differs across ranks; the collectives below will hang");
            Util::Assert(INFO, TEST(dt_lo == dt_hi),
                "timestep differs across ranks; the collectives below will hang");
        }
#endif

        if (n_voids_total >= nuc.max_voids) return;

        const Set::Scalar p_thr     = nuc.p_thresh;
        const Set::Scalar p_ref     = nuc.p_ref;
        const Set::Scalar hm_thr    = nuc.hm_thresh;
        const Set::Scalar base_rate = nuc.base_rate;
        const Set::Scalar void_r    = nuc.void_radius;
        const Set::Scalar excl_r    = nuc.exclusion_radius;
        const Set::Scalar hm_alpha  = nuc.hm_alpha;
        const Set::Scalar hm_beta   = nuc.hm_beta;
        const Set::Scalar sig_spall = nuc.sigma_spall;   // use_stress==4
        const Set::Scalar g1_ref    = nuc.g1_ref;        // use_stress==4,5
        const Set::Scalar m_exp     = nuc.m_exp;         // use_stress==4
        const int         use_stress = nuc.use_stress;

        if (base_rate <= 0.0) return;

        const Set::Scalar dt_eff = this->timestep * (Set::Scalar)nuc.check_freq;
        if (dt_eff <= 0.0) return;

        //  the mt19937 stream is gone; see CellRand above.
        const std::uint64_t lseed = (std::uint64_t)nuc.seed;

        std::vector<std::pair<Set::Scalar, Set::Scalar>> candidates;

        int n_eligible_total  = 0;
        int n_eligible_stress = 0;
        int n_eligible_both   = 0;
        Set::Scalar max_p_tens_seen = -1.0e100;
        Set::Scalar max_P_accept    = 0.0;
        Set::Scalar max_rate_seen   = 0.0;

        int lev = finest_level;
        const Set::Scalar* DX      = geom[lev].CellSize();
        const Set::Scalar* prob_lo = geom[lev].ProbLo();
        const Set::Scalar* prob_hi = geom[lev].ProbHi();
        const Set::Scalar bnd      = nuc.bnd_excl_cells * DX[0];
        const Set::Scalar xlo      = prob_lo[0] + bnd;
        const Set::Scalar xhi      = prob_hi[0] - bnd;
        const Set::Scalar ylo      = prob_lo[1] + bnd;
        const Set::Scalar yhi      = prob_hi[1] - bnd;

        amrex::Box cdomain = geom[lev].Domain();
        const int cdlo_x = amrex::lbound(cdomain).x;
        const int cdhi_x = amrex::ubound(cdomain).x;
        pressure_mf[lev]->FillBoundary(geom[lev].periodicity());

        for (amrex::MFIter mfi(*d_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            auto const& d_arr   = d_mf[lev]->const_array(mfi);
            auto const& nuc_arr = nucleated_mf[lev]->const_array(mfi);
            auto const& p_arr   = pressure_mf[lev]->const_array(mfi);
            auto const& pr_arr  = pr_mf[lev]->array(mfi);

            const auto lo = amrex::lbound(bx);
            const auto hi = amrex::ubound(bx);

            for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                int k = 0;

                pr_arr(i,j,k) = 0.0;

                if (nuc_arr(i,j,k) > 0.5) continue;
                if (d_arr(i,j,k)   > 0.01) continue;

                Set::Scalar x = prob_lo[0] + (i + 0.5) * DX[0];
                Set::Scalar y = prob_lo[1] + (j + 0.5) * DX[1];
                if (x < xlo || x > xhi || y < ylo || y > yhi) continue;

                ++n_eligible_total;
                const int im1 = (i - 1 < cdlo_x) ? cdlo_x : i - 1;
                const int ip1 = (i + 1 > cdhi_x) ? cdhi_x : i + 1;
                Set::Scalar p_tens = -0.25 * (p_arr(im1,j,k)
                                            + 2.0 * p_arr(i,j,k)
                                            + p_arr(ip1,j,k));
                if (p_tens > max_p_tens_seen) max_p_tens_seen = p_tens;

                if (use_stress != 0 && p_tens <= p_thr) continue;
                ++n_eligible_stress;

                Set::Scalar hm = 1.0;
                if (use_heatmap)
                {
                    auto const& hm_arr = heatmap_mf[lev]->const_array(mfi);
                    Set::Scalar hmR = hm_arr(i,j,k,0);
                    Set::Scalar hmG = hm_arr(i,j,k,1);
                    Set::Scalar hmB = hm_arr(i,j,k,2);
                    hm = std::max({hmR, hmG, hmB});
                    if (hm <= hm_thr) continue;
                    hm = hm - hm_thr;
                    if (hm < 0.0) hm = 0.0;
                    if (hm > 1.0) hm = 1.0;
                }
                ++n_eligible_both;

                Set::Scalar rate_cell = 0.0;
                Set::Scalar prefactor = (use_heatmap && hm_beta > 0.0)
                                      ? std::pow(hm, hm_beta)
                                      : 1.0;

                if (use_stress == 0)
                {
                    rate_cell = base_rate * prefactor;
                }
                else if (use_stress == 1)
                {
                    rate_cell = base_rate * prefactor;
                }
                else if (use_stress == 2)
                {
                    Set::Scalar w_p = (p_tens - p_thr) / p_ref;
                    if (w_p <= 0.0) continue;
                    if (w_p > 1.0) w_p = 1.0;
                    rate_cell = base_rate * w_p * prefactor;
                }
                else if (use_stress == 3)
                {
                    Set::Scalar barrier_factor = 1.0 - hm_alpha * hm;
                    if (barrier_factor < 0.0) barrier_factor = 0.0;
                    Set::Scalar exponent = p_ref * barrier_factor / p_tens;
                    rate_cell = base_rate * prefactor * std::exp(-exponent);
                }
                else if (use_stress == 4)
                {

                    Set::Scalar dG1 = g1_ref * (1.0 - hm);
                    if (dG1 < 0.0) dG1 = 0.0;

                    Set::Scalar ratio = sig_spall / p_tens;
                    Set::Scalar dG2 = (m_exp == 1.0)
                                    ? ratio
                                    : std::pow(ratio, m_exp);

                    Set::Scalar exponent = dG1 + dG2;

                    if (exponent > 60.0) exponent = 60.0;

                    rate_cell = base_rate * std::exp(-exponent);
                }
                else if (use_stress == 5)
                {
                    Set::Scalar dG1 = g1_ref * (1.0 - hm);
                    if (dG1 < 0.0)  dG1 = 0.0;
                    if (dG1 > 60.0) dG1 = 60.0;
                    rate_cell = base_rate * std::exp(-dG1);
                }
                else
                {
                    rate_cell = base_rate * prefactor;
                }

                if (!std::isfinite(rate_cell) || rate_cell <= 0.0) continue;
                if (rate_cell > max_rate_seen) max_rate_seen = rate_cell;

                Set::Scalar P_accept = 1.0 - std::exp(-rate_cell * dt_eff);
                if (!std::isfinite(P_accept) || P_accept < 0.0) P_accept = 0.0;
                if (P_accept > 1.0) P_accept = 1.0;
                pr_arr(i,j,k) = P_accept;
                if (P_accept > max_P_accept) max_P_accept = P_accept;
                if (P_accept <= 0.0) continue;


                if (CellRand(lseed, i, j, k, a_step) < P_accept)
                    candidates.emplace_back(x, y);
            }
            }
        }

        {
            int local_n = (int)candidates.size();
            std::vector<int> counts(amrex::ParallelDescriptor::NProcs(), 0);
            counts[amrex::ParallelDescriptor::MyProc()] = local_n;
            amrex::ParallelDescriptor::ReduceIntSum(
                counts.data(), (int)counts.size());

            int total_n = 0;
            for (int c : counts) total_n += c;
            if (total_n == 0)
            {
                if (amrex::ParallelDescriptor::IOProcessor())
                {
                    Util::Message(INFO,
                        "[NUC]   no voids this check "
                        "(n_eligible_total=", n_eligible_total,
                        " pass_stress=", n_eligible_stress,
                        " pass_both=", n_eligible_both,
                        " max_p_tens=", max_p_tens_seen,
                        " max_rate=", max_rate_seen,
                        " max_P_accept=", max_P_accept,
                        " dt_eff=", dt_eff, ")");
                }
                return;
            }

            std::vector<Set::Scalar> xs(total_n, 0.0), ys(total_n, 0.0);
            int off = 0;
            for (int r = 0; r < amrex::ParallelDescriptor::MyProc(); ++r)
                off += counts[r];
            for (int q = 0; q < local_n; ++q) {
                xs[off + q] = candidates[q].first;
                ys[off + q] = candidates[q].second;
            }
            amrex::ParallelDescriptor::ReduceRealSum(xs.data(), total_n);
            amrex::ParallelDescriptor::ReduceRealSum(ys.data(), total_n);

            candidates.clear();
            for (int q = 0; q < total_n; ++q)
                candidates.emplace_back(xs[q], ys[q]);
        }

        std::vector<std::pair<Set::Scalar, Set::Scalar>> accepted;
        for (const auto& c : candidates)
        {
            bool too_close = false;
            for (const auto& a : accepted)
            {
                Set::Scalar dx = c.first  - a.first;
                Set::Scalar dy = c.second - a.second;
                if (dx*dx + dy*dy < excl_r * excl_r) { too_close = true; break; }
            }
            if (!too_close) accepted.push_back(c);
        }

        if (accepted.empty()) return;

        int budget = nuc.max_voids - n_voids_total;
        if (budget <= 0) return;
        if ((int)accepted.size() > budget)
            accepted.resize(budget);

        n_voids_total += (int)accepted.size();

        const Set::Scalar lvoid_r = void_r;
        const Set::Scalar leps    = (nuc.eps > 0.0) ? nuc.eps : 1.0e-12;
        const Set::Scalar lexcl_r = excl_r;

        for (int ilev = 0; ilev <= finest_level; ++ilev)
        {
            const Set::Scalar* DXl = geom[ilev].CellSize();
            const Set::Scalar* lol = geom[ilev].ProbLo();

            for (amrex::MFIter mfi(*d_mf[ilev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.tilebox();
                auto const& d         = d_mf[ilev]->array(mfi);
                auto const& nuc_arr_w = nucleated_mf[ilev]->array(mfi);

                for (const auto& p : accepted)
                {
                    const Set::Scalar cx = p.first;
                    const Set::Scalar cy = p.second;

                    const int ic = (int)std::floor((cx - lol[0]) / DXl[0]);
                    const int jc = (int)std::floor((cy - lol[1]) / DXl[1]);

                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                    {
                        Set::Scalar x = lol[0] + (i + 0.5) * DXl[0];
                        Set::Scalar y = lol[1] + (j + 0.5) * DXl[1];
                        Set::Scalar r = std::sqrt((x - cx)*(x - cx)
                                                + (y - cy)*(y - cy));

                        const bool is_core_cell = (i == ic && j == jc);

                        if (r > lvoid_r + 4.0 * leps && r >= lexcl_r
                            && !is_core_cell)
                            return;

                        Set::Scalar m =
                            0.5 * (1.0 + std::erf((lvoid_r - r) / leps));

                        if (m > 0.5 || is_core_cell)
                        {
                            d(i,j,k)         = 1.0;
                            nuc_arr_w(i,j,k) = 1.0;
                            return;
                        }

                        if (m > 0.0)
                        {
                            Set::Scalar dnew = d(i,j,k) + m * (1.0 - d(i,j,k));
                            if (dnew > 1.0) dnew = 1.0;
                            d(i,j,k) = dnew;
                        }

                        if (r < lexcl_r && nuc_arr_w(i,j,k) < 0.5)
                            nuc_arr_w(i,j,k) = 0.5;
                    });
                }
            }
            FillDHalo(ilev, a_time);   // ### BC FIX
            nucleated_mf[ilev]->FillBoundary(geom[ilev].periodicity());
        }

        if (amrex::ParallelDescriptor::IOProcessor())
            Util::Message(INFO, "LagrangianShockCPDamage: nucleated ",
                           (int)accepted.size(), " voids (total: ",
                           n_voids_total, "/", nuc.max_voids,
                           ", hard d=1 core, smooth-stamp r_v=", lvoid_r,
                           " eps=", leps, ")");
    }

    void UpdateDamageDiagnostics(int lev)
    {
        if (dbg_bypass) return;   // ### BISECT SWITCH
        Set::Scalar local_max_d   = 0.0;
        Set::Scalar local_max_p   = -1.0e100;
        Set::Scalar local_max_psi = 0.0;

        for (amrex::MFIter mfi(*d_mf[lev], false); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            auto const& d_arr   = d_mf[lev]->const_array(mfi);
            auto const& p_arr   = p_drv_mf[lev]->const_array(mfi);
            auto const& psi_arr = psi_plus_mf[lev]->const_array(mfi);

            amrex::Loop(bx, [&](int i, int j, int k)
            {
                if (d_arr(i,j,k)   > local_max_d)   local_max_d   = d_arr(i,j,k);
                if (p_arr(i,j,k)   > local_max_p)   local_max_p   = p_arr(i,j,k);
                if (psi_arr(i,j,k) > local_max_psi) local_max_psi = psi_arr(i,j,k);
            });
        }

        //
        amrex::ParallelDescriptor::ReduceRealMax(local_max_d);
        amrex::ParallelDescriptor::ReduceRealMax(local_max_p);
        amrex::ParallelDescriptor::ReduceRealMax(local_max_psi);

        diag_max_d        = local_max_d;
        diag_max_p_tens   = local_max_p;
        diag_max_psi_plus = local_max_psi;
    }

    Set::Field<Set::Scalar> d_mf;
    Set::Field<Set::Scalar> d_old_mf;
    Set::Field<Set::Scalar> d_node_mf;   
    Set::Field<Set::Scalar> p_drv_mf;
    Set::Field<Set::Scalar> pr_mf;
    Set::Field<Set::Scalar> nucleated_mf;
    Set::Field<Set::Scalar> heatmap_mf;
    
    Set::Field<Set::Scalar> ddw_diag_mf;
    Set::Field<Set::Scalar> psi_plus_mf;
    Set::Field<Set::Scalar> hist_mf;
    Set::Field<Set::Scalar> edot_p_eq_mf;
    Set::Field<Set::Scalar> triaxiality_mf;

    IC::IC<Set::Scalar>* ic_heatmap = nullptr;
    bool use_heatmap = false;

    BC::BC<Set::Scalar>* bc_d              = nullptr;
    BC::BC<Set::Scalar>* bc_nothing_scalar = nullptr;

    
    int dbg_bypass = 0;

    Set::Scalar m_eta_ref_threshold = 1e-2;

    Set::Scalar C11_pristine = 0.0;
    Set::Scalar C12_pristine = 0.0;
    Set::Scalar C44_pristine = 0.0;

    struct {

        DamageModel model           = DamageModel::DoubleWell;
        Set::Scalar M               = 0.0;
        Set::Scalar kappa           = 0.0;

        Set::Scalar cp              = 0.0;

        Set::Scalar barrier_height  = 0.0;

        Set::Scalar Gc              = 0.0;
        Set::Scalar ell             = 0.0;
        Set::Scalar cw              = 8.0 / 3.0;
        Set::Scalar small           = 1e-3;

        Set::Scalar d_void          = 0.95;

        
        Set::Scalar void_res        = 0.02;

       
        Set::Scalar d_seed_floor    = 0.55;

        Set::Scalar hist_seed_factor = 0.0;

        Set::Scalar p_tens_floor    = 0.0;
    } dmg;


    struct {

        Set::Scalar M         = 0.0;

        bool        propagate = false;

        Set::Scalar T_cap     = 3.0;
    } grow;

    struct {
        int         enable            = 0;
        int         use_stress        = 1;
        Set::Scalar base_rate         = 0.0;
        Set::Scalar p_thresh          = 0.0;
        Set::Scalar p_ref             = 1.0e9;
        Set::Scalar hm_alpha          = 0.0;
        Set::Scalar hm_beta           = 0.0;
        Set::Scalar hm_thresh         = 0.0;
        Set::Scalar sigma_spall       = 0.0;   
        Set::Scalar g1_ref            = 0.0;   // dG1 barrier height (modes 4,5)
        Set::Scalar m_exp             = 1.0;   // overstress exponent m   (mode 4)
        Set::Scalar void_radius       = 0.0;
        Set::Scalar eps               = 0.0;
        Set::Scalar exclusion_radius  = 0.0;
        int         check_freq        = 10;
        int         bnd_excl_cells    = 4;
        int         max_voids         = 10;
        unsigned    seed              = 12345;
    } nuc;



    int n_voids_total = 0;

    struct {
        int         enabled = 0;
        Set::Scalar x       = 0.0;
        Set::Scalar y       = 0.0;
        Set::Scalar radius  = 0.0;
    } seed;

    Set::Scalar diag_max_d      = 0.0;
    Set::Scalar diag_max_p_tens = 0.0;
    Set::Scalar diag_max_psi_plus = 0.0;

    int         xt_enable         = 0;
    std::string xt_filename       = "lagshock_xt.dat";
    Set::Scalar xt_dt             = 1.0e-8;
    int         xt_jprobe         = -1;      // -1 => domain centerline
    Set::Scalar next_xt_time      = 0.0;
    std::vector<XTReduce>    xt_modes{XTReduce::Row};
    std::vector<std::string> xt_files;
    std::vector<char>        xt_header_written;

    using Base::pressure_mf;
    using Base::pressure_half_mf;
    using Base::density_mf;
    using Base::density_half_mf;
    using Base::energy_mf;
    using Base::energy_half_mf;
    using Base::soundspeed_mf;
    using Base::soundspeed_half_mf;
    using Base::artvisc2d_mf;
    using Base::artvisc2d_half_mf;
    using Base::artvisc_mf;
    using Base::cell_stress_mf;
    using Base::cell_stress_raw_mf;
    using Base::stress_scale;
    using Base::m_type;
    using Base::finest_level;
    using Base::geom;
    using Base::model_mf;
    using Base::timestep;
    using Base::disp_mf;
    using Base::disp_old_mf;
    using Base::vel_mf;
    using Base::vel_old_mf;
    using Base::stress_mf;
    using Base::strain_mf;
    using Base::rhs_mf;
    using Base::bc;
    using Base::eos;
    using Base::hydro_rho0;
    using Base::av_CQ;
    using Base::av_CL;
    using Base::hydro_cfl;
    using Base::diag_fs_vel;
    using Base::diag_max_pressure;
    using Base::diag_min_pressure;
    using Base::next_diag_time;
    using Base::diag_dt;
    using Base::check_sync;
};

}

#endif