#include "IC/IC.H"
#include "IC/SodShock.H"
#include "Set/Set.H"
#include "Util/Util.H"
#include "Util/ScimitarX_Util.H"

namespace IC {

SodShock::SodShock(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name,
                   const Util::ScimitarX_Util::getVariableIndex& precomputed_indices)
    : IC(_geom), variable_indices(&precomputed_indices), requires_variable_indices(true) {
    initialize(pp, name);
}

SodShock::SodShock(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name)
    : IC(_geom), variable_indices(nullptr), requires_variable_indices(false) {
    initialize(pp, name);
}

void SodShock::Add(const int& lev, Set::Field<Set::Scalar>& a_phi, Set::Scalar) {
    int ncomp = a_phi[0]->nComp();

    for (amrex::MFIter mfi(*a_phi[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box bx = mfi.growntilebox();
        amrex::IndexType type = a_phi[lev]->ixType();

        amrex::Array4<Set::Scalar> const& phi = a_phi[lev]->array(mfi);

        int dens_idx = requires_variable_indices ? variable_indices->DENS : -1;
        int uvel_idx = requires_variable_indices ? variable_indices->UVEL : -1;
#if AMREX_SPACEDIM >= 2
        int vvel_idx = requires_variable_indices ? variable_indices->VVEL : -1;
#endif
#if AMREX_SPACEDIM == 3
        int wvel_idx = requires_variable_indices ? variable_indices->WVEL : -1;
#endif
        int ie_idx = requires_variable_indices ? variable_indices->IE : -1;

        for (int n = 0; n < ncomp; ++n) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Set::Vector x = Set::Position(i, j, k, geom[lev], type);

                if (mf_name == "ic.pvec.sodshock") {
                    if (x(0) < shock_xpos) {
                        if (n == dens_idx) phi(i, j, k, n) = rho_left;
                        else if (n == uvel_idx) phi(i, j, k, n) = u_left;
#if AMREX_SPACEDIM >= 2
                        else if (n == vvel_idx) phi(i, j, k, n) = v_left;
#endif
#if AMREX_SPACEDIM == 3
                        else if (n == wvel_idx) phi(i, j, k, n) = w_left;
#endif
                        else if (n == ie_idx) phi(i, j, k, n) = 0.0;
                    } else {
                        if (n == dens_idx) phi(i, j, k, n) = rho_right;
                        else if (n == uvel_idx) phi(i, j, k, n) = u_right;
#if AMREX_SPACEDIM >= 2
                        else if (n == vvel_idx) phi(i, j, k, n) = v_right;
#endif
#if AMREX_SPACEDIM == 3
                        else if (n == wvel_idx) phi(i, j, k, n) = w_right;
#endif
                        else if (n == ie_idx) phi(i, j, k, n) = 0.0;
                    }
                } else if (mf_name == "ic.pressure.sodshock") {
                    phi(i, j, k, n) = (x(0) < shock_xpos) ? p_left : p_right;
                } else {
                    Util::Abort(INFO, "Unknown MultiFab name: " + mf_name);
                }
            });
        }
    }
}

void SodShock::initialize(IO::ParmParse& pp, const std::string& name) {
    // Extract type from the passed name
    std::string type = "sodshock";
    pp.query(name.c_str(), type);
    mf_name = name; // Keep full name for accurate matching

    Util::Message(INFO, "DEBUG: Initializing SodShock with mf_name = " + mf_name);

    if (type != "sodshock") {
        Util::Abort(INFO, "Unknown type in input: " + type);
    }

    // Parse variables based on mf_name
    if (mf_name == "ic.pvec.sodshock") {
        pp.query((name + ".left.DENS").c_str(), rho_left);
        pp.query((name + ".right.DENS").c_str(), rho_right);
        pp.query((name + ".left.UVEL").c_str(), u_left);
        pp.query((name + ".right.UVEL").c_str(), u_right);
#if AMREX_SPACEDIM >= 2
        pp.query((name + ".left.VVEL").c_str(), v_left);
        pp.query((name + ".right.VVEL").c_str(), v_right);
#endif
#if AMREX_SPACEDIM == 3
        pp.query((name + ".left.WVEL").c_str(), w_left);
        pp.query((name + ".right.WVEL").c_str(), w_right);
#endif

        // Debug parsed values
        Util::Message(INFO, "DEBUG: Parsed ic.pvec.sodshock values:");
        Util::Message(INFO, "  rho_left = " + std::to_string(rho_left));
        Util::Message(INFO, "  rho_right = " + std::to_string(rho_right));
        Util::Message(INFO, "  u_left = " + std::to_string(u_left) + ", u_right = " + std::to_string(u_right));
#if AMREX_SPACEDIM >= 2
        Util::Message(INFO, "  v_left = " + std::to_string(v_left) + ", v_right = " + std::to_string(v_right));
#endif
#if AMREX_SPACEDIM == 3
        Util::Message(INFO, "  w_left = " + std::to_string(w_left) + ", w_right = " + std::to_string(w_right));
#endif
    } else if (mf_name == "ic.pressure.sodshock") {
        pp.query((name + ".left").c_str(), p_left);
        pp.query((name + ".right").c_str(), p_right);

        // Debug parsed values
        Util::Message(INFO, "DEBUG: Parsed ic.pressure.sodshock values:");
        Util::Message(INFO, "  p_left = " + std::to_string(p_left));
        Util::Message(INFO, "  p_right = " + std::to_string(p_right));
    }

    // Parse shock positions
    pp.query("shock.xpos", shock_xpos);
    pp.query("shock.ypos", shock_ypos);
    pp.query("shock.zpos", shock_zpos);

    // Debug shock positions
    Util::Message(INFO, "DEBUG: Shock position: xpos = " + std::to_string(shock_xpos) +
                         ", ypos = " + std::to_string(shock_ypos) +
                         ", zpos = " + std::to_string(shock_zpos));
}

} // namespace IC


