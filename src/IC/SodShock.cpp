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

                if (mf_name == "ic.pvec") {
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
                } else if (mf_name == "ic.pressure") {
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
    mf_name = name.substr(0, name.find(".sodshock"));

    if (type != "sodshock") {
        Util::Abort(INFO, "Unknown type in input: " + type);
    }

    // Parse variables based on mf_name
    if (mf_name == "PVec") {
        pp.query((name + ".left.density").c_str(), rho_left);
        pp.query((name + ".right.density").c_str(), rho_right);
        pp.query((name + ".left.uvel").c_str(), u_left);
        pp.query((name + ".right.uvel").c_str(), u_right);
        pp.query((name + ".left.vvel").c_str(), v_left);
        pp.query((name + ".right.vvel").c_str(), v_right);
        pp.query((name + ".left.wvel").c_str(), w_left);
        pp.query((name + ".right.wvel").c_str(), w_right);
    } else if (mf_name == "Pressure") {
        pp.query((name + ".left").c_str(), p_left);
        pp.query((name + ".right").c_str(), p_right);
    }

    // Parse shock positions
    pp.query("shock.xpos", shock_xpos);
    pp.query("shock.ypos", shock_ypos);
    pp.query("shock.zpos", shock_zpos);
}

} // namespace IC


