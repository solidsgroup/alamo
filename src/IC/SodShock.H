 #ifndef IC_SODSHOCK_H_
 #define IC_SODSHOCK_H_

#include "Set/Set.H"
#include "IC/IC.H"
using namespace std;
#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include "mpi.h"
#include "IO/ParmParse.H"
#include "Util/ScimitarX_SetIndex.H"

namespace IC {
class SodShock : public IC {
public:
    // Constructor for qvec requiring variable indices
    SodShock(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name,
             const Util::ScimitarX_SetIndex::VariableIndicesResult& precomputed_indices)
        : IC(_geom), variable_indices(&precomputed_indices), requires_variable_indices(true) {
        initialize(pp, name);
    }

    // Constructor for pressure (no variable indices required)
    SodShock(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name)
        : IC(_geom), variable_indices(nullptr), requires_variable_indices(false) {
        initialize(pp, name);
    }

    void Define() {}

    void Add(const int& lev, Set::Field<Set::Scalar>& a_phi, Set::Scalar) {
        
        int ncomp = a_phi[0]->nComp();   

        for (amrex::MFIter mfi(*a_phi[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            amrex::Box bx = mfi.growntilebox();
            amrex::IndexType type = a_phi[lev]->ixType();

            amrex::Array4<Set::Scalar> const& phi = a_phi[lev]->array(mfi);

            int dens_idx = requires_variable_indices ? variable_indices->DENS : -1;
            int uvel_idx = requires_variable_indices ? variable_indices->UVEL : -1;
#if AMREX_SPACEDIM >=2
            int vvel_idx = requires_variable_indices ? variable_indices->VVEL : -1;
#endif
#if AMREX_SPACEDIM ==3
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
#if AMREX_SPACEDIM >=2
                    else if (n == vvel_idx) phi(i, j, k, n) = v_left;
#endif
#if AMREX_SPACEDIM ==3
                    else if (n == wvel_idx) phi(i, j, k, n) = w_left;
#endif
                    else if (n == ie_idx) phi(i, j, k, n) = 0.0;
                } else {
                    if (n == dens_idx) phi(i, j, k, n) = rho_right;
                    else if (n == uvel_idx) phi(i, j, k, n) = u_right;
#if AMREX_SPACEDIM >=2
                    else if (n == vvel_idx) phi(i, j, k, n) = v_right;
#endif
#if AMREX_SPACEDIM ==3
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

private:

    void initialize(IO::ParmParse& pp, const std::string& name) {
        // Extract type from the passed name
        std::string type = "sodshock";
        pp.query(name.c_str(), type); // Use c_str() for compatibility with ParmParse::query
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

    // MultiFab type name
    std::string mf_name;

    // Variable indices (only required for qvec)
    const Util::ScimitarX_SetIndex::VariableIndicesResult* variable_indices;
    bool requires_variable_indices;

    // QVec initialization variables
    Set::Scalar rho_left = 1.0;  // Left density
    Set::Scalar rho_right = 0.125; // Right density
    Set::Scalar u_left = 0.0;     // Left velocity (x-direction)
    Set::Scalar u_right = 0.0;    // Right velocity (x-direction)
    Set::Scalar v_left = 0.0;     // Left velocity (y-direction)
    Set::Scalar v_right = 0.0;    // Right velocity (y-direction)
    Set::Scalar w_left = 0.0;     // Left velocity (z-direction)
    Set::Scalar w_right = 0.0;    // Right velocity (z-direction)

    // Pressure initialization variables
    Set::Scalar p_left = 1.0;     // Left pressure
    Set::Scalar p_right = 0.1;    // Right pressure

    // Shock position
    Set::Scalar shock_xpos = 0.5; // Shock position in x-direction
    Set::Scalar shock_ypos = 0.0; // Shock position in y-direction
    Set::Scalar shock_zpos = 0.0; // Shock position in z-direction
};
}


#endif
