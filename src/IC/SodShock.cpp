#include "IC/IC.H"
#include "IC/SodShock.H"
#include "Set/Set.H"
#include "Util/Util.H"
#include "Util/ScimitarX_Util.H"
#include "Model/Fluid/Fluid.H"

namespace IC {

using namespace Model::Fluid;

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
        amrex::Box bx       = mfi.growntilebox();
        amrex::IndexType type = a_phi[lev]->ixType();
        amrex::Array4<Set::Scalar> const& phi = a_phi[lev]->array(mfi);

        // Retrieve variable indices.
        int dens_idx = requires_variable_indices ? variable_indices->DENS : -1;
        int uvel_idx = requires_variable_indices ? variable_indices->UVEL : -1;
#if AMREX_SPACEDIM >= 2
        int vvel_idx = requires_variable_indices ? variable_indices->VVEL : -1;
#endif
#if AMREX_SPACEDIM == 3
        int wvel_idx = requires_variable_indices ? variable_indices->WVEL : -1;
#endif
        int ie_idx = requires_variable_indices ? variable_indices->IE : -1;

        Model::Fluid::Fluid fluid_model;
        Util::ScimitarX_Util::Debug debug;

        const Set::Scalar xmin = geom[lev].ProbLo()[0];
        const Set::Scalar xmax = geom[lev].ProbHi()[0];

        for (int n = 0; n < ncomp; ++n) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Set::Vector x = Set::Position(i, j, k, geom[lev], type);
                Set::Scalar gamma = 1.4;
                
                const std::size_t num_zones = shock_positions.size() + 1;
                
                for (std::size_t idx = 0; idx < num_zones; ++idx) {
                    Set::Scalar lower = (idx == 0) ? xmin : shock_positions[idx - 1];
                    Set::Scalar upper = (idx == shock_positions.size()) ? xmax : shock_positions[idx];
                    
                    if (x(0) >= lower && x(0) < upper) {
                        if (mf_name == "ic.pvec.sodshock") {
                            if (n == dens_idx) {
                                phi(i, j, k, n) = density_zone[idx];
                               
                            }
                            else if (n == uvel_idx)
                                phi(i, j, k, n) = uvel_zone[idx];
#if AMREX_SPACEDIM >= 2
                            else if (n == vvel_idx)
                                phi(i, j, k, n) = vvel_zone[idx];
#endif
#if AMREX_SPACEDIM == 3
                            else if (n == wvel_idx)
                                phi(i, j, k, n) = wvel_zone[idx];
#endif
                            else if (n == ie_idx) {
                                Set::Scalar tmp_density = density_zone[idx];
                                Set::Scalar tmp_pressure = pressure_zone[idx];

                                //printf("phi(%d,%d,%d,%d) = %e at x = %f (zone %zu)\n", i, j, k, n, phi(i, j, k, n), x(0), idx);
                                phi(i, j, k, n) = fluid_model.ComputeInternalEnergyFromDensityAndPressure(
                                    tmp_density, tmp_pressure, gamma);

                                if (phi(i, j, k, n) == 0.0)
                                    debug.DebugComputeInternalEnergyFromDensityAndPressure(
                                        i, j, k, lev, tmp_density, tmp_pressure,
                                        gamma, true, "print", "Zone " + std::to_string(idx));
                            }
                        } else if (mf_name == "ic.pressure.sodshock") {
                            phi(i, j, k, n) = pressure_zone[idx];
                        }
                        break;
                    }
                }
            });
        }
    }
}

void SodShock::initialize(IO::ParmParse& pp, const std::string& name) {
    std::string type = "sodshock";
    pp.query(name.c_str(), type);
    mf_name = name;

    Util::Message(INFO, "DEBUG: Initializing SodShock with mf_name = " + mf_name);

    if (type != "sodshock") {
        Util::Abort(INFO, "Unknown type in input: " + type);
    }
    
    if (mf_name == "ic.pvec.sodshock") {
        pp.query("ic.number_shock_positions", number_of_shocks);
        pp.queryarr("ic.pvec.sodshock.positions", shock_positions);
        pp.queryarr("ic.pvec.sodshock.density", density_zone);
        pp.queryarr("ic.pvec.sodshock.uvel", uvel_zone);
#if AMREX_SPACEDIM >= 2
        pp.queryarr("ic.pvec.sodshock.vvel", vvel_zone);
#endif
#if AMREX_SPACEDIM == 3
        pp.queryarr("ic.pvec.sodshock.wvel", wvel_zone);
#endif
        // Prevent segfault 
        pp.queryarr("ic.pvec.sodshock.pressure", pressure_zone);

        // Ensure all zone arrays have the correct size
        size_t num_zones = shock_positions.size() + 1;
        if (uvel_zone.size() != num_zones ||
#if AMREX_SPACEDIM >= 2
            vvel_zone.size() != num_zones ||
#endif
#if AMREX_SPACEDIM == 3
            wvel_zone.size() != num_zones ||
#endif
            density_zone.size() != num_zones ) {
            Util::Abort(INFO, "Zone arrays have incorrect size for ic.pvec.sodshock");
        }

        Util::Message(INFO, "DEBUG: Parsed ic.pvec.sodshock values:");
        Set::Scalar xmin = geom[0].ProbLo()[0];
        Set::Scalar xmax = geom[0].ProbHi()[0];
        for (size_t i = 0; i < num_zones; ++i) {
            Set::Scalar lower = (i == 0) ? xmin : shock_positions[i - 1];
            Set::Scalar upper = (i == shock_positions.size()) ? xmax : shock_positions[i];
            Util::Message(INFO, "  Zone[" + std::to_string(i) + "]: from " + std::to_string(lower) + " to " + std::to_string(upper));
            Util::Message(INFO, "    Density = " + std::to_string(density_zone[i]));
            Util::Message(INFO, "    u = " + std::to_string(uvel_zone[i]));
#if AMREX_SPACEDIM >= 2
            Util::Message(INFO, "    v = " + std::to_string(vvel_zone[i]));
#endif
#if AMREX_SPACEDIM == 3
            Util::Message(INFO, "    w = " + std::to_string(wvel_zone[i]));
#endif
            
        }
    } else if (mf_name == "ic.pressure.sodshock") {
        pp.query("ic.number_shock_positions", number_of_shocks);
        pp.queryarr("ic.pvec.sodshock.positions", shock_positions);
        pp.queryarr("ic.pressure.sodshock.pressure", pressure_zone);

        // Ensure pressure_zone has the correct size
        size_t num_zones = shock_positions.size() + 1;
        if (pressure_zone.size() != num_zones) {
            Util::Abort(INFO, "Pressure zone array has incorrect size for ic.pressure.sodshock");
        } 


        Util::Message(INFO, "DEBUG: Parsed ic.pressure.sodshock values:");
        Set::Scalar xmin = geom[0].ProbLo()[0];
        Set::Scalar xmax = geom[0].ProbHi()[0];
        for (size_t i = 0; i < num_zones; ++i) {
            Set::Scalar lower = (i == 0) ? xmin : shock_positions[i - 1];
            Set::Scalar upper = (i == shock_positions.size()) ? xmax : shock_positions[i];
            Util::Message(INFO, "  Zone[" + std::to_string(i) + "]: from " + std::to_string(lower) + " to " + std::to_string(upper));
            Util::Message(INFO, "    Pressure = " + std::to_string(pressure_zone[i]));
        }
    }
}

} // end namespace IC