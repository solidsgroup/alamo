#include "IC/IC.H"
#include "IC/Shock.H"
#include "Set/Set.H"
#include "Util/Util.H"
#include "Util/ScimitarX_Util.H"
#include "Numeric/IntegratorVariableAccessLayer.H"
#include "Model/Fluid/Fluid.H"

namespace IC {

using namespace Model::Fluid;

Shock::Shock(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name,
                    const Numeric::GenericVariableAccessor::VariableIndices& precomputed_indices)
    : IC(_geom), variable_indices(&precomputed_indices), requires_variable_indices(true) {
    initialize(pp, name);
}

Shock::Shock(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name)
    : IC(_geom), variable_indices(nullptr), requires_variable_indices(false) {
    initialize(pp, name);
}

void Shock::Add(const int& lev, Set::Field<Set::Scalar>& a_phi, Set::Scalar) {
    int ncomp = a_phi[0]->nComp();

    for (amrex::MFIter mfi(*a_phi[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box bx = mfi.growntilebox();
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

        Set::Scalar gamma = 1.4;

       //Use the pre-computed direction index for efficiency
       int dir_idx = direction_index;
       const std::size_t num_zones = zone_lower_bounds.size();

        for (int n = 0; n < ncomp; ++n) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Set::Vector x = Set::Position(i, j, k, geom[lev], type);

                                
                for (std::size_t idx = 0; idx < num_zones; ++idx) {
                    
                    if (x(dir_idx) >= zone_lower_bounds[idx] && x(dir_idx) < zone_upper_bounds[idx]) {

                        if (mf_name == "ic.shock.pvec") {
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

                                phi(i, j, k, n) = fluid_model.ComputeInternalEnergyFromDensityAndPressure(
                                    tmp_density, tmp_pressure, gamma);

                                if (phi(i, j, k, n) == 0.0)
                                    debug.DebugComputeInternalEnergyFromDensityAndPressure(
                                        i, j, k, lev, tmp_density, tmp_pressure,
                                        gamma, true, "print", "Zone " + std::to_string(idx));
                            }
                        } else if (mf_name == "ic.shock.pressure") {
                            phi(i, j, k, n) = pressure_zone[idx];
                        }
                        break;
                    }
                }
            });
        }
    }
}

void Shock::initialize(IO::ParmParse& pp, const std::string& name) {
    std::string type = "shock";
    pp.query(name.c_str(), type);
    mf_name = name;

    Util::Message(INFO, "Initializing SodShock with mf_name = " + mf_name);

    if (type != "shock") {
        Util::Abort(INFO, "Unknown type in input: " + type);
    }
    
    // Parse common parameters for both PVec and Pressure
    number_of_shocks = 0; // Default to zero
    pp.query("ic.shock.count", number_of_shocks);

    // Handle shock position - may be empty if number_of_shocks = 0
    if (number_of_shocks > 0) {
       pp.queryarr("ic.shock.positions", shock_positions);

     // Validate that we have the currect number of positions
     if (shock_positions.size() != static_cast<size_t>(number_of_shocks)) {
        
        Util::Warning(INFO, "Mismatch between ic.shock.count (" + std::to_string(number_of_shocks) +
                      ") and shock_positions size (" + std::to_string(shock_positions.size()) + ")");
        // Use the actual size from the positions array
        number_of_shocks = static_cast<int>(shock_positions.size());
     }
    }   

    // Get domain bounds for all directions
    std::vector<Set::Scalar> domain_min(AMREX_SPACEDIM);
    std::vector<Set::Scalar> domain_max(AMREX_SPACEDIM);

    for (int d = 0; d < AMREX_SPACEDIM; d++) {
        domain_min[d] = geom[0].ProbLo()[d];
        domain_max[d] = geom[0].ProbHi()[d];
    }

    //Initialize shock directions (default to x-direction)
    shock_direction = "xdir";
    direction_index = 0;

    //Parse direction information if provided
    std::string dir;
    if (pp.query("ic.shock.direction", dir)) {
        
        if (dir == "xdir") { 
            direction_index = 0; 
        }
#if AMREX_SPACEDIM >= 2        
        else if (dir == "ydir") {
            direction_index = 1;
        }
#endif
#if AMREX_SPACEDIM == 3        
        else if (dir == "zdir") {
            direction_index = 2;
        }
#endif
        else {
            Util::Warning(INFO, "Invalid direction '" + dir + "', using xdir");
            shock_direction = "xdir";
            direction_index = 0;            
        }    
    }

    // Validate that the direction is valid for the current dimension
    if (direction_index >= AMREX_SPACEDIM) {
        Util::Warning(INFO, "Direction '" + shock_direction + 
                     "' is invalid for current dimensionality, using xdir");
        shock_direction = "xdir";
        direction_index = 0;
    }

    // Pre-compute zone boundaries based on shock positions and the single direction
    // When number_of_shocks = 0, we have 1 zone covering the entire domain
    const std::size_t num_zones = (number_of_shocks > 0) ? shock_positions.size() + 1 : 1;

    zone_lower_bounds.resize(num_zones);
    zone_upper_bounds.resize(num_zones);

    if (number_of_shocks == 0) {
        // Special case: no shocks means one zone covering the entire domain
        zone_lower_bounds[0] = domain_min[direction_index];
        zone_upper_bounds[0] = domain_max[direction_index];
    } else {
        // Normal case: multiple zones defined by shock positions
        for (std::size_t idx = 0; idx < num_zones; ++idx) {
            zone_lower_bounds[idx] = (idx == 0) ? domain_min[direction_index] : shock_positions[idx-1];
            zone_upper_bounds[idx] = (idx == shock_positions.size()) ? domain_max[direction_index] : shock_positions[idx];
        }
    }

    if (mf_name.find("pvec") != std::string::npos ) {
        pp.queryarr("ic.shock.pvec.density", density_zone);
        pp.queryarr("ic.shock.pvec.uvel", uvel_zone);
#if AMREX_SPACEDIM >= 2
        pp.queryarr("ic.shock.pvec.vvel", vvel_zone);
#endif
#if AMREX_SPACEDIM == 3
        pp.queryarr("ic.shock.pvec.wvel", wvel_zone);
#endif
        // Prevent segfault 
        pp.queryarr("ic.shock.pressure", pressure_zone);

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
            Util::Abort(INFO, "Zone arrays have incorrect size for ic.shock.pvec");
        }

        Util::Message(INFO, "DEBUG: Parsed ic.shock.pvec values:");
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
    } else if (mf_name == "ic.shock.pressure") {
        pp.queryarr("ic.shock.pressure", pressure_zone);

        // Ensure pressure_zone has the correct size
        size_t num_zones = shock_positions.size() + 1;
        if (pressure_zone.size() != num_zones) {
            Util::Abort(INFO, "Pressure zone array has incorrect size for ic.shock.pressure");
        } 


        Util::Message(INFO, "DEBUG: Parsed ic.shock.pressure values:");
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

}  // namespace IC
