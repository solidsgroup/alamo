#include "IC/IC.H"
#include "IC/Shock.H"
#include "Set/Set.H"
#include "Util/Util.H"
#include "Util/ScimitarX_Util.H"
#include "Numeric/IntegratorVariableAccessLayer.H"
#include "Model/Fluid/Fluid.H"
#include "IO/ParmParse.H"
#include "AMReX_Parser.H"

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

// Main Add method that delegates to the appropriate implementation
void Shock::Add(const int& lev, Set::Field<Set::Scalar>& a_phi, Set::Scalar time) {
    if (use_expressions) {
        AddExpression(lev, a_phi, time);
    } else {
        AddConstant(lev, a_phi, time);
    }
}

void Shock::AddConstant(const int& lev, Set::Field<Set::Scalar>& a_phi, Set::Scalar) {
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

// Implementation for expression-based initialization
void Shock::AddExpression(const int& lev, Set::Field<Set::Scalar>& a_phi, Set::Scalar time) {
    int ncomp = a_phi[0]->nComp();

    for (amrex::MFIter mfi(*a_phi[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box bx = mfi.growntilebox();
        amrex::IndexType type = a_phi[lev]->ixType();

        amrex::Array4<Set::Scalar> const& phi = a_phi[lev]->array(mfi);

        // Retrieve variable indices
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
        Set::Scalar gamma = 1.4;

        // Use the pre-computed direction index for efficiency
        int dir_idx = direction_index;
        const std::size_t num_zones = zone_lower_bounds.size();

        for (int n = 0; n < ncomp; ++n) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Set::Vector x = Set::Position(i, j, k, geom[lev], type);
                
                for (std::size_t idx = 0; idx < num_zones; ++idx) {
                    if (x(dir_idx) >= zone_lower_bounds[idx] && x(dir_idx) < zone_upper_bounds[idx]) {
                        if (mf_name == "ic.shock.pvec") {
                            if (n == dens_idx) {
                                #if AMREX_SPACEDIM == 1
                                phi(i, j, k, n) = density_func[idx](x(0), 0.0, 0.0, time);
                                #elif AMREX_SPACEDIM == 2
                                phi(i, j, k, n) = density_func[idx](x(0), x(1), 0.0, time);
                                #elif AMREX_SPACEDIM == 3
                                phi(i, j, k, n) = density_func[idx](x(0), x(1), x(2), time);
                                #endif
                            }
                            else if (n == uvel_idx) {
                                #if AMREX_SPACEDIM == 1
                                phi(i, j, k, n) = uvel_func[idx](x(0), 0.0, 0.0, time);
                                #elif AMREX_SPACEDIM == 2
                                phi(i, j, k, n) = uvel_func[idx](x(0), x(1), 0.0, time);
                                #elif AMREX_SPACEDIM == 3
                                phi(i, j, k, n) = uvel_func[idx](x(0), x(1), x(2), time);
                                #endif
                            }
#if AMREX_SPACEDIM >= 2
                            else if (n == vvel_idx) {
                                #if AMREX_SPACEDIM == 2
                                phi(i, j, k, n) = vvel_func[idx](x(0), x(1), 0.0, time);
                                #elif AMREX_SPACEDIM == 3
                                phi(i, j, k, n) = vvel_func[idx](x(0), x(1), x(2), time);
                                #endif
                            }
#endif
#if AMREX_SPACEDIM == 3
                            else if (n == wvel_idx) {
                                phi(i, j, k, n) = wvel_func[idx](x(0), x(1), x(2), time);
                            }
#endif
                            else if (n == ie_idx) {
                                // For internal energy, compute from pressure and density
                                #if AMREX_SPACEDIM == 1
                                Set::Scalar tmp_density = density_func[idx](x(0), 0.0, 0.0, time);
                                Set::Scalar tmp_pressure = pressure_func[idx](x(0), 0.0, 0.0, time);
                                #elif AMREX_SPACEDIM == 2
                                Set::Scalar tmp_density = density_func[idx](x(0), x(1), 0.0, time);
                                Set::Scalar tmp_pressure = pressure_func[idx](x(0), x(1), 0.0, time);
                                #elif AMREX_SPACEDIM == 3
                                Set::Scalar tmp_density = density_func[idx](x(0), x(1), x(2), time);
                                Set::Scalar tmp_pressure = pressure_func[idx](x(0), x(1), x(2), time);
                                #endif

                                phi(i, j, k, n) = fluid_model.ComputeInternalEnergyFromDensityAndPressure(
                                    tmp_density, tmp_pressure, gamma);
                            }
                        } else if (mf_name == "ic.shock.pressure") {
                            #if AMREX_SPACEDIM == 1
                            phi(i, j, k, n) = pressure_func[idx](x(0), 0.0, 0.0, time);
                            #elif AMREX_SPACEDIM == 2
                            phi(i, j, k, n) = pressure_func[idx](x(0), x(1), 0.0, time);
                            #elif AMREX_SPACEDIM == 3
                            phi(i, j, k, n) = pressure_func[idx](x(0), x(1), x(2), time);
                            #endif
                        }
                        break;
                    }
                }
            });
        }
    }
}

// Helper method to compile a parser from an expression string
void Shock::SetupParser(std::vector<amrex::Parser>& parser,
                        std::vector<amrex::ParserExecutor<4>>& func,
                        const std::vector<std::string>& expr_list,
                        const std::set<std::pair<std::string, Set::Scalar>>& constants) {
    
    parser.clear();
    func.clear();
    
    for (const auto& expr : expr_list) {
        parser.push_back(amrex::Parser(expr));
        
        // Add constants
        for (const auto& [name, value] : constants) {
            parser.back().setConstant(name, value);
        }
        
        parser.back().registerVariables({"x", "y", "z", "t"});
        func.push_back(parser.back().compile<4>());
    }
}

void Shock::initialize(IO::ParmParse& pp, const std::string& name) {
    std::string type = "shock";
    mf_name = name;

    Util::Message(INFO, "Initializing Shock IC with mf_name = " + mf_name);
    
    // Parse common parameters for both PVec and Pressure
    number_of_shocks = 0; // Default to zero
    pp.query("ic.shock.count", number_of_shocks); // Number of Shocks

    // Handle shock position - may be empty if number_of_shocks = 0
    if (number_of_shocks > 0) {
        pp.queryarr("ic.shock.positions", shock_positions); // Positions of Shocks

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

    // Check if expressions should be used
    pp.query("ic.shock.use_expressions", use_expressions); // use_espressions bool

    // Collect parser constants upfront
    std::set<std::pair<std::string, Set::Scalar>> constants;
    if (use_expressions) {

        std::string prefix = "ic.shock.expressions.constant";
        std::set<std::string> entries = pp.getEntries(prefix);

        for (const auto& entry : entries) {
            std::string fullname = entry;
            Set::Scalar val = NAN;
            pp.query(fullname.c_str(), val);
            
            std::size_t lastDot = fullname.find_last_of('.');
            std::string constant_name = fullname.substr(lastDot + 1);
            
            constants.insert({constant_name, val});
            Util::Message(INFO, "Added constant: " + constant_name + " = " + std::to_string(val));
        }
    }

    if (mf_name.find("pvec") != std::string::npos) {
        if (use_expressions) {
            // Read expression arrays directly like constant arrays
            pp.queryarr("ic.shock.expressions.pvec.density", density_expr); // density expressions for all zones
            pp.queryarr("ic.shock.expressions.pvec.uvel", uvel_expr); // U-Velocity expressions for all zones
#if AMREX_SPACEDIM >= 2
            pp.queryarr("ic.shock.expressions.pvec.vvel", vvel_expr); // V-velocity expressions for all zones
#endif
#if AMREX_SPACEDIM == 3
            pp.queryarr("ic.shock.expressions.pvec.wvel", wvel_expr); // W-velocity expressions for all zones
#endif
            pp.queryarr("ic.shock.expressions.pressure", pressure_expr); // pressure expressions for all zones

            // Validate array sizes
            if (density_expr.size() != num_zones ||
                uvel_expr.size() != num_zones ||
#if AMREX_SPACEDIM >= 2
                vvel_expr.size() != num_zones ||
#endif
#if AMREX_SPACEDIM == 3
                wvel_expr.size() != num_zones ||
#endif
                pressure_expr.size() != num_zones) {
                Util::Abort(INFO, "Expression arrays have incorrect size for ic.shock.pvec expressions");
            }

            // Set up parsers and executors
            SetupParser(density_parser, density_func, density_expr, constants);
            SetupParser(uvel_parser, uvel_func, uvel_expr, constants);
#if AMREX_SPACEDIM >= 2
            SetupParser(vvel_parser, vvel_func, vvel_expr, constants);
#endif
#if AMREX_SPACEDIM == 3
            SetupParser(wvel_parser, wvel_func, wvel_expr, constants);
#endif
            SetupParser(pressure_parser, pressure_func, pressure_expr, constants);

            // Log the expressions for debugging
            Util::Message(INFO, "DEBUG: Parsed ic.shock.pvec expressions:");
            for (size_t i = 0; i < num_zones; ++i) {
                Set::Scalar lower = (i == 0) ? domain_min[direction_index] : shock_positions[i - 1];
                Set::Scalar upper = (i == shock_positions.size()) ? domain_max[direction_index] : shock_positions[i];
                Util::Message(INFO, "Zone[" + std::to_string(i) + "]: from " + std::to_string(lower) + " to " + std::to_string(upper));
                Util::Message(INFO, "  Density expression = " + density_expr[i]);
                Util::Message(INFO, "  u expression = " + uvel_expr[i]);
#if AMREX_SPACEDIM >= 2
                Util::Message(INFO, "  v expression = " + vvel_expr[i]);
#endif
#if AMREX_SPACEDIM == 3
                Util::Message(INFO, "  w expression = " + wvel_expr[i]);
#endif
                Util::Message(INFO, "  Pressure expression = " + pressure_expr[i]);
            }
        } else {
            // Original implementation - read arrays of constant values
            pp.queryarr("ic.shock.pvec.density", density_zone);  // density values for all zones
            pp.queryarr("ic.shock.pvec.uvel", uvel_zone);  // U-velocity values for all zones
#if AMREX_SPACEDIM >= 2
            pp.queryarr("ic.shock.pvec.vvel", vvel_zone); // V-velocity values for all zones
#endif
#if AMREX_SPACEDIM == 3
            pp.queryarr("ic.shock.pvec.wvel", wvel_zone); // W-velocity values for all zones
#endif
            pp.queryarr("ic.shock.pressure", pressure_zone); // Pressure values for all zones

            // Ensure all zone arrays have the correct size
            if (density_zone.size() != num_zones ||
                uvel_zone.size() != num_zones ||
#if AMREX_SPACEDIM >= 2
                vvel_zone.size() != num_zones ||
#endif
#if AMREX_SPACEDIM == 3
                wvel_zone.size() != num_zones ||
#endif
                pressure_zone.size() != num_zones) {
                Util::Abort(INFO, "Zone arrays have incorrect size for ic.shock.pvec");
            }

            Util::Message(INFO, "DEBUG: Parsed ic.shock.pvec values:");
            for (size_t i = 0; i < num_zones; ++i) {
                Set::Scalar lower = (i == 0) ? domain_min[direction_index] : shock_positions[i - 1];
                Set::Scalar upper = (i == shock_positions.size()) ? domain_max[direction_index] : shock_positions[i];
                Util::Message(INFO, "  Zone[" + std::to_string(i) + "]: from " + std::to_string(lower) + " to " + std::to_string(upper));
                Util::Message(INFO, "    Density = " + std::to_string(density_zone[i]));
                Util::Message(INFO, "    u = " + std::to_string(uvel_zone[i]));
#if AMREX_SPACEDIM >= 2
                Util::Message(INFO, "    v = " + std::to_string(vvel_zone[i]));
#endif
#if AMREX_SPACEDIM == 3
                Util::Message(INFO, "    w = " + std::to_string(wvel_zone[i]));
#endif
                Util::Message(INFO, "    Pressure = " + std::to_string(pressure_zone[i]));
            }
        }
    } else if (mf_name == "ic.shock.pressure") {
        if (use_expressions) {
            // Read expression array directly like constant array
            pp.queryarr("ic.shock.expressions.pressure", pressure_expr);

            // Validate array size
            if (pressure_expr.size() != num_zones) {
                Util::Abort(INFO, "Pressure expression array has incorrect size for ic.shock.pressure");
            }

            // Set up parser and executor
            SetupParser(pressure_parser, pressure_func, pressure_expr, constants);

            // Log the expressions for debugging
            Util::Message(INFO, "DEBUG: Parsed ic.shock.pressure expressions:");
            for (size_t i = 0; i < num_zones; ++i) {
                Set::Scalar lower = (i == 0) ? domain_min[direction_index] : shock_positions[i - 1];
                Set::Scalar upper = (i == shock_positions.size()) ? domain_max[direction_index] : shock_positions[i];
                Util::Message(INFO, "  Zone[" + std::to_string(i) + "]: from " + std::to_string(lower) + " to " + std::to_string(upper));
                Util::Message(INFO, "    Pressure expression = " + pressure_expr[i]);
            }
        } else {
            // Original implementation - read array of constant values
            pp.queryarr("ic.shock.pressure", pressure_zone);

            // Ensure pressure_zone has the correct size
            if (pressure_zone.size() != num_zones) {
                Util::Abort(INFO, "Pressure zone array has incorrect size for ic.shock.pressure");
            } 

            Util::Message(INFO, "DEBUG: Parsed ic.shock.pressure values:");
            for (size_t i = 0; i < num_zones; ++i) {
                Set::Scalar lower = (i == 0) ? domain_min[direction_index] : shock_positions[i - 1];
                Set::Scalar upper = (i == shock_positions.size()) ? domain_max[direction_index] : shock_positions[i];
                Util::Message(INFO, "  Zone[" + std::to_string(i) + "]: from " + std::to_string(lower) + " to " + std::to_string(upper));
                Util::Message(INFO, "    Pressure = " + std::to_string(pressure_zone[i]));
            }
        }
    }

}

}  // namespace IC
