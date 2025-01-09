#ifndef UTIL_SCIMITARX_UTIL_H
#define UTIL_SCIMITARX_UTIL_H

#include <string>
#include <iostream>
#include <vector>
#include <array>
#include "Util/Util.H"
#include "Set/Set.H"

namespace Util
{
class ScimitarX_Util
{
public:
    ScimitarX_Util () {};

    struct Constants {
        static constexpr double zeroval = 0.0;
        static constexpr double one = 1.0;
        static constexpr double two = 2.0;
        static constexpr double three = 3.0;
        static constexpr double four = 4.0;
        static constexpr double five = 5.0;
        static constexpr double twothird = 0.6666666666666666667;
        static constexpr double half = 0.5;
        static constexpr double third = twothird * half;
        static constexpr double sixth = twothird * half * half;
        static constexpr double eighth = 0.125;
        static constexpr double sixteenth = 0.0625;
        static constexpr double twelveth = 0.08333333333;
        static constexpr double fourth = 0.25;
        static constexpr double tenth = 0.1;
        static constexpr double threetenth = 0.3;
        static constexpr double threefourth = 0.75;
        static constexpr double PI = 3.14159265358979323846264338327950288;
        static constexpr double TINYVAL = 1e-10;
        static constexpr double HUGEVAL = 1e+10;
    };

    // NOTE: When modifying solver types, ensure the corresponding variables are added here to enable dynamic indexing.
    struct getVariableIndex {
        std::map<int, int> variableIndexMap; // Mapping from variable enums to indices
        int NVAR_MAX = 0;                    // Total number of variables
        int DENS = -1;                       // Index for Density
        int UVEL = -1;                       // Index for U-velocity
        int VVEL = -1;                       // Index for V-velocity
        int WVEL = -1;                       // Index for W-velocity
        int IE = -1;                         // Index for Internal Energy
    };

    // Nested Debug class containing debug APIs
    class Debug {
    
    private:
        static amrex::IntVect target_location;  // Store target cell index (i, j, k)
    
    public:
        // 1. Function to set the target debug location
        static void SetTargetDebugLocation(double x, double y, double z, const amrex::Geometry& geom) {
        // Create a RealVect using AMREX_D_DECL to handle different dimensions
        amrex::Array<amrex::Real, AMREX_SPACEDIM> pos = {AMREX_D_DECL(x, y, z)};
    
        // Convert to cell index (i, j, k)
        target_location = geom.CellIndex(pos.data());  // Use .getVect() to pass the internal pointer

        Util::Message(INFO, "Target debug location set to: (" +
                           std::to_string(target_location[0]) + ", " +
                           std::to_string(target_location[1]) + ", " +
                           (AMREX_SPACEDIM == 3 ? std::to_string(target_location[2]) : "N/A") + ")");
        }

        static void SetTargetDebugLocationIndices(int i, int j, int k) {
        // Set the target location directly using grid indices (i, j, k)
        target_location = amrex::IntVect(AMREX_D_DECL(i, j, k));  // Create IntVect with (i, j, k)

        Util::Message(INFO, "Target debug location set using indices: (" +
                           std::to_string(target_location[0]) + ", " +
                           std::to_string(target_location[1]) + ", " +
                           (AMREX_SPACEDIM == 3 ? std::to_string(target_location[2]) : "N/A") + ")");
        }

        static bool IsTargetLocation(int i, int j, int k) {
            if (AMREX_SPACEDIM >= 1 && target_location[0] != i) return false;
            if (AMREX_SPACEDIM >= 2 && target_location[1] != j) return false;
            if (AMREX_SPACEDIM == 3 && target_location[2] != k) return false;
        return true;
        }

        // Inline debug function that checks the current cell's value

        static void DebugValuesIfTarget(int i, int j, int k, const Set::Scalar& value, const std::string& field_name,
                                const std::string& context, bool abort_if_nan, int component = -1) {
    
            if (IsTargetLocation(i, j, k)) {
            // Build debug message
            std::string debug_message = "Debug at (" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) +
                                    ") in context '" + context + "' for field '" + field_name + "'";
        
            // Add component information if applicable
            if (component >= 0) {
            debug_message += " (Component: " + std::to_string(component) + ")";
            }

            debug_message += " = " + std::to_string(value);

            // Print debug message
            Util::Message(INFO, debug_message);

                // Check for NaN or invalid value
                if (std::isnan(value) || std::isinf(value)) {
                    std::string error_message = "Invalid value detected at (" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) +
                                        ") for field '" + field_name + "' in context '" + context + "'. Value = " + std::to_string(value);
                if (component >= 0) {
                    error_message += " (Component: " + std::to_string(component) + ")";
                }

                        if (abort_if_nan) {
                            Util::Abort(INFO, error_message);
                        } else{
                            Util::Warning(INFO, error_message);
                        }       
                    }
                }
        }           

    static void CheckNaNAndAbort(int i, int j, int k, const Set::Scalar& value, const std::string& field_name, const std::string& context, bool abort_if_nan = true) {
        // Check if the value is NaN or Inf
        if (std::isnan(value) || std::isinf(value)) {
        std::string error_message = "Invalid value detected at (" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) +
                                    ") for field '" + field_name + "' in context '" + context + "'. Value = " + std::to_string(value);

        if (abort_if_nan) {
            Util::Abort(INFO, error_message);  // Abort with the error message
        } else {
            Util::Warning(INFO, error_message);  // Print a warning instead of aborting
        }
        }   
    }

        // API to debug internal energy computation
        static void DebugComputeInternalEnergyFromDensityAndPressure(int i, int j, int k, int lev,
                                                                     Set::Scalar pressure, Set::Scalar density, Set::Scalar gamma,
                                                                     bool abort_on_checkpoint, const std::string& mode = "check",
                                                                     const std::string& context = "Unknown");
        static void DebugEulerForwardStep(int i, int j, int k,
                                  Set::Scalar xflux_right, Set::Scalar xflux_left,
                                  Set::Scalar yflux_top, Set::Scalar yflux_bottom,
                                  Set::Scalar zflux_front, Set::Scalar zflux_back,
                                  Set::Scalar q_old, Set::Scalar dt,
                                  const std::string& context, bool abort_if_nan = true) {

            if (IsTargetLocation(i, j, k)) {
                // Compute the flux divergence
                Set::Scalar flux_div_x = (xflux_right - xflux_left);
                Set::Scalar flux_div_y = (yflux_top - yflux_bottom);
                Set::Scalar flux_div_z = (AMREX_SPACEDIM == 3) ? (zflux_front - zflux_back) : 0.0;
            
                Set::Scalar total_flux_div = flux_div_x + flux_div_y + flux_div_z;
            
                // Compute the new value
                Set::Scalar new_value = q_old - dt * total_flux_div;
            
                // Print the debug message
                Util::Message(INFO, "EulerForward Debug - Context: '" + context + "'");
                Util::Message(INFO, "Index: (" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + ")");
                Util::Message(INFO, "  xflux_right = " + std::to_string(xflux_right) + ", xflux_left = " + std::to_string(xflux_left));
                Util::Message(INFO, "  yflux_top = " + std::to_string(yflux_top) + ", yflux_bottom = " + std::to_string(yflux_bottom));
                if (AMREX_SPACEDIM == 3) {
                    Util::Message(INFO, "  zflux_front = " + std::to_string(zflux_front) + ", zflux_back = " + std::to_string(zflux_back));
                }
                Util::Message(INFO, "  q_old = " + std::to_string(q_old));
                Util::Message(INFO, "  dt = " + std::to_string(dt));
                Util::Message(INFO, "  Total Flux Divergence: " + std::to_string(total_flux_div));
                Util::Message(INFO, "  New Value: " + std::to_string(new_value));
            
                // Check for NaN or Inf and abort if required
                if (std::isnan(new_value) || std::isinf(new_value)) {
                    std::string error_message = "Invalid value detected during EulerForward at (" +
                                                std::to_string(i) + ", " + std::to_string(j) + ", " +
                                                std::to_string(k) + "). Context: '" + context +
                                                "', New Value = " + std::to_string(new_value);
            
                        if (abort_if_nan) {
                            Util::Abort(INFO, error_message);  // Abort if NaN/Inf
                        } else {
                            Util::Warning(INFO, error_message);  // Print warning instead of aborting
                        }
                    }
                }  
            }
            
            static void DebugLocalLaxFriedrichsFlux(int i, int j, int k,
                                                    Set::Scalar rho_L, Set::Scalar u_L, Set::Scalar v_L, Set::Scalar w_L, Set::Scalar e_L, Set::Scalar pressure_L,
                                                    Set::Scalar rho_R, Set::Scalar u_R, Set::Scalar v_R, Set::Scalar w_R, Set::Scalar e_R, Set::Scalar pressure_R,
                                                    Set::Scalar c_L, Set::Scalar c_R, Set::Scalar lambda_max,
                                                    Set::Scalar flux_L[], Set::Scalar flux_R[],
                                                    const std::string& context, bool abort_if_nan = true) {

            if (IsTargetLocation(i, j, k)) {
                // Print debug information for left and right states
                Util::Message(INFO, "Local Lax-Friedrichs Flux Debug - Context: '" + context + "'");
                Util::Message(INFO, "Index: (" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + ")");
                Util::Message(INFO, "  Left State (QL):");
                Util::Message(INFO, "    rho_L = " + std::to_string(rho_L) + ", u_L = " + std::to_string(u_L) + ", v_L = " + std::to_string(v_L) +
                                         ", w_L = " + std::to_string(w_L) + ", e_L = " + std::to_string(e_L) +
                                         ", pressure_L = " + std::to_string(pressure_L));
                Util::Message(INFO, "  Right State (QR):");
                Util::Message(INFO, "    rho_R = " + std::to_string(rho_R) + ", u_R = " + std::to_string(u_R) + ", v_R = " + std::to_string(v_R) +
                                         ", w_R = " + std::to_string(w_R) + ", e_R = " + std::to_string(e_R) +
                                         ", pressure_R = " + std::to_string(pressure_R));
                Util::Message(INFO, "  Wave Speeds:");
                Util::Message(INFO, "    c_L = " + std::to_string(c_L) + ", c_R = " + std::to_string(c_R) + ", lambda_max = " + std::to_string(lambda_max));
            
                // Print fluxes
                Util::Message(INFO, "  Flux Components (F_L, F_R):");
                for (int n = 0; n < 5; ++n) {
                    Util::Message(INFO, "    Flux component " + std::to_string(n) +
                                           ": F_L[" + std::to_string(n) + "] = " + std::to_string(flux_L[n]) +
                                           ", F_R[" + std::to_string(n) + "] = " + std::to_string(flux_R[n]));
                }
            
                // Check for NaN or Inf in lambda_max and fluxes
                if (std::isnan(lambda_max) || std::isinf(lambda_max)) {
                    std::string error_message = "Invalid wave speed detected at (" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) +
                                                ") in context '" + context + "'. lambda_max = " + std::to_string(lambda_max);
            
                    if (abort_if_nan) {
                        Util::Abort(INFO, error_message);  // Abort if NaN or Inf
                    } else {
                        Util::Warning(INFO, error_message);  // Print warning instead of aborting
                    }
                }
            
                // Check for NaN or Inf in the fluxes
                for (int n = 0; n < 5; ++n) {
                    if (std::isnan(flux_L[n]) || std::isinf(flux_L[n]) || std::isnan(flux_R[n]) || std::isinf(flux_R[n])) {
                        std::string error_message = "Invalid flux component detected at (" + std::to_string(i) + ", " + std::to_string(j) + ", " +
                                                    std::to_string(k) + ") in context '" + context + "'. Flux component " + std::to_string(n);
                            if (abort_if_nan) {
                                Util::Abort(INFO, error_message);
                            } else {
                                Util::Warning(INFO, error_message);
                            }
                        }
                    }
                } 
            }


    //Probe ghost cell given target boundary indices
    static void DebugTargetBoundaryAndGhostCell(int i, int j, int k, int lev, const amrex::MultiFab& mf, const amrex::Geometry& geom, const std::string& field_name, const std::string& context, bool abort_if_nan = true) {
        // Set target indices (boundary cell)
        SetTargetDebugLocationIndices(i, j, k);

        amrex::IntVect domain_lo = geom.Domain().smallEnd();  // Lower boundary
        amrex::IntVect domain_hi = geom.Domain().bigEnd();    // Upper boundary

        // Calculate ghost cell indices based on boundary
        int ghost_i = (i == domain_lo[0]) ? i - 1 : (i == domain_hi[0] ? i + 1 : i);
        int ghost_j = (j == domain_lo[1]) ? j - 1 : (j == domain_hi[1] ? j + 1 : j);
        int ghost_k = (k == domain_lo[2]) ? k - 1 : (k == domain_hi[2] ? k + 1 : k);

        // Check if this ghost cell is the target
        if (IsTargetLocation(ghost_i, ghost_j, ghost_k)) {
            auto const& arr = mf.array(lev);
            amrex::Real value = arr(ghost_i, ghost_j, ghost_k);

            // Debug message
            Util::Message(INFO, "Debug (Ghost Cell) for boundary cell (" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) +
                                   ") at ghost location (" + std::to_string(ghost_i) + ", " + std::to_string(ghost_j) + ", " +
                                   (AMREX_SPACEDIM == 3 ? std::to_string(ghost_k) : "N/A") + ") for field '" + field_name +
                                   "' in context '" + context + "'. Value: " + std::to_string(value));

            // Check for NaN/Inf values
            if (std::isnan(value) || std::isinf(value)) {
                std::string error_message = "Invalid value detected at ghost cell (" + std::to_string(ghost_i) + ", " + std::to_string(ghost_j) +
                                            ", " + (AMREX_SPACEDIM == 3 ? std::to_string(ghost_k) : "N/A") + ") for field '" +
                                            field_name + "' in context '" + context + "'. Value = " + std::to_string(value);
                if (abort_if_nan) {
                    Util::Abort(INFO, error_message);
                } else {
                    Util::Warning(INFO, error_message);
                }
            }
        }
    }

    //Sweep entire boundary direction and debug variable values
    static void DebugSweepBoundaryInDirection(const amrex::MultiFab& mf, const amrex::Geometry& geom, int lev, const std::string& field_name, const std::string& context, const std::string& direction, bool abort_if_nan = true) {
        const amrex::Box& domain = geom.Domain();  // Computational domain
        int num_components = mf.nComp();           // Number of components

        bool sweep_lo = (direction == "lo");
        bool sweep_hi = (direction == "hi");

        for (amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.growntilebox();
            auto const& arr = mf.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                bool is_boundary_x = (sweep_lo && i == domain.smallEnd(0)) || (sweep_hi && i == domain.bigEnd(0));
                bool is_boundary_y = (sweep_lo && j == domain.smallEnd(1)) || (sweep_hi && j == domain.bigEnd(1));
                bool is_boundary_z = AMREX_SPACEDIM == 3 && ((sweep_lo && k == domain.smallEnd(2)) || (sweep_hi && k == domain.bigEnd(2)));

                if (is_boundary_x || is_boundary_y || is_boundary_z) {
                    if (IsTargetLocation(i, j, k)) {
                        for (int n = 0; n < num_components; ++n) {
                            amrex::Real value = arr(i, j, k, n);
                            Util::Message(INFO, "Debug (Boundary Sweep) at (" + std::to_string(i) + ", " + std::to_string(j) + ", " +
                                                   (AMREX_SPACEDIM == 3 ? std::to_string(k) : "N/A") + ") Component " +
                                                   std::to_string(n) + ": " + std::to_string(value));

                            if (std::isnan(value) || std::isinf(value)) {
                                std::string error_message = "Invalid value detected at (" + std::to_string(i) + ", " + std::to_string(j) + ", " +
                                                            (AMREX_SPACEDIM == 3 ? std::to_string(k) : "N/A") + ") for field '" +
                                                            field_name + "' in context '" + context + "'. Value = " + std::to_string(value);
                                if (abort_if_nan) {
                                    Util::Abort(INFO, error_message);
                                } else {
                                    Util::Warning(INFO, error_message);
                                }
                            }
                        }
                    }
                }
            });
        }

        Util::Message(INFO, "Completed boundary sweep in direction '" + direction + "' for field: " + field_name);
    }            

// Set the target indices for boundary cell
//Util::ScimitarX_Util::Debug::SetDebugTargetLocationIndices(0, 0, 0);

// Debug ghost cell corresponding to (0, 0, 0)
//Util::ScimitarX_Util::Debug::DebugTargetBoundaryAndGhostCell(0, 0, 0, lev, *QVec_mf[lev], geom[lev], "QVec", "Boundary Ghost Cell Debug", true);
// Set a new target location
//Util::ScimitarX_Util::Debug::SetTargetIndices(0, domain_lo[1], domain_lo[2]);

// Sweep and debug the lower boundary
//Util::ScimitarX_Util::Debug::DebugSweepBoundaryInDirection(*QVec_mf[lev], geom[lev], lev, "QVec", "Boundary Sweep Debug", "lo", true);

    };


};
}

#endif // UTIL_SCIMITARX_UTIL_H

