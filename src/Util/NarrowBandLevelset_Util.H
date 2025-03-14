#ifndef UTIL_NarrowbandLevelset_UTIL_H
#define UTIL_NarrowbandLevelset_UTIL_H

#include <string>
#include <iostream>
#include <vector>
#include <array>
#include "Util/Util.H"
#include "Set/Set.H"

namespace Util
{
class NarrowBandLevelset_Util
{
public:
    NarrowBandLevelset_Util () {};

    // Nested Debug class containing debug APIs
    class Debug {
    
    private:
        static amrex::IntVect target_location;  // Store target cell index (i, j, k)
    
    public:
        // 1. Function to set the target debug location
        static void SetTargetDebugLocation(double x, double y, [[maybe_unused]] double z, const amrex::Geometry& geom) {
        // Create a RealVect using AMREX_D_DECL to handle different dimensions
        amrex::Array<amrex::Real, AMREX_SPACEDIM> pos = {AMREX_D_DECL(x, y, z)};
    
        // Convert to cell index (i, j, k)
        target_location = geom.CellIndex(pos.data());  // Use .getVect() to pass the internal pointer

        Util::Message(INFO, "Target debug location set to: (" +
                            std::to_string(target_location[0]) + ", " +
                            std::to_string(target_location[1]) + ", " +
                            (AMREX_SPACEDIM == 3 ? std::to_string(target_location[2]) : "N/A") + ")");
        }

        static void SetTargetDebugLocationIndices(int i, int j, [[maybe_unused]] int k) {
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
        template <typename T>
        static void DebugValuesIfTarget(int i, int j, int k, const T& value, const std::string& field_name,
                                        const std::string& context, bool abort_if_nan, int component = -1) {
            if (IsTargetLocation(i, j, k)) {
                // Build debug message
                std::string debug_message = "Debug at (" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) +
                                            ") in context '" + context + "' for field '" + field_name + "'";
        
                // Append the value
                debug_message += " = " + std::to_string(value);
        
                // Print debug message
                Util::Message(INFO, debug_message);
                
                // Print neighbor values
                //printf("i: %d, j: %d, phi: %f\n", i, j, phi);
                //printf("phi(i-1,j+1): %f, phi(i,j+1): %f, phi(i+1,j+1): %f\n", ls_arr(i-1,j+1,k), ls_arr(i,j+1,k), ls_arr(i+1,j+1,k));
                //printf("phi(i-1,j): %f, phi(i+1,j): %f\n", ls_arr(i-1,j,k), ls_arr(i+1,j,k));
                //printf("phi(i-1,j-1): %f, phi(i,j-1): %f, phi(i+1,j-1): %f\n", ls_arr(i-1,j-1,k), ls_arr(i,j-1,k), ls_arr(i+1,j-1,k)); 
        
                // Check for NaN or invalid value (only for floating-point types)
                if constexpr (std::is_floating_point<T>::value) {
                    if (std::isnan(value) || std::isinf(value)) {
                        std::string error_message = "Invalid value detected at (" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) +
                                                    ") for field '" + field_name + "' in context '" + context + "'. Value = " + std::to_string(value); 
                        if (abort_if_nan) {
                            Util::Abort(INFO, error_message);
                        } else {
                            Util::Warning(INFO, error_message);
                        }
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
            
            //static void Debug   
            
            
            /*static void DebugLocalLaxFriedrichsFlux(int i, int j, int k, Set::Scalar number_of_components,
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
                for (int n = 0; n < number_of_components; ++n) {
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
                for (int n = 0; n < number_of_components; ++n) {
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
            }*/


static void DebugGhostCellsAndBoundary(const amrex::MultiFab& mf, const amrex::Geometry& geom, int /*lev*/, const std::string& field_name, const std::string& context, bool abort_if_nan = true) {
    const amrex::Box& domain = geom.Domain();  // Physical domain without ghost cells

    for (amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.growntilebox();  // Tile box including ghost cells
        auto const& arr = mf.array(mfi);  // Access the MultiFab array

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            bool is_ghost = (i < domain.smallEnd(0) || i > domain.bigEnd(0) ||
                            j < domain.smallEnd(1) || j > domain.bigEnd(1) ||
                            (AMREX_SPACEDIM == 3 && (k < domain.smallEnd(2) || k > domain.bigEnd(2))));

            bool is_boundary = !is_ghost &&  // Must be inside the domain
                                (i == domain.smallEnd(0) || i == domain.bigEnd(0) ||
                                j == domain.smallEnd(1) || j == domain.bigEnd(1) ||
                                (AMREX_SPACEDIM == 3 && (k == domain.smallEnd(2) || k == domain.bigEnd(2))));

            amrex::Real value = arr(i, j, k);

            if (is_ghost) {
                Util::Message(INFO, "Debug (Ghost Cell) at (" + std::to_string(i) + ", " + std::to_string(j) + ", " +
                                        (AMREX_SPACEDIM == 3 ? std::to_string(k) : "N/A") + ") for field '" +
                                        field_name + "' in context '" + context + "'. Value: " + std::to_string(value));
            }
            else if (is_boundary) {
                Util::Message(INFO, "Debug (Boundary Cell) at (" + std::to_string(i) + ", " + std::to_string(j) + ", " +
                                        (AMREX_SPACEDIM == 3 ? std::to_string(k) : "N/A") + ") for field '" +
                                        field_name + "' in context '" + context + "'. Value: " + std::to_string(value));
            }

            // Check for NaN or Inf values
            if (std::isnan(value) || std::isinf(value)) {
                std::string error_message = "Invalid value detected at (" + std::to_string(i) + ", " +
                                            std::to_string(j) + ", " + (AMREX_SPACEDIM == 3 ? std::to_string(k) : "N/A") +
                                            ") for field '" + field_name + "' in context '" + context +
                                            "'. Value = " + std::to_string(value);
                if (abort_if_nan) {
                    Util::Abort(INFO, error_message);
                } else {
                    Util::Warning(INFO, error_message);
                }
            }
        });
    }

    Util::Message(INFO, "Completed debug of ghost and boundary cells for field: " + field_name);
}

static void DebugIndicesAndBounds(int i, int j, int k,
                                const int index[3], const int left_index[3], const int right_index[3],
                                const int lower_bounds[3], const int upper_bounds[3],
                                const std::string& context) {
    std::ostringstream message;
    message << "=== Debug Indices and Bounds in Context: " << context << " ===\n";
    message << "Cell Center Indices (i, j, k): (" << i << ", " << j << ", " << k << ")\n";

    message << "Index: (" << index[0] << ", " << index[1] << ", " << index[2] << ")\n";
    message << "Left Index: (" << left_index[0] << ", " << left_index[1] << ", " << left_index[2] << ")\n";
    message << "Right Index: (" << right_index[0] << ", " << right_index[1] << ", " << right_index[2] << ")\n";

    message << "Lower Bounds: (" << lower_bounds[0] << ", " << lower_bounds[1] << ", " << lower_bounds[2] << ")\n";
    message << "Upper Bounds: (" << upper_bounds[0] << ", " << upper_bounds[1] << ", " << upper_bounds[2] << ")\n";

    message << "==========================================\n";

    Util::Message(INFO, message.str());  // Print debug message to the log
}

    };


};
}

#endif // UTIL_NarrowbandLevelset_UTIL_H
