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

    // NOTE: When modifying solver types, ensure the corresponding variables are added here to enable dynamic indexing.
    struct getVariableIndex {
        std::map<int, int> variableIndexMap; // Mapping from variable enums to indices
        int NVAR_MAX = 0;                    // Total number of variables
        int DENS = -1;                       // Index for Density
        int UVEL = -1;                       // Index for U-velocity
        int VVEL = -1;                       // Index for V-velocity
        int WVEL = -1;                       // Index for W-velocity
        int IE = -1;                         // Index for Internal Energy
    }; // End of Debug class

    // Nested Debug class containing debug APIs
    class Debug {
        
    public:

        static amrex::IntVect& GetTargetLocation() {
                static amrex::IntVect target_location = amrex::IntVect::TheZeroVector();
                return target_location;
        }

        // 1. Function to set the target debug location
        static void SetTargetDebugLocation(double x, double y, [[maybe_unused]] double z, const amrex::Geometry& geom, bool enableDebug = true) {
            if (!enableDebug) return;
            
            // Create a RealVect using AMREX_D_DECL to handle different dimensions
            amrex::Array<amrex::Real, AMREX_SPACEDIM> pos = {AMREX_D_DECL(x, y, z)};
        
            // Convert to cell index (i, j, k)
            GetTargetLocation() = geom.CellIndex(pos.data());  // Use .getVect() to pass the internal pointer

            Util::Message(INFO, "Target debug location set to: (" +
                                std::to_string(GetTargetLocation()[0]) + ", " +
                                std::to_string(GetTargetLocation()[1]) + ", " +
                                (AMREX_SPACEDIM == 3 ? std::to_string(GetTargetLocation()[2]) : "N/A") + ")");
        }

        static void SetTargetDebugLocationIndices(int i, int j, [[maybe_unused]] int k, bool enableDebug = true) {
            if (!enableDebug) return;
            
            // Set the target location directly using grid indices (i, j, k)
            GetTargetLocation() = amrex::IntVect(AMREX_D_DECL(i, j, k));  // Create IntVect with (i, j, k)

            Util::Message(INFO, "Target debug location set using indices: (" +
                                std::to_string(GetTargetLocation()[0]) + ", " +
                                std::to_string(GetTargetLocation()[1]) + ", " +
                                (AMREX_SPACEDIM == 3 ? std::to_string(GetTargetLocation()[2]) : "N/A") + ")");
        }

        static bool IsTargetLocation(int i, int j, int k) {
            if (AMREX_SPACEDIM >= 1 && GetTargetLocation()[0] != i) return false;
            if (AMREX_SPACEDIM >= 2 && GetTargetLocation()[1] != j) return false;
            if (AMREX_SPACEDIM == 3 && GetTargetLocation()[2] != k) return false;
            return true;
        }

        // Inline debug function that checks the current cell's value
        template <typename T>
        static void DebugValuesIfTarget(int i, int j, int k, const T& value, const std::string& field_name,
                                        const std::string& context, bool abort_if_nan = false, int component = -1,
                                        bool enableDebug = false) {
            if (!enableDebug || !IsTargetLocation(i, j, k)) return;
                                            
            // Build debug message
            std::string debug_message = "Debug at (" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) +
                                        ") in context '" + context + "' for field '" + field_name + "'";
    
            // Add component information if applicable
            if (component >= 0) {
                debug_message += " (Component: " + std::to_string(component) + ")";
            }
    
            // Append the value
            debug_message += " = " + std::to_string(value);
    
            // Print debug message
            Util::Message(INFO, debug_message);
    
            // Check for NaN or invalid value (only for floating-point types)
            if constexpr (std::is_floating_point<T>::value) {
                if (std::isnan(value) || std::isinf(value)) {
                    std::string error_message = "Invalid value detected at (" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) +
                                                ") for field '" + field_name + "' in context '" + context + "'. Value = " + std::to_string(value);
                    if (component >= 0) {
                        error_message += " (Component: " + std::to_string(component) + ")";
                    }
    
                    if (abort_if_nan) {
                        Util::Abort(INFO, error_message);
                    } else {
                        Util::Warning(INFO, error_message);
                    }
                }
            }
        }


        static void CheckNaNAndAbort(int i, int j, int k, const Set::Scalar& value, const std::string& field_name, 
                                    const std::string& context, bool abort_if_nan = true, bool enableDebug = false) {
            if (!enableDebug) return;
            
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
                                                                    const std::string& context = "Unknown", bool enableDebug = false) {
            if (!enableDebug || !IsTargetLocation(i, j, k)) return;
            
            std::ostringstream message;
            message << "Debugging ComputeInternalEnergyFromDensityAndPressure (" << context << ") at cell: (" 
                    << i << ", " << j << ", " << k << ") at level " << lev << "\n";
            message << "Context: " << context << "\n";
            message << "Density: " << density << "\n";
            message << "Pressure: " << pressure << "\n";
            message << "Gamma: " << gamma << "\n";

            if (mode == "print") {
                Util::Message(INFO, message.str());
            } else if (mode == "check") {
                if (density <= 0.0 || pressure < 0.0 || std::isnan(pressure) || std::isnan(density) || std::isnan(gamma)) {
                    message << "Error: Invalid values detected in context: " << context << "\n";
                    if (abort_on_checkpoint) {
                        Util::Abort(INFO, message.str(), "Checkpoint reached. Aborting.");
                    } else {
                        Util::Warning(INFO, message.str());
                    }
                } else {
                    Util::Message(INFO, "Values are valid: " + message.str());
                }
            } else {
                Util::Warning(INFO, "Unknown debug mode specified. Defaulting to 'print'.");
                Util::Message(INFO, message.str());
            }

            if (abort_on_checkpoint) {
                Util::Abort(INFO, message.str(), "Checkpoint reached. Aborting as requested.");
            }
        }
        
        static void DebugEulerForwardStep(int i, int j, int k,
                                Set::Scalar xflux_right, Set::Scalar xflux_left,
                                Set::Scalar yflux_top, Set::Scalar yflux_bottom,
                                Set::Scalar zflux_front, Set::Scalar zflux_back,
                                Set::Scalar q_old, Set::Scalar dt,
                                const std::string& context, bool abort_if_nan = true, bool enableDebug = false) {
            if (!enableDebug || !IsTargetLocation(i, j, k)) return;

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
            
        static void DebugLocalLaxFriedrichsFlux(int i, int j, int k, Set::Scalar number_of_components,
                                                Set::Scalar rho_L, Set::Scalar u_L, Set::Scalar v_L, Set::Scalar w_L, Set::Scalar e_L, Set::Scalar pressure_L,
                                                Set::Scalar rho_R, Set::Scalar u_R, Set::Scalar v_R, Set::Scalar w_R, Set::Scalar e_R, Set::Scalar pressure_R,
                                                Set::Scalar c_L, Set::Scalar c_R, Set::Scalar lambda_max,
                                                Set::Scalar flux_L[], Set::Scalar flux_R[],
                                                const std::string& context, bool abort_if_nan = true, bool enableDebug = false) {
            if (!enableDebug || !IsTargetLocation(i, j, k)) return;

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

        static void DebugGhostCellsAndBoundary(const amrex::MultiFab& mf, const amrex::Geometry& geom, int /*lev*/, 
                                             const std::string& field_name, const std::string& context, 
                                             bool abort_if_nan = true, bool enableDebug = false) {
            if (!enableDebug) return;
            
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

        static void CheckNegativePressureAndAbort(int i, int j, int k, const Set::Scalar& pressure, const std::string& field_name, 
                                                const std::string& context, bool abort_if_negative = true, bool enableDebug = false) {
            if (!enableDebug || !IsTargetLocation(i, j, k)) return;
            
            // Check if the pressure value is negative
            if (pressure < 0.0) {
                std::string error_message = "Negative pressure detected at (" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) +
                                            ") for field '" + field_name + "' in context '" + context + "'. Pressure = " + std::to_string(pressure);

                if (abort_if_negative) {
                    Util::Abort(INFO, error_message);  // Abort with the error message if pressure is negative
                } else {
                    Util::Warning(INFO, error_message);  // Print a warning instead of aborting
                }
            }
        }

        static void DebugIndicesAndBounds(int i, int j, int k,
                                        const int index[3], const int left_index[3], const int right_index[3],
                                        const int lower_bounds[3], const int upper_bounds[3],
                                        const std::string& context, bool enableDebug = false) {
            if (!enableDebug || !IsTargetLocation(i, j, k)) return;
            
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

        // Add these debug functions to the ScimitarX_Util::Debug class in ScimitarX_Util.H

        template <typename MatrixType>
        static void DebugMatrixIfTarget(int i, int j, int k, const MatrixType& matrix, 
                                       const std::string& name, const std::string& context, bool enableDebug = false) {
            if (!enableDebug || !IsTargetLocation(i, j, k)) return;
            
            // Print determinant
            Set::Scalar det = matrix.determinant();
            Util::Message(INFO, context + " - Matrix " + name + " at (" + 
                         std::to_string(i) + ", " + std::to_string(j) + ", " + 
                         std::to_string(k) + "):");
            Util::Message(INFO, "  Determinant: " + std::to_string(det));
            
            // Check if near-singular
            if (std::abs(det) < 1e-10) {
                Util::Warning(INFO, "Near-singular matrix detected in " + context);
            }
            
            // Check condition number (approximate) using singular values
            Eigen::JacobiSVD<MatrixType> svd(matrix);
            Set::Scalar cond = svd.singularValues()(0) / 
                             svd.singularValues()(svd.singularValues().size()-1);
            Util::Message(INFO, "  Condition number: " + std::to_string(cond));
            
            // Print matrix entries
            Util::Message(INFO, "  Matrix entries:");
            for (int r = 0; r < matrix.rows(); r++) {
                for (int c = 0; c < matrix.cols(); c++) {
                    Util::Message(INFO, "    [" + std::to_string(r) + "," + 
                                 std::to_string(c) + "]: " + std::to_string(matrix(r,c)));
                    CheckNaNAndAbort(i, j, k, matrix(r,c), 
                                   name + "[" + std::to_string(r) + "," + std::to_string(c) + "]", 
                                   context);
                }
            }
        }

        template <typename VectorType>
        static void DebugVectorIfTarget(int i, int j, int k, const VectorType& vector, 
                                       const std::string& name, const std::string& context,
                                       bool check_negative_values = false, bool enableDebug = false) {
            if (!enableDebug || !IsTargetLocation(i, j, k)) return;
            
            Util::Message(INFO, context + " - Vector " + name + " at (" + 
                         std::to_string(i) + ", " + std::to_string(j) + ", " + 
                         std::to_string(k) + "):");
            
            // Print vector entries
            Util::Message(INFO, "  Vector entries:");
            for (int r = 0; r < vector.size(); r++) {
                Util::Message(INFO, "    [" + std::to_string(r) + "]: " + 
                             std::to_string(vector(r)));
                
                // Check for NaN values
                CheckNaNAndAbort(i, j, k, vector(r), 
                               name + "[" + std::to_string(r) + "]", 
                               context);
                
                // Optionally check for negative values
                if (check_negative_values && vector(r) < 0.0) {
                    Util::Warning(INFO, "Negative value detected in " + name + 
                                "[" + std::to_string(r) + "] = " + 
                                std::to_string(vector(r)) + " at (" + 
                                std::to_string(i) + ", " + std::to_string(j) + ", " + 
                                std::to_string(k) + ") in context " + context);
                }
            }
        }

        static void CheckDivisionBySmallandAbort(int i, int j, int k, 
                                               const Set::Scalar& denominator, 
                                               const std::string& field_name, 
                                               const std::string& context, 
                                               bool abort_if_small = true,
                                               Set::Scalar threshold = 1e-10,
                                               bool enableDebug = false) {
            if (!enableDebug || !IsTargetLocation(i, j, k)) return;
            
            // Check if the denominator is too small
            if (std::abs(denominator) < threshold) {
                std::string error_message = "Division by small value detected at (" + 
                                          std::to_string(i) + ", " + std::to_string(j) + ", " + 
                                          std::to_string(k) + ") for field '" + field_name + 
                                          "' in context '" + context + "'. Value = " + 
                                          std::to_string(denominator);

                if (abort_if_small) {
                    Util::Abort(INFO, error_message);
                } else {
                    Util::Warning(INFO, error_message);
                }
            }
        }

        static void DebugCharacteristicTransformation(int i, int j, int k,
                                                    const Set::MultiVector& input_vars,
                                                    const Set::MultiMatrix& transform_matrix,
                                                    const Set::MultiVector& output_vars,
                                                    const std::string& transform_type,
                                                    const std::string& context,
                                                    bool enableDebug = false) {
            if (!enableDebug || !IsTargetLocation(i, j, k)) return;
            
            Util::Message(INFO, context + " - " + transform_type + " at (" + 
                         std::to_string(i) + ", " + std::to_string(j) + ", " + 
                         std::to_string(k) + "):");
            
            // Print input variables
            Util::Message(INFO, "  Input variables:");
            for (int r = 0; r < input_vars.size(); r++) {
                Util::Message(INFO, "    [" + std::to_string(r) + "]: " + 
                             std::to_string(input_vars(r)));
            }
            
            // Check transform matrix
            Set::Scalar det = transform_matrix.determinant();
            Util::Message(INFO, "  Transform matrix determinant: " + std::to_string(det));
            CheckDivisionBySmallandAbort(i, j, k, det, "matrix_determinant", context, false, 1e-10, enableDebug);
            
            // Print output variables
            Util::Message(INFO, "  Output variables:");
            for (int r = 0; r < output_vars.size(); r++) {
                Util::Message(INFO, "    [" + std::to_string(r) + "]: " + 
                             std::to_string(output_vars(r)));
                CheckNaNAndAbort(i, j, k, output_vars(r), 
                               "output_var_" + std::to_string(r), context, true, enableDebug);
            }
        }
    

        // Debug WENO smoothness indicators calculation
        template <typename StencilArray, typename IndicatorArray, typename OffsetArray, typename TermWeightArray>
        AMREX_GPU_DEVICE
        static void DebugWENOSmoothnessIndicators(
            int i, int j, int k,
            const StencilArray& stencil,
            const IndicatorArray& indicators,
            const std::vector<Set::MultiMatrix>& coeffMatrices,
            const TermWeightArray& termWeights,
            const OffsetArray& stencilOffsets,
            int numCandidateStencils,
            int totalPoints,
            int numDifferenceTerms,
            bool enableDebug = false
        )
        {
            if (!enableDebug || !IsTargetLocation(i, j, k)) return;
            
            // Print stencil values with dynamic size
            Message(INFO, "WENO Stencil Values at (" + 
                std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + "):");
                
            for (int idx = 0; idx < totalPoints; ++idx) {
                int rel_idx = idx - totalPoints/2;  // Relative position
                Message(INFO, "  Stencil[" + std::to_string(rel_idx) + "] = " + 
                       std::to_string(stencil[idx]));
            }
            
            // Print indicator details for each stencil
            for (int r = 0; r < numCandidateStencils; ++r) {
                Message(INFO, "WENO Stencil " + std::to_string(r) + " Smoothness: " + 
                       std::to_string(indicators[r]));
                
                // Check for potential issues
                if (indicators[r] < 1e-10) {
                    Message(INFO, "  WARNING: Very small smoothness indicator detected!");
                }
                else if (indicators[r] > 1e10) {
                    Message(INFO, "  WARNING: Very large smoothness indicator detected!");
                }
                
                // Print details about the terms used in this smoothness indicator
                if (numDifferenceTerms > 0) {
                    Message(INFO, "  Difference terms weights:");
                    for (int term = 0; term < numDifferenceTerms; ++term) {
                        Message(INFO, "    Term " + std::to_string(term) + 
                               " weight: " + std::to_string(termWeights[term]));
                    }
                }
            }
        }
        
        // Debug WENO weights calculation - generalized for any number of candidate stencils
        template <typename IndicatorArray, typename WeightArray, typename OptimalWeightArray>
        AMREX_GPU_DEVICE
        static void DebugWENOWeights(
            int i, int j, int k,
            const IndicatorArray& indicators,
            const WeightArray& weights,
            const OptimalWeightArray& optimalWeights,
            Set::Scalar epsilon,
            int power,
            int numCandidateStencils,
            bool enableDebug = false
        )
        {
            if (!enableDebug || !IsTargetLocation(i, j, k)) return;
            
            Message(INFO, "WENO Weights at (" + 
                std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + "):");
                
            Message(INFO, "  Parameters: ε = " + std::to_string(epsilon) + 
                   ", Power = " + std::to_string(power));
            
            // Print indicator, weight, and deviation information for each stencil
            Set::Scalar totalWeight = 0.0;
            for (int r = 0; r < numCandidateStencils; ++r) {
                Set::Scalar deviation = weights[r] - optimalWeights[r];
                Set::Scalar rel_deviation = (std::abs(optimalWeights[r]) > 1e-10) ? 
                                           (deviation / optimalWeights[r]) : deviation;
                
                Message(INFO, "  Stencil " + std::to_string(r) + 
                       ": β = " + std::to_string(indicators[r]) + 
                       ", ω = " + std::to_string(weights[r]) + 
                       ", Optimal ω = " + std::to_string(optimalWeights[r]) +
                       ", Deviation: " + std::to_string(rel_deviation * 100.0) + "%");
                
                if (std::abs(rel_deviation) > 1.0) {
                    Message(INFO, "    WARNING: Weight significantly deviates from optimal!");
                }
                
                totalWeight += weights[r];
            }
            
            // Check if weights sum to 1 (convexity property)
            Set::Scalar weightSumDeviation = std::abs(totalWeight - 1.0);
            if (weightSumDeviation > 1e-10) {
                Message(INFO, "  WARNING: Weights sum to " + std::to_string(totalWeight) + 
                       " (deviation from 1.0: " + std::to_string(weightSumDeviation) + ")");
            } else {
                Message(INFO, "  Weight sum check: OK (sum = " + std::to_string(totalWeight) + ")");
            }
        }
        
        // Debug WENO final reconstruction - generalized for any number of candidate fluxes
        template <typename FluxArray, typename WeightArray, typename OptimalWeightArray>
        AMREX_GPU_DEVICE
        static void DebugWENOReconstruction(
            int i, int j, int k,
            const FluxArray& candidateFluxes,
            const WeightArray& weights,
            const OptimalWeightArray& optimalWeights,
            Set::Scalar reconstructedValue,
            int numCandidateStencils,
            bool enableDebug = false
        )
        {
            if (!enableDebug || !IsTargetLocation(i, j, k)) return;
            
            // Calculate linear reconstruction for comparison
            Set::Scalar linearValue = 0.0;
            for (int r = 0; r < numCandidateStencils; ++r) {
                linearValue += optimalWeights[r] * candidateFluxes[r];
            }
            
            Message(INFO, "WENO Reconstruction at (" + 
                std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + "):");
                
            Message(INFO, "  WENO value: " + std::to_string(reconstructedValue));
            Message(INFO, "  Linear value: " + std::to_string(linearValue));
            
            // Diffusion metric
            Set::Scalar diff = reconstructedValue - linearValue;
            Set::Scalar rel_diff = (std::abs(linearValue) > 1.0e-10) ? 
                                  (diff / std::abs(linearValue)) : 0.0;
            
            Message(INFO, "  Numerical diffusion: " + std::to_string(diff) + 
                   " (" + std::to_string(rel_diff * 100.0) + "%)");
            
            if (std::abs(rel_diff) > 0.1) {
                Message(INFO, "  WARNING: Potential excessive numerical diffusion!");
            }
            
            // Print detailed contribution of each candidate flux
            Message(INFO, "  Flux contributions:");
            for (int r = 0; r < numCandidateStencils; ++r) {
                Message(INFO, "    Stencil " + std::to_string(r) + 
                       ": flux = " + std::to_string(candidateFluxes[r]) + 
                       ", weight = " + std::to_string(weights[r]) + 
                       ", contribution = " + std::to_string(weights[r] * candidateFluxes[r]));
            }
        }
        
        // Debug WENO candidate fluxes calculation - generalized for any stencil configuration
        template <typename StencilArray, typename FluxArray, typename CoeffArray, typename OffsetArray>
        AMREX_GPU_DEVICE 
        static void DebugWENOCandidateFluxes(
            int i, int j, int k,
            const StencilArray& stencil,
            const FluxArray& candidateFluxes,
            const CoeffArray& coefficients,
            const OffsetArray& stencilOffsets,
            int numCandidateStencils,
            int pointsPerStencil,
            int totalPoints,
            bool enableDebug = false
        )
        {
            if (!enableDebug || !IsTargetLocation(i, j, k)) return;
            
            Message(INFO, "WENO Candidate Fluxes at (" + 
                std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + "):");
            
            for (int r = 0; r < numCandidateStencils; ++r) {
                // Collect stencil points for this candidate stencil
                Message(INFO, "  Stencil " + std::to_string(r) + ":");
                
                // Print stencil points
                Message(INFO, "    Points:");
                for (int p = 0; p < pointsPerStencil; ++p) {
                    int stencil_idx = stencilOffsets(r, p) + totalPoints/2; // Assuming TotalPoints/2 center
                    Message(INFO, "      Point " + std::to_string(p) + 
                           " (offset " + std::to_string(stencilOffsets(r, p)) + "): " +
                           std::to_string(stencil[stencil_idx]));
                }
                
                // Print coefficients for this stencil if available
                Message(INFO, "    Coefficients:");
                for (int p = 0; p < pointsPerStencil; ++p) {
                    Message(INFO, "      Coeff " + std::to_string(p) + ": " +
                           std::to_string(coefficients(r, p)));
                }
                
                // Print the resulting candidate flux
                Message(INFO, "    Flux: " + std::to_string(candidateFluxes[r]));
            }
        }

    };     
};
}

#endif // UTIL_SCIMITARX_UTIL_H
