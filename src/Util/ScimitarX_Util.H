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


       // Debug function to check eigenvector matrix validity
       AMREX_GPU_DEVICE
       static void DebugEigenvectorMatrix(
           int i, int j, int k,
           const Set::MultiMatrix& L_n,
           const Set::MultiVector& avg_state,
           const std::string& context = "Left Eigenvector Computation",
           bool abort_if_nan = false,
           bool enableDebug = false)
       {
           if (!enableDebug || !IsTargetLocation(i, j, k)) return;
           
           // Check if matrix contains NaN
           bool has_nan = false;
           for (int r = 0; r < L_n.rows(); ++r) {
               for (int c = 0; c < L_n.cols(); ++c) {
                   if (std::isnan(L_n(r, c)) || std::isinf(L_n(r, c))) {
                       has_nan = true;
                       break;
                   }
               }
               if (has_nan) break;
           }
           
           if (has_nan || enableDebug) {
               Message(INFO, "Eigenvector Matrix Debug at (" + 
                   std::to_string(i) + ", " + std::to_string(j) + ", " + 
                   std::to_string(k) + ") - Context: " + context);
               
               // Print matrix details
               Message(INFO, "  Matrix dimensions: " + std::to_string(L_n.rows()) + 
                       "×" + std::to_string(L_n.cols()));
               
               // Print matrix norm
               Set::Scalar frobenius_norm = 0.0;
               for (int r = 0; r < L_n.rows(); ++r) {
                   for (int c = 0; c < L_n.cols(); ++c) {
                       frobenius_norm += L_n(r, c) * L_n(r, c);
                   }
               }
               frobenius_norm = std::sqrt(frobenius_norm);
               Message(INFO, "  Matrix Frobenius norm: " + std::to_string(frobenius_norm));
               
               // Print full matrix for detailed inspection
               for (int r = 0; r < L_n.rows(); ++r) {
                   std::string row_str = "  Row " + std::to_string(r) + ": ";
                   for (int c = 0; c < L_n.cols(); ++c) {
                       row_str += std::to_string(L_n(r, c)) + " ";
                   }
                   Message(INFO, row_str);
               }
               
               // Print average state
               Message(INFO, "  Average State:");
               for (int n = 0; n < avg_state.size(); ++n) {
                   Message(INFO, "    [" + std::to_string(n) + "]: " + 
                          std::to_string(avg_state(n)));
               }
               
               if (has_nan && abort_if_nan) {
                   Util::Abort(INFO, "NaN detected in eigenvector matrix in context: " + context);
               }
           }
       }
       
       // Debug function to check sound speed computation
       AMREX_GPU_DEVICE
       static void DebugSoundSpeedComputation(
           int i, int j, int k,
           Set::Scalar H, Set::Scalar Va2, Set::Scalar gamma,
           Set::Scalar a2, Set::Scalar a,
           const std::string& context = "Sound Speed Computation",
           bool abort_on_negative = true,
           bool enableDebug = false)
       {
           if (!enableDebug || !IsTargetLocation(i, j, k)) return;
           
           Set::Scalar raw_a2 = H - 0.5 * Va2;
           bool is_negative = (raw_a2 <= 0.0);
           
           if (is_negative || enableDebug) {
               Message(INFO, "Sound Speed Debug at (" + 
                   std::to_string(i) + ", " + std::to_string(j) + ", " + 
                   std::to_string(k) + ") - Context: " + context);
               
               Message(INFO, "  Total enthalpy (H): " + std::to_string(H));
               Message(INFO, "  Velocity magnitude squared (Va2): " + std::to_string(Va2));
               Message(INFO, "  H - 0.5*Va2: " + std::to_string(raw_a2));
               Message(INFO, "  Gamma: " + std::to_string(gamma));
               
               if (is_negative) {
                   Message(INFO, "  WARNING: Negative value detected for sound speed squared!");
                   
                   // Apply a correction
                   a2 = std::max(gamma * raw_a2, 1e-10);
                   a = std::sqrt(a2);
                   
                   Message(INFO, "  Corrected a²: " + std::to_string(a2));
                   Message(INFO, "  Corrected a: " + std::to_string(a));
                   
                   if (abort_on_negative) {
                       Util::Abort(INFO, "Negative sound speed squared detected in context: " + context);
                   }
               }
           }
       }
       
       // Debug function to verify transformation of stencils to characteristic space
       AMREX_GPU_DEVICE
       static void DebugCharacteristicTransformation(
           int i, int j, int k,
           const Set::MultiMatrix& stencil_matrix,
           const Set::MultiMatrix& char_stencil_matrix,
           const Set::MultiMatrix& L_n,
           int stencil_point,
           const std::string& context = "Characteristic Transformation",
           bool abort_if_nan = true,
           bool enableDebug = false)
       {
           if (!enableDebug || !IsTargetLocation(i, j, k)) return;
           
           // Check if characteristic stencil contains NaN
           bool has_nan = false;
           for (int r = 0; r < char_stencil_matrix.rows(); ++r) {
               for (int c = 0; c < char_stencil_matrix.cols(); ++c) {
                   if (std::isnan(char_stencil_matrix(r, c)) || std::isinf(char_stencil_matrix(r, c))) {
                       has_nan = true;
                       break;
                   }
               }
               if (has_nan) break;
           }
           
           if (has_nan || enableDebug) {
               Message(INFO, "Characteristic Transformation Debug at (" + 
                   std::to_string(i) + ", " + std::to_string(j) + ", " + 
                   std::to_string(k) + ") - Stencil Point: " + std::to_string(stencil_point) +
                   " - Context: " + context);
               
               // Print original stencil values for this point
               Message(INFO, "  Original stencil values:");
               for (int c = 0; c < stencil_matrix.cols(); ++c) {
                   Message(INFO, "    Component " + std::to_string(c) + ": " + 
                          std::to_string(stencil_matrix(stencil_point, c)));
               }
               
               // Print transformed values
               Message(INFO, "  Transformed characteristic values:");
               for (int c = 0; c < char_stencil_matrix.cols(); ++c) {
                   Message(INFO, "    Component " + std::to_string(c) + ": " + 
                          std::to_string(char_stencil_matrix(stencil_point, c)));
               }
               
               if (has_nan && abort_if_nan) {
                   Util::Abort(INFO, "NaN detected in characteristic transformation in context: " + context);
               }
           }
       }
       
       // Debug function to check split flux computation
       AMREX_GPU_DEVICE
       static void DebugSplitFluxes(
           int i, int j, int k,
           const Set::Scalar* FL, const Set::Scalar* FR,
           const Set::Scalar* QL_cons, const Set::Scalar* QR_cons,
           Set::Scalar lambda_max,
           const Set::Scalar* FluxPos, const Set::Scalar* FluxNeg,
           int num_components,
           const std::string& context = "Split Flux Computation",
           bool abort_if_nan = true,
           bool enableDebug = false)
       {
           if (!enableDebug || !IsTargetLocation(i, j, k)) return;
           
           // Check for NaN in split fluxes
           bool has_nan = false;
           for (int n = 0; n < num_components; ++n) {
               if (std::isnan(FluxPos[n]) || std::isinf(FluxPos[n]) || 
                   std::isnan(FluxNeg[n]) || std::isinf(FluxNeg[n])) {
                   has_nan = true;
                   break;
               }
           }
           
           if (has_nan || enableDebug) {
               Message(INFO, "Split Flux Debug at (" + 
                   std::to_string(i) + ", " + std::to_string(j) + ", " + 
                   std::to_string(k) + ") - Context: " + context);
               
               Message(INFO, "  Lambda Max: " + std::to_string(lambda_max));
               
               // Print original fluxes, conservative variables, and split fluxes
               for (int n = 0; n < num_components; ++n) {
                   Message(INFO, "  Component " + std::to_string(n) + ":");
                   Message(INFO, "    F_L: " + std::to_string(FL[n]) + 
                          ", F_R: " + std::to_string(FR[n]));
                   Message(INFO, "    U_L: " + std::to_string(QL_cons[n]) + 
                          ", U_R: " + std::to_string(QR_cons[n]));
                   Message(INFO, "    F+: " + std::to_string(FluxPos[n]) + 
                          ", F-: " + std::to_string(FluxNeg[n]));
                   
                   // Check flux formulation
                   Set::Scalar computed_pos = 0.5 * FL[n] + 0.5 * lambda_max * QL_cons[n];
                   Set::Scalar computed_neg = 0.5 * FR[n] - 0.5 * lambda_max * QR_cons[n];
                   
                   Message(INFO, "    Computed F+: " + std::to_string(computed_pos) + 
                          ", Computed F-: " + std::to_string(computed_neg));
                   
                   // Check for large discrepancies
                   if (std::abs(computed_pos - FluxPos[n]) > 1e-10 || 
                       std::abs(computed_neg - FluxNeg[n]) > 1e-10) {
                       Message(INFO, "    WARNING: Discrepancy in flux computation!");
                   }
               }
               
               if (has_nan && abort_if_nan) {
                   Util::Abort(INFO, "NaN detected in split fluxes in context: " + context);
               }
           }
       }
       
       // Debug function for average state check
       AMREX_GPU_DEVICE
       static void DebugAverageState(
           int i, int j, int k,
           const Set::MultiVector& WL, const Set::MultiVector& WR,
           const Set::MultiVector& Wavg,
           const std::string& context = "Roe Average Computation",
           bool abort_if_invalid = true,
           bool enableDebug = false)
       {
           if (!enableDebug || !IsTargetLocation(i, j, k)) return;
           
           // Check for non-physical values in average state
           bool invalid_state = false;
           if (Wavg(0) <= 0.0) {  // Negative density
               invalid_state = true;
           }
           
           // Check for NaN values
           bool has_nan = false;
           for (int n = 0; n < Wavg.size(); ++n) {
               if (std::isnan(Wavg(n)) || std::isinf(Wavg(n))) {
                   has_nan = true;
                   break;
               }
           }
           
           if (invalid_state || has_nan || enableDebug) {
               Message(INFO, "Average State Debug at (" + 
                   std::to_string(i) + ", " + std::to_string(j) + ", " + 
                   std::to_string(k) + ") - Context: " + context);
               
               // Print left state
               Message(INFO, "  Left State (WL):");
               for (int n = 0; n < WL.size(); ++n) {
                   Message(INFO, "    [" + std::to_string(n) + "]: " + 
                          std::to_string(WL(n)));
               }
               
               // Print right state
               Message(INFO, "  Right State (WR):");
               for (int n = 0; n < WR.size(); ++n) {
                   Message(INFO, "    [" + std::to_string(n) + "]: " + 
                          std::to_string(WR(n)));
               }
               
               // Print average state
               Message(INFO, "  Average State (Wavg):");
               for (int n = 0; n < Wavg.size(); ++n) {
                   Message(INFO, "    [" + std::to_string(n) + "]: " + 
                          std::to_string(Wavg(n)));
               }
               
               // Additional physics-based diagnostics
               if (Wavg.size() >= 5) {  // Assuming standard Euler variables
                   Set::Scalar rho = Wavg(0);
                   Set::Scalar u = Wavg(1);
                   Set::Scalar v = Wavg(2);
                   Set::Scalar w = (Wavg.size() > 3) ? Wavg(3) : 0.0;
                   Set::Scalar H = Wavg(4);
                   
                   Set::Scalar Va2 = u*u + v*v + w*w;
                   Set::Scalar raw_a2 = H - 0.5*Va2;
                   
                   Message(INFO, "  Physics diagnostics:");
                   Message(INFO, "    Density: " + std::to_string(rho));
                   Message(INFO, "    Velocity magnitude squared: " + std::to_string(Va2));
                   Message(INFO, "    Enthalpy: " + std::to_string(H));
                   Message(INFO, "    H - 0.5*Va2: " + std::to_string(raw_a2));
                   
                   if (raw_a2 <= 0.0) {
                       Message(INFO, "    WARNING: Non-physical state detected (negative sound speed squared)");
                   }
               }
               
               if ((invalid_state || has_nan) && abort_if_invalid) {
                   Util::Abort(INFO, "Invalid average state detected in context: " + context);
               }
           }
       }
       
       // Debug function for the full characteristic-based reconstruction process
       AMREX_GPU_DEVICE
       static void DebugCharacteristicBasedReconstruction(
           int i, int j, int k,
           const Set::MultiMatrix& stencil_matrix,
           const Set::MultiMatrix& char_stencil_matrix,
           const Set::MultiVector& reconstructed_values,
           int direction,
           const std::string& context = "Characteristic-based Reconstruction",
           bool abort_if_nan = true,
           bool enableDebug = false)
       {
           if (!enableDebug || !IsTargetLocation(i, j, k)) return;
           
           // Check for NaN in reconstructed values
           bool has_nan = false;
           for (int n = 0; n < reconstructed_values.size(); ++n) {
               if (std::isnan(reconstructed_values(n)) || std::isinf(reconstructed_values(n))) {
                   has_nan = true;
                   break;
               }
           }
           
           if (has_nan || enableDebug) {
               Message(INFO, "Characteristic Reconstruction Debug at (" + 
                   std::to_string(i) + ", " + std::to_string(j) + ", " + 
                   std::to_string(k) + ") - Direction: " + std::to_string(direction) +
                   " - Context: " + context);
               
               // Print input stencil summary
               Message(INFO, "  Input Stencil Summary:");
               Message(INFO, "    Stencil size: " + std::to_string(stencil_matrix.rows()) + 
                       "×" + std::to_string(stencil_matrix.cols()));
               
               // Print characteristic stencil summary
               Message(INFO, "  Characteristic Stencil Summary:");
               Message(INFO, "    Stencil size: " + std::to_string(char_stencil_matrix.rows()) + 
                       "×" + std::to_string(char_stencil_matrix.cols()));
               
               // Print reconstructed values
               Message(INFO, "  Reconstructed Values:");
               for (int n = 0; n < reconstructed_values.size(); ++n) {
                   Message(INFO, "    Component " + std::to_string(n) + ": " + 
                          std::to_string(reconstructed_values(n)));
               }
               
               if (has_nan && abort_if_nan) {
                   Util::Abort(INFO, "NaN detected in characteristic reconstruction in context: " + context);
               }
           }
       }
       
    };     
};
}

#endif // UTIL_SCIMITARX_UTIL_H
