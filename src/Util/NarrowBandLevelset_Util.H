#ifndef UTIL_NarrowBandLevelset_UTIL_H
#define UTIL_NarrowBandLevelset_UTIL_H

#include <string>
#include <iostream>
#include <vector>
#include <array>
#include "Util/Util.H"
#include "Util/ScimitarX_Util.H"
#include "Set/Set.H"

namespace Util
{
class NarrowBandLevelset_Util
{
public:
    NarrowBandLevelset_Util () {};
    
    class Debug {
    
    public:

    static void DebugLocalLaxFriedrichsFlux(int i, int j, int k, Set::Scalar number_of_components,
        Set::Scalar u_L, Set::Scalar v_L, Set::Scalar w_L,
        Set::Scalar u_R, Set::Scalar v_R, Set::Scalar w_R, 
        Set::Scalar c_L, Set::Scalar c_R, Set::Scalar lambda_max,
        Set::Scalar flux_L[], Set::Scalar flux_R[],
        const std::string& context, bool abort_if_nan = true) {

        if (Util::ScimitarX_Util::Debug::IsTargetLocation(i, j, k)) {
            // Print debug information for left and right states
            Util::Message(INFO, "Local Lax-Friedrichs Flux Debug - Context: '" + context + "'");
            Util::Message(INFO, "Index: (" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + ")");
            Util::Message(INFO, "  Left State (QL):");
            Util::Message(INFO, "    u_L = " + std::to_string(u_L) + ", v_L = " + std::to_string(v_L) +
                ", w_L = " + std::to_string(w_L));
            Util::Message(INFO, "  Right State (QR):");
            Util::Message(INFO, "    u_R = " + std::to_string(u_R) + ", v_R = " + std::to_string(v_R) +
                ", w_R = " + std::to_string(w_R));
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
    }

    static void CheckTileBoundaryNeighbors(){
    }
    
    static void printNeighbors2D(int i_min, int i_max, int j_min, int j_max, amrex::MultiFab& mf, const amrex::Geometry& geom, int /*lev*/, 
                                const std::string& field_name, const std::string& context){
        
        // Print the debug context
        printf("%s\n", context.c_str());
        
        for (amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.tilebox();  
            /*printf("Tile box: [x: %d, %d], [y: %d, %d]\n",
                bx.smallEnd(0), bx.bigEnd(0),
                bx.smallEnd(1), bx.bigEnd(1));
                printf("done\n");*/
                
            mf.FillBoundary();
            auto const& arr = mf.array(mfi);  // Access the MultiFab array

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                if (i >= i_min  && i <= i_max){
                    if (j >= j_min  && j <= j_max){
                        printf("i: %d, j: %d, %s(i,j): %f\n", i, j, field_name.c_str(), arr(i,j,k));
                        printf("%s(i-1,j+1): %f, %s(i,j+1): %f, %s(i+1,j+1): %f\n", field_name.c_str(), arr(i-1,j+1,k), 
                            field_name.c_str(), arr(i,j+1,k), field_name.c_str(), arr(i+1,j+1,k));
                        printf("%s(i-1,j): %f, %s(i+1,j): %f\n", field_name.c_str(), arr(i-1,j,k), field_name.c_str(), arr(i+1,j,k));
                        printf("%s(i-1,j-1): %f, %s(i,j-1): %f, %s(i+1,j-1): %f\n", field_name.c_str(), arr(i-1,j-1,k), 
                            field_name.c_str(), arr(i,j-1,k), field_name.c_str(), arr(i+1,j-1,k)); 
                    }
                }
            });
        }                       
    }

    }; // end class Debug
}; // end class NarrowBandLevelset_Util
} // end namespace util

#endif // UTIL_NarrowBandLevelset_UTIL_H
