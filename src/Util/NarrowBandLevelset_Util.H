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

    /*static void PrintZeroLS(const amrex::iMultiFab& imf, const std::string& context) {
        // Print the debug context
        printf("%s\n", context.c_str());
    
        for (amrex::MFIter mfi(imf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.tilebox();
            auto const& arr = imf.const_array(mfi);  // Use `const_array` since imf is `const`
    
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                printf("i: %d, j: %d, Zerols_imf(i,j,k): %d\n", i, j, k, arr(i, j, k));
            });
        }
    }*/   

    /*for (amrex::MFIter mfi(*Zerols_imf, false); mfi.isValid(); ++mfi) { 
        const amrex::Box& bx = mfi.validbox();
        auto const& zerols_arr = Zerols_imf->array(mfi);
        auto const& ls_arr = ls_mf[lev]->array(mfi);
    
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
            for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                printf("After Update: Zerols_imf(%d, %d, 0) = %d\n", i, j, zerols_arr(i, j, 0));
                printf("After Update: ls_mf(%d, %d, 0) = %f\n", i, j, ls_arr(i, j, 0, ls_id));
            }
        }
    }*/
    
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
