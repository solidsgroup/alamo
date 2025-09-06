#ifndef NUMERIC_SOURCETERM_H
#define NUMERIC_SOURCETERM_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include "Set/Set.H"
#include "Integrator/ScimitarX.H"
#include "Util/Util.H"
#include "Util/ScimitarX_Util.H"
#include "Model/Fluid/Fluid.H"
#include "Numeric/Stencil.H"
#include "Numeric/NumericTypes.H"
#include "Numeric/SolverCapabilities.H"
#include "Numeric/IntegratorVariableAccessLayer.H"
#include "Numeric/SymmetryPreservingRoeAveragingOperations.H"
#include "Integrator/ScimitarX.H"
#include "Numeric/FluxHandler.H"
#include <memory>
#include <cmath>
#include <algorithm>

//forward declaration for ScimitarX
namespace Integrator {
    class ScimitarX;
}

namespace Numeric {

AMREX_GPU_DEVICE
inline void ClampIndices(int indices[3], const int lower_bounds[3], const int upper_bounds[3]);

template <typename IntegratorDerivedClass>
class SourceTerm {
public:

Set::Scalar mu;
static void Parse(SourceTerm& value, IO::ParmParse& pp)
{
pp.query_default("mu", value.mu, 0.0);
}

void ComputeSourceTerm(int lev, IntegratorDerivedClass* solver) const {

    const int nghosts = solver->number_of_ghost_cells;

    IO::ParmParse pp;
    SourceTerm::Parse(*const_cast<SourceTerm*>(this), pp);

    const Set::Scalar dx = solver->geom[lev].CellSize()[0];  // Access the geometry at level `lev`
    const Set::Scalar dy = solver->geom[lev].CellSize()[1];
#if (AMREX_SPACEDIM == 3)
    const Set::Scalar dz = solver->geom[lev].CellSize()[2];
#endif	

    const Set::Scalar inU = solver->variableIndex.UVEL;
    const Set::Scalar inV = solver->variableIndex.VVEL;
#if (AMREX_SPACEDIM == 3)	
	const Set::Scalar inW = solver->variableIndex.WVEL;
#endif	

    // Loop through the grid
    for (amrex::MFIter mfi(*solver->SourceTermVec[lev], false); mfi.isValid(); ++mfi) {

        const amrex::Box& bx_with_ghosts = mfi.growntilebox(nghosts);
        //const amrex::Box& bx_with_ghosts = mfi.validbox();

        // Get arrays for pArr and shear stress output
        auto const& pArr = solver->PVec_mf.Patch(lev, mfi);  
        auto const& sourceTermArr = solver->SourceTermVec.Patch(lev, mfi);

        // Parallel loop over cells
        amrex::ParallelFor(bx_with_ghosts, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
             
            // Finding pos-neg diagonal stresses on faces
            // "Xpos" means the face center on the positive (right) side of the cell along x-rection, etc.
            // These are the shear stress values on the cell faces from the linearly interpolated velocity gradients.

            // tao_xx 

            Set::Scalar tao_xx_Xpos = 2.0*mu*((pArr(i+1, j, k, inU)-pArr(i, j, k, inU))/dx) 
                                        - 2.0/3.0*mu*((pArr(i, j+1, k, inV)+pArr(i+1, j+1, k, inV)-
					pArr(i, j-1, k, inV)-pArr(i+1, j-1, k, inV))/(4.0*dy) + (pArr(i+1, j, k, inU)-pArr(i, j, k, inU))/dx);
#if (AMREX_SPACEDIM == 3)		
            tao_xx_Xpos = tao_xx_Xpos -2.0/3.0*mu*(pArr(i, j, k+1, inW)+pArr(i+1, j, k+1, inW)-
                                       pArr(i, j, k-1, inW)-pArr(i+1, j, k-1, inW))/(4.0*dz) ;
#endif				
						
            Set::Scalar tao_xx_Xneg = 2.0*mu*((pArr(i, j, k, inU)-pArr(i-1, j, k, inU))/dx) 
                                        - 2.0/3.0*mu*((pArr(i, j+1, k, inV)+pArr(i-1, j+1, k, inV)-
                                        pArr(i, j-1, k, inV)-pArr(i-1, j-1, k, inV))/(4.0*dy) + (pArr(i, j, k, inU)-pArr(i-1, j, k, inU))/dx);
#if (AMREX_SPACEDIM == 3)		
            tao_xx_Xneg = tao_xx_Xneg -2.0/3.0*mu*(pArr(i, j, k+1, inW)+pArr(i-1, j, k+1, inW)-
                                        pArr(i, j, k-1, inW)-pArr(i-1, j, k-1, inW))/(4.0*dz) ;
#endif
            // tao_yy
			
            Set::Scalar tao_yy_Ypos = 2.0*mu*((pArr(i, j+1, k, inV)-pArr(i, j, k, inV))/dy) 
                                        - 2.0/3.0*mu*((pArr(i+1, j, k, inU)+pArr(i+1, j+1, k, inU)-
                                        pArr(i-1, j, k, inU)-pArr(i-1, j+1, k, inU))/(4.0*dx) + (pArr(i, j+1, k, inV)-pArr(i, j, k, inV))/dy);
#if (AMREX_SPACEDIM == 3)		
            tao_yy_Ypos = tao_yy_Ypos -2.0/3.0*mu*(pArr(i, j, k+1, inW)+pArr(i, j+1, k+1, inW)-
										pArr(i, j, k-1, inW)-pArr(i, j+1, k-1, inW))/(4.0*dz) ;
#endif										
            Set::Scalar tao_yy_Yneg = 2.0*mu*((pArr(i, j, k, inV)-pArr(i, j-1, k, inV))/dy) 
                                        - 2.0/3.0*mu*((pArr(i+1, j, k, inU)+pArr(i+1, j-1, k, inU)-
                                        pArr(i-1, j, k, inU)-pArr(i-1, j-1, k, inU))/(4.0*dx) + (pArr(i, j, k, inV)-pArr(i, j-1, k, inV))/dy);
#if (AMREX_SPACEDIM == 3)		
            tao_yy_Yneg = tao_yy_Yneg -2.0/3.0*mu*(pArr(i, j, k+1, inW)+pArr(i, j-1, k+1, inW)-
                                       pArr(i, j, k-1, inW)-pArr(i, j-1, k-1, inW))/(4.0*dz) ;
#endif	
            // tao_zz
			
#if (AMREX_SPACEDIM == 3)										
            Set::Scalar tao_zz_Zpos = 2.0*mu*((pArr(i, j, k+1, inW)-pArr(i, j, k, inW))/dz) 
                                         - 2.0/3.0*mu*((pArr(i+1, j, k, inU)+pArr(i+1, j, k+1, inU)-
                                         pArr(i-1, j, k, inU)-pArr(i-1, j, k+1, inU))/(4.0*dx) + 
                                         (pArr(i, j+1, k, inV)+pArr(i, j+1, k+1, inV)-
                                         pArr(i, j-1, k, inV)-pArr(i, j-1, k+1, inV))/(4.0*dy) +
                                         (pArr(i, j, k+1, inW)-pArr(i, j, k, inW))/dz);
										
            Set::Scalar tao_zz_Zneg = 2.0*mu*((pArr(i, j, k, inW)-pArr(i, j, k-1, inW))/dz) 
                                         - 2.0/3.0*mu*((pArr(i+1, j, k, inU)+pArr(i+1, j, k-1, inU)-
                                         pArr(i-1, j, k, inU)-pArr(i-1, j, k-1, inU))/(4.0*dx) + 
                                         (pArr(i, j+1, k, inV)+pArr(i, j+1, k-1, inV)-
                                         pArr(i, j-1, k, inV)-pArr(i, j-1, k-1, inV))/(4.0*dy) +
                                         (pArr(i, j, k, inW)-pArr(i, j, k-1, inW))/dz);
#endif
		    
            // Finding off-diagonal terms one after another:

            // tao_xy
            Set::Scalar tao_xy_Xpos = mu*((pArr(i, j+1, k, inU)+pArr(i+1, j+1, k, inU)-
                                        pArr(i, j-1, k, inU)-pArr(i+1, j-1, k, inU))/(4.0*dy) +
                                        (pArr(i+1, j, k, inV)-pArr(i, j, k, inV))/dx);
										
	    Set::Scalar tao_xy_Xneg = mu*((pArr(i, j+1, k, inU)+pArr(i-1, j+1, k, inU)-
                                        pArr(i, j-1, k, inU)-pArr(i-1, j-1, k, inU))/(4.0*dy) +
                                        (pArr(i, j, k, inV)-pArr(i-1, j, k, inV))/dx);
										
            Set::Scalar tao_xy_Ypos = mu*((pArr(i+1, j, k, inV)+pArr(i+1, j+1, k, inV)-
                                        pArr(i-1, j, k, inV)-pArr(i-1, j+1, k, inV))/(4.0*dx) +
                                        (pArr(i, j+1, k, inU)-pArr(i, j, k, inU))/dy);
            
            Set::Scalar tao_xy_Yneg = mu*((pArr(i+1, j, k, inV)+pArr(i+1, j-1, k, inV)-
                                         pArr(i-1, j, k, inV)-pArr(i-1, j-1, k, inV))/(4.0*dx) +
                                         (pArr(i, j, k, inU)-pArr(i, j-1, k, inU))/dy);
										
	    // tao_yz
#if (AMREX_SPACEDIM == 3)				
            Set::Scalar tao_yz_Ypos = mu*((pArr(i, j, k+1, inV)+pArr(i, j+1, k+1, inV)-
										pArr(i, j, k-1, inV)-pArr(i, j+1, k-1, inV))/(4.0*dz) +
										(pArr(i, j+1, k, inW)-pArr(i, j, k, inW))/dy);	

            Set::Scalar tao_yz_Yneg = mu*((pArr(i, j, k+1, inV)+pArr(i, j-1, k+1, inV)-
										pArr(i, j, k-1, inV)-pArr(i, j-1, k-1, inV))/(4.0*dz) +
										(pArr(i, j, k, inW)-pArr(i, j-1, k, inW))/dy);

            Set::Scalar tao_yz_Zpos = mu*((pArr(i, j+1, k, inW)+pArr(i, j+1, k+1, inW)-
										pArr(i, j-1, k, inW)-pArr(i, j-1, k+1, inW))/(4.0*dy) +
										(pArr(i, j, k+1, inV)-pArr(i, j, k, inV))/dz);
										
            Set::Scalar tao_yz_Zneg = mu*((pArr(i, j+1, k, inW)+pArr(i, j+1, k-1, inW)-
										pArr(i, j-1, k, inW)-pArr(i, j-1, k-1, inW))/(4.0*dy) +
										(pArr(i, j, k, inV)-pArr(i, j, k-1, inV))/dz);	   
#endif						

            // tao_xz
#if (AMREX_SPACEDIM == 3)				
            Set::Scalar tao_xz_Xpos = mu*((pArr(i, j, k+1, inU)+pArr(i+1, j, k+1, inU)-
										pArr(i, j, k-1, inU)-pArr(i+1, j, k-1, inU))/(4.0*dz) +
										(pArr(i+1, j, k, inW)-pArr(i, j, k, inW))/dx);

            Set::Scalar tao_xz_Xneg = mu*((pArr(i, j, k+1, inU)+pArr(i-1, j, k+1, inU)-
										pArr(i, j, k-1, inU)-pArr(i-1, j, k-1, inU))/(4.0*dz) +
										(pArr(i, j, k, inW)-pArr(i-1, j, k, inW))/dx);	

            Set::Scalar tao_xz_Zpos = mu*((pArr(i+1, j, k, inW)+pArr(i+1, j, k+1, inW)-
										pArr(i-1, j, k, inW)-pArr(i-1, j, k+1, inW))/(4.0*dx) +
										(pArr(i, j, k+1, inU)-pArr(i, j, k, inU))/dz);

            Set::Scalar tao_xz_Zneg = mu*((pArr(i+1, j, k, inW)+pArr(i+1, j, k-1, inW)-
										pArr(i-1, j, k, inW)-pArr(i-1, j, k-1, inW))/(4.0*dx) +
										(pArr(i, j, k, inU)-pArr(i, j, k-1, inU))/dz);										
#endif
			
            // Finding gradients and divergence of shear stresses using the above relations

            Set::Scalar tao_xx_x = (tao_xx_Xpos - tao_xx_Xneg)/dx;
            Set::Scalar tao_yy_y = (tao_yy_Ypos - tao_yy_Yneg)/dy;
            Set::Scalar tao_zz_z = 0.0;
#if (AMREX_SPACEDIM == 3)			
            tao_zz_z = 	(tao_zz_Zpos - tao_zz_Zneg)/dz;	
#endif
            // -----
            Set::Scalar tao_xy_x = (tao_xy_Xpos - tao_xy_Xneg)/dx;	
            Set::Scalar tao_xy_y = (tao_xy_Ypos - tao_xy_Yneg)/dy;
			// -----
            Set::Scalar tao_yz_y = 0.0;
            Set::Scalar tao_yz_z = 0.0;
#if (AMREX_SPACEDIM == 3)			
            tao_yz_y =  (tao_yz_Ypos - tao_yz_Yneg)/dy;
            tao_yz_z =  (tao_yz_Zpos - tao_yz_Zneg)/dz;
#endif
            // -----
            Set::Scalar tao_xz_x = 0.0;
            Set::Scalar tao_xz_z = 0.0;
#if (AMREX_SPACEDIM == 3)			
            tao_xz_x =  (tao_xz_Xpos - tao_xz_Xneg)/dx;
            tao_xz_z =  (tao_xz_Zpos - tao_xz_Zneg)/dz;
#endif
            
            // Assigning values to source term multiFab
			
            sourceTermArr(i,j,k,0) = 0.0;
			
            sourceTermArr(i,j,k,1) = tao_xx_x + tao_xy_y + tao_xz_z;

            sourceTermArr(i,j,k,2) = tao_xy_x + tao_yy_y + tao_yz_z;
		
            sourceTermArr(i,j,k,3) = tao_xz_x + tao_yz_y + tao_zz_z;

            // The final energy term of the source term vector

            // Divide the total energy term into four parts, namely, S_E1, S_E2, S_E3 and S_E4 (S_E4 coming from the thermal contribution)
			
            Set::Scalar S_E1 = (((pArr(i, j, k, inU)+pArr(i+1, j, k, inU))/2.0*tao_xx_Xpos +
                                            (pArr(i, j, k, inV)+pArr(i+1, j, k, inV))/2.0*tao_xy_Xpos)-
                                            ((pArr(i, j, k, inU)+pArr(i-1, j, k, inU))/2.0*tao_xx_Xneg +
                                            (pArr(i, j, k, inV)+pArr(i-1, j, k, inV))/2.0*tao_xy_Xneg))/dx;
#if (AMREX_SPACEDIM == 3)			
            S_E1 = S_E1 + ((pArr(i, j, k, inW)+pArr(i+1, j, k, inW))/2.0*tao_xz_Xpos - 
		                    (pArr(i, j, k, inW)+pArr(i-1, j, k, inW))/2.0*tao_xz_Xneg)/dx;
#endif

            // -----
            
            Set::Scalar S_E2 = (((pArr(i, j, k, inU)+pArr(i, j+1, k, inU))/2.0*tao_xy_Ypos +
			                     (pArr(i, j, k, inV)+pArr(i, j+1, k, inV))/2.0*tao_yy_Ypos)-
			                    ((pArr(i, j, k, inU)+pArr(i, j-1, k, inU))/2.0*tao_xy_Yneg +
			                     (pArr(i, j, k, inV)+pArr(i, j-1, k, inV))/2.0*tao_yy_Yneg))/dy;
#if (AMREX_SPACEDIM == 3)			
            S_E2 = S_E2 + ((pArr(i, j, k, inW)+pArr(i, j+1, k, inW))/2.0*tao_yz_Ypos - 
		                    (pArr(i, j, k, inW)+pArr(i, j-1, k, inW))/2.0*tao_yz_Yneg)/dy;
#endif			

            // -----
			
            Set::Scalar S_E3 = 0.0;
#if (AMREX_SPACEDIM == 3)			
            S_E3 = (((pArr(i, j, k, inU)+pArr(i, j, k+1, inU))/2.0*tao_xz_Zpos +
                                 (pArr(i, j, k, inV)+pArr(i, j, k+1, inV))/2.0*tao_yz_Zpos +
                                 (pArr(i, j, k, inW)+pArr(i, j, k+1, inW))/2.0*tao_zz_Zpos)-
			         ((pArr(i, j, k, inU)+pArr(i, j, k-1, inU))/2.0*tao_xz_Zneg +
			         (pArr(i, j, k, inV)+pArr(i, j, k-1, inV))/2.0*tao_yz_Zneg +
                                 (pArr(i, j, k, inW)+pArr(i, j, k-1, inW))/2.0*tao_zz_Zneg))/dz;
#endif

            // -----
			
            Set::Scalar S_E4 = 0.0;  // Thermal Diffusion Contribution hasn't been added yet
			
            sourceTermArr(i,j,k,4) = S_E1 + S_E2 + S_E3 + S_E4;

	
            }); // index loop end
		
      } // MFIter loop end
	
} // ComputeSourceTerm function end

}; // end of class

} // end of Namespace

#endif	
