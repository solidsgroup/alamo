#ifndef ELASTOPLASTICSOURCETERM_H
#define ELASTOPLASTICSOURCETERM_H

#include "AMReX.H"
#include "AMReX_REAL.H"
#include "AMReX_MultiFab.H"
#include "Set/Set.H"
#include "Numeric/Stencil.H"
#include "Util/Util.H"
#include "Util/ScimitarX_Util.H"
#include <memory>
#include <map>
#include <string>
#include <vector>

#include "Model/Fluid/ElastoPlastic.H"

#include "Integrator/ScimitarX.H"
#include "Numeric/NumericTypes.H"
#include "Numeric/SolverCapabilities.H"
#include "Numeric/IntegratorVariableAccessLayer.H"
#include "Numeric/SymmetryPreservingRoeAveragingOperations.H"
#include "Integrator/ScimitarX.H"
#include "Numeric/FluxHandler.H"
#include <cmath>
#include <algorithm>

//forward declaration for ScimitarX                                                                      
namespace Integrator {
    class ScimitarX;
}

namespace Numeric {


template <typename IntegratorDerivedClass>
class ElastoPlasticSourceTerm {
public:


void ComputeSourceTerm(int lev, IntegratorDerivedClass* solver) const {


    const Set::Scalar DX = solver->geom[lev].CellSize()[0];  // Access the geometry at level `lev`       
    const Set::Scalar DY = solver->geom[lev].CellSize()[1];                                              
#if (AMREX_SPACEDIM == 3)                                                                                
    const Set::Scalar DZ = solver->geom[lev].CellSize()[2];
#endif
    const int idx_DENS = solver->variableIndex.DENS; 
    const int idx_U = solver->variableIndex.UVEL; 
    const int idx_V = solver->variableIndex.VVEL; 
#if (AMREX_SPACEDIM == 3) 
    const int idx_W = solver->variableIndex.WVEL; 
#endif 
    const int idx_IE = solver->variableIndex.IE;
    const int idx_SXX = solver->variableIndex.SXX; 
    const int idx_SXY = solver->variableIndex.SXY; 
    const int idx_SYY = solver->variableIndex.SYY; 
#if AMREX_SPACEDIM == 3
    const int idx_SXZ = solver->variableIndex.SXZ;
    const int idx_SYZ = solver->variableIndex.SYZ;
    const int idx_SZZ = solver->variableIndex.SZZ;
#endif 
    const int idx_EPS = solver->variableIndex.EPSBAR; 
    const int idx_IEE = solver->variableIndex.IE_ELASTIC;


    for (amrex::MFIter mfi(*solver->SourceTermVec[lev], false); mfi.isValid(); ++mfi) {
	
	const amrex::Box& bx = mfi.validbox();

        auto const& p_arr = solver->PVec_mf.Patch(lev, mfi);
	auto const& sourceTermArr = solver->SourceTermVec.Patch(lev, mfi);

	amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

	    const Set::Scalar rho = p_arr(i,j,k, idx_DENS);
	    const Set::Scalar uvel = p_arr(i,j,k, idx_U);
	    const Set::Scalar vvel = p_arr(i,j,k, idx_V);
#if AMREX_SPACEDIM == 3
	    const Set::Scalar wvel = p_arr(i,j,k, idx_W);
#endif 
	    const Set::Scalar sxx = p_arr(i,j,k, idx_SXX); 
	    const Set::Scalar sxy = p_arr(i,j,k, idx_SXY); 
	    const Set::Scalar syy = p_arr(i,j,k, idx_SYY);
#if AMREX_SPACEDIM == 3
	    const Set::Scalar sxz = p_arr(i,j,k, idx_SXZ);
	    const Set::Scalar syz = p_arr(i,j,k, idx_SYZ);
	    const Set::Scalar szz = p_arr(i,j,k, idx_SZZ);
#endif
	    //compute the derivative using the stencil struct first order derivatives function?
	    Set::Scalar dudx = Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(p_arr, i, j, k, idx_U, solver->geom[lev].CellSize());
	    Set::Scalar dudy = Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(p_arr, i, j, k, idx_U, solver->geom[lev].CellSize());
	    Set::Scalar dvdx = Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(p_arr, i, j, k, idx_V, solver->geom[lev].CellSize());
	    Set::Scalar dvdy = Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(p_arr, i, j, k, idx_V, solver->geom[lev].CellSize());
	    
#if AMREX_SPACEDIM == 3
	    Set::Scalar dudz = Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(p_arr, i, j, k, idx_U, solver->geom[lev].CellSize());
	    Set::Scalar dvdz = Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(p_arr, i, j, k, idx_V, solver->geom[lev].CellSize());
	    Set::Scalar dwdz = Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(p_arr, i, j, k, idx_W, solver->geom[lev].CellSize());
	    Set::Scalar dwdx = Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(p_arr, i, j, k, idx_W, solver->geom[lev].CellSize());
	    Set::Scalar dwdy = Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(p_arr, i, j, k, idx_W, solver->geom[lev].CellSize());
#endif

	    Set::Scalar dSxxdx = Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(p_arr, i, j, k, idx_SXX, solver->geom[lev].CellSize());
	    Set::Scalar dSyydy = Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(p_arr, i, j, k, idx_SYY, solver->geom[lev].CellSize());
	    Set::Scalar dSxydx = Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(p_arr, i, j, k, idx_SXY, solver->geom[lev].CellSize());
	    Set::Scalar dSxydy = Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(p_arr, i, j, k, idx_SXY, solver->geom[lev].CellSize());
	    
#if AMREX_SPACEDIM == 3
	    Set::Scalar dSxzdz = Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(p_arr, i, j, k, idx_SXZ, solver->geom[lev].CellSize());
	    Set::Scalar dSyzdz = Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(p_arr, i, j, k, idx_SYZ, solver->geom[lev].CellSize());
	    Set::Scalar dSzzdz = Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(p_arr, i, j, k, idx_SZZ, solver->geom[lev].CellSize());
	    Set::Scalar dSxzdx = Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(p_arr, i, j, k, idx_SXZ, solver->geom[lev].CellSize());
	    Set::Scalar dSyzdy = Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(p_arr, i, j, k, idx_SYZ, solver->geom[lev].CellSize());
#endif
/*
	    Set::Scalar duSxxdx = uvel*dSxxdx + sxx*dudx;
	    Set::Scalar dvSxydx = vvel*dSxydx + sxy*dvdx;
	    Set::Scalar duSxydy = uvel*dSxydy + sxy*dudy;
	    Set::Scalar dvSyydy = vvel*dSyydy + syy*dvdy;
*/	    

	    Set::Scalar duSxxdx = (p_arr(i+1,j,k, idx_U)*p_arr(i+1,j,k, idx_SXX) - p_arr(i-1,j,k, idx_U)*p_arr(i-1,j,k, idx_SXX)) / (2.0 * DX);
	    Set::Scalar dvSxydx = (p_arr(i+1,j,k, idx_V)*p_arr(i+1,j,k, idx_SXY) - p_arr(i-1,j,k, idx_V)*p_arr(i-1,j,k, idx_SXY)) / (2.0 * DX);
	    Set::Scalar duSxydy = (p_arr(i,j+1,k, idx_U)*p_arr(i,j+1,k, idx_SXY) - p_arr(i,j-1,k, idx_U)*p_arr(i,j-1,k, idx_SXY)) / (2.0 * DY);
	    Set::Scalar dvSyydy = (p_arr(i,j+1,k, idx_V)*p_arr(i,j+1,k, idx_SYY) - p_arr(i,j-1,k, idx_V)*p_arr(i,j-1,k, idx_SYY)) / (2.0 * DY);
	    
#if AMREX_SPACEDIM == 3
	    Set::Scalar dwSxzdx = (p_arr(i+1,j,k, idx_W)*p_arr(i+1,j,k, idx_SXZ) - p_arr(i-1,j,k, idx_W)*p_arr(i-1,j,k, idx_SXZ)) / (2.0 * DX);
	    Set::Scalar dwSyzdy = (p_arr(i+1,j,k, idx_W)*p_arr(i+1,j,k, idx_SYZ) - p_arr(i-1,j,k, idx_W)*p_arr(i-1,j,k, idx_SYZ)) / (2.0 * DY);
	    Set::Scalar duSxzdz = (p_arr(i+1,j,k, idx_U)*p_arr(i+1,j,k, idx_SXZ) - p_arr(i-1,j,k, idx_U)*p_arr(i-1,j,k, idx_SXZ)) / (2.0 * DZ);
	    Set::Scalar dvSyzdz = (p_arr(i+1,j,k, idx_V)*p_arr(i+1,j,k, idx_SYZ) - p_arr(i-1,j,k, idx_V)*p_arr(i-1,j,k, idx_SYZ)) / (2.0 * DZ);
	    Set::Scalar dwSzzdz = (p_arr(i+1,j,k, idx_W)*p_arr(i+1,j,k, idx_SZZ) - p_arr(i-1,j,k, idx_W)*p_arr(i-1,j,k, idx_SZZ)) / (2.0 * DZ);
#endif

	    //Spin tensor
	    Set::Scalar omega_xy = 0.5*(dudy - dvdx);
	    Set::Scalar omega_yx = -1.0*omega_xy;
	    Set::Scalar sgn_omega_xy = std::copysign(1.0, omega_xy); 
	    omega_xy = sgn_omega_xy * std::floor(std::abs(omega_xy * 10.0e16)) / 10.0e16;
	    omega_yx = -1.0*omega_xy;
#if AMREX_SPACEDIM == 3
	    Set::Scalar omega_yz = 0.5*(dvdz - dwdy);
	    Set::Scalar omega_zy = -1.0*omega_yz;
	    Set::Scalar omega_xz = 0.5*(dudz - dwdx);
	    Set::Scalar omega_zx = -1.0*omega_xz;
#endif
	    //Deviatoric Strain rate tensor
	    Set::Scalar divergenceV = dudx + dvdy;
#if AMREX_SPACEDIM == 3
	    Set::Scalar divergenceV = dudx + dvdy + dwdz;
#endif
	    Set::Scalar trace_strain = divergenceV/3.0; 
	    Set::Scalar edotxx = dudx - trace_strain;
	    Set::Scalar edotxy = 0.5*(dudy + dvdx);
	    Set::Scalar edotyx = edotxy;
	    Set::Scalar edotyy = dvdy - trace_strain;
#if AMREX_SPACEDIM == 3
	    Set::Scalar edotxz = 0.5*(dudz + dwdx);
	    Set::Scalar edotyz = 0.5*(dvdz + dwdy);
	    Set::Scalar edotzz = -(edotxx + edotyy);
#endif
	    Set::Scalar sgn_edotxx = std::copysign(1.0, edotxx);
	    edotxx = sgn_edotxx * std::floor(std::abs(edotxx * 10.0e16)) / 10.0e16;

	    Set::Scalar sgn_edotxy = std::copysign(1.0, edotxy);
            edotxy = sgn_edotxy * std::floor(std::abs(edotxy * 10.0e16)) / 10.0e16;

	    Set::Scalar sgn_edotyy = std::copysign(1.0, edotyy);
            edotyy = sgn_edotyy * std::floor(std::abs(edotyy * 10.0e16)) / 10.0e16;

	    edotyx = edotxy;

	    //change to a vector array
	    Set::Scalar Source_continuity = 0.0;
	    Set::Scalar Source_Umom = dSxxdx + dSxydy;
	    Set::Scalar Source_Vmom = dSxydx + dSyydy;
#if AMREX_SPACEDIM == 3
	    Set::Scalar Source_Umom = dSxxdx + dSxydy + dSzzdz;
            Set::Scalar Source_Vmom = dSxydx + dSyydy + dSyzdz; 
	    Set::Scalar Source_Wmom = dSxzdx + dSyzdy + dSzzdz;
#endif
	    Set::Scalar Source_energy = duSxxdx + dvSxydx + duSxydy + dvSyydy;

#if AMREX_SPACEDIM == 3
	    Set::Scalar Source_energy = duSxxdx + dvSxydx + duSxydy + dvSyydy
                    + dwSxzdx + dwSyzdy + duSxzdz + dvSyzdz + dwSzzdz;
#endif

	    Set::Scalar shear_mod = 8.53e10; //4.5e10; //Copper

	    Set::Scalar Source_sxx = 2.0*rho*((omega_xy*sxy) + (shear_mod*edotxx));
	    Set::Scalar Source_sxy = rho*omega_xy*(syy - sxx) + 2.0*rho*shear_mod*edotxy;
	    Set::Scalar Source_syy = 2.0*rho*((omega_yx*sxy) + (shear_mod*edotyy));
#if AMREX_SPACEDIM == 3
	    Set::Scalar Source_sxx = 2.0*rho*(omega_xy*sxy + shear_mod*edotxx + omega_xz*sxz);
	    Set::Scalar Source_sxy = rho*omega_xy*(syy - sxx) + 2.0*rho*shear_mod*edotxy + rho*omega_xz*szy - omega_zy*sxz;
	    Set::Scalar Source_syy = 2.0*rho*(omega_yx*sxy + shear_mod*edotyy + omega_yz*syz);
	    Set::Scalar Source_sxz = rho*omega_xz*(szz - sxx) + 2.0*rho*shear_mod*edotxz + rho*omega_xy*syz - omega_yz*sxy //3d flag
	    Set::Scalar Source_syz = rho*omega_yz*(szz - syy) + 2.0*rho*shear_mod*edotyz + rho*omega_yx*sxz - omega_xz*sxy //3d flag
	    Set::Scalar Source_szz = 2.0*rho*(omega_zx*sxz + shear_mod*edotzz + omega_zy*syz);//3d flag
#endif
	    Set::Scalar Source_ie_elastic = 0.0; //for now 
	    Set::Scalar Source_epsbar = 0.0; //for now and no rk3 for epsbar
	    
            //source array

	    sourceTermArr(i,j,k,idx_DENS) = Source_continuity;
            sourceTermArr(i,j,k,idx_U) = Source_Umom;
	    sourceTermArr(i,j,k,idx_V) = Source_Vmom;
            sourceTermArr(i,j,k,idx_IE) = Source_energy;
            sourceTermArr(i,j,k,idx_SXX) = Source_sxx;
	    sourceTermArr(i,j,k,idx_SYY) = Source_syy;
	    sourceTermArr(i,j,k,idx_SXY) = Source_sxy;
            sourceTermArr(i,j,k,idx_EPS) = Source_epsbar;
	    sourceTermArr(i,j,k,idx_IEE) = Source_ie_elastic;
#if AMREX_SPACEDIM == 3
	    sourceTermArr(i,j,k,idx_W) = Source_Wmom;
	    sourceTermArr(i,j,k,idx_SXZ) = Source_sxz;
            sourceTermArr(i,j,k,idx_SYZ) = Source_syz;
	    sourceTermArr(i,j,k,idx_SZZ) = Source_szz;
#endif

	}); //amrex
      } //for loop
  } //ComputeSourceTerms



// Start Radial Return here
void RadialReturn(int lev, Set::Scalar dt, IntegratorDerivedClass* solver ) const {

    const int idx_DENS = solver->variableIndex.DENS; 
    const int idx_SXX = solver->variableIndex.SXX; 
    const int idx_SXY = solver->variableIndex.SXY; 
    const int idx_SYY = solver->variableIndex.SYY; 
#if AMREX_SPACEDIM == 3
    const int idx_SXZ = solver->variableIndex.SXZ;
    const int idx_SYZ = solver->variableIndex.SYZ;
    const int idx_SZZ = solver->variableIndex.SZZ;
#endif 
    const int idx_EPS = solver->variableIndex.EPSBAR; 


    for (amrex::MFIter mfi(*solver->QVec_mf[lev], false); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();

	auto const& q_arr = solver->QVec_mf.Patch(lev, mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            Set::Scalar rho = q_arr(i, j, k, idx_DENS);

            Set::Scalar sxx = q_arr(i, j, k, idx_SXX)/rho;
            Set::Scalar sxy = q_arr(i, j, k, idx_SXY)/rho;
            Set::Scalar syy = q_arr(i, j, k, idx_SYY)/rho;
#if AMREX_SPACEDIM == 3
            Set::Scalar sxz = q_arr(i, j, k, idx_SXZ)/rho;
            Set::Scalar syz = q_arr(i, j, k, idx_SYZ)/rho;
            Set::Scalar szz = q_arr(i, j, k, idx_SZZ)/rho;
#endif

            Set::Scalar epsbar = q_arr(i, j, k, idx_EPS);
            //add epsdot
            //add temperature
            Set::Scalar Shear_mod = 8.53e10; //4.5e10; //Copper
            Set::Scalar Cv0 = 1; //Needs its own function

            Set::Scalar Sigma_Y = 9.79e8; //9.0e7; //Set for constant now. Needs a function to account for other types.
            Set::Scalar DsigmaDeps = 0; //Set for constant now. Needs a function to account for other types.

            Set::Scalar sum = sxx*sxx + syy*syy + sxx*syy + sxy*sxy;
#if (AMREX_SPACEDIM == 3)
            Set::Scalar sum = sxx*sxx + syy*syy + sxx*syy + sxy*sxy + sxz*sxz + syz*syz;
#endif
            Set::Scalar Sqrt_trace = std::sqrt(2.0 * sum);
            Set::Scalar Consistency_f = Sqrt_trace - (Sigma_Y * std::sqrt(2.0 / 3.0));

	    if (Consistency_f < 1.0e-10) {
		    return;
	    }


	    int count = 0;
	    Set::Scalar sxx_trial = sxx;
            Set::Scalar syy_trial = syy;
            Set::Scalar sxy_trial = sxy;
#if (AMREX_SPACEDIM == 3)
            Set::Scalar sxz_trial = sxz;
            Set::Scalar syz_trial = syz;
	    Set::Scalar szz_trial = szz;
#endif
            Set::Scalar epsbar_previous_timestep = epsbar;

            Set::Scalar Numerator = Sqrt_trace - (Sigma_Y * std::sqrt(2.0 / 3.0));
            Set::Scalar Denominator = 2.0 * Shear_mod * (1.0 + (DsigmaDeps / (3.0 * Shear_mod)));
            Set::Scalar LambdaOld = Numerator/Denominator;
            Set::Scalar Residue = 1e6;
	    Set::Scalar Lambdanew = LambdaOld; //0;
/*
            while (std::abs(Residue) > 1.0e-15) {
                    //No Assumption of Linear Hardening ~ Linearized Consistency
                    Set::Scalar epsbar0 = epsbar_previous_timestep + (std::sqrt(2.0 / 3.0) * LambdaOld);
                    Set::Scalar epsdot0 = std::sqrt(2.0 / 3.0) * LambdaOld / dt;
                    Sigma_Y = 9.0e7;
                    DsigmaDeps = 0; //Slope Stress Strain Curve. Still set for constant yield strength
                    Set::Scalar DsigmaDepsdot = 0; //Rate Sensitivity Stress Strain Curve. Set for constant
                    sum = sxx_trial*sxx_trial + syy_trial*syy_trial + sxx_trial*syy_trial + sxy_trial*sxy_trial;
#if (AMREX_SPACEDIM == 3)
                    sum = sxx_trial*sxx_trial + syy_trial*syy_trial + sxx_trial*syy_trial
                            + sxy_trial*sxy_trial + sxz_trial*sxz_trial + syz_trial*syz_trial;
#endif
                    Sqrt_trace = std::sqrt(2.0 * sum);
                    Numerator = Sqrt_trace - (Sigma_Y * std::sqrt(2.0 / 3.0));
                    Set::Scalar EffectiveSlope = DsigmaDeps + DsigmaDepsdot / dt; //define delt
                    Denominator = 2.0 * (Shear_mod + (EffectiveSlope / 3.0));
                    Set::Scalar ObjFunc = LambdaOld - (Numerator / Denominator);
                    Set::Scalar DervtObjFunc = 1.0 + ((2.0 / 3.0) * (EffectiveSlope / Denominator));
                    Residue = -1.0 * ObjFunc / DervtObjFunc;
		    Lambdanew = LambdaOld + Residue;
                    LambdaOld = Lambdanew;
                    count = count + 1;

                    if (count > 200) {
                            Lambdanew = std::abs(Lambdanew);
                            Residue = 1e-20;
                    }
            }
*/
	    Set::Scalar Nxx = sxx_trial / Sqrt_trace;
            Set::Scalar Nyy = syy_trial / Sqrt_trace;
            Set::Scalar Nxy = sxy_trial / Sqrt_trace;
#if (AMREX_SPACEDIM == 3)
            Set::Scalar Nxz = sxz_trial / Sqrt_trace;
            Set::Scalar Nyz = syz_trial / Sqrt_trace;
	    Set::Scalar Nzz = szz_trial / Sqrt_trace;
#endif

            sxx = sxx_trial - (2.0 * Shear_mod * Lambdanew * Nxx);
            syy = syy_trial - (2.0 * Shear_mod * Lambdanew * Nyy);
            sxy = sxy_trial - (2.0 * Shear_mod * Lambdanew * Nxy);
#if (AMREX_SPACEDIM == 3)
            sxz = sxz_trial - (2.0 * Shear_mod * Lambdanew * Nxz);
            syz = syz_trial - (2.0 * Shear_mod * Lambdanew * Nyz);
	    szz = szz_trial - (2.0 * Shear_mod * Lambdanew * Nzz);
#endif
            Set::Scalar epsbar0 = epsbar_previous_timestep + std::sqrt(2 / 3) * Lambdanew;
            Set::Scalar epsdot0 = std::sqrt(2 / 3) * Lambdanew / dt;

            q_arr(i, j, k, idx_EPS) = epsbar0;

	    q_arr(i, j, k, idx_SXX) = rho * sxx;
	    q_arr(i, j, k, idx_SYY) = rho * syy;
	    q_arr(i, j, k, idx_SXY) = rho * sxy;
#if (AMREX_SPACEDIM == 3)
            q_arr(i, j, k, idx_SXZ) = rho * sxz;
            q_arr(i, j, k, idx_SYZ) = rho * syz;
	    q_arr(i, j, k, idx_SZZ) = rho * szz;
#endif

        });
    }

    solver->QVec_mf[lev]->FillBoundary();

	} //RR



}; // ElastoPlastic

}; // Numeric


#endif
	
