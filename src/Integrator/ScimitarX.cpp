#include "IO/ParmParse.H"
#include "Integrator/ScimitarX.H"
#include "BC/BC.H"
#include "BC/Nothing.H"
#include "BC/Constant.H"
#include "IC/SodShock.H"
#include "Model/Fluid/Fluid.H"
#include "Numeric/Stencil.H"
#include "Numeric/TimeStepper.H"
#include "Numeric/FluxHandler.H"
#include "Util/Util.H"
#include "Util/ScimitarX_Util.H"

namespace Integrator
{

// Define the static member variable
Util::ScimitarX_Util::getVariableIndex ScimitarX::variableIndex;


ScimitarX::ScimitarX(IO::ParmParse& pp):ScimitarX()
{

    fluxHandler = new Numeric::FluxHandler<ScimitarX>();
    timeStepper = new Numeric::TimeStepper<ScimitarX>();

    pp.queryclass(*this);
}

void
ScimitarX::Parse(ScimitarX& value, IO::ParmParse& pp)
{
    BL_PROFILE("Integrator::ScimitarX::Parse()");
    {
        ScimitarX::SolverTypeManager& setIndex = ScimitarX::SolverTypeManager::getInstance();

        std::string solverTypeStr;

        if (pp.query("SolverType", solverTypeStr)) 
        {
            auto it = ScimitarX::stringToSolverType.find(solverTypeStr);
            if (it != ScimitarX::stringToSolverType.end()) {
                value.solverType = it->second;

                ScimitarX::variableIndex = setIndex.computeAndAssignVariableIndices(value.solverType);
                value.number_of_components = ScimitarX::variableIndex.NVAR_MAX;
               
                std::cout << "DEBUG: ScimitarX::variableIndex.NVAR_MAX = " 
                << ScimitarX::variableIndex.NVAR_MAX << std::endl;

                std::cout << "DEBUG: Variable indices:" << std::endl;
                for (const auto& [variable, index] : ScimitarX::variableIndex.variableIndexMap) {
                std::cout << "  Variable: " << variable << ", Index: " << index << std::endl;
                }

                // Call the setupBoundaryConditions function
               /* IO::ParmParse bc_pp = setupPVecBoundaryConditions(pp, ScimitarX::variableIndex);

                std::string result;
                if (bc_pp.query("bc.pvec.type.xlo", result)) {
                std::cout << "DEBUG: Queried type.xlo from bc_pp: " << result << std::endl;
                } else {
                std::cerr << "ERROR: Missing bc.pvec.type.xlo in bc_pp" << std::endl;
                }

                if (bc_pp.query("bc.pvec.val.xlo", result)) {
                std::cout << "DEBUG: Queried val.xlo from bc_pp: " << result << std::endl;
                } else {
                std::cerr << "ERROR: Missing bc.pvec.val.xlo in bc_pp" << std::endl;
                } */  

                value.bc_PVec = new BC::Constant(value.number_of_components, pp, "bc.pvec");
                value.bc_Pressure = new BC::Constant(1, pp, "bc.pressure");
            } else {
                Util::Abort(__FILE__, __func__, __LINE__, "Invalid SolverType: " + solverTypeStr);
            }
        }
    }
    
    // Register New Fabs
    {
        value.RegisterNewFab(value.QVec_mf, &value.bc_nothing, value.number_of_components, value.number_of_ghost_cells, "QVec", false);
        value.RegisterNewFab(value.QVec_old_mf, &value.bc_nothing, value.number_of_components, value.number_of_ghost_cells, "QVec_old", false);

        value.RegisterFaceFab<0>(value.XFlux_mf, &value.bc_nothing, value.number_of_components, value.number_of_ghost_cells, "xflux", false);
#if AMREX_SPACEDIM >= 2
        value.RegisterFaceFab<1>(value.YFlux_mf, &value.bc_nothing, value.number_of_components, value.number_of_ghost_cells, "yflux", false);
#endif
#if AMREX_SPACEDIM == 3
        value.RegisterFaceFab<2>(value.ZFlux_mf, &value.bc_nothing, value.number_of_components, value.number_of_ghost_cells, "zflux", false);
#endif
//        value.RegisterNewFab(value.Source_mf, &value.bc_nothing, value.number_of_components, value.number_of_ghost_cells, "SourceVec", false);
        value.RegisterNewFab(value.PVec_mf, value.bc_PVec, value.number_of_components, value.number_of_ghost_cells, "PrimitiveVec", true); 
        value.RegisterNewFab(value.Pressure_mf, value.bc_Pressure, 1, value.number_of_ghost_cells, "Pressure", true);

    }
    // Initial Conditions
    {
        std::string type = "constant";
        pp.query("ic.pvec.type", type);

        if (type == "sodshock") {
            value.ic_PVec = new IC::SodShock(value.geom, pp, "ic.pvec.sodshock", ScimitarX::variableIndex);
        } else {
            Util::Abort(__FILE__, __func__, __LINE__, "Invalid ic.pvec.type: " + type);
        }

        pp.query("ic.pressure.type", type);
        if (type == "sodshock") {
            value.ic_Pressure = new IC::SodShock(value.geom, pp, "ic.pressure.sodshock");
        } else {
            Util::Abort(__FILE__, __func__, __LINE__, "Invalid ic.pressure.type: " + type);
        }
    } 

    { 
    // Access the maps from the FeatureMaps singleton
    auto& fluxReconstructionMap = getFeatureMaps().getFluxReconstructionMap();
    auto& fluxSchemeMap = getFeatureMaps().getFluxSchemeMap();
    auto& timeSteppingSchemeMap = getFeatureMaps().getTimeSteppingSchemeMap();

    // Flux Reconstruction parsing
    std::string fluxReconstructionStr;
    if (pp.query("FluxReconstruction", fluxReconstructionStr)) {
        auto it = fluxReconstructionMap.find(fluxReconstructionStr);
        if (it != fluxReconstructionMap.end()) {
            value.reconstruction_method = it->second;
        } else {
            Util::Abort(__FILE__, __func__, __LINE__, "Invalid FluxReconstruction value: " + fluxReconstructionStr);
        }
    }

    // Flux Scheme parsing
    std::string fluxSchemeStr;
    if (pp.query("FluxScheme", fluxSchemeStr)) {
        auto it = fluxSchemeMap.find(fluxSchemeStr);
        if (it != fluxSchemeMap.end()) {
            value.flux_scheme = it->second;
        } else {
            Util::Abort(__FILE__, __func__, __LINE__, "Invalid FluxScheme value: " + fluxSchemeStr);
        }
    }

    // Time-Stepping Scheme parsing
    std::string timeSteppingSchemeStr;
    if (pp.query("TimeSteppingScheme", timeSteppingSchemeStr)) {
        auto it = timeSteppingSchemeMap.find(timeSteppingSchemeStr);
        if (it != timeSteppingSchemeMap.end()) {
            value.temporal_scheme = it->second;
        } else {
            Util::Abort(__FILE__, __func__, __LINE__, "Invalid TimeSteppingScheme value: " + timeSteppingSchemeStr);
        }
    }

    }
    // Read CFL number and initial time step
    pp.query_required("cflNumber", value.cflNumber);
    
}

void ScimitarX::Initialize(int lev)
{
    ic_PVec->Initialize(lev, PVec_mf);
    ic_Pressure->Initialize(lev, Pressure_mf);

    ScimitarX::ComputeConservedVariables<SolverType::SolveCompressibleEuler>(lev);
    std::swap(*QVec_old_mf[lev], *QVec_mf[lev]); 
}


void ScimitarX::TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags, Set::Scalar /*time*/, int /*ngrow*/)
{
    const Set::Scalar* DX = geom[lev].CellSize();
    Set::Scalar dr = sqrt(AMREX_D_TERM(DX[0] * DX[0], +DX[1] * DX[1], +DX[2] * DX[2]));
    Set::Scalar refinement_threshold = 10.0; // Set refinement threshold to 10

    //for (amrex::MFIter mfi(*Pressure_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    for (amrex::MFIter mfi(*Pressure_mf[lev], false); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<char> const& tags = a_tags.array(mfi);
        amrex::Array4<Set::Scalar> const& pressure = (*Pressure_mf[lev]).array(mfi);

        amrex::ParallelFor(bx, [=](int i, int j, int k) {
            Set::Vector grad = Numeric::Gradient(pressure, i, j, k, 0, DX);
            Set::Scalar grad_magnitude = grad.lpNorm<2>();

            if (grad_magnitude * dr > refinement_threshold) {
                tags(i, j, k) = amrex::TagBox::SET;
            }
        });
    }

        Util::Message(INFO, "Refinement threshold set to", refinement_threshold);
}


void ScimitarX::TimeStepBegin(Set::Scalar time, int lev) {

    Set::Scalar Coarsest_dt = 0.0;
    Coarsest_dt = ComputeAndSetNewTimeStep();  // Compute dt based on global `minDt`
    Util::Message(INFO, "Starting time step at time: ", Coarsest_dt);
}


void ScimitarX::TimeStepComplete(Set::Scalar time, int lev) {
    
}

void ScimitarX::Regrid(int lev, Set::Scalar time) {

}

void ScimitarX::ApplyBoundaryConditions(int lev, Set::Scalar time) {
        
Util::Message(INFO, "Number of Components " + std::to_string(number_of_components));
        
       // ApplyPatch(lev, time, QVec_mf, *QVec_mf[lev], *bc_PVec, 0);    
       // ApplyPatch(lev, time, QVec_old_mf, *QVec_old_mf[lev], *bc_PVec, 0);    
        Integrator::ApplyPatch(lev, time, PVec_mf, *PVec_mf[lev], *bc_PVec, 0);        
        Integrator:: ApplyPatch(lev, time, Pressure_mf, *Pressure_mf[lev], *bc_Pressure, 0); 
        Integrator::ApplyPatch(lev, time, QVec_mf, *QVec_mf[lev], bc_nothing, 0);        
        Integrator::ApplyPatch(lev, time, XFlux_mf, *XFlux_mf[lev],bc_nothing, 0);        
        Integrator::ApplyPatch(lev, time, YFlux_mf, *YFlux_mf[lev], bc_nothing, 0);
#if AMREX_SPACEDIM == 3        
        Integrator::ApplyPatch(lev, time, ZFlux_mf, *ZFlux_mf[lev], bc_nothing, 0);        
#endif
        // Fill ghost cells before computing fluxes
        //bc_PVec->FillBoundary(*PVec_mf[lev], 0, number_of_components, time); 
        //bc_Pressure->FillBoundary(*Pressure_mf[lev], 0, number_of_components, time); 

     // PVec_mf[lev]->FillBoundary();
     // Pressure_mf[lev]->FillBoundary();   
}

Set::Scalar ScimitarX::ComputeAndSetNewTimeStep() {
    // Compute the minimum time step over the entire domain using GetTimeStep
    Set::Scalar coarsest_dt = GetTimeStep();  // GetTimeStep already accounts for the CFL number

    // Start with the finest level time step
    //Set::Scalar coarsest_dt = finest_dt;

    // Adjust the time step for the coarsest level by multiplying back the refinement ratios
    /*for (int lev = finest_level; lev > 0; --lev) {
        // `refRatio(lev - 1)` returns an IntVect. Assume isotropic refinement for simplicity.
        int refinement_factor = refRatio(lev - 1)[0];  // Assuming refinement is isotropic (same value in all directions)

        coarsest_dt *= refinement_factor;  // Scale the time step conservatively for refinement
    }*/

    // Set the coarsest-level time step for all levels
    Integrator::SetTimestep(coarsest_dt);

    Util::Message(INFO, "New coarsest-level timestep set: ",coarsest_dt);
    return coarsest_dt;
}

// Function to compute the time step size based on CFL condition
Set::Scalar ScimitarX::GetTimeStep() {
    Set::Scalar minDt = std::numeric_limits<Set::Scalar>::max();  // Start with a large value       

    for (int lev = 0; lev <= maxLevel(); ++lev) {  // Use maxLevel() from the base class
        const Set::Scalar* dx = geom[lev].CellSize();  // Access the geometry at level `lev`

        //for (amrex::MFIter mfi(*PVec_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        for (amrex::MFIter mfi(*PVec_mf[lev], false); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.tilebox();  // Iterate over tiles in the multifab
            auto const& pArr = PVec_mf.Patch(lev, mfi);
            auto const& pressure = Pressure_mf.Patch(lev, mfi);

            Set::Scalar minDt_local = std::numeric_limits<Set::Scalar>::max();  // Thread-local minDt

            amrex::ParallelFor(bx, [=, &minDt_local](int i, int j, int k) noexcept {
                Set::Scalar rho = pArr(i, j, k, variableIndex.DENS);
                Set::Scalar u = pArr(i, j, k, variableIndex.UVEL);
                Set::Scalar v = pArr(i, j, k, variableIndex.VVEL);
#if (AMREX_SPACEDIM == 3)
                Set::Scalar w = pArr(i, j, k, variableIndex.WVEL);
#endif
                Set::Scalar ie = pArr(i, j, k, variableIndex.IE);
                Set::Scalar gamma = 1.4;
                Set::Scalar p = std::max(pressure(i, j, k),1e-6);
                Set::Scalar c = Model::Fluid::Fluid().ComputeWaveSpeed(rho, p, gamma);

                // Compute the maximum characteristic speed
                Set::Scalar maxSpeed = std::abs(u) + c;
#if (AMREX_SPACEDIM >= 2)
                maxSpeed = std::max(maxSpeed, std::abs(v) + c);
#endif
#if (AMREX_SPACEDIM == 3)
                maxSpeed = std::max(maxSpeed, std::abs(w) + c);
#endif

                // Compute local timestep for this cell
                Set::Scalar dtLocal = dx[0] / maxSpeed;
#if (AMREX_SPACEDIM >= 2)
                dtLocal = std::min(dtLocal, dx[1] / maxSpeed);
#endif
#if (AMREX_SPACEDIM == 3)
                dtLocal = std::min(dtLocal, dx[2] / maxSpeed);
#endif

                // Track the local minimum
                if (dtLocal < minDt_local) {
                    minDt_local = dtLocal;
                }
            });

            // Update the global minDt
            minDt = std::min(minDt, minDt_local);
        }
    }

    // Reduce across processes to find the global minimum timestep
    amrex::ParallelDescriptor::ReduceRealMin(minDt);
    return ScimitarX::cflNumber * minDt;  // Return CFL-adjusted time step for the finest level
}


IO::ParmParse ScimitarX::setupPVecBoundaryConditions(IO::ParmParse& pp, const Util::ScimitarX_Util::getVariableIndex& variableIndex)
{
    int n_components = variableIndex.NVAR_MAX;
    std::vector<std::string> component_names(n_components);

    for (const auto& [variable, index] : variableIndex.variableIndexMap) {
        component_names[index] = std::to_string(variable);
    }

    IO::ParmParse bc_pp;

    // Define sides based on AMREX_SPACEDIM
#if AMREX_SPACEDIM == 1
    const std::vector<std::string> sides = {"xlo", "xhi"};
#elif AMREX_SPACEDIM == 2
    const std::vector<std::string> sides = {"xlo", "xhi", "ylo", "yhi"};
#elif AMREX_SPACEDIM == 3
    const std::vector<std::string> sides = {"xlo", "xhi", "ylo", "yhi", "zlo", "zhi"};
#endif

    for (const std::string& side : sides) {
        for (int i = 0; i < n_components; ++i) {
            std::string type_key = "bc.pvec." + component_names[i] + ".type." + side;
            std::string val_key = "bc.pvec." + component_names[i] + ".val." + side;

            std::string type_value = "neumann"; // Default
            std::string value_str = "0.0";     // Default

            pp.query_default(type_key.c_str(), type_value, "neumann");
            pp.query_default(val_key.c_str(), value_str, "0.0");

            // Add individual entries directly to bc_pp
            std::string bc_type_key = "bc.pvec.type." + side;
            std::string bc_val_key = "bc.pvec.val." + side;

            bc_pp.addarr(bc_type_key.c_str(), {type_value});
            bc_pp.addarr(bc_val_key.c_str(), {value_str});
        }
    }

    return bc_pp;
}


} // namespace Integrator

