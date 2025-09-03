// Standard library includes
#include <iostream>
#include <fstream>
#include <iomanip>

// Alamo Includes
#include "IO/ParmParse.H"
#include "NarrowBandLevelset.H"

// IC
#include "IC/Constant.H"
#include "IC/Expression.H"
#include "IC/LS/Sphere.H"
#include "IC/LS/Zalesak.H"
#include "IC/LS/Expression.H"

// BC
#include "BC/Constant.H"
#include "BC/Expression.H"

// Numeric 
#include "Numeric/NarrowBandFluxHandler.H"
#include "Numeric/TimeStepper.H"

namespace Integrator
{

// Define constructor functions
// Empty constructor
NarrowBandLevelset::NarrowBandLevelset(int a_nghost) 
    : Integrator(), 
    number_of_ghost_cells(a_nghost) 
{}
    
// Constructor that triggers Parse
NarrowBandLevelset::NarrowBandLevelset(IO::ParmParse& pp) : NarrowBandLevelset() // Call default constructor
{
    fluxHandler = std::make_shared<Numeric::FluxHandler<NarrowBandLevelset>>();
    timeStepper = std::make_shared<Numeric::TimeStepper<NarrowBandLevelset>>();

    pp.queryclass(*this); // Call the static Parse function
}

// Define Parse function
void NarrowBandLevelset::Parse(NarrowBandLevelset& value, IO::ParmParse& pp){
    // Read CFL number
    pp.query_required("cflNumber", value.cflNumber);

    // Define levelset variables from input
    pp.query_default("levelset.number_of_levelsets", value.number_of_levelsets, 1);

    {// Resize the LSData structure to number of amr levels
    const int num_amr_levels = value.max_level + 1;
    value.LSData.objects.resize(num_amr_levels);
    value.LSData.levelset_object_map.resize(num_amr_levels);
    }

    // Get number of ghosts
    value.number_of_ghost_cells = value.LSData.GhostCells();

    // Get Reinitalization Order
    pp.query_validate("levelset.ReinitOrd", value.ReinitOrder, {1, 2});

    {// Define initial conditions for levelset field
    // Define the initial condition for the levelset
    pp.select_default<IC::LS::Sphere,IC::LS::Zalesak,IC::LS::Expression>("ls.ic",value.ic_ls,value.geom);

    // Define levelset bc
    pp.select_default<BC::Constant,BC::Expression>("ls.bc",value.bc_ls,value.number_of_levelsets);

    // Define levelset multifab objects - only old and new domain levelset fields
    value.RegisterNewFab(value.LS_mf, value.bc_ls, value.number_of_levelsets, value.number_of_ghost_cells, "LS", true);
    value.RegisterNewFab(value.LSold_mf, value.bc_ls, value.number_of_levelsets, value.number_of_ghost_cells, "LSold", true);
    }

    {// Store narrowband Tube properties in mask multifab
    value.RegisterNewFab(value.Tube_mf, value.bc_ls, value.number_of_levelsets, value.number_of_ghost_cells, "Tube", true);
    value.RegisterNewFab(value.Zero_mf, &value.bc_nothing, value.number_of_levelsets, value.number_of_ghost_cells, "Zerols", true);
    value.RegisterNewFab(value.CPT_mf, value.bc_ls, value.number_of_levelsets, value.number_of_ghost_cells, "CPT", true);
    value.RegisterNewFab(value.SignLS_mf, &value.bc_nothing, value.number_of_levelsets, value.number_of_ghost_cells, "SignPhi", true);
    }

    {// Define Geometry MultiFabs
    value.RegisterNewFab(value.NormalX_mf, value.bc_ls, value.number_of_levelsets, value.number_of_ghost_cells, "NormalX", true);
#if AMREX_SPACEDIM >= 2
    value.RegisterNewFab(value.NormalY_mf, value.bc_ls, value.number_of_levelsets, value.number_of_ghost_cells, "NormalY", true);
#endif
#if AMREX_SPACEDIM == 3
    value.RegisterNewFab(value.NormalZ_mf, value.bc_ls, value.number_of_levelsets, value.number_of_ghost_cells, "NormalZ", true);
#endif
    value.RegisterNewFab(value.Curvature_mf, &value.bc_nothing, value.number_of_levelsets, value.number_of_ghost_cells, "Curvature", true);   
    }
    
    {// Initialize velocity field. Will be removed once integrated with ScimitarX integrator
    pp.select_default<IC::Constant,IC::Expression>("velocity.ic",value.ic_velocity,value.geom); // velocity initial condtion
    pp.select_default<BC::Constant,BC::Expression>("velocity.bc",value.bc_velocity,AMREX_SPACEDIM); // velocty boundary condition
    value.RegisterNewFab(value.LSvel_mf, value.bc_velocity, AMREX_SPACEDIM, value.number_of_ghost_cells, "LSvel", true);
    
    // Select method to update velocity field based on validation method
    pp.query_validate("velocity.method",value.velocity_method,{"constant"});
    }

    {// Access the maps from the FeatureMaps singleton
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

    {// Set LSData to reflect FluxReconstruction, FluxMethod, and Temporal schemes
    // Reconstruction
    switch (value.reconstruction_method) {
        using FR = FluxReconstruction;
        case FR::FirstOrder:
            value.LSData.reconstruction_method = LevelSetData::FluxReconstruction::FirstOrder;
            value.fluxHandler->SetReconstruction(std::make_shared<Numeric::FirstOrderReconstruction<NarrowBandLevelset>>());
            break;

        case FR::ThirdOrderENO:
            Util::Abort(__FILE__, __func__, __LINE__, "ThirdOrderENO is not implemented yet.");
            break;

        case FR::FifthOrderWENO:
            value.LSData.reconstruction_method = LevelSetData::FluxReconstruction::FifthOrderWENO;
            value.fluxHandler->SetReconstruction(std::make_shared<Numeric::FifthOrderWENOReconstruction<NarrowBandLevelset>>());
            break;
    }

    // FluxMethod
    value.fluxHandler->SetFluxMethod(std::make_shared<Numeric::NarrowbandFluxComputation<NarrowBandLevelset>>());

    // Temporal
    if (value.temporal_scheme == TimeSteppingScheme::ForwardEuler) {
        value.timeStepper->SetTimeSteppingScheme(std::make_shared<Numeric::EulerForwardScheme<NarrowBandLevelset>>());
    } else {
        Util::Abort(__FILE__, __func__, __LINE__, "RK3 Timestepper not implemented yet.");
    }
    }
}

// Define required override functions
void NarrowBandLevelset::Initialize(int lev){
    // Initialize levelset data structure
    LSData.geom = geom[lev];
    LSData.flow_domain = geom[lev].Domain();
    LSData.ComputeTubeWidths();
    const int nghost = LSData.GhostCells();

    if (auto* expr = dynamic_cast<IC::LS::Expression*>(ic_ls)) {
        Tube_mf[lev]->setVal(NarrowBandTubeType::OUTSIDETUBE);
        expr->SetTubeField(&Tube_mf); 
        expr->SetZeroField(&Zero_mf);
        expr->SetCPTField(&CPT_mf);
        expr->SetLSoldField(&LSold_mf);          
        expr->SetLSData(&LSData);  
    }
    
    ic_ls->Initialize(lev, LS_mf);

    // Define velocity_mf containing velocity vector
    ic_velocity->Initialize(lev, LSvel_mf);

    // Initialize sign to -2
    SignLS_mf[lev]->setVal(-2, 0, number_of_levelsets, number_of_ghost_cells);

    // Finish initializing levelset map, geometry, and velocity information
    const int num_objects = LSData.num_objs;
    const int ngrow = LSData.GhostCells();

    // Ensure levelset_object_map has enough levels
    if (LSData.levelset_object_map.size() <= lev) {
        LSData.levelset_object_map.resize(lev + 1);
    }
    
    for (int obj = 0; obj < num_objects; obj++) {
        auto& object = LSData.objects[lev][obj];

        // Populate the ls_id → object index map
        LSData.levelset_object_map[lev][object.ls_id].push_back(obj);

        /*// Get geometries
        object.ComputeGeometryQuantities();
        
        // Copy for plotting
#if AMREX_SPACEDIM == 1
        object.CopyGeometryToFlowFields(*NormalX_mf[lev]);
#elif AMREX_SPACEDIM == 2
        object.CopyGeometryToFlowFields(*NormalX_mf[lev], *NormalY_mf[lev]);
#else
        object.CopyGeometryToFlowFields(*NormalX_mf[lev], *NormalY_mf[lev], *NormalZ_mf[lev]);
#endif*/

        // Get the ReinitInterval & Order
        object.ReinitInterval = object.FlowData->GetReinitInterval();
        object.ReinitOrder = ReinitOrder;
    }

    // Get the proper timestep
    ComputeAndSetNewTimeStep();
}

void NarrowBandLevelset::TimeStepBegin(Set::Scalar time, int lev) {

}

void NarrowBandLevelset::TimeStepComplete(Set::Scalar time, int lev) { 
    ComputeAndSetNewTimeStep();
}

void NarrowBandLevelset::ComputeAndSetNewTimeStep() {
    // Compute the minimum time step using the CFL condition
    Set::Scalar finest_dt = GetTimeStep();  // GetTimeStep() already applies the CFL number

    // Start with the finest-level time step
    Set::Scalar coarsest_dt = finest_dt;

    // Adjust time step for coarser levels based on refinement ratios
    for (int lev = finest_level; lev > 0; --lev) {
        int refinement_factor = refRatio(lev - 1)[0];  // Assume isotropic refinement
        coarsest_dt *= refinement_factor;  // Scale conservatively for refinement
    }

    // Set the time step for all levels
    Integrator::SetTimestep(coarsest_dt);
}

Set::Scalar NarrowBandLevelset::GetTimeStep() {
    Set::Scalar minDt = std::numeric_limits<Set::Scalar>::max();  // Start with a large value
    
    // Loop through all levelsets
    for (int ils=0; ils < 1; ils++){
        for (int lev = 0; lev <= maxLevel(); ++lev) {  // Iterate over AMR levels
            const Set::Scalar* dx = geom[lev].CellSize();  // Access the geometry at level `lev`
            const Set::Scalar min_DX = *std::min_element(dx, dx + AMREX_SPACEDIM);
    
            for (amrex::MFIter mfi(*LSvel_mf[lev], false); mfi.isValid(); ++mfi) {
                const amrex::Box& bx = mfi.tilebox();  // Iterate over tiles in the multifab
                auto const& velocity_arr = LSvel_mf.Patch(lev, mfi);  // Velocity field
    
                Set::Scalar minDt_local = std::numeric_limits<Set::Scalar>::max();  // Thread-local minDt
    
                amrex::ParallelFor(bx, [=, &minDt_local](int i, int j, int k) noexcept {
                    // Extract velocity components
                    Set::Scalar u = velocity_arr(i, j, k, 0);
    #if (AMREX_SPACEDIM >= 2)
                    Set::Scalar v = velocity_arr(i, j, k, 1);
    #endif
    #if (AMREX_SPACEDIM == 3)
                    Set::Scalar w = velocity_arr(i, j, k, 2);
    #endif
    
                    // Compute the maximum velocity magnitude
                    Set::Scalar maxSpeed = std::abs(u);
    #if (AMREX_SPACEDIM >= 2)
                    maxSpeed = std::max(maxSpeed, std::abs(v));
    #endif
    #if (AMREX_SPACEDIM == 3)
                    maxSpeed = std::max(maxSpeed, std::abs(w));
    #endif
    
                    // Compute local CFL time step restriction
                    if (maxSpeed > 1e-8) {  // Avoid division by zero
                        Set::Scalar dtLocal = dx[0] / maxSpeed;
    #if (AMREX_SPACEDIM >= 2)
                        dtLocal = std::min(dtLocal, dx[1] / maxSpeed);
    #endif
    #if (AMREX_SPACEDIM == 3)
                        dtLocal = std::min(dtLocal, dx[2] / maxSpeed);
    #endif
                        minDt_local = std::min(minDt_local, dtLocal);
                    }
                    else {
                        minDt_local = min_DX;
                    }
                });
    
                // Update the global minDt
                minDt = std::min(minDt, minDt_local);
            }
        }
    }

    // Reduce across MPI processes to find the global minimum timestep
    amrex::ParallelDescriptor::ReduceRealMin(minDt);
    
    return cflNumber * minDt;  // Return CFL-adjusted time step
}

void NarrowBandLevelset::Advance(int lev, Set::Scalar time, Set::Scalar dt) 
{
    // Increase timestep
    current_timestep ++;

    // Copy the old LS field to have LS values
    const int num_objs = LSData.num_objs;
    for (int obj=0; obj<num_objs; obj++) {
        auto& object = LSData.objects[lev][obj];
        object.CopyLSToLSOld();
    }

    // Advect
    Advect(lev, time, dt);

    SignLS_mf[lev]->setVal(-2.0);

    // Update Tube/geometries
    for (int obj=0; obj<num_objs; obj++) {
        auto& object = LSData.objects[lev][obj];

        // Reinitialize
        if (object.Reinitialize || current_timestep % object.ReinitInterval == 0) {
            if (object.ReinitOrder == 1) {
                object.ReinitializeLSFirstOrder(*LS_mf[lev], *SignLS_mf[lev]);
            } else {
                object.ReinitializeLSSecondOrder(*LS_mf[lev], *SignLS_mf[lev]);
            }
        }

        // Update the Narrowband cells
        object.UpdateNarrowbandTubeandMapping(*LS_mf[lev], *LSold_mf[lev], *Tube_mf[lev], *Zero_mf[lev], *CPT_mf[lev]);

        /*// Perform narrowband property computations
        object.ComputeGeometryQuantities();

        // Copy back for plotting
#if AMREX_SPACEDIM == 1
        object.CopyGeometryToFlowFields(*NormalX_mf[lev]);
#elif AMREX_SPACEDIM == 2
        object.CopyGeometryToFlowFields(*NormalX_mf[lev], *NormalY_mf[lev]);
#else
        object.CopyGeometryToFlowFields(*NormalX_mf[lev], *NormalY_mf[lev], *NormalZ_mf[lev]);
#endif*/
    }
}

void NarrowBandLevelset::Advect(int lev, Set::Scalar time, Set::Scalar dt){
    // Update the interface velocity
    UpdateInterfaceVelocity(lev);

    // Compute the flux
    switch (temporal_scheme) {
        case TimeSteppingScheme::ForwardEuler: {
            int numStages = timeStepper->GetNumberOfStages();
            // One-stage loop for Forward Euler
            for (int stage = 0; stage < numStages; ++stage) {

                // 2. Perform flux reconstruction and compute fluxes in all directions
                fluxHandler->ConstructFluxes(lev, this);

                // 3. Compute sub-step using the chosen time-stepping scheme
                timeStepper->ComputeSubStep(lev, dt, stage, this);
            }
            break;
        }
    } 

    // Ensure proper LS update and update CPT flags for reinitialization
    const int num_objs = LSData.num_objs;
    for (int obj=0; obj<num_objs; obj++) {
        auto& object = LSData.objects[lev][obj];
        object.CheckZeroCrossings(*LS_mf[lev], *LSold_mf[lev]);
        object.UpdateCPTFlagAfterAdvection();
    }
}

void NarrowBandLevelset::UpdateInterfaceVelocity(int lev){
    // Update flow domain velocity
    for (amrex::MFIter mfi(*LSvel_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto const& vel_arr = LSvel_mf.Patch(lev,mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            for (int d=0; d < AMREX_SPACEDIM; ++d){
                vel_arr(i,j,k,d) = vel_arr(i,j,k,d);
            }
        });
    }

    // Update object velocities
    const int num_objs = LSData.num_objs;

    for (int obj=0; obj < num_objs; obj++){
        auto& object = LSData.objects[lev][obj];

        object.LoadVelField(*LSvel_mf[lev]);
    }
}

void NarrowBandLevelset::TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow){
}

void NarrowBandLevelset::Regrid(int lev, Set::Scalar time) {
}

} // namespace Integrator
