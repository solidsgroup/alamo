// Standard library includes
#include <iostream>
#include <fstream>
#include <iomanip>

// Alamo Includes
#include "IO/ParmParse.H"
#include "NarrowBandLevelset.H"

// BC
#include "BC/BC.H"
#include "BC/Constant.H"
#include "BC/Expression.H"

// IC
#include "IC/Constant.H"
#include "IC/LS/Sphere.H"
#include "IC/LS/Zalesak.H"

// Numerc
#include "Numeric/NarrowBandFluxHandler.H"
#include "Numeric/TimeStepper.H"
//#include "Util/NarrowBandLevelset_Util.H"

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
    {// Define initial and boundary conditions
    
    // Query the IC assuming either LS::Sphere or LS::Zalesak
    pp.select_default<IC::LS::Sphere,IC::LS::Zalesak>("ic.ls",value.ic_ls,value.geom);
    
    // Assume Neumann BC for levelset field
    value.bc_ls = new BC::Constant(value.number_of_components, pp, "bc.ls");
    }
    
    {// Define levelset multifab objects - only old and new domain levelset fields
    value.RegisterNewFab(value.ls_mf, value.bc_ls, value.number_of_components, value.number_of_ghost_cells, "LS", true);
    value.RegisterNewFab(value.ls_old_mf, value.bc_ls, value.number_of_components, value.number_of_ghost_cells, "LS_old", false);
    }

    // Define levelset data structure
    value.level_sets.resize(value.number_of_components);
    
    {// Initialize constant velocity. Will be removed once integrated with ScimitarX integrator
    // Define constant velocity vector
    value.ic_velocity = new IC::Constant(value.geom, pp, "ic.velocity");
    value.bc_velocity = new BC::Constant(AMREX_SPACEDIM, pp, "bc.velocity");
    //pp_queryarr("ic.velocity.value", value.constant_velocity);

    // Define velocity, normal, and curvature multifabs
    for (int ils=0; ils < value.number_of_components; ils++){
        value.RegisterNewFab(value.level_sets[ils].velocity_mf, value.bc_velocity, AMREX_SPACEDIM, value.number_of_ghost_cells, "velocity", true); 
        value.RegisterNewFab(value.level_sets[ils].normal_mf, value.bc_velocity, AMREX_SPACEDIM, value.number_of_ghost_cells, "normal", true);
        value.RegisterNewFab(value.level_sets[ils].curvature_mf, &value.bc_nothing, 1, value.number_of_ghost_cells, "curvature", true);   
    }
    }

    {// Register face-centered flux fields
    value.RegisterFaceFab<0>(value.XFlux_mf, &value.bc_nothing, value.number_of_components, value.number_of_ghost_cells, "xflux", true);
    #if AMREX_SPACEDIM >= 2
    value.RegisterFaceFab<1>(value.YFlux_mf, &value.bc_nothing, value.number_of_components, value.number_of_ghost_cells, "yflux", true);
    #endif
    #if AMREX_SPACEDIM == 3
    value.RegisterFaceFab<2>(value.ZFlux_mf, &value.bc_nothing, value.number_of_components, value.number_of_ghost_cells, "zflux", true);
    #endif
    }

    {// Store Newfab to plot Zerols and Tube imfs
    value.RegisterNewFab(value.Zerols_mf, &value.bc_nothing, 1, 0, "ZEROLS", true);
    value.RegisterNewFab(value.cpt_mf, &value.bc_nothing, 1, 0, "CPT", true);
    value.RegisterNewFab(value.Tube_mf, &value.bc_nothing, value.number_of_components, 0, "NB_Tube", true);
    value.RegisterNewFab(value.BA_mf, &value.bc_nothing, value.number_of_components, 0, "NB_Boxes", true);
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

    // Read CFL number and initial time step
    pp.query_required("cflNumber", value.cflNumber);
}

// Define required override functions
void NarrowBandLevelset::Initialize(int lev){
    // Initialize levelset fields
    ic_ls->Initialize(lev, ls_old_mf);
    ic_ls->Initialize(lev, ls_mf);

    // Initialize the zero levelset
    Zerols_imf.reset(new amrex::iMultiFab(ls_mf[lev]->boxArray(), ls_mf[lev]->DistributionMap(), 1, number_of_ghost_cells));
    Zerols_imf->setVal(-1, number_of_ghost_cells);

    // Initialize the cpt flags
    cpt_imf.reset(new amrex::iMultiFab(ls_mf[lev]->boxArray(), ls_mf[lev]->DistributionMap(), 1, number_of_ghost_cells));
    cpt_imf->setVal(-1, number_of_ghost_cells);

    // Initialize the narrowband IntVects
    BA_imf.reset(new amrex::iMultiFab(ls_mf[lev]->boxArray(), ls_mf[lev]->DistributionMap(), 1, number_of_ghost_cells));
    BA_imf->setVal(-1, number_of_ghost_cells);

    // Loop through all levelsets
    for (int ils=0; ils < number_of_components; ils++){
        // Get structure id number
        level_sets[ils].id = ils;

        // Apply boundary conditions for ls ghost cells
        Integrator::ApplyPatch(lev, 0, ls_mf, *ls_mf[lev], *bc_ls, ils);

        // Initialize zerols
        InitializeCPT(lev, ils);

        // Initialize the narrowband tube, box array, and distribution mapping
        auto& tube_imf = level_sets[ils].Tube_imf;
        tube_imf.reset(new amrex::iMultiFab(ls_mf[lev]->boxArray(), ls_mf[lev]->DistributionMap(), 1, number_of_ghost_cells));
        tube_imf->setVal(NarrowBandTubeType::InnerTube, number_of_ghost_cells); // need ghost cells for flux computation
        level_sets[ils].narrowband_ba = ls_mf[lev]->boxArray();
        level_sets[ils].narrowband_dm = ls_mf[lev]->DistributionMap();

        // Define initial narrowband boxarray and distribution mapping
        //amrex::Print() << "Initializing tube" << std::endl;
        UpdateNarrowbandTubeandMapping(lev, ils);

        // Save Tube_imf to Tube_imf
        CopyTubeIMFtoMF(lev, ils);

        // Define velocity_mf containing velocity vector
        ic_velocity->Initialize(lev, level_sets[ils].velocity_mf);

        // Define Geometric quantities - must be after narrowband to update narrowband boxes!
        level_sets[ils].normal_mf[lev]->setVal(0.0, number_of_ghost_cells);
        level_sets[ils].curvature_mf[lev]->setVal(0.0, number_of_ghost_cells);
        //ComputeGeometryQuantities(lev, ils);
    } 

    // Save Zerols_imf to Zerols_mf
    CopyZerolsAndCPTIMFtoMF(lev);

    // Get the proper timestep
    ComputeAndSetNewTimeStep();
    
    // After Initialization, check flux ix Types
    //amrex::Print() << "XFlux ix Type: " << XFlux_mf[lev]->boxArray().ixType() << std::endl;
    //amrex::Print() << "YFlux ix Type: " << YFlux_mf[lev]->boxArray().ixType() << std::endl;
}

void NarrowBandLevelset::InitializeCPT(int lev, int ls_id){
    for (amrex::MFIter mfi(*cpt_imf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box valid_bx = mfi.validbox();

        const auto& ls_arr = ls_mf.Patch(lev, mfi);
        const auto& cpt_arr = cpt_imf->array(mfi);

        amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const Set::Scalar LS = ls_arr(i, j, k, ls_id);
            if (LS < 0.0) cpt_arr(i, j, k) = ls_id;
        });
    }
}

void NarrowBandLevelset::ZeroTheZeros(amrex::iMultiFab& zerols, amrex::iMultiFab& cpt, int ls_id){
    // Loop through Zerols_imf
    for (amrex::MFIter mfi(zerols, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        // Define valid box
        const amrex::Box ghost_bx = mfi.growntilebox(number_of_ghost_cells);

        // Define arrays
        const auto& zerols_arr = zerols.array(mfi);
        const auto& cpt_arr = cpt.array(mfi);

        amrex::ParallelFor(ghost_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (zerols_arr(i, j, k) == ls_id) zerols_arr(i, j, k) = -1;
            if (cpt_arr(i, j, k) == ls_id) cpt_arr(i, j, k) = -1;
        });
    }
}

void NarrowBandLevelset::CopyTubeIMFtoMF(int lev, int ls_id) {
    auto const& tube_imf = level_sets[ls_id].Tube_imf;

    for (amrex::MFIter mfi(*tube_imf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& valid_box = mfi.validbox();
        auto const& Tube_imf_arr = tube_imf->array(mfi);
        auto const& Tube_mf_arr = Tube_mf.Patch(lev, mfi);

        // Convert and copy values
        amrex::ParallelFor(valid_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Tube_mf_arr(i,j,k) = static_cast<amrex::Real>(Tube_imf_arr(i,j,k));  // Ensure type conversion
        });
    }
}

void NarrowBandLevelset::CopyZerolsAndCPTIMFtoMF(int lev){
    for (amrex::MFIter mfi(*Zerols_imf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& valid_box = mfi.validbox();
        auto const& zerols_imf_arr = Zerols_imf->const_array(mfi);  // FIXED ACCESS
        auto const& zerols_mf_arr = Zerols_mf.Patch(lev, mfi);

        auto const& cpt_imf_arr = cpt_imf->const_array(mfi);  // FIXED ACCESS
        auto const& cpt_mf_arr = cpt_mf.Patch(lev, mfi);

        auto const& ba_imf_arr = BA_imf->const_array(mfi);  // FIXED ACCESS
        auto const& ba_mf_arr = BA_mf.Patch(lev, mfi);

        // Scan through the box to find narrow band region
        amrex::ParallelFor(valid_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Convert integer to double
            zerols_mf_arr(i,j,k) = static_cast<Set::Scalar>(zerols_imf_arr(i,j,k));
            cpt_mf_arr(i,j,k) = static_cast<Set::Scalar>(cpt_imf_arr(i,j,k));
            ba_mf_arr(i,j,k) = static_cast<Set::Scalar>(ba_imf_arr(i,j,k));

        });
    }
}

// Utility to convert IntVect list to a simplified BoxArray
amrex::BoxArray MakeBoxArrayFromIntVects(const amrex::Vector<amrex::IntVect>& ivects, bool simplify = true) {
    amrex::BoxList bl;
    for (const auto& iv : ivects) {
        bl.push_back(amrex::Box(iv, iv));
    }
    if (simplify) bl.simplify();
    return amrex::BoxArray(bl);
}

// UpdateNarrowband: replicated Fortran LSTubeInfo logic
void NarrowBandLevelset::UpdateNarrowband(int lev, int ls_id) {
    // Define geometry constants
    //const amrex::Box& domain_box = geom[lev].Domain();
    const Set::Scalar* DX = geom[lev].CellSize();
    const Set::Scalar min_DX = *std::min_element(DX, DX + AMREX_SPACEDIM);
    const Set::Scalar inner_tube_width = inner_narrow_band_width * min_DX;
    const Set::Scalar tube_width = narrow_band_width * min_DX;

    // Define narrowband constant values
    const int Interface      = NarrowBandTubeType::Interface;
    const int InnerTube      = NarrowBandTubeType::InnerTube;
    const int OuterTube      = NarrowBandTubeType::OuterTube;
    const int InnerEdge      = NarrowBandTubeType::InnerEdge;
    const int OuterEdge      = NarrowBandTubeType::OuterEdge;
    const int OutsideBandNeg = NarrowBandTubeType::OutsideNarrowBandNeg;
    const int OutsideBandPos = NarrowBandTubeType::OutsideNarrowBandPos;

    // Define an IntVect to store indices of cells within narrowband
    //amrex::Vector<amrex::IntVect> narrowband_cells_host;

    // Save level_set data structure for reference
    auto& ls_data = level_sets[ls_id];
    const amrex::IntVect debug(AMREX_D_DECL(12, 12, 0));

    // Define MultiFab properties of old band to create temporary (i)MulitFab objects
    const amrex::BoxArray& ba = ls_data.narrowband_ba;
    const amrex::DistributionMapping& dm = ls_data.narrowband_dm;

    // Initialize zerols iMultifab to store interface cells
    amrex::iMultiFab zerols(ba, dm, 1, number_of_ghost_cells);
    zerols.ParallelCopy(*Zerols_imf, 0, 0, 1, number_of_ghost_cells, number_of_ghost_cells);

    // Initialize CPT iMultiFab for object tagging
    amrex::iMultiFab cpt(ba, dm, 1, number_of_ghost_cells);
    cpt.ParallelCopy(*cpt_imf, 0, 0, 1, number_of_ghost_cells, number_of_ghost_cells);

    // Initialize temporary narrowband iMultifab to store flags
    amrex::iMultiFab narrowband(ba, dm, 1, number_of_ghost_cells);
    narrowband.setVal(InnerTube, number_of_ghost_cells);

    // Initialize temporary levelset MultiFab 
    amrex::MultiFab ls(ba, dm, 1, number_of_ghost_cells);
    ls.ParallelCopy(*ls_mf[lev], ls_id, 0, 1, number_of_ghost_cells, number_of_ghost_cells);

    // Reset zerols, cpt, and narrowband imfs
    ZeroTheZeros(zerols, cpt, ls_id);

    // Fill all temporary ghost cells
    zerols.FillBoundary();
    cpt.FillBoundary();
    narrowband.FillBoundary();
    ls.FillBoundary();

    // Perform initial loop through narrowband cells to identify interface and tube cells
    for (amrex::MFIter mfi(ls, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        // Define valid and ghost boxes
        const amrex::Box& ghost_bx = mfi.growntilebox(number_of_ghost_cells);
        //const amrex::Box& valid_bx = mfi.validbox();

        // Define array views
        const auto& zero_arr = zerols.array(mfi);
        const auto& nb_arr = narrowband.array(mfi);
        const auto& cpt_arr = cpt.array(mfi);
        const auto& ls_arr = ls.array(mfi); 

        // Loop through ghost box to update narrowband ghost cells
        // Since there are no automatic nuemann conditions
        ParallelFor(ghost_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const amrex::IntVect coord(AMREX_D_DECL(i, j, k)); 
            //if (coord == debug) amrex::Print() << "Debugging point " << debug << " in narrowband tube update" << std::endl;

            // Clamp levelset value and define abs|ls|
            Set::Scalar LS = ls_arr(coord);
            LS = std::clamp(LS, -tube_width, tube_width);
            const Set::Scalar abs_LS = std::abs(LS);
            ls_arr(coord) = LS;
            //if (coord == debug) amrex::Print() << "LS: " << ls_arr(coord) << std::endl;

            // Assign cpt flag
            if (LS < 0.0) cpt_arr(i, j, k) = ls_id;

            // Check if cell is within narrowband
            if (abs_LS < tube_width) {
                // If not within inner_tube_width, skip Interface check
                // and set to Outer/Inner Tube
                if (abs_LS >= inner_tube_width) {
                    nb_arr(i, j, k) = LS > 0 ? OuterTube : InnerTube;
                    return;
                }

                // Cell is within inner_tube_width. Check neighbors for zero crossings
                // and set to Interface if true
                for (int d = 0; d < Neighbors::num_neighbors; ++d) {
                    int ni = i + Neighbors::offsets[d][0];
                    int nj = j + Neighbors::offsets[d][1];
                    int nk = k + Neighbors::offsets[d][2];

                    // Define nbr index
                    amrex::IntVect nbr(AMREX_D_DECL(ni, nj, nk));

                    // Ensure neighbor is within ghost range
                    if (!ghost_bx.contains(nbr)) continue;

                    // Define nbr Tube value
                    int nb_nbr = nb_arr(nbr);

                    // Ensure nbr is within Tube
                    if (std::abs(nb_nbr) == OutsideBandPos) continue;

                    // Get the levelset value of the nbr
                    Set::Scalar LSnbr = ls_arr(nbr);

                    // If LS * LSnbr is negative, mark as Interface and break loop
                    if (LS * LSnbr <= 0.0){
                        // Mark zerols with ls_id and cell to Interface
                        zero_arr(i, j, k) = ls_id;
                        nb_arr(i, j, k) = Interface;
                        break;
                    }
                }

                // If cell is not an Interface cell, mark as Outer/Inner Tube
                if (nb_arr(i, j, k) != Interface) {
                    nb_arr(i, j, k) = LS > 0 ? OuterTube : InnerTube;
                }
            } 
            // If cell is not within Tube, mark as outside
            else {
                nb_arr(i, j, k) = LS > 0 ? OutsideBandPos : OutsideBandNeg;
            }

            //if (coord == debug) amrex::Print() << "NB: " << nb_arr(coord) << std::endl;
        });
    }

    // Update ghost cells for new narrowband tube
    narrowband.FillBoundary();

    // Perform second loop to pad the narrowband one cell for stencil construction
    for (amrex::MFIter mfi(narrowband, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        // Get ghost and valid boxes
        const amrex::Box& ghost_bx = mfi.growntilebox(number_of_ghost_cells);

        // Define array view
        const auto& nb_arr = narrowband.array(mfi);

        ParallelFor(ghost_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // Get coordinate
            const amrex::IntVect coord(AMREX_D_DECL(i, j, k));

            // Get the cell narrowband value
            const int nb = nb_arr(coord);
            const int abs_nb = std::abs(nb);

            // Skip over Interface and outside tube cells
            if (abs_nb != OuterTube) return;

            // Skip over ghost cells
            //if (!domain_box.contains(coord)) return;

            // nbr loop
            for (int d = 0; d < Neighbors::num_neighbors; ++d) {
                int ni = i + Neighbors::offsets[d][0];
                int nj = j + Neighbors::offsets[d][1];
                int nk = k + Neighbors::offsets[d][2];
                
                // Define nbr coords
                amrex::IntVect nbr(AMREX_D_DECL(ni, nj, nk));
            
                // Skip if outside ghost range
                if (!ghost_bx.contains(nbr)) continue;

                // Get neighbor Tube value
                const int nb_nbr = nb_arr(nbr);

                // Skip if nbr is within Tube
                if (std::abs(nb_nbr) <= OuterTube) continue;

                // If nbr is Outside Tube, promote to narrowband
                if (nb_nbr == OutsideBandPos) nb_arr(nbr) = OuterEdge;
                if (nb_nbr == OutsideBandNeg) nb_arr(nbr) = InnerEdge;
            }
            //if (coord == debug) amrex::Print() << "NB after edge loop: " << nb_arr(coord) << std::endl;
        });

        /*// Loop through valid_bx and store narrowband cell IntVects
        const auto& lo = valid_bx.smallEnd();
        const auto& hi = valid_bx.bigEnd();

        // Loop through valid box to save the indices of narrowband cells
        // depending on AMREX_SPACEDIM - WILL GPU LATER
        #if (AMREX_SPACEDIM == 1)   
        for (int i = lo[0]; i <= hi[0]; ++i) {
            if (std::abs(nb_arr(i, 0, 0)) < OutsideBandPos) {
                narrowband_cells_host.emplace_back(amrex::IntVect(AMREX_D_DECL(i, 0, 0)));
            }
        }
        #elif (AMREX_SPACEDIM == 2)
        for (int j = lo[1]; j <= hi[1]; ++j) {
            for (int i = lo[0]; i <= hi[0]; ++i) {
                if (std::abs(nb_arr(i, j, 0)) < OutsideBandPos) {
                    narrowband_cells_host.emplace_back(amrex::IntVect(AMREX_D_DECL(i, j, 0)));
                }
            }
        }
        #elif (AMREX_SPACEDIM == 3)
        for (int k = lo[2]; k <= hi[2]; ++k) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
                for (int i = lo[0]; i <= hi[0]; ++i) {
                    if (std::abs(nb_arr(i, j, k)) < OutsideBandPos) {
                        narrowband_cells_host.emplace_back(amrex::IntVect(AMREX_D_DECL(i, j, k)));
                    }
                }
            }
        }
        #endif*/
    }

    /*// Compute Box Array and Distribution Mapping using IntVects
    if (!narrowband_cells_host.empty()){
        ls_data.has_narrowband = true;
        ls_data.narrowband_ba = MakeBoxArrayFromIntVects(narrowband_cells_host);
        ls_data.narrowband_dm = amrex::DistributionMapping(ls_data.narrowband_ba);
    }
    else{
        ls_data.has_narrowband = false;
        ls_data.narrowband_ba = ls_mf[lev]->boxArray(); // Make these empty?
        ls_data.narrowband_dm = ls_mf[lev]->distributionMap;
    }

    // Update narrowband box flags
    const amrex::BoxArray& full_domain_ba = ls_mf[lev]->boxArray();
    amrex::Vector<int> flags(full_domain_ba.size(), 0);

    for (const auto& iv : narrowband_cells_host) {
        for (int i = 0; i < full_domain_ba.size(); ++i) {
            if (full_domain_ba[i].contains(iv)) {
                flags[i] = 1;
                break;
            }
        }
    }

    ls_data.narrowband_flags = std::move(flags);*/

    // Fill narrowband ghost cells after editing
    narrowband.FillBoundary();

    // Create Box Array from narrowband mask
    amrex::BoxList narrowband_boxes;

    for (amrex::MFIter mfi(narrowband, false); mfi.isValid(); ++mfi) {
        const auto& mask_arr = narrowband.const_array(mfi);
        const amrex::Box& bx = mfi.validbox();

        bool has_owned_cell = false;
        for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
            const int abs_nb = std::abs(mask_arr(iv));
            if (abs_nb < OutsideBandPos) {
                has_owned_cell = true;
                break;
            }
        }

        if (has_owned_cell) {
            narrowband_boxes.push_back(bx);
        }
    }

    // Create box array if list is not empty
    if ((narrowband_boxes.isNotEmpty())){
        // Create and store BoxArray and Distrbution Mapping
        amrex::BoxArray narrowband_ba(narrowband_boxes);
        narrowband_ba.coarsen(1);  // optional
        narrowband_ba.refine(1);

        ls_data.has_narrowband = true;
        ls_data.narrowband_ba = narrowband_ba;
        ls_data.narrowband_dm = amrex::DistributionMapping(narrowband_ba);
    }
    else{
        ls_data.has_narrowband = false;
        ls_data.narrowband_ba = ls_mf[lev]->boxArray();
        ls_data.narrowband_dm = ls_mf[lev]->distributionMap;
    }

    const amrex::BoxArray& full_ba = ls_mf[lev]->boxArray();  // Full domain boxes
    const int nboxes = full_ba.size();
    
    amrex::Vector<int> narrowband_flags(nboxes, 0);
    
    // Convert BoxList (narrowband_boxes) to BoxArray
    amrex::BoxArray narrowband_ba(narrowband_boxes);
    
    // Scratch space to hold intersection results
    std::vector<std::pair<int, amrex::Box>> isects;
    
    // Loop over each narrowband box and find which full domain boxes intersect
    for (int nbx = 0; nbx < narrowband_ba.size(); ++nbx) {
        const amrex::Box& nb = narrowband_ba[nbx];
    
        isects.clear();
        full_ba.intersections(nb, isects);
    
        for (const auto& [i, _] : isects) {
            narrowband_flags[i] = 1;  // Mark intersecting box
        }
    }
    
    // Store result in the level set structure
    ls_data.narrowband_flags = std::move(narrowband_flags);

    // Fill all ghosts before copying back
    zerols.FillBoundary();
    cpt.FillBoundary();
    narrowband.FillBoundary();
    ls.FillBoundary();

    // Copy back to full domain (i)MultiFabs 
    Zerols_imf->ParallelCopy(zerols, 0, 0, 1, number_of_ghost_cells, number_of_ghost_cells);
    cpt_imf->ParallelCopy(cpt, 0, 0, 1, number_of_ghost_cells, number_of_ghost_cells);
    level_sets[ls_id].Tube_imf->ParallelCopy(narrowband, 0, 0, 1, number_of_ghost_cells, number_of_ghost_cells);
    ls_mf[lev]->ParallelCopy(ls, 0, ls_id, 1, number_of_ghost_cells, number_of_ghost_cells);
}

/*void NarrowBandLevelset::ComputeNarrowBandBoxList(int lev, int ls_id) {
    auto& ls_data = level_sets[ls_id];
    auto& tube_imf = ls_data.Tube_imf;

    amrex::BoxList narrow_band_boxes;
    amrex::Gpu::DeviceVector<int> narrow_band_flags(tube_imf->local_size(), 0);
    auto* d_narrow_band_flags = narrow_band_flags.data();

    amrex::Vector<amrex::Box> temp_boxes;
    int box_index = 0;

    for (amrex::MFIter mfi(*tube_imf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& valid_box = mfi.tilebox();
        auto const& nb_arr = tube_imf->array(mfi);
        int* flag_ptr = d_narrow_band_flags + box_index;

        // Define OutsideNarrowband value
        const int OutsideNarrowband = NarrowBandTubeType::OutsideNarrowBandPos;

#ifdef AMREX_USE_GPU
        // GPU path — search entire tile in parallel
        amrex::ParallelFor(valid_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (nb_arr(i, j, k) != OutsideNarrowband)) {
                amrex::Gpu::Atomic::Or(flag_ptr, 1);
            }
        });
#else
        // CPU path — loop based on AMREX_SPACEDIM
        const auto& lo = valid_box.smallEnd();
        const auto& hi = valid_box.bigEnd();

#if (AMREX_SPACEDIM == 1)
        for (int i = lo[0]; i <= hi[0]; ++i) {
            if (std::abs(nb_arr(i, 0, 0)) < OutsideNarrowband) {
                *flag_ptr = 1;
                goto done;
            }
        }
#elif (AMREX_SPACEDIM == 2)
        for (int j = lo[1]; j <= hi[1]; ++j) {
            for (int i = lo[0]; i <= hi[0]; ++i) {
                if (std::abs(nb_arr(i, j, 0)) < OutsideNarrowband) {
                    *flag_ptr = 1;
                    goto done;
                }
            }
        }
#elif (AMREX_SPACEDIM == 3)
        for (int k = lo[2]; k <= hi[2]; ++k) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
                for (int i = lo[0]; i <= hi[0]; ++i) {
                    if (std::abs(nb_arr(i, j, k)) < OutsideNarrowband) {
                        *flag_ptr = 1;
                        goto done;
                    }
                }
            }
        }
#endif
done:;
#endif
        temp_boxes.push_back(valid_box);
        ++box_index;
    }

    amrex::Gpu::Device::synchronize();

    std::vector<int> h_flags(narrow_band_flags.size());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, narrow_band_flags.begin(), narrow_band_flags.end(), h_flags.begin());

    for (int i = 0; i < temp_boxes.size(); ++i) {
        if (h_flags[i] != 0) {
            narrow_band_boxes.push_back(temp_boxes[i]);
        }
    }

    ls_data.narrowband_boxes = narrow_band_boxes;
}

void NarrowBandLevelset::ComputeNarrowBandMapping(int lev, int ls_id){
    auto& ls_data = level_sets[ls_id];

    // Check if there are any narrowband boxes
    amrex::BoxList narrowband_boxes = ls_data.narrowband_boxes;
    if (narrowband_boxes.size()> 0) {
        ls_data.has_narrowband = true;  // Flag indicating valid narrowband
        amrex::BoxArray ba(narrowband_boxes);
        ls_data.narrowband_ba = ba;
        ls_data.narrowband_dm = amrex::DistributionMapping(ba);
    } else {
        ls_data.has_narrowband = false; // No narrowband region found
    }
}

void NarrowBandLevelset::UpdateNarrowBandFlags(int lev, int ls_id) {
    const amrex::BoxArray& ba = ls_mf[lev]->boxArray();  // Full domain boxes
    const int nboxes = ba.size();

    amrex::Vector<int> narrowband_flags(nboxes, 0);

    // Convert BoxList (narrowband_boxes) to BoxArray
    amrex::BoxArray narrowband_ba(level_sets[ls_id].narrowband_boxes);

    // Scratch space to hold intersection results
    std::vector<std::pair<int, amrex::Box>> isects;

    // Loop over each narrowband box and find which full domain boxes intersect
    for (int nbx = 0; nbx < narrowband_ba.size(); ++nbx) {
        const amrex::Box& nb = narrowband_ba[nbx];

        isects.clear();
        ba.intersections(nb, isects);

        for (const auto& [i, _] : isects) {
            narrowband_flags[i] = 1;  // Mark intersecting box
        }
    }

    // Store result in the level set structure
    level_sets[ls_id].narrowband_flags = std::move(narrowband_flags);
}

// Utility to convert IntVect list to a simplified BoxArray
amrex::BoxArray MakeBoxArrayFromIntVects(const amrex::Vector<amrex::IntVect>& ivects, bool simplify = true) {
    amrex::BoxList bl;
    for (const auto& iv : ivects) {
        bl.push_back(amrex::Box(iv, iv));
    }
    if (simplify) bl.simplify();
    return amrex::BoxArray(bl);
}

// Create BoxArray/DM from narrowband cells
void NarrowBandLevelset::ComputeNarrowBandMapping(int lev, int ls_id) {
    auto& ls_data = level_sets[ls_id];
    if (!ls_data.narrowband_cells.empty()) {
        ls_data.has_narrowband = true;
        ls_data.narrowband_ba = MakeBoxArrayFromIntVects(ls_data.narrowband_cells);
        ls_data.narrowband_dm = amrex::DistributionMapping(ls_data.narrowband_ba);
    } else {
        ls_data.has_narrowband = false;
    }
}

// Mark full-domain BoxArray entries that overlap narrowband cells
void NarrowBandLevelset::UpdateNarrowBandFlags(int lev, int ls_id) {
    const amrex::BoxArray& ba = ls_mf[lev]->boxArray();
    amrex::Vector<int> flags(ba.size(), 0);

    for (const auto& iv : level_sets[ls_id].narrowband_cells) {
        for (int i = 0; i < ba.size(); ++i) {
            if (ba[i].contains(iv)) {
                flags[i] = 1;
                break;
            }
        }
    }

    level_sets[ls_id].narrowband_flags = std::move(flags);
}*/

void NarrowBandLevelset::UpdateNarrowbandTubeandMapping(int lev, int ls_id){
    // Tag 0 levelset gridpoints - zero will be done outside loop
    UpdateNarrowband(lev, ls_id);

    // Define narrowband box list 
    //ComputeNarrowBandBoxList(lev, ls_id);

    // Get the distribution mapping for processors
    //ComputeNarrowBandMapping(lev, ls_id);

    // Update flags
    //UpdateNarrowBandFlags(lev, ls_id);
}

// CHECK VALID BOX TILING HERE - MAY ADD GHOST BOX AND APPLY TO GETSTENCIL
void NarrowBandLevelset::ComputeNormal(int lev, int ls_id) {
    // Define DX for gradient
    const Set::Scalar* DX = geom[lev].CellSize();

    // Define OutsideNarrowband constant
    const int OutsideNarrowband = NarrowBandTubeType::OutsideNarrowBandPos;

    // Create a temporary narrowband MultiFab
    amrex::BoxArray narrowband_ba = level_sets[ls_id].narrowband_ba;
    amrex::DistributionMapping narrowband_dm = level_sets[ls_id].narrowband_dm;
    amrex::iMultiFab narrowband(narrowband_ba, narrowband_dm, 1, number_of_ghost_cells);  
    amrex::MultiFab narrowband_phi(narrowband_ba, narrowband_dm, 1, number_of_ghost_cells);  // Level set values
    amrex::MultiFab narrowband_normal(narrowband_ba, narrowband_dm, AMREX_SPACEDIM, number_of_ghost_cells); // Normal field

    // Initialize narrowband_normal
    narrowband_normal.setVal(0.0);

    // Copy narrowband and level set values from full level set to narrowband (including ghost cells)
    narrowband.ParallelCopy(*level_sets[ls_id].Tube_imf, 0, 0, 1, number_of_ghost_cells, number_of_ghost_cells);
    narrowband_phi.ParallelCopy(*ls_mf[lev], ls_id, 0, 1, number_of_ghost_cells, number_of_ghost_cells);

    // Fill ghost cells
    narrowband_phi.FillBoundary();

    // Iterate over only the narrowband cells
    for (amrex::MFIter mfi(narrowband_phi, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& valid_bx = mfi.validbox();
        const amrex::Box& ghost_bx = mfi.growntilebox(number_of_ghost_cells);

        auto const& phi_arr = narrowband_phi.array(mfi);
        auto const& normal = narrowband_normal.array(mfi);
        auto const& nb = narrowband.array(mfi);

        amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            if (std::abs(nb(i,j,k)) == OutsideNarrowband) return; // Skip non-narrowband cells
 
            // Utilize Numeric::Gradient for simplified computations
            std::array<Numeric::StencilType, AMREX_SPACEDIM> stencil = Numeric::GetStencil(i, j, k, valid_bx);
            Set::Vector gradient = Numeric::Gradient(phi_arr, i, j, k, 0, DX, stencil);
            Set::Scalar grad_mag = gradient.norm();

            for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
                normal(i, j, k, dim) = (grad_mag > 0.0) ? (gradient[dim] / grad_mag) : 0.0;
            }
        });
    }

    // Fill boundary before copying back
    narrowband_normal.FillBoundary();

    // Copy computed normals back to the full normal MultiFab (including ghost cells)
    level_sets[ls_id].normal_mf[lev]->ParallelCopy(narrowband_normal, 0, 0, AMREX_SPACEDIM, number_of_ghost_cells, number_of_ghost_cells);
}

// CHECK VALID BOX TILING HERE - MAY ADD GHOST BOX AND APPLY TO GETSTENCIL
void NarrowBandLevelset::ComputeCurvature(int lev, int ls_id){
    // Define DX for gradient
    const Set::Scalar* DX = geom[lev].CellSize();

    // Define OutsideNarrowband constant
    const int OutsideNarrowband = NarrowBandTubeType::OutsideNarrowBandPos;

    // Create a temporary narrowband MultiFab
    amrex::BoxArray narrowband_ba = level_sets[ls_id].narrowband_ba;
    amrex::DistributionMapping narrowband_dm = level_sets[ls_id].narrowband_dm;
    amrex::iMultiFab narrowband(narrowband_ba, narrowband_dm, 1, number_of_ghost_cells);
    amrex::MultiFab narrowband_normal(narrowband_ba, narrowband_dm, AMREX_SPACEDIM, number_of_ghost_cells); 
    amrex::MultiFab narrowband_curve(narrowband_ba, narrowband_dm, 1, number_of_ghost_cells); 

    // Initialize narrowband curvature
    narrowband_curve.setVal(0.0);

    // Copy narrowband and normal
    narrowband.ParallelCopy(*level_sets[ls_id].Tube_imf, 0, 0, 1, number_of_ghost_cells, number_of_ghost_cells);
    narrowband_normal.ParallelCopy(*level_sets[ls_id].normal_mf[lev], 0, 0, AMREX_SPACEDIM, number_of_ghost_cells, number_of_ghost_cells);
    
    // Fill ghost cells
    narrowband_normal.FillBoundary();
    
    for (amrex::MFIter mfi(narrowband_curve, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& valid_bx = mfi.tilebox();
        auto const& normal = narrowband_normal.array(mfi); 
        auto const& curvature = narrowband_curve.array(mfi);
        auto const& nb = narrowband.array(mfi);

        amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            // skip if not within narowband
            if (std::abs(nb(i,j,k)) == OutsideNarrowband) return;
            
            // Get diagonals of matrix
            Set::Scalar sum = 0;
            std::array<Numeric::StencilType, AMREX_SPACEDIM> stencil = Numeric::GetStencil(i, j, k, valid_bx);
            for (int dim=0; dim < AMREX_SPACEDIM; dim++){
                // Utilize Numeric::Gradient for simplified computations
                Set::Vector gradient = Numeric::Gradient(normal, i, j, k, dim, DX, stencil);
                sum += gradient[dim];
            }
            curvature(i, j, k) = sum;
        });
    }  

    // Fill boundary before copying back
    narrowband_curve.FillBoundary();

    // Copy computed normals back to the full normal MultiFab (including ghost cells)
    level_sets[ls_id].curvature_mf[lev]->ParallelCopy(narrowband_curve, 0, 0, 1, number_of_ghost_cells, number_of_ghost_cells);
}

void NarrowBandLevelset::ComputeGeometryQuantities(int lev, int ls_id){
    ComputeNormal(lev, ls_id);
    ComputeCurvature(lev, ls_id);  
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
    for (int ils=0; ils < number_of_components; ils++){
        for (int lev = 0; lev <= maxLevel(); ++lev) {  // Iterate over AMR levels
            const Set::Scalar* dx = geom[lev].CellSize();  // Access the geometry at level `lev`
            const Set::Scalar min_DX = *std::min_element(dx, dx + AMREX_SPACEDIM);
    
            for (amrex::MFIter mfi(*level_sets[ils].velocity_mf[lev], false); mfi.isValid(); ++mfi) {
                const amrex::Box& bx = mfi.tilebox();  // Iterate over tiles in the multifab
                auto const& velocity_arr = level_sets[ils].velocity_mf.Patch(lev, mfi);  // Velocity field
    
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
    // Update current timestep (for reinitialize)
    current_timestep ++;

    // Advect
    Advect(lev, time, dt);

    // Update the narrowband
    for (int ls_id = 0; ls_id < number_of_components; ls_id++){
        // Apply Boundary conditions after advection
        Integrator::ApplyPatch(lev, time, ls_mf, *ls_mf[lev], *bc_ls, ls_id);

        // Reinitialize
        //amrex::Print() << "Reinitializing LS" << std::endl;
        Reinitialize(lev, ls_id);

        // Apply Boundary conditions after reinitialization
        Integrator::ApplyPatch(lev, time, ls_mf, *ls_mf[lev], *bc_ls, ls_id);
        
        // Update narrowband info
        //amrex::Print() << "Updating tube after Reinit" << std::endl;
        UpdateNarrowbandTubeandMapping(lev, ls_id);

        // Apply Boundary conditions after tube update
        Integrator::ApplyPatch(lev, time, ls_mf, *ls_mf[lev], *bc_ls, ls_id);

        // Compute geometries
        ComputeGeometryQuantities(lev, ls_id);

        // Save Tube_imf to Tube_imf
        CopyTubeIMFtoMF(lev, ls_id);
    }

    // Copy the zero
    CopyZerolsAndCPTIMFtoMF(lev);
}

void NarrowBandLevelset::Advect(int lev, Set::Scalar time, Set::Scalar dt){
    /*for (amrex::MFIter mfi(*ls_mf[lev]); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        const amrex::FArrayBox& ls_fab     = (*ls_mf[lev])[mfi];
        const amrex::FArrayBox& ls_old_fab = (*ls_old_mf[lev])[mfi];
    
        if (bx.contains(amrex::IntVect(12, 12))) {
            amrex::Print() << "Before swap:\n";
            amrex::Print() << "LS(12,12): " << ls_fab(amrex::IntVect(12, 12), 0) << "\n";
            amrex::Print() << "LS_OLD(12,12): " << ls_old_fab(amrex::IntVect(12, 12), 0) << "\n";
        }
    }*/
    
    // Perform the swap
    //std::swap(*ls_mf[lev], *ls_old_mf[lev]);
    amrex::MultiFab::Copy(*ls_old_mf[lev], *ls_mf[lev], 0, 0, 1, ls_mf[lev]->nGrow());
    
    /*for (amrex::MFIter mfi(*ls_mf[lev]); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        const amrex::FArrayBox& ls_fab     = (*ls_mf[lev])[mfi];
        const amrex::FArrayBox& ls_old_fab = (*ls_old_mf[lev])[mfi];
    
        if (bx.contains(amrex::IntVect(12, 12))) {
            amrex::Print() << "After swap:\n";
            amrex::Print() << "LS(12,12): " << ls_fab(amrex::IntVect(12, 12), 0) << "\n";
            amrex::Print() << "LS_OLD(12,12): " << ls_old_fab(amrex::IntVect(12, 12), 0) << "\n";
        }
    }*/

    // Update the interface velocity
    UpdateInterfaceVelocity(lev);

    // Compute the flux
    switch (temporal_scheme) {
        case TimeSteppingScheme::ForwardEuler: {
            
            fluxHandler->SetReconstruction(std::make_shared<Numeric::FirstOrderReconstruction<NarrowBandLevelset>>());
            fluxHandler->SetFluxMethod(std::make_shared<Numeric::LocalLaxFriedrichsMethod<NarrowBandLevelset>>());

            timeStepper->SetTimeSteppingScheme(std::make_shared<Numeric::EulerForwardScheme<NarrowBandLevelset>>());

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
}

void NarrowBandLevelset::UpdateInterfaceVelocity(int lev){
    for (int ils=0; ils < number_of_components; ils++){
        for (amrex::MFIter mfi(*level_sets[ils].velocity_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.tilebox();
            auto const& vel_arr = level_sets[ils].velocity_mf.Patch(lev,mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                for (int dim=0; dim < AMREX_SPACEDIM; ++dim){
                    vel_arr(i,j,k,dim) = vel_arr(i,j,k,dim);
                }
            });
        }
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Set::Scalar SafeLSAccess(const amrex::Array4<const Set::Scalar>& arr,
                        const amrex::IntVect& coord,
                        const amrex::IntVect& offset,
                        const amrex::Box& domain)
{
    const amrex::IntVect nbr = coord + offset;
    return domain.contains(nbr) ? arr(nbr) : arr(coord); // Neumann BC
}

void NarrowBandLevelset::Reinitialize(int lev, int ls_id) {
    // Define reinitialization constants
    const Set::Scalar reinit_tolerance = 1e-3;
    const int max_iterations = 50;
    const amrex::IntVect debug(AMREX_D_DECL(12, 12, 0));

    // Define geometry constants
    const amrex::Box& domain_box = geom[lev].Domain();
    const Set::Scalar* DX    = geom[lev].CellSize();
    const Set::Scalar min_DX = *std::min_element(DX, DX + AMREX_SPACEDIM);
    const Set::Scalar tau    = cflReinit * min_DX;

    // Define Narrowband constants
    const Set::Scalar INNERTUBE = inner_narrow_band_width * min_DX; // used for testing sign changes
    const int Interface         = NarrowBandTubeType::Interface;
    const int OutsideNarrowband = NarrowBandTubeType::OutsideNarrowBandPos;
    
    // Get narrowband geometry
    amrex::BoxArray ba            = level_sets[ls_id].narrowband_ba;
    amrex::DistributionMapping dm = level_sets[ls_id].narrowband_dm;
    
    // Allocate temporary (i)multifabs
    amrex::iMultiFab zerols(ba, dm, 1, number_of_ghost_cells); // Used to update zero after advection
    amrex::iMultiFab Narrowband(ba, dm, 1, number_of_ghost_cells);
    amrex::MultiFab LS(ba, dm, 1, number_of_ghost_cells);
    amrex::MultiFab LS_old(ba, dm, 1, number_of_ghost_cells);
    amrex::MultiFab error_mf(ba, dm, 1, 0); // Used to compute Linf norm

    // Initialize (i)multifabs
    // NOTE FOR FUTURE: WILL NEED TO FIND A WAY TO PROPERLY COPY Zerols_imf SO THAT
    // OTHER ls_ids ARE NOT OVERWRITTEN WHEN COPIED BACK. SAME FOR ANY OTHER FUNCTION 
    // THAT USES Zerols_imf
    zerols.setVal(-1, number_of_ghost_cells);
    error_mf.setVal(0.0);
    Narrowband.ParallelCopy(*level_sets[ls_id].Tube_imf, 0, 0, 1, number_of_ghost_cells, number_of_ghost_cells);
    LS.ParallelCopy(*ls_mf[lev], ls_id, 0, 1, number_of_ghost_cells, number_of_ghost_cells);
    LS_old.ParallelCopy(*ls_mf[lev], ls_id, 0, 1, number_of_ghost_cells, number_of_ghost_cells);

    // Fill ghost cells
    Narrowband.FillBoundary();
    zerols.FillBoundary();
    LS.FillBoundary();

    // Redefine the interface cells using exact logic as Fortran
    for (amrex::MFIter mfi(Narrowband, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& ghost_bx = mfi.growntilebox(number_of_ghost_cells); // Potentially use valid box?
        const auto& ls_arr = LS.const_array(mfi);
        const auto& zerols_arr = zerols.array(mfi);

        amrex::ParallelFor(ghost_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // Get current coord and LS value
            const amrex::IntVect coord(AMREX_D_DECL(i, j, k));
            const Set::Scalar LS = ls_arr(coord);
            /*if (coord == debug){
                amrex::Print() << "Debugging point " << debug << " in reinitialize" << std::endl;
                amrex::Print() << "LS in zero loop: " << LS << std::endl;
            }*/

            for (int d = 0; d < Neighbors::num_neighbors; ++d) { 
                int ni = i + Neighbors::offsets[d][0];
                int nj = j + Neighbors::offsets[d][1];
                int nk = k + Neighbors::offsets[d][2];

                // Define nbr coord and skip if outside ghost range
                amrex::IntVect nbr(AMREX_D_DECL(ni, nj, nk));
                if (!ghost_bx.contains(nbr)) continue;

                // Get LS for neighbor to compare signs
                const Set::Scalar LSnbr = ls_arr(nbr);

                if (LS * LSnbr < 0.0 && std::max(std::abs(LS), std::abs(LSnbr)) < INNERTUBE) {
                    zerols_arr(coord) = ls_id;
                    zerols_arr(nbr)   = ls_id;
                }
            }
            //if (coord == debug) amrex::Print() << "LS after zero loop: " << LS << std::endl;
        });
    }

    // Fill zerols ghost cells
    zerols.FillBoundary();
    
    // Perform first order PDE Reinitialization scheme
    for (int iter = 0; iter < max_iterations; iter++){
        //amrex::Print() << "Iteration: " << iter << std::endl;
        // Switch LS with LS_old
        std::swap(LS, LS_old);

        // Fill ghost cells for LS_old
        LS_old.FillBoundary();
        zerols.FillBoundary(); // May remove later 
        Narrowband.FillBoundary(); // May remove later

        // Perform main PDE loop
        for (amrex::MFIter mfi(LS, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& valid_bx = mfi.validbox();

            // Get arrays
            const auto& ls_arr = LS.array(mfi);
            const auto& ls_old_arr = LS_old.array(mfi);
            const auto& nb_arr = Narrowband.const_array(mfi);
            const auto& zerols_arr = zerols.const_array(mfi);
            const auto& error_arr = error_mf.array(mfi);

            amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::IntVect coord(AMREX_D_DECL(i, j, k));
                
                // Skip interface and outerband cells
                const int nb_val = std::abs(nb_arr(coord));
                const int zero_cell = zerols_arr(coord);
                if (nb_val == OutsideNarrowband || zero_cell == ls_id) return;

                // Define sign of ls for upwinding
                const Set::Scalar ls = ls_old_arr(coord);
                //if (coord == debug) amrex::Print() << "LS at beginning of iteration " << iter << " is " << ls << std::endl;
                const Set::Scalar sign_ls = (ls > 0.0 ? 1.0 : -1.0);

                // First-order upwind differences
                // USE Numeric::Stencil HERE?
                Set::Scalar phi_xm = SafeLSAccess(ls_old_arr, coord, amrex::IntVect(AMREX_D_DECL(-1,0,0)), domain_box);
                Set::Scalar phi_xp = SafeLSAccess(ls_old_arr, coord, amrex::IntVect(AMREX_D_DECL(1,0,0)), domain_box);
                Set::Scalar dxm = (ls_old_arr(coord) - phi_xm) / DX[0];
                Set::Scalar dxp = (phi_xp - ls_old_arr(coord)) / DX[0];
   
                /*Set::Scalar dxm = (ls_old_arr(coord) - ls_old_arr(i - 1, j, k)) / DX[0]; // backward
                Set::Scalar dxp = (ls_old_arr(i + 1, j, k) - ls_old_arr(coord)) / DX[0]; // forward*/
                Set::Scalar gx = std::max(
                    std::pow(std::max(sign_ls * dxm, 0.0), 2),
                    std::pow(std::min(sign_ls * dxp, 0.0), 2)
                );

                Set::Scalar gy = 0.0;
                #if AMREX_SPACEDIM >= 2
                Set::Scalar phi_ym = SafeLSAccess(ls_old_arr, coord, amrex::IntVect(AMREX_D_DECL(0,-1,0)), domain_box);
                Set::Scalar phi_yp = SafeLSAccess(ls_old_arr, coord, amrex::IntVect(AMREX_D_DECL(0,+1,0)), domain_box);
                Set::Scalar dym = (ls_old_arr(coord) - phi_ym) / DX[1];
                Set::Scalar dyp = (phi_yp - ls_old_arr(coord)) / DX[1];

                /*Set::Scalar dym = (ls_old_arr(coord) - ls_old_arr(i, j - 1, k)) / DX[1];
                Set::Scalar dyp = (ls_old_arr(i, j + 1, k) - ls_old_arr(coord)) / DX[1];*/
                gy = std::max(
                    std::pow(std::max(sign_ls * dym, 0.0), 2),
                    std::pow(std::min(sign_ls * dyp, 0.0), 2)
                );
                #endif

                Set::Scalar gz = 0.0;
                #if AMREX_SPACEDIM == 3
                Set::Scalar dzm = (ls_old_arr(i, j, k) - ls_old_arr(i, j, k - 1)) / DX[2];
                Set::Scalar dzp = (ls_old_arr(i, j, k + 1) - ls_old_arr(i, j, k)) / DX[2];
                gz = std::max(
                    std::pow(std::max(sign_ls * dzm, 0.0), 2),
                    std::pow(std::min(sign_ls * dzp, 0.0), 2)
                );
                #endif

                // Compute gradient and update
                Set::Scalar grad_phi = std::sqrt(gx + gy + gz);
                ls_arr(coord) = ls - tau * sign_ls * (grad_phi - 1.0);
                //if (coord == debug) amrex::Print() << "LS at end of iteration " << iter << " is " <<  ls_arr(coord) << std::endl;

                // Store error in error_mf
                error_arr(i, j, k) = std::abs(ls_arr(i, j, k) - ls);
            });
        }

        // Break if converged
        Set::Scalar max_error = error_mf.norm0();
        if (max_error < reinit_tolerance) break;
    }

    // Fill levelset and zero levelset neighbor values across processors
    zerols.FillBoundary();
    LS.FillBoundary();
    
    // Copy mapped values back to full domain
    Zerols_imf->ParallelCopy(zerols, 0, 0, 1, number_of_ghost_cells, number_of_ghost_cells);
    ls_mf[lev]->ParallelCopy(LS, 0, ls_id, 1, number_of_ghost_cells, number_of_ghost_cells);
}

void NarrowBandLevelset::TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow){
}

void NarrowBandLevelset::Regrid(int lev, Set::Scalar time) {
}

} // namespace Integrator
