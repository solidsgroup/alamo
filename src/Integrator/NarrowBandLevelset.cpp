// Standard library includes
#include <iostream>
#include <fstream>
#include <iomanip>

// Alamo Includes
#include "IO/ParmParse.H"
#include "NarrowBandLevelset.H"

// IC
#include "IC/Constant.H"
#include "IC/LS/Sphere.H"
#include "IC/LS/Zalesak.H"
#include "IC/Expression.H"

// BC
#include "BC/Constant.H"
#include "BC/Expression.H"

// Numerc
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

    // Define number of components
    pp_query_default("levelset.number_of_levelsets", value.number_of_components, 1);

    {// Define initial and boundary conditions
    
    // Query the IC assuming either LS::Sphere or LS::Zalesak
    pp.select_default<IC::LS::Sphere,IC::LS::Zalesak,IC::Expression>("ls.ic",value.ic_ls,value.geom);
    
    // Assume Neumann BC for levelset field
    pp.select_default<BC::Constant,BC::Expression>("ls.bc",value.bc_ls,value.number_of_components);
    }
    
    {// Define levelset multifab objects - only old and new domain levelset fields
    value.RegisterNewFab(value.ls_mf, value.bc_ls, value.number_of_components, value.number_of_ghost_cells, "LS", true);
    value.RegisterNewFab(value.ls_old_mf, value.bc_ls, value.number_of_components, value.number_of_ghost_cells, "LS_old", false);
    }

    // Resize levelset_data structure to be [amr_max_level+1][number_of_components]
    const int max_lev = value.max_level + 1;
    value.level_sets.resize(max_lev);
    for (int level = 0; level < max_lev; level++){
        value.level_sets[level].resize(value.number_of_components);
    }
    //value.level_sets.resize(value.number_of_components);
    
    {// Initialize velocity field. Will be removed once integrated with ScimitarX integrator
    pp.select_default<IC::Constant,IC::Expression>("velocity.ic",value.ic_velocity,value.geom); // velocity initial condtion
    pp.select_default<BC::Constant,BC::Expression>("velocity.bc",value.bc_velocity,AMREX_SPACEDIM); // velocty boundary condition
    value.RegisterNewFab(value.velocity_mf, value.bc_velocity, AMREX_SPACEDIM, value.number_of_ghost_cells, "velocity", true);
    
    // Select method to update velocity field based on validation method
    pp.query_validate("velocity.method",value.velocity_method,{"constant"});
    }

    {// Define normal and curvature multifabs
    for (int lev = 0; lev < max_lev; lev++){
        for (int ils = 0; ils < value.number_of_components; ils++){
            std::string ils_str = std::to_string(ils);
            std::string fab_name = "normal_" + ils_str + "_";
            value.RegisterNewFab(value.level_sets[lev][ils].normal_mf, value.bc_velocity, AMREX_SPACEDIM, value.number_of_ghost_cells, fab_name, true);
        }
    }
    /*for (size_t ils=0; ils<value.number_of_components; ils++){
        value.RegisterNewFab(value.level_sets[ils].normal_mf, value.bc_velocity, AMREX_SPACEDIM, value.number_of_ghost_cells, "normal", true);
    }*/
    value.RegisterNewFab(value.curvature_mf, &value.bc_nothing, value.number_of_components, 0, "curvature", true);   
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
}

// Define required override functions
void NarrowBandLevelset::Initialize(int lev){
    // Initialize levelset fields
    ic_ls->Initialize(lev, ls_old_mf);
    ic_ls->Initialize(lev, ls_mf);
    ClampPhi(lev);

    // Define the ba/dm for pointers
    amrex::BoxArray ba = ls_mf[lev]->boxArray();
    amrex::DistributionMapping dm = ls_mf[lev]->DistributionMap();

    // Initialize cpt
    cpt_mf[lev]->setVal(-1);

    // Initialize zerols
    Zerols_mf[lev]->setVal(-1);

    // Initialize tube mf 
    Tube_mf[lev]->setVal(NarrowBandTubeType::InsideTube);

    // Define velocity_mf containing velocity vector
    ic_velocity->Initialize(lev, velocity_mf);

    // Initialize structure
    for (int ils = 0; ils < number_of_components; ils++){
        // Save level_set data structure for reference
        auto& ls_data = level_sets[lev][ils];

        // Get structure id number
        ls_data.id = ils;

        // Initialize boolean
        ls_data.has_narrowband = true;
        ls_data.reinitialize = false;
        ls_data.reconstruct_tube = true;

        // Initialize the box array and distribution mapping
        ls_data.narrowband_ba = ba;
        ls_data.narrowband_dm = dm;

        // Initialize ptrs
        // Levelset fields
        ls_data.Tube_ls = std::make_unique<amrex::MultiFab>(ba, dm, 1, number_of_ghost_cells);
        auto& Tube_ls = *ls_data.Tube_ls;
        Tube_ls.ParallelCopy(*ls_mf[lev], ils, 0, 1, number_of_ghost_cells, number_of_ghost_cells);

        ls_data.Tube_ls_old = std::make_unique<amrex::MultiFab>(ba, dm, 1, number_of_ghost_cells);
        auto& Tube_ls_old = *ls_data.Tube_ls_old;
        Tube_ls_old.ParallelCopy(*ls_old_mf[lev], ils, 0, 1, number_of_ghost_cells, number_of_ghost_cells);
        
        // Define initial narrowband boxarray and distribution mapping
        UpdateNarrowbandTubeandMapping(lev, ils);

        // Get normal/curvature
        ComputeGeometryQuantities(lev, ils);
    } 

    // Get the proper timestep
    ComputeAndSetNewTimeStep();
}

void NarrowBandLevelset::ClampPhi(int lev) {
    // Define geometry constants
    const Set::Scalar* DX = geom[lev].CellSize();
    const Set::Scalar min_DX = *std::min_element(DX, DX + AMREX_SPACEDIM);
    const Set::Scalar tube_width = narrow_band_width * min_DX;

    for (amrex::MFIter mfi(*ls_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& ghost_box = mfi.growntilebox(number_of_ghost_cells);
        auto const& ls_arr = ls_mf.Patch(lev, mfi);
        auto const& ls_old_arr = ls_old_mf.Patch(lev, mfi);

        // Convert and copy values
        amrex::ParallelFor(ghost_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            for (int ils=0; ils < number_of_components; ils++){
                // Get coord
                amrex::IntVect coord(AMREX_D_DECL(i,j,k));

                // Define initial LS value
                Set::Scalar LS = ls_arr(coord, ils);
                Set::Scalar LS_old = ls_old_arr(coord, ils);

                // Clamp
                LS = std::clamp(LS, -tube_width, tube_width);
                LS_old = std::clamp(LS_old, -tube_width, tube_width);

                // Save
                ls_arr(coord, ils) = LS;
                ls_old_arr(coord, ils) = LS_old;
            }
        });
    }
}

void FillTubeNeumannBC(
    amrex::MultiFab& phi,
    const amrex::Box& tube_domain,
    const amrex::Box& physical_domain,
    int nghost)
{
    // Grow the tube_geom by nghost and clip to physical domain
    amrex::Box grown_tube = tube_domain;
    grown_tube.grow(nghost);
    grown_tube &= physical_domain; // ensure stays in bounds

    for (amrex::MFIter mfi(phi); mfi.isValid(); ++mfi)
    {
        auto const& phi_arr = phi.array(mfi);
        const amrex::Box& gbx = mfi.growntilebox(); // includes ghosts

        amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            amrex::IntVect iv(AMREX_D_DECL(i, j, k));

            if (!tube_domain.contains(iv))
            {
                amrex::IntVect clamped = iv;

                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    if (iv[d] < tube_domain.smallEnd(d)) {
                        clamped[d] = tube_domain.smallEnd(d);
                    } else if (iv[d] > tube_domain.bigEnd(d)) {
                        clamped[d] = tube_domain.bigEnd(d);
                    }
                }

                if (gbx.contains(clamped)){
                    phi_arr(i, j, k) = phi_arr(clamped, 0);
                }
            }
        });
    }
}

void NarrowBandLevelset::BuildNarrowbandMask(
    amrex::MultiFab& ls,
    amrex::iMultiFab& narrowband,
    amrex::iMultiFab& zerols,
    amrex::iMultiFab& cpt,
    const int ls_id,
    const Set::Scalar tube_width,
    const Set::Scalar inner_tube_width)
{
    // Define constants
    constexpr int Interface   = NarrowBandTubeType::Interface;
    constexpr int InsideTube  = NarrowBandTubeType::InsideTube;
    constexpr int OutsideTube = NarrowBandTubeType::OutsideTube;

    for (amrex::MFIter mfi(ls, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& ghost_bx = mfi.growntilebox(number_of_ghost_cells);

        const auto& ls_arr = ls.array(mfi);
        const auto& nb_arr = narrowband.array(mfi);
        const auto& zero_array = zerols.array(mfi);
        const auto& cpt_arr = cpt.array(mfi);

        ParallelFor(ghost_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const amrex::IntVect coord(AMREX_D_DECL(i, j, k));

            Set::Scalar LS = std::clamp(ls_arr(coord, 0), -tube_width, tube_width);
            Set::Scalar absLS = amrex::Math::abs(LS);
            ls_arr(coord) = LS;

            int in_tube = (absLS < tube_width) ? 1 : 0;
            int in_inner = (absLS < inner_tube_width) ? 1 : 0;

            int nb_val = NarrowBandTubeType::OutsideTube;

            // Set CPT flag (inside zero levelset)
            if (LS < 0.0) cpt_arr(coord, 0) = ls_id;

            if (in_tube == 1) {
                nb_val = NarrowBandTubeType::InsideTube;  // default within tube

                if (in_inner == 1) {
                    // Check for zero-crossing
                    int is_interface = 0;
                    for (int d = 0; d < Neighbors::num_neighbors; ++d) {
                        const int ni = i + Neighbors::offsets[d][0];
                        const int nj = j + Neighbors::offsets[d][1];
                        const int nk = k + Neighbors::offsets[d][2];
                        const amrex::IntVect nbr(AMREX_D_DECL(ni, nj, nk));

                        if (!ghost_bx.contains(nbr)) continue;

                        const int nb_nbr = nb_arr(nbr, 0);
                        const bool nbr_in_tube = (amrex::Math::abs(nb_nbr) != OutsideTube);
                        const Set::Scalar LSnbr = ls_arr(nbr, 0);
                        const bool sign_change = (LS * LSnbr <= 0.0);

                        is_interface = is_interface || (nbr_in_tube && sign_change);
                    }

                    if (is_interface) {
                        nb_val = Interface;
                        zero_array(coord) = ls_id;
                    }
                }
            }

            nb_arr(coord) = nb_val;
        });
    }

    // Fill boundary
    narrowband.FillBoundary();
    ls.FillBoundary();
}

std::pair<amrex::BoxArray, amrex::DistributionMapping>
NarrowBandLevelset::BuildTileBoxes(const amrex::iMultiFab& narrowband_mask,
                                    const amrex::Box& physical_domain,
                                    const int number_of_ghost_cells)
{
    // Parameters
    constexpr int min_tile_size = 8;
    constexpr int max_tile_size = 32;
    constexpr int interface_flag = NarrowBandTubeType::Interface; // value for Interface
    constexpr int inside_flag = NarrowBandTubeType::InsideTube;    // value for InsideTube

    // Step 1: Collect all narrowband cell indices and compute bounds
    amrex::Vector<amrex::IntVect> narrowband_cells;
    amrex::IntVect min_iv(AMREX_D_DECL(INT_MAX, INT_MAX, INT_MAX));
    amrex::IntVect max_iv(AMREX_D_DECL(INT_MIN, INT_MIN, INT_MIN));

    for (amrex::MFIter mfi(narrowband_mask); mfi.isValid(); ++mfi) {
        const auto& mask_arr = narrowband_mask.const_array(mfi);
        const amrex::Box& bx = mfi.validbox();

        AMREX_LOOP_3D(bx, i, j, k, {
            int val = std::abs(mask_arr(i, j, k));
            if (val <= inside_flag) {
                amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                narrowband_cells.push_back(iv);

                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    min_iv[d] = std::min(min_iv[d], iv[d]);
                    max_iv[d] = std::max(max_iv[d], iv[d]);
                }
            }
        });
    }

    // Check if narrowband is empty
    if (narrowband_cells.empty()) {
        return {{}, {}};
    }

    for (const auto& iv : narrowband_cells) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            min_iv[d] = std::min(min_iv[d], iv[d]);
            max_iv[d] = std::max(max_iv[d], iv[d]);
        }
    }

    // clip  geom 
    amrex::Box Tube_domain(min_iv, max_iv);
    Tube_domain &= physical_domain; 

    // ---------------------------------------------------------------------
    // Step 2: Adaptively select tile_size and build tile boxes
    // ---------------------------------------------------------------------
    std::map<amrex::IntVect, amrex::Box> tile_boxes;
    amrex::IntVect tile_size;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        tile_size[d] = std::max(8, std::min(32, Tube_domain.length(d) / 8));
    }

    // Create tile boxes from narrowband cells
    for (const auto& iv : narrowband_cells)
    {
        // Compute tile index
        amrex::IntVect tile_idx = iv;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            tile_idx[d] /= tile_size[d];
        }

        // Anchor tile box
        amrex::IntVect small = tile_idx * tile_size;
        amrex::IntVect big = small + tile_size - 1;
        amrex::Box tile_box(small, big);

        // clip
        tile_box.grow(number_of_ghost_cells);
        tile_box &= Tube_domain;

        if (tile_box.ok()) {
            auto it = tile_boxes.find(tile_idx);
            if (it == tile_boxes.end()) {
                tile_boxes[tile_idx] = tile_box;
            } else {
                it->second = amrex::minBox(it->second, tile_box);
            }
        }
    } 

    // Convert to BoxList and build BoxArray
    amrex::BoxList bl;
    for (const auto& pair : tile_boxes) {
        bl.push_back(pair.second);
    }
    bl.simplify();

    amrex::BoxArray tube_ba(bl);
    amrex::DistributionMapping tube_dm(tube_ba);

    return {tube_ba, tube_dm};
}

void NarrowBandLevelset::UpdateNarrowbandTubeandMapping(int lev, int ls_id) {
        // Define geometry constants
    const Set::Scalar* DX = geom[lev].CellSize();
    const Set::Scalar min_DX = *std::min_element(DX, DX + AMREX_SPACEDIM);
    const Set::Scalar inner_tube_width = inner_narrow_band_width * min_DX;
    const Set::Scalar tube_width = narrow_band_width * min_DX;

    // Save level_set data structure for reference
    auto& ls_data = level_sets[lev][ls_id];

    // Define MultiFab properties of old band to create temporary (i)MulitFab objects
    const amrex::BoxArray& ba = ls_data.narrowband_ba;
    const amrex::DistributionMapping& dm = ls_data.narrowband_dm;

    // Create temp multifab objects for parallel copy before static casting to int 
    // to match alamo Set::Scalar types
    amrex::MultiFab tmp_cpt_mf(ba, dm, 1, number_of_ghost_cells);
    amrex::MultiFab tmp_Zerols_mf(ba, dm, 1, number_of_ghost_cells);
    amrex::MultiFab tmp_Tube_mf(ba, dm, 1, number_of_ghost_cells);
    tmp_cpt_mf.ParallelCopy(*cpt_mf[lev], 0, 0, 1, 0, 0);
    tmp_Zerols_mf.ParallelCopy(*Zerols_mf[lev], 0, 0, 1, 0, 0);
    tmp_Tube_mf.ParallelCopy(*Tube_mf[lev], ls_id, 0, 1, 0, 0);
    
    // Create narrowband imultifab variables
    amrex::iMultiFab cpt(ba, dm, 1, number_of_ghost_cells);
    amrex::iMultiFab zerols(ba, dm, 1, number_of_ghost_cells);
    amrex::iMultiFab narrowband(ba, dm, 1, number_of_ghost_cells);

    // Run parallel for to static cast Set::Scalar to int
    for (amrex::MFIter mfi(tmp_Tube_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& valid_box = mfi.validbox();

        auto const& CPT_mf_arr  = tmp_cpt_mf.const_array(mfi);
        auto const& Zero_mf_arr = tmp_Zerols_mf.const_array(mfi);
        auto const& Tube_mf_arr = tmp_Tube_mf.const_array(mfi);

        auto const& CPT_imf_arr  = cpt.array(mfi);
        auto const& Zero_imf_arr = zerols.array(mfi);
        auto const& Tube_imf_arr = narrowband.array(mfi);

        // Convert and copy values
        amrex::ParallelFor(valid_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            CPT_imf_arr(i,j,k)  = static_cast<int>(CPT_mf_arr(i,j,k)); 
            Zero_imf_arr(i,j,k) = static_cast<int>(Zero_mf_arr(i,j,k)); 
            Tube_imf_arr(i,j,k) = static_cast<int>(Tube_mf_arr(i,j,k));
            if(CPT_imf_arr(i,j,k) == ls_id) CPT_imf_arr(i,j,k) = -1;
            if(Zero_imf_arr(i,j,k) == ls_id) Zero_imf_arr(i,j,k) = -1;
        });
    }

    // Create temporary levelset variables
    auto& ls = *ls_data.Tube_ls;
    auto& ls_old = *ls_data.Tube_ls_old;

    // Update narrowband mask
    BuildNarrowbandMask(ls, narrowband, zerols, cpt, ls_id, tube_width, inner_tube_width);

    // Copy back ls
    ls_mf[lev]->ParallelCopy(ls, 0, ls_id, 1, 0, 0);

    // Update narrowband boxes
    amrex::Box physical_domain = geom[lev].Domain();
    amrex::BoxArray Tube_ba = ls_data.narrowband_ba;
    if (ls_data.reconstruct_tube){
        const auto [Tube_ba, Tube_dm] = BuildTileBoxes(narrowband, physical_domain, number_of_ghost_cells);
        ls_data.reconstruct_tube = false;
    }
    
    // Update ls_data structure
    bool has_narrowband = true;
    if (!Tube_ba.ok()) has_narrowband = false;
    ls_data.has_narrowband = has_narrowband;

    if (!has_narrowband) {
        // reset
        ls_data.Tube_domain = {};
        ls_data.narrowband_ba.clear();
        ls_data.narrowband_dm = {};

        ls_data.Tube_zerols.reset();
        ls_data.Tube.reset();
        ls_data.Tube_ls.reset();
        ls_data.Tube_ls_old.reset();
        ls_data.Tube_velocity.reset();

    } else{
        auto Tube        = std::make_unique<amrex::iMultiFab>(ba, dm, 1, number_of_ghost_cells);
        auto Tube_cpt    = std::make_unique<amrex::iMultiFab>(ba, dm, 1, number_of_ghost_cells);
        auto Tube_zerols = std::make_unique<amrex::iMultiFab>(ba, dm, 1, number_of_ghost_cells);
        auto Tube_ls     = std::make_unique<amrex::MultiFab>(ba, dm, 1, number_of_ghost_cells);
        auto Tube_ls_old = std::make_unique<amrex::MultiFab>(ba, dm, 1, number_of_ghost_cells);
        auto velocity    = std::make_unique<amrex::MultiFab>(ba, dm, AMREX_SPACEDIM, number_of_ghost_cells);

        Tube->ParallelCopy(narrowband, 0, 0, 1, 0, 0);
        Tube_cpt->ParallelCopy(cpt, 0, 0, 1, 0, 0);
        Tube_zerols->ParallelCopy(zerols, 0, 0, 1, 0, 0);

        Tube_ls->ParallelCopy(ls, 0, 0, 1, number_of_ghost_cells, number_of_ghost_cells);
        Tube_ls->FillBoundary();
        FillTubeNeumannBC(*Tube_ls, Tube_ba.minimalBox(), physical_domain, number_of_ghost_cells);
        Tube_ls->FillBoundary();

        Tube_ls_old->ParallelCopy(ls_old, 0, 0, 1, 0, 0);

        // Update ls_data
        ls_data.Tube_domain = Tube->boxArray().minimalBox();
        ls_data.narrowband_ba = ba;
        ls_data.narrowband_dm = dm;

        ls_data.Tube_zerols = std::move(Tube_zerols);
        ls_data.Tube = std::move(Tube);
        ls_data.Tube_ls = std::move(Tube_ls);
        ls_data.Tube_ls_old = std::move(Tube_ls_old);
        ls_data.Tube_velocity = std::move(velocity);
    }

    // Copy iMFs to MF
    for (amrex::MFIter mfi(tmp_Tube_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& valid_box = mfi.validbox();

        auto const& CPT_mf_arr  = tmp_cpt_mf.array(mfi);
        auto const& Zero_mf_arr = tmp_Zerols_mf.array(mfi);
        auto const& Tube_mf_arr = tmp_Tube_mf.array(mfi);

        auto const& CPT_imf_arr  = cpt.const_array(mfi);
        auto const& Zero_imf_arr = zerols.const_array(mfi);
        auto const& Tube_imf_arr = narrowband.const_array(mfi);

        // Convert and copy values
        amrex::ParallelFor(valid_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            CPT_mf_arr(i,j,k)  = static_cast<Set::Scalar>(CPT_imf_arr(i,j,k)); 
            Zero_mf_arr(i,j,k) = static_cast<Set::Scalar>(Zero_imf_arr(i,j,k)); 
            Tube_mf_arr(i,j,k) = static_cast<Set::Scalar>(Tube_imf_arr(i,j,k));
        });
    }

    // Copy iMFs to full domain
    cpt_mf[lev]->ParallelCopy(tmp_cpt_mf, 0, 0, 1, 0, 0);
    Zerols_mf[lev]->ParallelCopy(tmp_Zerols_mf, 0, 0, 1, 0, 0);
    Tube_mf[lev]->ParallelCopy(tmp_Tube_mf, 0, ls_id, 1, 0, 0);
}

/*// UpdateNarrowband: replicated Fortran LSTubeInfo logic
void NarrowBandLevelset::UpdateNarrowbandTubeandMapping(int lev, int ls_id) {
    // Define geometry constants
    const Set::Scalar* DX = geom[lev].CellSize();
    const Set::Scalar min_DX = *std::min_element(DX, DX + AMREX_SPACEDIM);
    const Set::Scalar inner_tube_width = inner_narrow_band_width * min_DX;
    const Set::Scalar tube_width = narrow_band_width * min_DX;

    // Define narrowband constant values
    const int Interface   = NarrowBandTubeType::Interface;
    const int InsideTube  = NarrowBandTubeType::InsideTube;
    const int OutsideTube = NarrowBandTubeType::OutsideTube;

    // Save level_set data structure for reference
    auto& ls_data = level_sets[ls_id];

    // Define MultiFab properties of old band to create temporary (i)MulitFab objects
    const amrex::BoxArray& ba = ls_data.narrowband_ba;
    const amrex::DistributionMapping& dm = ls_data.narrowband_dm;

    // Create temp multifab objects for parallel copy before static casting to int 
    // to match alamo Set::Scalar types
    amrex::MultiFab tmp_cpt_mf(ba, dm, 1, number_of_ghost_cells);
    amrex::MultiFab tmp_Zerols_mf(ba, dm, 1, number_of_ghost_cells);
    amrex::MultiFab tmp_Tube_mf(ba, dm, 1, number_of_ghost_cells);
    tmp_cpt_mf.ParallelCopy(*cpt_mf[lev], 0, 0, 1, 0, 0);
    tmp_Zerols_mf.ParallelCopy(*Zerols_mf[lev], 0, 0, 1, 0, 0);
    tmp_Tube_mf.ParallelCopy(*Tube_mf[lev], ls_id, 0, 1, 0, 0);
    
    // Create narrowband imultifab variables
    amrex::iMultiFab cpt(ba, dm, 1, number_of_ghost_cells);
    amrex::iMultiFab zerols(ba, dm, 1, number_of_ghost_cells);
    amrex::iMultiFab narrowband(ba, dm, 1, number_of_ghost_cells);

    // Run parallel for to static cast Set::Scalar to int
    for (amrex::MFIter mfi(tmp_Tube_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& valid_box = mfi.validbox();

        auto const& CPT_mf_arr  = tmp_cpt_mf.const_array(mfi);
        auto const& Zero_mf_arr = tmp_Zerols_mf.const_array(mfi);
        auto const& Tube_mf_arr = tmp_Tube_mf.const_array(mfi);

        auto const& CPT_imf_arr  = cpt.array(mfi);
        auto const& Zero_imf_arr = zerols.array(mfi);
        auto const& Tube_imf_arr = narrowband.array(mfi);

        // Convert and copy values
        amrex::ParallelFor(valid_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            CPT_imf_arr(i,j,k)  = static_cast<int>(CPT_mf_arr(i,j,k)); 
            Zero_imf_arr(i,j,k) = static_cast<int>(Zero_mf_arr(i,j,k)); 
            Tube_imf_arr(i,j,k) = static_cast<int>(Tube_mf_arr(i,j,k));
        });
    }

    // Create temporary levelset variables
    auto const& temp_ls = *ls_data.Tube_ls;
    auto const& temp_ls_old = *ls_data.Tube_ls_old;
    amrex::MultiFab ls(ba,dm, 1, number_of_ghost_cells);
    amrex::MultiFab ls_old(ba,dm, 1, number_of_ghost_cells);
    ls.ParallelCopy(temp_ls, 0, 0, 1, 0, 0);
    ls_old.ParallelCopy(temp_ls_old, 0, 0, 1, 0, 0);

    // Perform initial loop through narrowband cells to identify interface and tube cells
    for (amrex::MFIter mfi(ls, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        // Define valid and ghost box to update all ghost cells
        const amrex::Box& ghost_bx = mfi.growntilebox(number_of_ghost_cells);

        // Define array views
        const auto& cpt_arr = cpt.array(mfi);
        const auto& zero_arr = zerols.array(mfi);
        const auto& nb_arr = narrowband.array(mfi);
        const auto& ls_arr = ls.array(mfi); 

        // Loop through ghost box
        ParallelFor(ghost_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const amrex::IntVect coord(AMREX_D_DECL(i, j, k)); 

            // Clamp levelset value and define abs|ls|
            Set::Scalar LS = ls_arr(coord);
            LS = std::clamp(LS, -tube_width, tube_width);
            const Set::Scalar abs_LS = std::abs(LS);
            ls_arr(coord) = LS;

            // Assign cpt flag
            if (LS < 0.0) cpt_arr(coord) = ls_id;

            // Check if cell is within narrowband
            if (abs_LS < tube_width) {
                // If not within inner_tube_width, skip Interface check
                // and set to InsideTube
                if (abs_LS >= inner_tube_width) {
                    nb_arr(i, j, k) = InsideTube;
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
                    if (std::abs(nb_nbr) == OutsideTube) continue;

                    // Get the levelset value of the nbr
                    Set::Scalar LSnbr = ls_arr(nbr);
            
                    // If LS * LSnbr is negative, mark as Interface and break loop
                    if (LS * LSnbr <= 0.0){
                        // Mark zerols with ls_id and cell to Interface
                        zero_arr(coord) = ls_id;
                        nb_arr(coord) = Interface;
                        break;
                    }
                }

                // If cell is not an Interface cell, mark as InsideTube
                if (nb_arr(coord) != Interface) nb_arr(coord) = InsideTube;
            } 
            // If cell is not within Tube, mark as OutsideTube
            else {
                nb_arr(coord) = OutsideTube;
            }
        });
    }

    // Fill Boundary for narrowband
    narrowband.FillBoundary();
    ls.FillBoundary();

    // ---------------------------------------------------------------------
    // Step 1: Compute bounding box (Tube_geom) of the narrowband region
    // ---------------------------------------------------------------------
    // Collect narrowband cell locations
    amrex::Vector<amrex::IntVect> narrowband_cells;

    for (amrex::MFIter mfi(narrowband); mfi.isValid(); ++mfi) {
        const auto& mask_arr = narrowband.const_array(mfi);
        const amrex::Box& bx = mfi.validbox();

        AMREX_LOOP_3D(bx, i, j, k, {
            if (std::abs(mask_arr(i, j, k)) <= InsideTube) {
                narrowband_cells.push_back(amrex::IntVect(AMREX_D_DECL(i, j, k)));
            }
        });
    }

    // Check if narrowband cells exist and save flag
    bool has_narrowband = false;
    if (!narrowband_cells.empty()) has_narrowband = true;
    ls_data.has_narrowband = has_narrowband;
    if (has_narrowband){
        amrex::IntVect min_iv(AMREX_D_DECL(INT_MAX, INT_MAX, INT_MAX));
        amrex::IntVect max_iv(AMREX_D_DECL(INT_MIN, INT_MIN, INT_MIN));

        for (const auto& iv : narrowband_cells) {
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                min_iv[d] = std::min(min_iv[d], iv[d]);
                max_iv[d] = std::max(max_iv[d], iv[d]);
            }
        }

        // clip  geom 
        amrex::Box physical_domain = geom[lev].Domain();
        amrex::Box Tube_domain(min_iv, max_iv);
        Tube_domain &= physical_domain; 

        // ---------------------------------------------------------------------
        // Step 2: Adaptively select tile_size and build tile boxes
        // ---------------------------------------------------------------------
        std::map<amrex::IntVect, amrex::Box> tile_boxes;
        amrex::IntVect tile_size;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            tile_size[d] = std::max(8, std::min(32, Tube_domain.length(d) / 8));
        }

        // Create tile boxes from narrowband cells
        for (const auto& iv : narrowband_cells)
        {
            // Compute tile index
            amrex::IntVect tile_idx = iv;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                tile_idx[d] /= tile_size[d];
            }

            // Anchor tile box
            amrex::IntVect small = tile_idx * tile_size;
            amrex::IntVect big = small + tile_size - 1;
            amrex::Box tile_box(small, big);

            // clip
            tile_box.grow(number_of_ghost_cells);
            tile_box &= Tube_domain;

            if (tile_box.ok()) {
                auto it = tile_boxes.find(tile_idx);
                if (it == tile_boxes.end()) {
                    tile_boxes[tile_idx] = tile_box;
                } else {
                    it->second = amrex::minBox(it->second, tile_box);
                }
            }
        } 

        // Convert to BoxList and build BoxArray
        amrex::BoxList bl;
        for (const auto& pair : tile_boxes) {
            bl.push_back(pair.second);
        }
        bl.simplify();

        amrex::BoxArray tube_ba(bl);
        amrex::DistributionMapping tube_dm(tube_ba);

        // Store variables to ls_data
        ls_data.Tube_domain = Tube_domain;
        ls_data.narrowband_ba = tube_ba;
        ls_data.narrowband_dm = tube_dm;

        // Update ptrs
        // Zerols
        ls_data.Tube_zerols = std::make_unique<amrex::iMultiFab>(ba, dm, 1, number_of_ghost_cells);
        auto& Tube_zerols = *ls_data.Tube_zerols;
        Tube_zerols.ParallelCopy(zerols, 0, 0, 1, number_of_ghost_cells, number_of_ghost_cells);

        // Narrowband
        ls_data.Tube = std::make_unique<amrex::iMultiFab>(ba, dm, 1, number_of_ghost_cells);
        auto& Tube = *ls_data.Tube;
        Tube.ParallelCopy(narrowband, 0, 0, 1, number_of_ghost_cells, number_of_ghost_cells);

        // Levelset fields
        ls_data.Tube_ls = std::make_unique<amrex::MultiFab>(ba, dm, 1, number_of_ghost_cells);
        auto& Tube_ls = *ls_data.Tube_ls;
        Tube_ls.ParallelCopy(ls, 0, 0, 1, number_of_ghost_cells, number_of_ghost_cells);
        // Apply Nuemann BCs to ls field
        FillTubeNeumannBC(Tube_ls, Tube_domain, physical_domain, number_of_ghost_cells);
        Tube_ls.FillBoundary();

        ls_data.Tube_ls_old = std::make_unique<amrex::MultiFab>(ba, dm, 1, number_of_ghost_cells);
        auto& Tube_ls_old = *ls_data.Tube_ls_old;
        Tube_ls_old.ParallelCopy(ls_old, 0, 0, 1, number_of_ghost_cells, number_of_ghost_cells);

        // Velocity
        ls_data.Tube_velocity = std::make_unique<amrex::MultiFab>(ba, dm, AMREX_SPACEDIM, number_of_ghost_cells);
    } else {
        // Create an empty BoxArray and DistributionMapping
        ls_data.Tube_domain = amrex::Box();
        ls_data.narrowband_ba.clear();
        ls_data.narrowband_dm = amrex::DistributionMapping();  // Default constructed empty mapping

        // Reset smart pointers
        ls_data.Tube_zerols.reset();
        ls_data.Tube.reset();
        ls_data.Tube_ls.reset();
        ls_data.Tube_ls_old.reset();
        ls_data.Tube_velocity.reset();
    }

    // Copy iMFs to MF
    for (amrex::MFIter mfi(tmp_Tube_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& valid_box = mfi.validbox();

        auto const& CPT_mf_arr  = tmp_cpt_mf.array(mfi);
        auto const& Zero_mf_arr = tmp_Zerols_mf.array(mfi);
        auto const& Tube_mf_arr = tmp_Tube_mf.array(mfi);

        auto const& CPT_imf_arr  = cpt.const_array(mfi);
        auto const& Zero_imf_arr = zerols.const_array(mfi);
        auto const& Tube_imf_arr = narrowband.const_array(mfi);

        // Convert and copy values
        amrex::ParallelFor(valid_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            CPT_mf_arr(i,j,k)  = static_cast<Set::Scalar>(CPT_imf_arr(i,j,k)); 
            Zero_mf_arr(i,j,k) = static_cast<Set::Scalar>(Zero_imf_arr(i,j,k)); 
            Tube_mf_arr(i,j,k) = static_cast<Set::Scalar>(Tube_imf_arr(i,j,k));
        });
    }

    // Copy iMFs to full domain
    cpt_mf[lev]->ParallelCopy(tmp_cpt_mf, 0, 0, 1, 0, 0);
    Zerols_mf[lev]->ParallelCopy(tmp_Zerols_mf, 0, 0, 1, 0, 0);
    Tube_mf[lev]->ParallelCopy(tmp_Tube_mf, 0, ls_id, 1, 0, 0);

    // Copy LS to full domain
    auto& final_ls = *ls_data.Tube_ls;
    ls_mf[lev]->ParallelCopy(final_ls, 0, ls_id, 1, 0, 0);
}*/

void NarrowBandLevelset::ComputeGeometryQuantities(int lev, int ls_id) {
    // Define DX for gradient
    const Set::Scalar* DX = geom[lev].CellSize();
    amrex::Box physical_domain = geom[lev].Domain();
    //const amrex::IntVect debug(AMREX_D_DECL(14,49,0));

    // Define OutsideNarrowband constant
    const int InsideTube = NarrowBandTubeType::InsideTube;

    // Save level_set data structure for reference
    auto& ls_data = level_sets[lev][ls_id];

    // Get narrowband variables from data structure
    const amrex::BoxArray& ba = ls_data.narrowband_ba;
    const amrex::DistributionMapping dm = ls_data.narrowband_dm;
    auto const& narrowband = *ls_data.Tube;
    auto& ls = *ls_data.Tube_ls;
    auto& normal_mf = *ls_data.normal_mf[lev];

    // Define temporary normal
    amrex::MultiFab normal(ba, dm, AMREX_SPACEDIM, 1);
    normal.ParallelCopy(normal_mf, 0, 0, AMREX_SPACEDIM, 1, 1);

    // Fill ghosts
    ls.FillBoundary();
    normal.FillBoundary();

    // Iterate over only the narrowband cells
    for (amrex::MFIter mfi(ls, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& valid_bx = mfi.validbox();

        auto const& nb_arr = narrowband.const_array(mfi);
        auto const& ls_arr = ls.const_array(mfi);
        auto const& normal_arr = normal.array(mfi);
        
        amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            amrex::IntVect coord(AMREX_D_DECL(i,j,k));
            if (std::abs(nb_arr(coord, 0)) > InsideTube) return; // Skip non-narrowband cells

            // Utilize Numeric::Gradient for simplified computations
            std::array<Numeric::StencilType, AMREX_SPACEDIM> stencil = Numeric::GetStencil(i, j, k, physical_domain);
            Set::Vector gradient = Numeric::Gradient(ls_arr, i, j, k, 0, DX, stencil);

            Set::Scalar grad_mag = gradient.norm();

            //if (debug==coord) amrex::Print() << ls_arr(i-1, j, k) << std::endl;

            for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
                normal_arr(i, j, k, dim) = (grad_mag > 0.0) ? (gradient[dim] / grad_mag) : 0.0;
                //if (debug==coord) amrex::Print() << "normal" << dim << ": " << normal_arr(debug, dim) << std::endl;
            }
        });
    }  

    // Fill ghost cells
    normal.FillBoundary();

    // Copy normal back for plotting
    normal_mf.ParallelCopy(normal, 0, 0, AMREX_SPACEDIM, 0, 0);

    // Get temp variable for curvature
    amrex::MultiFab curve(ba, dm, 1, 0); 
    curve.ParallelCopy(*curvature_mf[lev], ls_id, 0, 1, 0, 0);

    for (amrex::MFIter mfi(normal, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& valid_bx = mfi.validbox();

        auto const& nb_arr     = narrowband.const_array(mfi);
        auto const& normal_arr = normal.const_array(mfi); 
        auto const& curvature  = curve.array(mfi);

        amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            // skip if not within narowband
            amrex::IntVect coord(AMREX_D_DECL(i,j,k));
            if (std::abs(nb_arr(coord, 0)) > InsideTube) return; // Skip non-narrowband cells
            
            // Get curvature as divergence of normal
            std::array<Numeric::StencilType, AMREX_SPACEDIM> stencil = Numeric::GetStencil(i, j, k, physical_domain);
            curvature(i, j, k) = Numeric::Divergence(normal_arr, i, j, k, 0, DX, stencil);
        });
    }  

    // Copy computed normals back to the full normal MultiFab 
    curvature_mf[lev]->ParallelCopy(curve, 0, ls_id, 1, 0, 0);
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
    
            for (amrex::MFIter mfi(*velocity_mf[lev], false); mfi.isValid(); ++mfi) {
                const amrex::Box& bx = mfi.tilebox();  // Iterate over tiles in the multifab
                auto const& velocity_arr = velocity_mf.Patch(lev, mfi);  // Velocity field
    
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
    /*for (amrex::MFIter mfi(*ls_mf[lev]); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        const amrex::FArrayBox& ls_fab     = (*ls_mf[lev])[mfi];
        const amrex::FArrayBox& ls_old_fab = (*ls_old_mf[lev])[mfi];
    
        if (bx.contains(amrex::IntVect(20, 49))) {
            amrex::Print() << "Before swap:\n";
            amrex::Print() << "LS(20,49): " << ls_fab(amrex::IntVect(20, 49), 0) << "\n";
            amrex::Print() << "LS_OLD(20,49): " << ls_old_fab(amrex::IntVect(20, 49), 0) << "\n";
        }
    }*/

    std::swap(*ls_mf[lev], *ls_old_mf[lev]);

    /*for (amrex::MFIter mfi(*ls_mf[lev]); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        const amrex::FArrayBox& ls_fab     = (*ls_mf[lev])[mfi];
        const amrex::FArrayBox& ls_old_fab = (*ls_old_mf[lev])[mfi];
    
        if (bx.contains(amrex::IntVect(20, 49))) {
            amrex::Print() << "After swap:\n";
            amrex::Print() << "LS(20,49): " << ls_fab(amrex::IntVect(20, 49), 0) << "\n";
            amrex::Print() << "LS_OLD(20,49): " << ls_old_fab(amrex::IntVect(20, 49), 0) << "\n";
        }
    }*/
    
    // Advect
    Advect(lev, time, dt);

    /*for (amrex::MFIter mfi(*ls_mf[lev]); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        const amrex::FArrayBox& ls_fab     = (*ls_mf[lev])[mfi];
    
        if (bx.contains(amrex::IntVect(46, 46))) {
            amrex::Print() << "After Advect:\n";
            amrex::Print() << "LS(46,46): " << ls_fab(amrex::IntVect(46, 46), 0) << "\n";
        }
    }*/

    /*for (amrex::MFIter mfi(*ls_mf[lev]); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        const auto& arr = ls_mf.Patch(lev, mfi);

        amrex::ParallelFor(bx, ls_mf[lev]->nComp(), [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
            if (amrex::isnan(arr(i,j,k,n))) {
                printf("NaN at (%d, %d, %d), component %d\n", i, j, k, n);
            }
        });
    }*/

    // Apply Boundary conditions after advection
    //Integrator::ApplyPatch(lev, time, ls_mf, *ls_mf[lev], *bc_ls, 0);

    // Reinitialize every 5 timesteps
    if (current_timestep % 100 == 0){
        for (int ls_id =0; ls_id < number_of_components; ls_id++){
            Reinitialize(lev, ls_id);
            // Activate flag to update box
        }

        // Apply Boundary conditions after reinitialization
        Integrator::ApplyPatch(lev, time, ls_mf, *ls_mf[lev], *bc_ls, 0);
    }

    // Update the narrowband
    for (int ls_id = 0; ls_id < number_of_components; ls_id++){
        // Update narrowband info
        UpdateNarrowbandTubeandMapping(lev, ls_id);

        // Compute geometries
        //ComputeGeometryQuantities(lev, ls_id);

        // Save Tube_imf to Tube_imf
        //CopyTubeIMFtoMF(lev, ls_id);
    }

    /*for (amrex::MFIter mfi(*ls_mf[lev]); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        const amrex::FArrayBox& ls_fab     = (*ls_mf[lev])[mfi];
    
        if (bx.contains(amrex::IntVect(46, 46))) {
            amrex::Print() << "After Tube Update:\n";
            amrex::Print() << "LS(46,46): " << ls_fab(amrex::IntVect(46, 46), 0) << "\n";
        }
    }

    // Copy the zero
    //CopyZerolsAndCPTIMFtoMF(lev);

    for (amrex::MFIter mfi(*ls_mf[lev]); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        const amrex::FArrayBox& ls_fab     = (*ls_mf[lev])[mfi];
        const amrex::FArrayBox& ls_old_fab = (*ls_old_mf[lev])[mfi];
    
        if (bx.contains(amrex::IntVect(20, 49))) {
            amrex::Print() << "At end of Advance:\n";
            amrex::Print() << "LS(20,49): " << ls_fab(amrex::IntVect(20, 49), 0) << "\n";
            amrex::Print() << "LS_OLD(20,49): " << ls_old_fab(amrex::IntVect(20, 49), 0) << "\n";
        }
    }*/
}

void NarrowBandLevelset::Advect(int lev, Set::Scalar time, Set::Scalar dt){
    // Perform the swap
    //amrex::MultiFab::Copy(*ls_old_mf[lev], *ls_mf[lev], 0, 0, 1, ls_mf[lev]->nGrow());

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
    for (int ils = 0; ils < number_of_components; ils++){
        auto const& ls_data = level_sets[lev][ils];
        auto& velocity = *ls_data.Tube_velocity;
        velocity.ParallelCopy(*velocity_mf[lev], 0, 0, AMREX_SPACEDIM, number_of_ghost_cells, number_of_ghost_cells);
    }
}

/*void NarrowBandLevelset::UpdateInterfaceVelocity(int lev){
    for (int ils=0; ils < number_of_components; ils++){
        for (amrex::MFIter mfi(*velocity_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.tilebox();
            auto const& vel_arr = velocity_mf.Patch(lev,mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                for (int dim=0; dim < AMREX_SPACEDIM; ++dim){
                    vel_arr(i,j,k,dim) = vel_arr(i,j,k,dim);
                }
            });
        }
    }
}*/

void NarrowBandLevelset::Reinitialize(int lev, int ls_id) {
    // Define reinitialization constants
    const Set::Scalar reinit_tolerance = 1e-3;
    const int max_iterations = 50;

    // Define geometry constants
    const amrex::Box& domain_box = geom[lev].Domain();
    const Set::Scalar* DX    = geom[lev].CellSize();
    const Set::Scalar min_DX = *std::min_element(DX, DX + AMREX_SPACEDIM);
    const Set::Scalar tau    = cflReinit * min_DX;

    // Define Narrowband constants
    const Set::Scalar INNERTUBE = inner_narrow_band_width * min_DX; // used for testing sign changes
    //const int Interface       = NarrowBandTubeType::Interface;
    const int OutsideTube = NarrowBandTubeType::OutsideTube;
    
    // Get narrowband geometry
    amrex::BoxArray ba            = level_sets[lev][ls_id].narrowband_ba;
    amrex::DistributionMapping dm = level_sets[lev][ls_id].narrowband_dm;
    
    // Allocate temporary (i)multifabs
    amrex::iMultiFab zerols(ba, dm, 1, 1); // Used to update zero after advection
    amrex::iMultiFab Narrowband(ba, dm, 1, 0);
    amrex::MultiFab LS(ba, dm, 1, 1);
    amrex::MultiFab LS_old(ba, dm, 1, 1);
    amrex::MultiFab error_mf(ba, dm, 1, 0); // Used to compute Linf norm

    // Initialize (i)multifabs
    // NOTE FOR FUTURE: WILL NEED TO FIND A WAY TO PROPERLY COPY Zerols_imf SO THAT
    // OTHER ls_ids ARE NOT OVERWRITTEN WHEN COPIED BACK. SAME FOR ANY OTHER FUNCTION 
    // THAT USES Zerols_imf
    //zerols.setVal(-1);
    //zerols.ParallelCopy(*Zerols_imf, 0, 0, 1, 0, 0);
    //ZeroTheZeros(ls_id, zerols, false, nullptr);
    Narrowband.ParallelCopy(*level_sets[lev][ls_id].Tube, 0, 0, 1, 0, 0);
    LS.ParallelCopy(*ls_mf[lev], ls_id, 0, 1, 1, 1);
    LS_old.ParallelCopy(*ls_mf[lev], ls_id, 0, 1, 1, 1);
    error_mf.setVal(0.0);

    // Fill ghost cells
    zerols.FillBoundary();
    LS.FillBoundary();

    // Redefine the interface cells using exact logic as Fortran
    for (amrex::MFIter mfi(Narrowband, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& ghost_bx = mfi.growntilebox(1); 
        const auto& ls_arr = LS.const_array(mfi);
        const auto& zerols_arr = zerols.array(mfi);

        amrex::ParallelFor(ghost_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // Get current coord and LS value
            const amrex::IntVect coord(AMREX_D_DECL(i, j, k));
            const Set::Scalar LS = ls_arr(coord);

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
        });
    }

    // Fill zerols ghost cells
    zerols.FillBoundary();
    
    // Perform first order PDE Reinitialization scheme
    for (int iter = 0; iter < max_iterations; iter++){
        //amrex::Print() << "iter: " << iter << std::endl;
        
        // Switch LS with LS_old
        std::swap(LS, LS_old);

        // Fill ghost cells for LS_old
        LS_old.FillBoundary();

        // Perform main PDE loop
        for (amrex::MFIter mfi(LS, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Define upper and lower bounds of valid_box for stencil boundaries
            const amrex::Box& valid_bx = mfi.validbox();
            const amrex::Dim3 lo = amrex::lbound(valid_bx);
            const amrex::Dim3 hi = amrex::ubound(valid_bx);

            // Get arrays
            const auto& ls_arr = LS.array(mfi);
            const auto& ls_old_arr = LS_old.const_array(mfi);
            const auto& nb_arr = Narrowband.const_array(mfi);
            const auto& zerols_arr = zerols.const_array(mfi);
            const auto& error_arr = error_mf.array(mfi);

            amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::IntVect coord(AMREX_D_DECL(i, j, k));
                
                // Skip interface and outerband cells
                const int abs_nb = std::abs(nb_arr(coord));
                const int zero_cell = zerols_arr(coord);
                if (abs_nb == OutsideTube || zero_cell == ls_id) return;

                // Define sign of ls for upwinding
                const Set::Scalar ls = ls_old_arr(coord);
                const Set::Scalar sign_ls = (ls >= 0.0 ? 1.0 : -1.0);

                // Use first order forward/backward difference stencils
                auto stencil_lo = Numeric::GetStencil(i, j, k, valid_bx);
                auto stencil_hi = stencil_lo;

                stencil_lo[0] = Numeric::StencilType::Lo;
                stencil_hi[0] = Numeric::StencilType::Hi;

                // Compute first order differences
                // x-direction
                Set::Scalar dxm = (i > lo.x)
                    ? Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(ls_old_arr, i, j, k, 0, DX, stencil_lo)
                    : 0.0;

                Set::Scalar dxp = (i < hi.x)
                    ? Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(ls_old_arr, i, j, k, 0, DX, stencil_hi)
                    : 0.0;

                Set::Scalar gx = std::max(
                    std::pow(std::max(sign_ls * dxm, 0.0), 2),
                    std::pow(std::min(sign_ls * dxp, 0.0), 2)
                );

                // y-direction
                Set::Scalar gy = 0.0;
                #if AMREX_SPACEDIM >= 2
                stencil_lo[1] = Numeric::StencilType::Lo;
                stencil_hi[1] = Numeric::StencilType::Hi;

                Set::Scalar dym = (j > lo.y)
                    ? Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(ls_old_arr, i, j, k, 0, DX, stencil_lo)
                    : 0.0;

                Set::Scalar dyp = (j < hi.y)
                    ? Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(ls_old_arr, i, j, k, 0, DX, stencil_hi)
                    : 0.0;

                gy = std::max(
                    std::pow(std::max(sign_ls * dym, 0.0), 2),
                    std::pow(std::min(sign_ls * dyp, 0.0), 2)
                );
                #endif

                // z-direction
                Set::Scalar gz = 0.0;
                #if AMREX_SPACEDIM == 3
                stencil_lo[2] = Numeric::StencilType::Lo;
                stencil_hi[2] = Numeric::StencilType::Hi;

                Set::Scalar dzm = (k > lo.z)
                    ? Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(ls_old_arr, i, j, k, 0, DX, stencil_lo)
                    : 0.0;

                Set::Scalar dzp = (k < hi.z)
                    ? Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(ls_old_arr, i, j, k, 0, DX, stencil_hi)
                    : 0.0;

                gz = std::max(
                    std::pow(std::max(sign_ls * dzm, 0.0), 2),
                    std::pow(std::min(sign_ls * dzp, 0.0), 2)
                ); 
                #endif

                // Compute gradient and update
                Set::Scalar grad_phi = std::sqrt(gx + gy + gz);
                ls_arr(coord) = ls - tau * sign_ls * (grad_phi - 1.0);

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
    //Zerols_imf->ParallelCopy(zerols, 0, 0, 1, 0, 0);
    ls_mf[lev]->ParallelCopy(LS, 0, ls_id, 1, 1, 1);
}

void NarrowBandLevelset::TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow){
}

void NarrowBandLevelset::Regrid(int lev, Set::Scalar time) {
}

} // namespace Integrator
