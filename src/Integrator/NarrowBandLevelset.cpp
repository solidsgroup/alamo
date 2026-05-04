#include "NarrowBandLevelset.H"
#include "IC/Constant.H"
#include "IC/Expression.H"
#include "BC/Constant.H"
#include "BC/Nothing.H"
#include "BC/Expression.H"
#include "BC/LevelSetNeumann.H"
#include "Numeric/NBStencil.H"
#include "Numeric/NBFluxHandler.H"
#include "Util/NB_Util.H"

#include <AMReX_Reduce.H>

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Scalar to_scalar(int v)
{
    return static_cast<Set::Scalar>(v);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real upwind_flux_component(
    amrex::Real dE, amrex::Real dW,
#if (AMREX_SPACEDIM >= 2)
    amrex::Real dN, amrex::Real dS,
#endif
#if (AMREX_SPACEDIM == 3)
    amrex::Real dT, amrex::Real dB,
#endif
    amrex::Real Sx, 
#if (AMREX_SPACEDIM >= 2)
    amrex::Real Sy,
#endif
#if (AMREX_SPACEDIM == 3)
    amrex::Real Sz,
#endif
    amrex::Real absSx, 
#if (AMREX_SPACEDIM >= 2)
    amrex::Real absSy,
#endif
#if (AMREX_SPACEDIM == 3)
    amrex::Real absSz,
#endif
    const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>& dx)
{
    amrex::Real flux =
        (Sx - absSx)*dE/dx[0] +
        (Sx + absSx)*dW/dx[0];

#if (AMREX_SPACEDIM >= 2)
    flux += 
        (Sy - absSy)*dN/dx[1] +
        (Sy + absSy)*dS/dx[1];
#endif

#if (AMREX_SPACEDIM == 3)
    flux +=
        (Sz - absSz)*dT/dx[2] +
        (Sz + absSz)*dB/dx[2];
#endif

    return flux;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real enforce_signed_band(
    amrex::Real phi,
    int cpt,
    amrex::Real outer_cutoff)
{
    const amrex::Real target_sign = (cpt > 0) ? -1.0 : 1.0;

    // wrong sign OR outside band → clamp
    if (phi * target_sign <= 0.0 || amrex::Math::abs(phi) >= outer_cutoff)
        return target_sign * outer_cutoff;

    return phi;
}

// =====================================
// Level set object
// =====================================
struct LevelSetObject
{
    int obj_ID; // Global across all ls fields

    bool active = true;
    bool moving = true;
    bool needs_reinit = false;

    amrex::Real volume = 0.0;
    amrex::Real centroid[AMREX_SPACEDIM] = {0.0};
};

// -------------------------------
// Level set field (AMR-ready)
// -------------------------------
struct LevelSetField
{
    // Define Field data
    Set::Field<Set::Scalar> LS;
    Set::Field<Set::Scalar> LSold;

    Set::Field<Set::Scalar> LSvel;
    Set::Field<Set::Scalar> LSnormal;

    Set::Field<Set::Scalar> LSkappa;

    Set::Field<Set::Scalar> TUBE;

    // Store collection of objects in this level set field
    amrex::Vector<LevelSetObject> objects; // Number of objects per field (global)

    amrex::Vector<amrex::Vector<char>> tile_has_nb; // [AMR lev][Num tiles]

    // Store IC/BC for level set field
    std::unique_ptr<IC::IC<Set::Scalar>> ic_LS = nullptr;
    std::unique_ptr<BC::BC<Set::Scalar>> bc_LS = nullptr;

    std::unique_ptr<BC::BC<Set::Scalar>> bc_LSnormal = nullptr;

    // Store bandwidth constants
    const Set::Scalar HALFINNERTUBE = 4.0;
    const Set::Scalar TUBECUTOFF    = 6.0;

    // Save boolean to force reinitialization of entire level set field
    bool force_reinit_all = false;

    // Helper function for connected component analysis (Placeholder)
    int GetNumObjects() {
        return 1;
    }

    // Helper to resize tile data
    void resizeNBTiles(int finest_lev) {
        tile_has_nb.resize(finest_lev + 1);
    }

    // Helper function to initialize OBJ_ID (CTP) flags
    // ASSUMES SINGLE OBJECT - WILL CHANGE IN FUTURE
    void initializeObjID(int lev, const amrex::Geometry& geom) {
        using Tidx = NarrowBandData::TubeIdx;

        const Set::Scalar* DX = geom.CellSize();
        const Set::Scalar min_dx = *std::min_element(DX, DX + AMREX_SPACEDIM);

        const Set::Scalar outer_cutoff = TUBECUTOFF * min_dx;

        for (const auto& obj : objects) {
            const Set::Scalar CPT = static_cast<Set::Scalar>(obj.obj_ID);

            for (amrex::MFIter mfi(*TUBE[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const amrex::Box& tilebox = mfi.tilebox();

                Set::Patch<Set::Scalar> ls_arr = LS.Patch(lev, mfi);
                Set::Patch<Set::Scalar> tube_arr = TUBE.Patch(lev, mfi);

                amrex::ParallelFor(tilebox, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Set::Scalar& phi = ls_arr(i,j,k);
                    phi = amrex::Clamp(phi, -outer_cutoff, outer_cutoff);

                    tube_arr(i,j,k,Tidx::OBJ_ID) = phi < 0.0 ? CPT : Set::Scalar(0.0);
                });
            }
        }

        LS[lev]->FillBoundary();
        TUBE[lev]->FillBoundary();
    }

    // Helper function to update narrowband tube
    void updateTube(int lev, const amrex::Geometry& geom)
    {
        using TT   = NarrowBandData::TubeType;
        using Tidx = NarrowBandData::TubeIdx;

        // --- Grid spacing ---
        const Set::Scalar* DX = geom.CellSize();
        Set::Scalar min_dx = *std::min_element(DX, DX + AMREX_SPACEDIM);

        const Set::Scalar inner_cutoff = HALFINNERTUBE * min_dx;
        const Set::Scalar outer_cutoff = TUBECUTOFF    * min_dx;

        // --- Domain bounds ---
        const amrex::Box domain = geom.Domain();
        const auto dlo = domain.smallEnd();
        const auto dhi = domain.bigEnd();

        // Reset tile size and index
        auto& tiles = tile_has_nb[lev];
        
        const int ntiles = TUBE[lev]->local_size();
        tiles.assign(ntiles, 0);

        int tile_idx = 0;
        
        for (amrex::MFIter mfi(*TUBE[lev], amrex::TilingIfNotGPU());
            mfi.isValid(); ++mfi, ++tile_idx)
        {
            const amrex::Box& tilebox = mfi.tilebox();

            Set::Patch<const Set::Scalar> phi_arr  = LS.Patch(lev, mfi);
            Set::Patch<      Set::Scalar> tube_arr = TUBE.Patch(lev, mfi);

            // -----------------------------
            // Fused compute + reduction
            // -----------------------------
            amrex::ReduceOps<amrex::ReduceOpMax> reduce_op;
            amrex::ReduceData<int> reduce_data(reduce_op);

            using ReduceTuple = amrex::GpuTuple<int>;

            reduce_op.eval(tilebox, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                const Set::Scalar phi_c = phi_arr(i,j,k);
                const Set::Scalar absphi = amrex::Math::abs(phi_c);

                bool inside_tube = absphi < outer_cutoff;
                bool inner_tube  = absphi < inner_cutoff;

                // --- neighbors ---
                Set::Scalar xm = (i != dlo[0]) ? phi_arr(i-1,j,k) : phi_c;
                Set::Scalar xp = (i != dhi[0]) ? phi_arr(i+1,j,k) : phi_c;

    #if AMREX_SPACEDIM >= 2
                Set::Scalar ym = (j != dlo[1]) ? phi_arr(i,j-1,k) : phi_c;
                Set::Scalar yp = (j != dhi[1]) ? phi_arr(i,j+1,k) : phi_c;
    #endif
    #if AMREX_SPACEDIM == 3
                Set::Scalar zm = (k != dlo[2]) ? phi_arr(i,j,k-1) : phi_c;
                Set::Scalar zp = (k != dhi[2]) ? phi_arr(i,j,k+1) : phi_c;
    #endif

                // --- edge detection ---
                bool is_edge = false;
                if (!inside_tube) {
                    is_edge =
                        (amrex::Math::abs(xm) < outer_cutoff) ||
                        (amrex::Math::abs(xp) < outer_cutoff);
    #if AMREX_SPACEDIM >= 2
                    is_edge |=
                        (amrex::Math::abs(ym) < outer_cutoff) ||
                        (amrex::Math::abs(yp) < outer_cutoff);
    #endif
    #if AMREX_SPACEDIM == 3
                    is_edge |=
                        (amrex::Math::abs(zm) < outer_cutoff) ||
                        (amrex::Math::abs(zp) < outer_cutoff);
    #endif
                }

                // --- interface detection ---
                int iface_cpt = 0;
                const int cpt1 = static_cast<int>(tube_arr(i,j,k,Tidx::OBJ_ID));

                if (inner_tube) {
                    const bool pos  = (phi_c > 0.0);
                    const bool zero = (phi_c == 0.0);

                    auto check_nbr = [&](Set::Scalar phi_n, int nbr_cpt)
                    {
                        // Exit early if cpt flag != 0
                        if (iface_cpt != 0) return; 

                        // Exit early if nbr has same sign
                        if (!zero && (phi_n > 0.0) == pos) return;

                        iface_cpt = amrex::max(cpt1, nbr_cpt);
                    };

                    check_nbr(xm, tube_arr(i-1,j,k,Tidx::OBJ_ID));
                    check_nbr(xp, tube_arr(i+1,j,k,Tidx::OBJ_ID));
    #if AMREX_SPACEDIM >= 2
                    check_nbr(ym, tube_arr(i,j-1,k,Tidx::OBJ_ID));
                    check_nbr(yp, tube_arr(i,j+1,k,Tidx::OBJ_ID));
    #endif
    #if AMREX_SPACEDIM == 3
                    check_nbr(zm, tube_arr(i,j,k-1,Tidx::OBJ_ID));
                    check_nbr(zp, tube_arr(i,j,k+1,Tidx::OBJ_ID));
    #endif
                }

                // --- tube classification ---
                int type = TT::OUTSIDETUBE;
                if (is_edge)      type = TT::EDGEPOINT;
                if (inside_tube)  type = TT::INSIDETUBE;
                if (inner_tube)   type = TT::INNERBAND;
                if (iface_cpt != 0) type = TT::INTERFACE;

                tube_arr(i,j,k,Tidx::TUBE_TYPE) = type;

                // --- narrowband mask ---
                int mask = (type != TT::OUTSIDETUBE);

                tube_arr(i,j,k,Tidx::NB_MASK) = mask ? 1.0 : 0.0;

                // --- interface id ---
                tube_arr(i,j,k,Tidx::ZERO_ID) = Set::Scalar(iface_cpt);

                // --- reduction contribution ---
                return {mask};
            });

            int tile_flag = amrex::get<0>(reduce_data.value());
            tiles[tile_idx] = static_cast<char>(tile_flag);
        }

        TUBE[lev]->FillBoundary();
    }

    void extendVelocityPDE(int lev, const amrex::Geometry& geom) {
        using Tidx = NarrowBandData::TubeIdx;

        const Set::Scalar dtau = 0.4;

        // Create tmp MultiFab for velold
        amrex::MultiFab vel_old;

        vel_old.define(
            LSvel[lev]->boxArray(),
            LSvel[lev]->DistributionMap(),
            LSvel[lev]->nComp(),
            1
        );

        // Get variables
        const auto& dx = geom.CellSizeArray();

        const auto& tiles = tile_has_nb[lev];

        // Perform 25 timesteps (10 is sufficient for convergence)
        for (int i = 0; i < 25; ++i) {
            // Copy new velocity to old velocity
            amrex::MultiFab::Copy(
                vel_old,
                *LSvel[lev],
                0,
                0,
                AMREX_SPACEDIM,
                1
            );

            // Fill ghosts
            vel_old.FillBoundary(geom.periodicity());

            int tile_idx = 0;

            // Perform tile loop
            for (amrex::MFIter mfi(*LSvel[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi, ++tile_idx) {
                // Skip non-narrowband tiles
                if (!tiles[tile_idx]) continue;

                const auto& tilebox = mfi.tilebox();

                Set::Patch<const Set::Scalar> phi_arr = LS.Patch(lev, mfi);
                Set::Patch<const Set::Scalar> tube_arr = TUBE.Patch(lev, mfi);
                Set::Patch<const Set::Scalar> norm_arr = LSnormal.Patch(lev, mfi);
                Set::Patch<      Set::Scalar> vel_arr  = LSvel.Patch(lev, mfi);
                Set::Patch<      Set::Scalar> vel_old_arr = vel_old.array(mfi);

                amrex::ParallelFor(tilebox, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    // Skip non-narrowband cells
                    if (tube_arr(i,j,k,Tidx::NB_MASK) == 0.0) return;

                    const Set::Scalar H = 
                        (tube_arr(i,j,k,Tidx::ZERO_ID) != 0.0) ? 0.0 : 1.0;

                    const Set::Scalar sgn = 
                        phi_arr(i,j,k) >= 0.0 ? 1.0 : -1.0;

                    const Set::Scalar Sx = sgn * norm_arr(i,j,k,0);
                    const Set::Scalar absSx = amrex::Math::abs(Sx);
    #if AMREX_SPACEDIM >=2
                    const Set::Scalar Sy = sgn * norm_arr(i,j,k,1);
                    const Set::Scalar absSy = amrex::Math::abs(Sy);
    #endif
    #if AMREX_SPACEDIM == 3
                    const Set::Scalar Sz = sgn * norm_arr(i,j,k,2);
                    const Set::Scalar absSz = amrex::Math::abs(Sz);
    #endif

                    // X-direction flux
                    {
                        const int comp = 0;
                        const Set::Scalar u0 = vel_old_arr(i,j,k,comp);

                        const Set::Scalar DE = vel_old_arr(i+1,j,k,comp) - u0;
                        const Set::Scalar DW = u0 - vel_old_arr(i-1,j,k,comp);
    #if AMREX_SPACEDIM >= 2
                        const Set::Scalar DN = vel_old_arr(i,j+1,k,comp) - u0;
                        const Set::Scalar DS = u0 - vel_old_arr(i,j-1,k,comp);
    #endif
    #if AMREX_SPACEDIM == 3
                        const Set::Scalar DT = vel_old_arr(i,j,k+1,comp) - u0;
                        const Set::Scalar DB = u0 - vel_old_arr(i,j,k-1,comp);
    #endif
                        const Set::Scalar flux = upwind_flux_component(
                            DE, DW,
    #if AMREX_SPACEDIM >= 2
                            DN, DS,
    #endif
    #if AMREX_SPACEDIM == 3
                            DT, DB,
    #endif
                            Sx,
    #if AMREX_SPACEDIM >= 2
                            Sy,
    #endif
    #if AMREX_SPACEDIM == 3
                            Sz,
    #endif
                            absSx,
    #if AMREX_SPACEDIM >= 2
                            absSy,
    #endif
    #if AMREX_SPACEDIM == 3
                            absSz,
    #endif
                            dx
                        );

                    vel_arr(i,j,k,comp) = u0 - 0.5*H*dtau*flux;
                    }

    #if AMREX_SPACEDIM >= 2
                    // Y-direction flux
                    {
                        const int comp = 1;
                        const Set::Scalar v0 = vel_old_arr(i,j,k,comp);

                        const Set::Scalar DE = vel_old_arr(i+1,j,k,comp) - v0;
                        const Set::Scalar DW = v0 - vel_old_arr(i-1,j,k,comp);
                        const Set::Scalar DN = vel_old_arr(i,j+1,k,comp) - v0;
                        const Set::Scalar DS = v0 - vel_old_arr(i,j-1,k,comp);
    #if AMREX_SPACEDIM == 3
                        const Set::Scalar DT = vel_old_arr(i,j,k+1,comp) - v0;
                        const Set::Scalar DB = v0 - vel_old_arr(i,j,k-1,comp);
    #endif
                        const Set::Scalar flux = upwind_flux_component(
                            DE, DW,
                            DN, DS,
    #if AMREX_SPACEDIM == 3
                            DT, DB,
    #endif
                            Sx,
                            Sy,
    #if AMREX_SPACEDIM == 3
                            Sz,
    #endif
                            absSx,
                            absSy,
    #if AMREX_SPACEDIM == 3
                            absSz,
    #endif
                            dx
                        );

                    vel_arr(i,j,k,comp) = v0 - 0.5*H*dtau*flux;
                    }
    #endif

    #if AMREX_SPACEDIM == 3
                     // Z-direction flux
                    {
                        const int comp = 2;
                        const Set::Scalar w0 = vel_old_arr(i,j,k,comp);

                        const Set::Scalar DE = vel_old_arr(i+1,j,k,comp) - w0;
                        const Set::Scalar DW = w0 - vel_old_arr(i-1,j,k,comp);
                        const Set::Scalar DN = vel_old_arr(i,j+1,k,comp) - w0;
                        const Set::Scalar DS = w0 - vel_old_arr(i,j-1,k,comp);
                        const Set::Scalar DT = vel_old_arr(i,j,k+1,comp) - w0;
                        const Set::Scalar DB = w0 - vel_old_arr(i,j,k-1,comp);

                        const Set::Scalar flux = upwind_flux_component(
                            DE, DW,
                            DN, DS,
                            DT, DB,
                            Sx,
                            Sy,
                            Sz,
                            absSx,
                            absSy,
                            absSz,
                            dx
                        );

                    vel_arr(i,j,k,comp) = w0 - 0.5*H*dtau*flux;
                    }    
    #endif               
                });
            }
        }
    }

    void compute_normal(int lev, const amrex::Geometry& geom) {
        // Set normal to 0.0
        LSnormal[lev]->setVal(0.0);

        // Get dx per dimension outside loop
        auto const dx = geom.CellSizeArray();

        // Get tileboxes
        const auto& tiles = tile_has_nb[lev];

        int tile_idx = 0;

        for (amrex::MFIter mfi(*LSnormal[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi, ++tile_idx)
        {
            // Skip non-narrowband tiles
            if (!tiles[tile_idx]) continue;

            const amrex::Box& tilebox = mfi.tilebox();

            Set::Patch<const Set::Scalar> phi_arr  = LS.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> tube_arr = TUBE.Patch(lev, mfi);
            Set::Patch<      Set::Scalar> norm_arr = LSnormal.Patch(lev, mfi);

            amrex::ParallelFor(tilebox, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                // Skip non-narrowband cells
                const int bandval = static_cast<int>(tube_arr(i, j, k, NarrowBandData::TubeIdx::NB_MASK));
                if (bandval == 0) return;

                // Get the normal
                Set::Vector normal = NBStencil::normal(i,j,k,phi_arr,tube_arr,dx);
                for (size_t d = 0; d < AMREX_SPACEDIM; ++d) {
                    norm_arr(i,j,k,d) = normal[d];
                } 
            });
        }

        // Fill ghosts
        LSnormal[lev]->FillBoundary();
    }

    void compute_curvature(int lev, const amrex::Geometry& geom) {
        // Set curvature to 0.0
        LSkappa[lev]->setVal(0.0);

        // Get dx per dimension outside loop
        auto const dx = geom.CellSizeArray();

        // Get tileboxes
        const auto& tiles = tile_has_nb[lev];

        int tile_idx = 0;

        for (amrex::MFIter mfi(*LSkappa[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi, ++tile_idx)
        {
            // Skip non-narrowband tiles
            if (!tiles[tile_idx]) continue;

            const amrex::Box& tilebox = mfi.tilebox();

            Set::Patch<const Set::Scalar> tube_arr  = TUBE.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> norm_arr  = LSnormal.Patch(lev, mfi);
            Set::Patch<      Set::Scalar> kappa_arr = LSkappa.Patch(lev, mfi);

            amrex::ParallelFor(tilebox, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                // Skip non-narrowband cells
                const int bandval = static_cast<int>(tube_arr(i, j, k, NarrowBandData::TubeIdx::NB_MASK));
                if (bandval == 0) return;

                // Get the curvature
                kappa_arr(i,j,k) = NBStencil::curvature(
                    i,j,k,norm_arr,tube_arr,dx
                );
            });
        }

        LSkappa[lev]->FillBoundary();
    }

    void compute_geometry (int lev, const amrex::Geometry& geom) {
        compute_normal(lev, geom);
        compute_curvature(lev, geom);
    }

    void postfixLSAdvect(int lev, const amrex::Geometry& geom) {
        using Tidx = NarrowBandData::TubeIdx;

        // Get tiles
        const auto& tiles = tile_has_nb[lev];

        int tile_idx = 0;

        for (amrex::MFIter mfi(*LSkappa[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi, ++tile_idx) {
            // Skip non-narrowband tiles
            if (!tiles[tile_idx]) continue;

            const amrex::Box& tilebox = mfi.tilebox();

            Set::Patch<const Set::Scalar> Tube_arr  = TUBE.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> LSold_arr = LSold.Patch(lev, mfi);
            Set::Patch<      Set::Scalar> LS_arr    = LS.Patch(lev, mfi);

            amrex::ParallelFor(tilebox, [=] AMREX_GPU_DEVICE (int i, int j, int k) 
            {
                // Skip non-narrowband cells
                const int band_val = static_cast<int>(Tube_arr(i,j,k, Tidx::NB_MASK));
                if (band_val == 0) return;

                // Skip cells where phi is positive
                Set::Scalar& phi = LS_arr(i,j,k);
                const Set::Scalar phiold = LSold_arr(i,j,k);
                const int zero = static_cast<int>(Tube_arr(i,j,k,Tidx::ZERO_ID));

                if (phi > 0.0) return;

                const bool false_crossing = (phi * phiold < 0.0) && (zero == 0);

                phi = false_crossing ? phiold : phi;
            });
        }

        LS[lev]->FillBoundary();
    }

    void updateCPT(int lev) {
        using Tidx = NarrowBandData::TubeIdx;

        // Get tiles
        const auto& tiles = tile_has_nb[lev];

        int tile_idx = 0;

        for (amrex::MFIter mfi(*LSkappa[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi, ++tile_idx) {
            // Skip non-narrowband tiles
            if (!tiles[tile_idx]) continue;

            const amrex::Box& tilebox = mfi.tilebox();

            Set::Patch<      Set::Scalar> Tube_arr  = TUBE.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> LS_arr = LS.Patch(lev, mfi);

            amrex::ParallelFor(tilebox, [=] AMREX_GPU_DEVICE (int i, int j, int k) 
            {
                // Skip non-narrowband cells
                const int band_val = static_cast<int>(Tube_arr(i,j,k, Tidx::NB_MASK));
                if (band_val == 0) return;

                // Define variables
                const Set::Scalar phi  = LS_arr(i,j,k);
                const Set::Scalar zero = Tube_arr(i,j,k,Tidx::ZERO_ID);
                const Set::Scalar cptold = Tube_arr(i,j,k,Tidx::OBJ_ID);
                Set::Scalar& cptnew = Tube_arr(i,j,k,Tidx::OBJ_ID);

                cptnew = 
                    (phi <= 0.0) ? zero :
                    (cptold == zero) ? 0.0 : 
                                    cptold;
            }); 
        }    
        
        TUBE[lev]->FillBoundary(Tidx::OBJ_ID, 1);
    }

    void reinitialize(int lev, const amrex::Geometry& geom) 
    {
        using Tidx = NarrowBandData::TubeIdx;
        using XDir = NBStencil::XDir;
    #if AMREX_SPACEDIM >= 2
        using YDir = NBStencil::YDir;
    #endif
    #if AMREX_SPACEDIM == 3
        using ZDir = NBStencil::ZDir;
    #endif

        // --- Grid spacing ---
        const Set::Scalar* DX = geom.CellSize();
        Set::Scalar min_dx = *std::min_element(DX, DX + AMREX_SPACEDIM);
        auto const dx = geom.CellSizeArray();

        const Set::Scalar inner_cutoff = HALFINNERTUBE * min_dx;
        const Set::Scalar outer_cutoff = TUBECUTOFF    * min_dx;

        // --- Domain bounds ---
        const amrex::Box domain = geom.Domain();
        const auto dlo = domain.smallEnd();
        const auto dhi = domain.bigEnd();

        constexpr int NEVER = -1;
        constexpr int NO    =  0;
        constexpr int YES   =  1;

        // Allocate reinit/reinit_old
        amrex::MultiFab LSreinit, LSreinitold;
        amrex::iMultiFab compflag;

        LSreinit.define(
            LS[lev]->boxArray(),
            LS[lev]->DistributionMap(),
            1,
            2 // Eno requires 2 ghosts
        );

        LSreinitold.define(
            LS[lev]->boxArray(),
            LS[lev]->DistributionMap(),
            1,
            2 
        );

        compflag.define(
            LS[lev]->boxArray(),
            LS[lev]->DistributionMap(),
            1,
            2 
        );

        // Initialize variables
        LSreinit.setVal(10.0*outer_cutoff);
        compflag.setVal(NO);

        // Get tiles
        const auto& tiles = tile_has_nb[lev];
        int tile_idx = 0;

        for (amrex::MFIter mfi(LSreinit, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi, ++tile_idx) {
            // Skip non-narrowband tiles
            if (!tiles[tile_idx]) continue;

            const amrex::Box& tilebox = mfi.tilebox();

            Set::Patch<const Set::Scalar> tube_arr = TUBE.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> ls_arr = LS.Patch(lev, mfi);
            Set::Patch<      Set::Scalar> reinit_arr = LSreinit.array(mfi);
            const auto&                   compflag_arr = compflag.array(mfi);

            amrex::ParallelFor(tilebox, [=] AMREX_GPU_DEVICE (int i, int j, int k) 
            {
                // Skip non-narrowband cells
                const int mask = static_cast<int>(tube_arr(i,j,k,Tidx::NB_MASK));

                if (!mask) return;

                // If narrowband, set LSreinit to LS with proper sign
                const Set::Scalar phi = ls_arr(i,j,k);
                const Set::Scalar absphi = amrex::Math::abs(phi);
                const Set::Scalar cpt = tube_arr(i,j,k,Tidx::OBJ_ID);
                const Set::Scalar zero = tube_arr(i,j,k,Tidx::ZERO_ID);

                reinit_arr(i,j,k) = (cpt > 0.0) ? -absphi : phi;

                // --- neighbors ---
                Set::Scalar xm  = (i != dlo[0]) ? ls_arr(i-1,j,k) : phi;
                Set::Scalar x0m = (i != dlo[0]) ? tube_arr(i-1,j,k,Tidx::ZERO_ID) : 0.0;

                Set::Scalar xp  = (i != dhi[0]) ? ls_arr(i+1,j,k) : phi;
                Set::Scalar x0p = (i != dhi[0]) ? tube_arr(i+1,j,k,Tidx::ZERO_ID) : 0.0;
    #if AMREX_SPACEDIM >= 2
                Set::Scalar ym  = (j != dlo[1]) ? ls_arr(i,j-1,k) : phi;
                Set::Scalar y0m = (j != dlo[1]) ? tube_arr(i,j-1,k,Tidx::ZERO_ID) : 0.0;

                Set::Scalar yp  = (j != dhi[1]) ? ls_arr(i,j+1,k) : phi;
                Set::Scalar y0p = (j != dhi[1]) ? tube_arr(i,j+1,k,Tidx::ZERO_ID) : 0.0;
    #endif
    #if AMREX_SPACEDIM == 3
                Set::Scalar zm  = (k != dlo[2]) ? ls_arr(i,j,k-1) : phi;
                Set::Scalar z0m = (k != dlo[2]) ? tube_arr(i,j,k-1,Tidx::ZERO_ID) : 0.0;

                Set::Scalar zp  = (k != dhi[2]) ? ls_arr(i,j,k+1) : phi;
                Set::Scalar z0p = (k != dhi[2]) ? tube_arr(i,j,k+1,Tidx::ZERO_ID) : 0.0;
    #endif

                bool mark_never = false;

                auto check_nbr = [=] AMREX_GPU_DEVICE
                    (Set::Scalar phi_nbr, Set::Scalar zero_nbr, bool& mark_never_local)
                {
                    if (mark_never_local) return;

                    const bool is_zero_pair = (zero > 0.0) || (zero_nbr > 0.0);

                    if (is_zero_pair &&
                        phi * phi_nbr < 0.0 &&
                        amrex::max(phi, phi_nbr) < inner_cutoff)
                    {
                        mark_never_local = true;
                    }
                };

                check_nbr(xm, x0m, mark_never);
                check_nbr(xp, x0p, mark_never);
    #if AMREX_SPACEDIM >= 2
                check_nbr(ym, y0m, mark_never);
                check_nbr(yp, y0p, mark_never);
    #endif
    #if AMREX_SPACEDIM == 3
                check_nbr(zm, z0m, mark_never);
                check_nbr(zp, z0p, mark_never);
    #endif

                if (mark_never) compflag_arr(i,j,k) = NEVER;
            });
        }

        // Fill boundary
        compflag.FillBoundary();
        LSreinit.FillBoundary();

        // Perform tube expansion
        tile_idx = 0;

        for (amrex::MFIter mfi(LSreinit, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi, ++tile_idx) {
            // Skip non-narrowband tiles
            if (!tiles[tile_idx]) continue;

            const amrex::Box& tilebox = mfi.tilebox();

            Set::Patch<const Set::Scalar> tube_arr = TUBE.Patch(lev, mfi);
            Set::Patch<      Set::Scalar> reinit_arr = LSreinit.array(mfi);
            const auto&                   compflag_arr = compflag.array(mfi);

            amrex::ParallelFor(tilebox, [=] AMREX_GPU_DEVICE (int i, int j, int k) 
            {
                // Skip non-narrowband cells
                const int mask = static_cast<int>(tube_arr(i,j,k,Tidx::NB_MASK));

                if (!mask) return;

                // Skip if already marked
                if (compflag_arr(i,j,k) != NO) return;

                // Inherit sign from neighbor
                bool found = false;
                amrex::Real sgn = 0.0;

                auto try_nbr = [=] AMREX_GPU_DEVICE
                    (int ni,int nj,int nk, bool& found_local, amrex::Real& sgn_local)
                {
                    if (found_local) return;

                    if (compflag_arr(ni,nj,nk) != NEVER)
                    {
                        sgn_local = amrex::Math::copysign(1.0, reinit_arr(ni,nj,nk));
                        found_local = true;
                    }
                };

                try_nbr(i-1,j,k,found,sgn); try_nbr(i+1,j,k,found,sgn); 
    #if AMREX_SPACEDIM >= 2
                try_nbr(i,j-1,k,found,sgn); try_nbr(i,j+1,k,found,sgn); 
    #endif
    #if AMREX_SPACEDIM == 3
                try_nbr(i,j,k-1,found,sgn); try_nbr(i,j,k+1,found,sgn); 
    #endif
                if (found)
                {
                    reinit_arr(i,j,k) = amrex::Math::copysign(reinit_arr(i,j,k), sgn);
                    compflag_arr(i,j,k) = YES;
                }
            });
        }

        // Fill ghosts
        LSreinit.FillBoundary();
        compflag.FillBoundary();

        // Copy LSreinitold
        amrex::MultiFab::Copy(
            LSreinitold,
            LSreinit,
            0,
            0,
            LSreinit.nComp(),
            LSreinit.nGrow()
        );

        // Define variables for timestep loop
        amrex::MultiFab error;
        error.define(
            LSreinit.boxArray(),
            LSreinit.DistributionMap(),
            1,
            0
        );

        error.setVal(0.0);

        const Set::Scalar toler_error_reinit = 1e-3;
        const int LSREINIT_MAX_ITERS = 100;
        const Set::Scalar cflReinit = 0.4;
        const Set::Scalar DTreinit = cflReinit * min_dx;

        Set::Scalar error_max = 2.0 * toler_error_reinit;
        int icount = 0;
        bool CHECK_ERROR = false;

        // Begin artifical timestep loop
        while (error_max > toler_error_reinit && icount <= LSREINIT_MAX_ITERS) {
            // Increase iteration count
            ++ icount;

            // Check error
            CHECK_ERROR = (icount % 5 == 0) ? true : false;

            // Copy old reinit value and fill ghosts
            amrex::MultiFab::Copy(
                LSreinitold,
                LSreinit,
                0,
                0,
                1,
                2
            );
            LSreinitold.FillBoundary();

            tile_idx = 0;

            for (amrex::MFIter mfi(LSreinit, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi, ++tile_idx) {
                // Skip non-narrowband tiles
                if (!tiles[tile_idx]) continue;

                const amrex::Box& tilebox = mfi.tilebox();

                Set::Patch<const Set::Scalar> LS_arr = LS.Patch(lev, mfi);
                Set::Patch<const Set::Scalar> tube_arr = TUBE.Patch(lev, mfi);
                Set::Patch<      Set::Scalar> reinit_arr = LSreinit.array(mfi);
                Set::Patch<      Set::Scalar> reinitold_arr = LSreinitold.array(mfi);
                const auto&                   compflag_arr = compflag.const_array(mfi);

                const auto&                   error_arr = error.array(mfi);

                amrex::ParallelFor(tilebox, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    // Skip non-narrowband cells
                    const int mask = static_cast<int>(tube_arr(i,j,k,Tidx::NB_MASK));

                    if (!mask) return;

                    // Skip NO/NEVER cells
                    if (compflag_arr(i,j,k) <= NO) return;

                    // Compute ENO fluxes
                    const Set::Scalar phi = LS_arr(i,j,k);
                    const Set::Scalar sgn = amrex::Math::copysign(1.0, phi);

                    Set::Scalar dfx = NB::FluxHandler::eno2_flux<XDir>(
                        i,j,k,sgn,reinitold_arr,tube_arr,dx
                    );
                    Set::Scalar dfy = 0.0;
                    Set::Scalar dfz = 0.0;
    #if AMREX_SPACEDIM >= 2
                    dfy = NB::FluxHandler::eno2_flux<YDir>(
                        i,j,k,sgn,reinitold_arr,tube_arr,dx
                    );
    #endif
    #if AMREX_SPACEDIM == 3
                    dfz = NB::FluxHandler::eno2_flux<ZDir>(
                        i,j,k,sgn,reinitold_arr,tube_arr,dx
                    );
    #endif

                    Set::Scalar nablaG = std::sqrt(dfx*dfx + dfy*dfy + dfx*dfz);

                    const Set::Scalar& reold = reinitold_arr(i,j,k);
                    Set::Scalar& re = reinit_arr(i,j,k);

                    re = (
                        reold - DTreinit*sgn*(nablaG - 1.0)
                    );
                    re = amrex::Math::copysign(re, sgn);

                    // Check error
                    if (re <= outer_cutoff && compflag_arr(i,j,k) > NO && CHECK_ERROR) {
                        error_arr(i,j,k) = amrex::Math::abs(nablaG - 1.0);
                    }
                });
            }

            // Find max error
            error_max = 10.0;
            if (CHECK_ERROR) {
                error_max = error.max(0);
                error.setVal(0.0);
            }
        }

        // Copy updated band values and enforce clamping
        tile_idx = 0;

        for (amrex::MFIter mfi(LSreinit, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi, ++tile_idx) {
            // Skip non-narrowband tiles
            if (!tiles[tile_idx]) continue;

            const amrex::Box& tilebox = mfi.tilebox();

            Set::Patch<      Set::Scalar> LS_arr = LS.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> tube_arr = TUBE.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> reinit_arr = LSreinit.const_array(mfi);

            amrex::ParallelFor(tilebox, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                // Skip non-narrowband cells
                const int mask = static_cast<int>(tube_arr(i,j,k,Tidx::NB_MASK));

                if (!mask) return;

                // Copy LSreinit to LS
                Set::Scalar& phi = LS_arr(i,j,k);
                phi = reinit_arr(i,j,k);

                // Apply clamping
                const int cpt = static_cast<int>(tube_arr(i,j,k,Tidx::OBJ_ID));
                phi = enforce_signed_band(phi, cpt, outer_cutoff);
            });
        }
    }
};

// ======================================
// Simulation Domain (AMR FIXED)
// ======================================
struct SimulationDomain
{
    // Store reference to integrator
    Integrator::Integrator* integrator = nullptr;

    // Save ghost cell width (3 for 5th order WENO)
    int nghost = 3;

    // Store collection of level set fields
    amrex::Vector<std::unique_ptr<LevelSetField>> fields;

    // Store velocity field (will copy from flow domain in future)
    Set::Field<Set::Scalar> Flowvel;

    std::unique_ptr<IC::IC<Set::Scalar>> ic_Vel;
    std::unique_ptr<BC::BC<Set::Scalar>> bc_Vel;

    // Define function to attach integrator to retrieve geometry information
    void Define(Integrator::Integrator* _integrator)
    {
        integrator = _integrator;
    }

    // Define helper functions to retrieve geometry information from intergrator
    const amrex::Geometry& Geom(int lev) const
    {
        return integrator->Geom(lev);
    }

    const amrex::BoxArray& BA(int lev) const
    {
        return integrator->boxArray(lev);
    }

    const amrex::DistributionMapping& DM(int lev) const
    {
        return integrator->DistributionMap(lev);
    }

    int FinestLevel() const
    {
        return integrator->finestLevel();
    }
};

// Define Integrator namespace
namespace Integrator
{
    // Define Constructor/Destructor
    NarrowBandLevelset::NarrowBandLevelset() = default;
    NarrowBandLevelset::NarrowBandLevelset(IO::ParmParse& pp)
        : Integrator()
    {
        Parse(*this, pp);
    }

    NarrowBandLevelset::~NarrowBandLevelset() = default;

    // Define Parser
    void NarrowBandLevelset::Parse(
        NarrowBandLevelset& value,
        IO::ParmParse& pp)
    {
        // -------------------------
        // Allocate domain
        // -------------------------
        value.domain_ = std::make_unique<SimulationDomain>();

        // Attach integrator
        value.domain_->Define(&value);

        // -------------------------
        // Initialize velocity field (Will pull from flow domain in future)
        // -------------------------
        IC::IC<Set::Scalar>* ic_vel_tmp = nullptr;
        BC::BC<Set::Scalar>* bc_vel_tmp = nullptr;

        pp.select_default<IC::Constant, IC::Expression>(
            "velocity.ic",
            ic_vel_tmp,
            value.geom
        );

        value.domain_->ic_Vel =
            std::unique_ptr<IC::IC<Set::Scalar>>(ic_vel_tmp);

        pp.select_default<BC::Constant, BC::Expression>(
            "velocity.bc",
            bc_vel_tmp,
            AMREX_SPACEDIM
        );
        value.domain_->bc_Vel = std::unique_ptr<BC::BC<Set::Scalar>>(bc_vel_tmp);

        value.AddField<Set::Scalar, Set::Hypercube::Cell>(
            value.domain_->Flowvel,          // Set::Field<Set::Vector>
            nullptr,                         // BC::BC<Set::Vector>*
            AMREX_SPACEDIM,             
            value.domain_->nghost,         
            /* name    */ "Flowvel",         // plotfile base name
            /* writeout*/ true,             // include in plotfiles
            /* evolving*/ false              // evolves in time

        );
        
        // -------------------------
        // Number of level sets
        // -------------------------
        int nls = 1;
        pp_query_default("ls.number_of_levelsets", nls, 1);

        // -------------------------
        // Allocate fields
        // -------------------------
        value.domain_->fields.resize(nls);

        for (int i = 0; i < nls; ++i)
        {
            value.domain_->fields[i] = std::make_unique<LevelSetField>();
            LevelSetField& ls = *value.domain_->fields[i];

            std::string prefix = "ls" + std::to_string(i);
            IO::ParmParse pp_ls(prefix.c_str());

            ls.ic_LS = nullptr;
            ls.bc_LS = nullptr;

            IC::IC<Set::Scalar>* ic_tmp = nullptr;
            BC::BC<Set::Scalar>* bc_tmp = nullptr;

            // -------------------------
            // IC (AMR FIXED)
            // -------------------------
            pp_ls.select_default<IC::Constant, IC::Expression>(
                "ic",
                ic_tmp,
                value.geom
            );
            ls.ic_LS = std::unique_ptr<IC::IC<Set::Scalar>>(ic_tmp);

            // -------------------------
            // BC
            // -------------------------
            pp.select_default<BC::Constant, BC::Expression>(
                "ls.bc",
                bc_tmp,
                1
            );
            ls.bc_LS = std::unique_ptr<BC::BC<Set::Scalar>>(bc_tmp);
            // ls.bc_LS = std::make_unique<BC::LevelSetNeumann>();

            pp.select_default<BC::Constant, BC::Expression>(
                "normal.bc",
                bc_tmp,
                AMREX_SPACEDIM
            );
            ls.bc_LSnormal = std::unique_ptr<BC::BC<Set::Scalar>>(bc_tmp);

            // -------------------------
            // Register fields (AMR handled internally)
            // -------------------------
            value.RegisterOneField(ls, prefix, value.domain_->bc_Vel);
        }
    }

    void NarrowBandLevelset::RegisterOneField(
        LevelSetField& ls,
        const std::string& prefix,
        const std::unique_ptr<BC::BC<Set::Scalar>>& bc_vel)
    {
        // ====================================================
        // Evolving level set fields (PDE state)
        // ====================================================

        RegisterNewFab(
            ls.LS,
            ls.bc_LS.get(),
            1,
            domain_->nghost,
            prefix + "_LS",
            true
        );

        RegisterNewFab(
            ls.LSold,
            ls.bc_LS.get(),
            1,
            domain_->nghost,
            prefix + "_LSold",
            true
        );

        // ====================================================
        // Derived / diagnostic fields (non-evolving, no BCs)
        // ====================================================

        // Level-set velocity (vector)
        AddField<Set::Scalar, Set::Hypercube::Cell>(
            ls.LSvel,
            bc_vel.get(),         
            AMREX_SPACEDIM,     
            domain_->nghost,  
            prefix + "_LSvel",
            true,             // plot
            false             // NOT evolving
        );

        // Level-set normal (vector)
        AddField<Set::Scalar, Set::Hypercube::Cell>(
            ls.LSnormal,
            ls.bc_LSnormal.get(),
            AMREX_SPACEDIM,
            domain_->nghost,
            prefix + "_LSnormal",
            true,
            false
        );

        // Curvature (scalar)
        AddField<Set::Scalar, Set::Hypercube::Cell>(
            ls.LSkappa,
            ls.bc_LS.get(),
            1,
            domain_->nghost,
            prefix + "_LSkappa",
            true,
            false
        );

        // Tube data (int flags packed as scalars)
        AddField<Set::Scalar, Set::Hypercube::Cell>(
            ls.TUBE,
            nullptr,
            NarrowBandData::TubeIdx::NCOMP,                
            domain_->nghost,
            prefix + "_TUBE",
            true,
            false
        );
    }


    // Placeholder override functions
    void NarrowBandLevelset::Initialize(int lev) {
        // Initialize global velocity field (will use flow domain in future)
        domain_->ic_Vel->Initialize(lev, domain_->Flowvel);

        // Loop over all level set fields and initialize
        int global_obj_id = 0;
        for (auto& ls_ptr : domain_->fields)
        {
            // Get the pointer to the current level set field
            LevelSetField& ls = *ls_ptr;
            const auto& geom = domain_->Geom(lev);

            // If first level, size all vector data
            if (lev == 0) {
                const int finest_lev = domain_->FinestLevel();

                // Resize tile vector based on num AMR levels
                ls.resizeNBTiles(finest_lev);

                // Get number of objects and resize
                const int num_obj = ls.GetNumObjects();
                ls.objects.resize(num_obj);

                // Loop through objects and assign global ID
                for (int i = 0; i < num_obj; ++i) {
                    LevelSetObject& obj = ls.objects[i];

                    obj.obj_ID = ++global_obj_id;
                }
            }

            // Initialize LS from IC
            ls.ic_LS->Initialize(lev, ls.LS);

            // Initialize CPT
            ls.initializeObjID(lev, geom);

            // Update tube
            ls.updateTube(lev, geom);

            // Copy LS to LSold
            amrex::MultiFab::Copy(*ls.LSold[lev],
                *ls.LS[lev],
                0,
                0,
                ls.LS[lev]->nComp(),
                ls.LS[lev]->nGrow()
            );

            // Copy domain velocity to level set velocity
            amrex::MultiFab::Copy(
                *ls.LSvel[lev],
                *domain_->Flowvel[lev],
                0, 
                0, 
                ls.LSvel[lev]->nComp(), 
                ls.LSvel[lev]->nGrow()
            );

            // Compute geometry
            ls.compute_geometry(lev, geom);
        }
    }

    void NarrowBandLevelset::TimeStepBegin(Set::Scalar time, int iter) {
        // Compute dt based on max flow velocity/CFL condition
        return;
    }

    void NarrowBandLevelset::TimeStepComplete(Set::Scalar time, int lev) {
        return;
    }

    void NarrowBandLevelset::Advance(int lev, Set::Scalar time, Set::Scalar dt) {
        // For each levelset field in domain_
            // Copy LS to LSold

            // Compute extended velocities

            // Perform Advection

            // Perform Advection levelset fix

            // Update CPT flag

            // if step % 5 == 0 reinitialize

            // Update Tube

            // Compute geometries (normals/curvature)
    }

    void NarrowBandLevelset::TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) {
        return;
    }

    void NarrowBandLevelset::Regrid(int lev, Set::Scalar time) {
        return;
    }

}
