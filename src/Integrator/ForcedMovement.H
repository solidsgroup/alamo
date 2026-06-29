//
//

#ifndef INTEGRATOR_FORCEDMOVEMENT_H // Include guards
#define INTEGRATOR_FORCEDMOVEMENT_H // 

// Standard library includes
#include <iostream>
#include <fstream>
#include <iomanip>

// AMReX Includes
#include "AMReX.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_ParmParse.H"

// Alamo Includes
#include "IO/ParmParse.H"
#include "Integrator/Integrator.H"
#include "BC/Constant.H"
#include "IC/IC.H"
#include "IC/Sphere.H"
#include "IC/Constant.H"
#include "IC/Expression.H"
#include "IC/BMP.H"
#include "IC/PNG.H"
#include "IC/PSRead.H"
#include "Numeric/Stencil.H"
#include "IC/Random.H"
#include "Unit/Unit.H"
#include "IC/PointList.H"

namespace Integrator
{

/// This is the definition of the ForcedMovement class
class ForcedMovement : virtual public Integrator
{
public:
    static constexpr const char* name = "forcedmovement";

    // Empty constructor
    ForcedMovement() : Integrator()
    {}

    // Constructor that triggers parse
    ForcedMovement(IO::ParmParse& pp) : ForcedMovement()
    {
        Parse(*this, pp);
    }

    ~ForcedMovement()
    {
        delete ic;
        delete bc;
    }

    // The Parse function initializes the ForcedMovement object using
    // a parser, pp. 
    // Note that this is a static function, which means it does not have
    // direct access to member variables. Instead, it initializes the variables
    // inside the argument, "value", and so all references to member items are
    // prefixed by "value."
    static void Parse(ForcedMovement& value, IO::ParmParse& pp)
    {
        // Criterion for mesh refinement [0.01]
        pp_query_default("refinement_threshold", value.refinement_threshold, 0.01);

        // Set the initial condition for the alpha field
        pp.select_default<IC::Sphere, IC::Constant, IC::Expression,IC::BMP, IC::PNG, IC::Random, IC::PSRead, IC::PointList>("alpha.ic", value.ic, value.geom);
        
        // Set the initial condition for the velocity field
        pp.select_default<IC::Sphere, IC::Constant, IC::Expression,IC::BMP, IC::PNG, IC::Random, IC::PSRead, IC::PointList>("velocity.ic", value.velocity_ic, value.geom);

        // Use a constant BC object for temperature
        // value.bc = new BC::Constant(1);
        // :ref:`BC::Constant` parameters 
        pp.select_default<BC::Constant>("alpha.bc", value.bc, 1);
        // pp.select_default<BC::Constant>("position_new_bc.bc", value.position_new_bc, 2);
        pp.select_default<BC::Constant>("velocity.bc", value.velocity_bc, 2);
        // pp.select_default<BC::Constant>("alpha.bc", alpha.bc, 1);

        // Register the temperature and old temperature fields.
        // alpha_mf and alpha_old_mf are defined near the bottom of this Header file.
        value.RegisterNewFab(value.alpha_mf, value.bc, value.number_of_components, value.number_of_ghost_cells, "alpha", true);
        value.RegisterNewFab(value.alpha_old_mf, value.bc, value.number_of_components, value.number_of_ghost_cells, "alpha_old", false);
        value.RegisterNewFab(value.velocity_mf, value.velocity_bc, 2, value.number_of_ghost_cells, "velocity_bc", true ,true, {"x","y"});
        // value.RegisterNewFab(value.position_new_mf, value.position_new_bc, 2, value.number_of_ghost_cells, "position_new", true ,true, {"x","y"});
    }

protected:

    // Use the ic object to initialize the order parameter
    void Initialize(int lev)
    {
        ic->Initialize(lev, alpha_mf);
        ic->Initialize(lev, alpha_old_mf);
        velocity_ic->Initialize(lev, velocity_mf);
    }

    // Integrate the forced movement to update the alpha field
    void Advance(int lev, Set::Scalar time, Set::Scalar dt)
    {   
        std::swap(*alpha_mf[lev], *alpha_old_mf[lev]);
        Util::RealFillBoundary(*(alpha_old_mf)[lev], geom[lev]);
        Util::RealFillBoundary(*(velocity_mf)[lev], geom[lev]);
        amrex::Real total_alpha_before = alpha_old_mf[lev]->sum();

        const Set::Scalar *DX = this->geom[lev].CellSize();
        // Evolve alpha

        for (amrex::MFIter mfi(*alpha_mf[lev], true); mfi.isValid(); ++mfi)
        {   

            amrex::Box bx;
            amrex::IndexType type = alpha_mf[lev]->ixType();
            if (type == amrex::IndexType::TheCellType())      bx = mfi.tilebox();
            else if (type == amrex::IndexType::TheNodeType()) bx = mfi.nodaltilebox();
            else Util::Abort(INFO, "Unkonwn index type");

            Set::Patch<Set::Scalar>       alpha_new = alpha_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> alpha     = alpha_old_mf.Patch(lev,mfi);

            Set::Patch<Set::Scalar> position_new = position_new_mf.Patch(lev,mfi);
            Set::Patch<Set::Scalar> velocity = velocity_mf.Patch(lev,mfi);
            // Set::Scalar vel_x = 0.01;
            // Set::Scalar vel_y = 0.01;

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                // Gather-form conservative remap:
                // For target cell (i,j,k), pull mass from each of its 3x3 neighbors j,
                // weighted by the fraction of neighbor j's displaced cell that overlaps
                // this target cell. Each neighbor's weights sum to 1 over the stencil
                // (tent function partition of unity), so this conserves alpha*V exactly
                // (on a uniform-volume stencil) rather than diffusing it.

                Set::Vector TargetCenter = Set::Position(i, j, k, geom[lev], type);

                Set::Scalar Q_new = 0.0; // accumulator for alpha*V gathered into this cell

                Set::Scalar weight_sum = 0.0;   // <-- track this explicitly
                Set::Scalar F_stencil[3][3];    // cache weights, avoid recomputing

                for (int dj = -1; dj <= 1; ++dj)
                for (int di = -1; di <= 1; ++di)
                {
                    int ii = i + di;
                    int jj = j + dj;

                    // Neighbor's own displacement (cell-centered velocity * dt)
                    Set::Vector u_n;
                    u_n(0) = velocity(ii,jj,k,0) * dt;
                    u_n(1) = velocity(ii,jj,k,1) * dt;

                    if (std::abs(u_n(0)) >= DX[0] || std::abs(u_n(1)) >= DX[1])
                    {
                        amrex::Print() << "WARNING: displacement exceeds cell width at ("
                                        << i << "," << j << "," << k << "): "
                                        << "u = (" << u_n(0) << ", " << u_n(1) << "), "
                                        << "DX = (" << DX[0] << ", " << DX[1] << ")\n";
                    }

                    // Neighbor's undeformed center
                    Set::Vector NeighborCenter = Set::Position(ii, jj, k, geom[lev], type);

                    // Displaced neighbor center
                    Set::Scalar xC = NeighborCenter(0) + u_n(0);
                    Set::Scalar yC = NeighborCenter(1) + u_n(1);

                    // Offset of displaced neighbor center from THIS target cell's center
                    Set::Scalar dx = xC - TargetCenter(0);
                    Set::Scalar dy = yC - TargetCenter(1);

                    // 1D tent overlap fractions (DX[0], DX[1] = cell widths)
                    Set::Scalar fx = std::max(0.0, DX[0] - std::abs(dx)) / DX[0];
                    Set::Scalar fy = std::max(0.0, DX[1] - std::abs(dy)) / DX[1];

                    Set::Scalar F = fx * fy;

                    F_stencil[dj+1][di+1] = F;
                    weight_sum += F;
                    }

                    if (weight_sum > 1.0e-12)
                    {
                        for (int dj = -1; dj <= 1; ++dj)
                        for (int di = -1; di <= 1; ++di)
                        {
                            int ii = i + di;
                            int jj = j + dj;
                            Q_new += alpha(ii,jj,k) * (F_stencil[dj+1][di+1] / weight_sum);
                        }
                    }
                    else
                    {
                        // Degenerate case: no overlap found at all (shouldn't happen if
                        // |u| < DX, but guard anyway)
                        Q_new = alpha(i,j,k);
                    }

                    alpha_new(i, j, k) = Q_new;
            });
        }

        amrex::Real total_alpha_after = alpha_mf[lev]->sum();
        amrex::Print() << "Q before = " << total_alpha_before
                    << "  Q after = " << total_alpha_after
                    << "  drift = " << (total_alpha_after - total_alpha_before) << "\n";
        Util::RealFillBoundary(*(alpha_mf)[lev], geom[lev]);
    }

    // Tag cells for mesh refinement based on temperature gradient
    void TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags, Set::Scalar /*time*/, int /*ngrow*/)
    {
        // Get cell dimensions as done above.
        const Set::Scalar* DX = geom[lev].CellSize();
        // Calculate the diagonal.
        Set::Scalar dr = sqrt(AMREX_D_TERM(DX[0] * DX[0], +DX[1] * DX[1], +DX[2] * DX[2]));

        // Iterate over the patches on this level
        for (amrex::MFIter mfi(*alpha_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Get the box and handles as done above.
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<char>         const& tags = a_tags.array(mfi);
            amrex::Array4<Set::Scalar>  const& alpha = (*alpha_mf[lev]).array(mfi);

            // Iterate over the grid as done above.
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                // Calculate the temperature gradient.
                Set::Vector grad = Numeric::Gradient(alpha, i, j, k, 0, DX);

                // Is the gradient * cell_size too big? If so, then
                // mark this cell as needing refinement.
                if (grad.lpNorm<2>() * dr > refinement_threshold)
                    tags(i, j, k) = amrex::TagBox::SET;
            });
        }
    }

protected:
    Set::Field<Set::Scalar> alpha_mf;         /// Temperature field variable (current timestep)
    Set::Field<Set::Scalar> alpha_old_mf;     /// Temperature field variable (previous timestep)
    IC::IC<Set::Scalar>* ic;                  /// Object used to initialize temperature field
    IC::IC<Set::Scalar>* velocity_ic;                  /// Object used to initialize velocity field

    Set::Field<Set::Scalar> &eta_mf = alpha_mf;
    Set::Field<Set::Scalar> &eta_old_mf = alpha_old_mf;

    Set::Field<Set::Vector> position_mf;
    Set::Field<Set::Scalar> position_new_mf; // Updated position of each point after move is applied
    Set::Field<Set::Scalar> velocity_mf; // Velocity of each point

private:
    int number_of_components = 1;            // Number of components
    int number_of_ghost_cells = 2;           // Number of ghost cells

    Set::Scalar refinement_threshold = NAN; // Criterion for cell refinement

    BC::BC<Set::Scalar>* bc;                 // Object used to update temp field boundary ghost cells
    BC::BC<Set::Scalar>* position_new_bc;
    BC::BC<Set::Scalar>* velocity_bc;
};
} // namespace Integrator
#endif
