#ifndef IC_RIEMANN2D_H_
#define IC_RIEMANN2D_H_

#include "IC/Shock.H"
#include "AMReX.H"
#include "AMReX_ParmParse.H"
#include "AMReX_Geometry.H"
#include "AMReX_MultiFab.H"
#include "AMReX_Array4.H"
#include "IC/IC.H"
#include "Set/Set.H"
#include "Util/Util.H"
#include "Util/ScimitarX_Util.H"
#include "Numeric/IntegratorVariableAccessLayer.H"
#include "Model/Fluid/Fluid.H"
#include <string>

namespace IC {

using namespace Model::Fluid;

class Riemann2D : public IC<Set::Scalar> {
public:
    // Constructor requiring variable indices
    Riemann2D(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name,
            const Numeric::GenericVariableAccessor::VariableIndices& precomputed_indices)
        : IC(_geom), mf_name(name), variable_indices(&precomputed_indices), requires_variable_indices(true) 
    {
        initialize(pp, name);
    }

    // Constructor for pressure (no variable indices required)
    Riemann2D(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name)
        : IC(_geom), mf_name(name), variable_indices(nullptr), requires_variable_indices(false)
    {
        initialize(pp, name);
    }

    // Implementation of Add function to handle 2D Riemann problems
    virtual void Add(const int& lev, Set::Field<Set::Scalar>& a_phi, Set::Scalar /*time*/) override {
        int ncomp = a_phi[0]->nComp();

        for (amrex::MFIter mfi(*a_phi[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            amrex::Box bx = mfi.growntilebox();
            amrex::IndexType type = a_phi[lev]->ixType();

            amrex::Array4<Set::Scalar> const& phi = a_phi[lev]->array(mfi);

            // Retrieve variable indices if required
            int dens_idx = requires_variable_indices ? variable_indices->DENS : -1;
            int uvel_idx = requires_variable_indices ? variable_indices->UVEL : -1;
#if AMREX_SPACEDIM >= 2
            int vvel_idx = requires_variable_indices ? variable_indices->VVEL : -1;
#endif
#if AMREX_SPACEDIM == 3
            int wvel_idx = requires_variable_indices ? variable_indices->WVEL : -1;
#endif
            int ie_idx = requires_variable_indices ? variable_indices->IE : -1;

            // For proper pressure-to-energy conversion
            Model::Fluid::Fluid fluid_model;

            for (int n = 0; n < ncomp; ++n) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    Set::Vector x = Set::Position(i, j, k, geom[lev], type);
                    Set::Scalar gamma = 1.4;

                    int quadrant;
                    if (x(0) < riemann_shock_xpos && x(1) < riemann_shock_ypos) {
                        quadrant = 3;  // Bottom left quadrant
                    } else if (x(0) >= riemann_shock_xpos && x(1) < riemann_shock_ypos) {
                        quadrant = 4;  // Bottom right quadrant
                    } else if (x(0) < riemann_shock_xpos && x(1) >= riemann_shock_ypos) {
                        quadrant = 2;  // Top left quadrant
                    } else {
                        quadrant = 1;  // Top right quadrant
                    }

                    if (mf_name.find("pvec") != std::string::npos) {
                        // Handle primitive variable vector
                        if (n == dens_idx) {
                            // Density assignment
                            switch (quadrant) {
                                case 1: phi(i, j, k, n) = q1_density; break;
                                case 2: phi(i, j, k, n) = q2_density; break;
                                case 3: phi(i, j, k, n) = q3_density; break;
                                case 4: phi(i, j, k, n) = q4_density; break;
                            }
                        } else if (n == uvel_idx) {
                            // X-velocity assignment
                            switch (quadrant) {
                                case 1: phi(i, j, k, n) = q1_velocity[0]; break;
                                case 2: phi(i, j, k, n) = q2_velocity[0]; break;
                                case 3: phi(i, j, k, n) = q3_velocity[0]; break;
                                case 4: phi(i, j, k, n) = q4_velocity[0]; break;
                            }
                        }
#if AMREX_SPACEDIM >= 2
                        else if (n == vvel_idx) {
                            // Y-velocity assignment
                            switch (quadrant) {
                                case 1: phi(i, j, k, n) = q1_velocity[1]; break;
                                case 2: phi(i, j, k, n) = q2_velocity[1]; break;
                                case 3: phi(i, j, k, n) = q3_velocity[1]; break;
                                case 4: phi(i, j, k, n) = q4_velocity[1]; break;
                            }
                        }
#endif
#if AMREX_SPACEDIM == 3
                        else if (n == wvel_idx) {
                            // Z-velocity assignment (usually 0 for 2D Riemann problems)
                            switch (quadrant) {
                                case 1: phi(i, j, k, n) = q1_velocity[2]; break;
                                case 2: phi(i, j, k, n) = q2_velocity[2]; break;
                                case 3: phi(i, j, k, n) = q3_velocity[2]; break;
                                case 4: phi(i, j, k, n) = q4_velocity[2]; break;
                            }
                        }
#endif
                        else if (n == ie_idx) {
                            // Internal energy assignment (computed from pressure)
                            Set::Scalar density = 0.0, pressure = 0.0;
                            
                            switch (quadrant) {
                                case 1:
                                    density = q1_density;
                                    pressure = q1_pressure;
                                    break;
                                case 2:
                                    density = q2_density;
                                    pressure = q2_pressure;
                                    break;
                                case 3:
                                    density = q3_density;
                                    pressure = q3_pressure;
                                    break;
                                case 4:
                                    density = q4_density;
                                    pressure = q4_pressure;
                                    break;
                            }
                            
                            phi(i, j, k, n) = fluid_model.ComputeInternalEnergyFromDensityAndPressure(
                                density, pressure, gamma);
                        }
                    } else if (mf_name.find("pressure") != std::string::npos) {
                        // Handle pressure values
                        switch (quadrant) {
                            case 1: phi(i, j, k, n) = q1_pressure; break;
                            case 2: phi(i, j, k, n) = q2_pressure; break;
                            case 3: phi(i, j, k, n) = q3_pressure; break;
                            case 4: phi(i, j, k, n) = q4_pressure; break;
                        }
                    }
                });
            }
        }
    }

private:
    void initialize(IO::ParmParse& pp, const std::string& name) {
        Util::Message(INFO, "Initializing Riemann2D with mf_name = " + name);

        if (name.find("pvec") != std::string::npos) {
            // Parse density values for each quadrant
            pp.query("ic.riemann2d.pvec.q1.DENS", q1_density);
            pp.query("ic.riemann2d.pvec.q2.DENS", q2_density);
            pp.query("ic.riemann2d.pvec.q3.DENS", q3_density);
            pp.query("ic.riemann2d.pvec.q4.DENS", q4_density);

            // Parse velocity values for each quadrant
            pp.query("ic.riemann2d.pvec.q1.UVEL", q1_velocity[0]);
            pp.query("ic.riemann2d.pvec.q2.UVEL", q2_velocity[0]);
            pp.query("ic.riemann2d.pvec.q3.UVEL", q3_velocity[0]);
            pp.query("ic.riemann2d.pvec.q4.UVEL", q4_velocity[0]);

#if AMREX_SPACEDIM >= 2
            pp.query("ic.riemann2d.pvec.q1.VVEL", q1_velocity[1]);
            pp.query("ic.riemann2d.pvec.q2.VVEL", q2_velocity[1]);
            pp.query("ic.riemann2d.pvec.q3.VVEL", q3_velocity[1]);
            pp.query("ic.riemann2d.pvec.q4.VVEL", q4_velocity[1]);
#endif

#if AMREX_SPACEDIM == 3
            // Default z-velocity to zero if not specified
            pp.query("ic.riemann2d.pvec.q1.WVEL", q1_velocity[2]);
            pp.query("ic.riemann2d.pvec.q2.WVEL", q2_velocity[2]);
            pp.query("ic.riemann2d.pvec.q3.WVEL", q3_velocity[2]);
            pp.query("ic.riemann2d.pvec.q4.WVEL", q4_velocity[2]);
#endif

            // Parse pressure values for internal energy calculation
            pp.query("ic.riemann2d.pressure.q1", q1_pressure);
            pp.query("ic.riemann2d.pressure.q2", q2_pressure);
            pp.query("ic.riemann2d.pressure.q3", q3_pressure);
            pp.query("ic.riemann2d.pressure.q4", q4_pressure);

            // Debug output
            Util::Message(INFO, "Parsed Riemann2D parameters for pvec:");
            Util::Message(INFO, "Quadrant 1 (top-right): ρ=" + std::to_string(q1_density) + 
                                " u=" + std::to_string(q1_velocity[0]) + 
                                " v=" + std::to_string(q1_velocity[1]) + 
                                " p=" + std::to_string(q1_pressure));
            Util::Message(INFO, "Quadrant 2 (top-left): ρ=" + std::to_string(q2_density) + 
                                " u=" + std::to_string(q2_velocity[0]) + 
                                " v=" + std::to_string(q2_velocity[1]) + 
                                " p=" + std::to_string(q2_pressure));
            Util::Message(INFO, "Quadrant 3 (bottom-left): ρ=" + std::to_string(q3_density) + 
                                " u=" + std::to_string(q3_velocity[0]) + 
                                " v=" + std::to_string(q3_velocity[1]) + 
                                " p=" + std::to_string(q3_pressure));
            Util::Message(INFO, "Quadrant 4 (bottom-right): ρ=" + std::to_string(q4_density) + 
                                " u=" + std::to_string(q4_velocity[0]) + 
                                " v=" + std::to_string(q4_velocity[1]) + 
                                " p=" + std::to_string(q4_pressure));
        } 
        else if (name.find("pressure") != std::string::npos) {
            // Parse pressure values only
            pp.query("ic.riemann2d.pressure.q1", q1_pressure);
            pp.query("ic.riemann2d.pressure.q2", q2_pressure);
            pp.query("ic.riemann2d.pressure.q3", q3_pressure);
            pp.query("ic.riemann2d.pressure.q4", q4_pressure);

            // Debug output
            Util::Message(INFO, "Parsed Riemann2D parameters for pressure:");
            Util::Message(INFO, "Quadrant 1 (top-right): p=" + std::to_string(q1_pressure));
            Util::Message(INFO, "Quadrant 2 (top-left): p=" + std::to_string(q2_pressure));
            Util::Message(INFO, "Quadrant 3 (bottom-left): p=" + std::to_string(q3_pressure));
            Util::Message(INFO, "Quadrant 4 (bottom-right): p=" + std::to_string(q4_pressure));
        }

        // Parse shock position
        pp.query("shock.xpos", riemann_shock_xpos);
        pp.query("shock.ypos", riemann_shock_ypos);

        Util::Message(INFO, "Riemann2D interface position: x=" + std::to_string(riemann_shock_xpos) + 
                            " y=" + std::to_string(riemann_shock_ypos));
    }

    // MultiFab name for determining parameter set and variable types
    std::string mf_name;
    
    // Variable indices (only required for pvec)
    const Numeric::GenericVariableAccessor::VariableIndices* variable_indices;
    bool requires_variable_indices;

    // Riemann2D specific shock positions
    Set::Scalar riemann_shock_xpos = 0.5;
    Set::Scalar riemann_shock_ypos = 0.5;

    // Quadrant-specific parameters with default values
    Set::Scalar q1_density = 1.0;
    Set::Scalar q2_density = 0.5;
    Set::Scalar q3_density = 0.25;
    Set::Scalar q4_density = 0.5;

    Set::Scalar q1_pressure = 1.0;
    Set::Scalar q2_pressure = 0.4;
    Set::Scalar q3_pressure = 0.2;
    Set::Scalar q4_pressure = 0.4;

    Set::Vector q1_velocity = Set::Vector::Zero();
    Set::Vector q2_velocity = Set::Vector::Zero();
    Set::Vector q3_velocity = Set::Vector::Zero();
    Set::Vector q4_velocity = Set::Vector::Zero();
};

} // namespace IC

#endif // IC_RIEMANN2D_H_
