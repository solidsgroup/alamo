#ifndef IC_RIEMANN2D_H_
#define IC_RIEMANN2D_H_

#include "IC/SodShock.H"
#include "AMReX.H"
#include "AMReX_ParmParse.H"
#include "AMReX_Geometry.H"
#include "AMReX_MultiFab.H"
#include "AMReX_Array4.H"
#include <string>
#include "IC/IC.H"
#include "IC/SodShock.H"
#include "Set/Set.H"
#include "Util/Util.H"
#include "Util/ScimitarX_Util.H"
#include "Model/Fluid/Fluid.H"


namespace IC {

using namespace Model::Fluid;

class Riemann2D : public SodShock {
public:
    // Constructor requiring variable indices
    Riemann2D(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name,
              const Util::ScimitarX_Util::getVariableIndex& precomputed_indices)
        : SodShock(_geom, pp, name, precomputed_indices) {
        initialize(pp, name);
    }

    // Constructor for pressure (no variable indices required)
    Riemann2D(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name)
        : SodShock(_geom, pp, name) {
        initialize(pp, name);
    }

    // Override Add function to handle 2D Riemann problems
    virtual void Add(const int& lev, Set::Field<Set::Scalar>& a_phi, Set::Scalar time) override {
        int ncomp = a_phi[0]->nComp();

        for (amrex::MFIter mfi(*a_phi[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            amrex::Box bx = mfi.growntilebox();
            amrex::IndexType type = a_phi[lev]->ixType();

            amrex::Array4<Set::Scalar> const& phi = a_phi[lev]->array(mfi);

            for (int n = 0; n < ncomp; ++n) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    Set::Vector x = Set::Position(i, j, k, geom[lev], type);

                    if (mf_name == "ic.pvec.riemann2d") {
                        if (x(0) < shock_xpos && x(1) < shock_ypos) {
                            assignQuadrantValues(phi, i, j, k, n, 3);
                        } else if (x(0) >= shock_xpos && x(1) < shock_ypos) {
                            assignQuadrantValues(phi, i, j, k, n, 4);
                        } else if (x(0) < shock_xpos && x(1) >= shock_ypos) {
                            assignQuadrantValues(phi, i, j, k, n, 1);
                        } else {
                            assignQuadrantValues(phi, i, j, k, n, 2);
                        }
                    } else if (mf_name == "ic.pressure.riemann2d") {
                        phi(i, j, k, n) = (x(0) < shock_xpos && x(1) < shock_ypos) ? q3_pressure :
                                          (x(0) >= shock_xpos && x(1) < shock_ypos) ? q4_pressure :
                                          (x(0) < shock_xpos && x(1) >= shock_ypos) ? q1_pressure : q2_pressure;
                    } else {
                        Util::Abort(INFO, "Unknown MultiFab name: " + mf_name);
                    }
                });
            }
        }
    }

private:
    void initialize(IO::ParmParse& pp, const std::string& name) {
        mf_name = name;

        if (mf_name == "ic.pvec.riemann2d") {
            pp.query("ic.pvec.q1.DENS", q1_density);
            pp.query("ic.pvec.q2.DENS", q2_density);
            pp.query("ic.pvec.q3.DENS", q3_density);
            pp.query("ic.pvec.q4.DENS", q4_density);

            pp.query("ic.pvec.q1.UVEL", q1_velocity[0]);
            pp.query("ic.pvec.q2.UVEL", q2_velocity[0]);
            pp.query("ic.pvec.q3.UVEL", q3_velocity[0]);
            pp.query("ic.pvec.q4.UVEL", q4_velocity[0]);

            pp.query("ic.pvec.q1.VVEL", q1_velocity[1]);
            pp.query("ic.pvec.q2.VVEL", q2_velocity[1]);
            pp.query("ic.pvec.q3.VVEL", q3_velocity[1]);
            pp.query("ic.pvec.q4.VVEL", q4_velocity[1]);

            pp.query("ic.pressure.q1", q1_pressure);
            pp.query("ic.pressure.q2", q2_pressure);
            pp.query("ic.pressure.q3", q3_pressure);
            pp.query("ic.pressure.q4", q4_pressure);
        } else if (mf_name == "ic.pressure.riemann2d") {
            pp.query("ic.pressure.q1", q1_pressure);
            pp.query("ic.pressure.q2", q2_pressure);
            pp.query("ic.pressure.q3", q3_pressure);
            pp.query("ic.pressure.q4", q4_pressure);
        }

        pp.query("shock.xpos", shock_xpos);
        pp.query("shock.ypos", shock_ypos);

        Util::Message(INFO, "Initialized Riemann2D with shock positions: x=" + std::to_string(shock_xpos) + ", y=" + std::to_string(shock_ypos));
    }

    void assignQuadrantValues(amrex::Array4<Set::Scalar> const& phi, int i, int j, int k, int n, int quadrant) const {

        Model::Fluid::Fluid fluid_model;   // Create an instance of the Fluid model
          
        Set::Scalar gamma = 1.4

        switch (quadrant) {
            case 1:
                phi(i, j, k, n) = (n == variable_indices->DENS) ? q1_density :
                                   (n == variable_indices->UVEL) ? q1_velocity[0] :
                                   (n == variable_indices->VVEL) ? q1_velocity[1] :
                                   (n == variable_indices->IE) ? fluid_model.ComputeInternalEnergyFromDensityAndPressure(q1_density, q1_pressure, gamma) : 0.0;
                break;
            case 2:
                phi(i, j, k, n) = (n == variable_indices->DENS) ? q2_density :
                                   (n == variable_indices->UVEL) ? q2_velocity[0] :
                                   (n == variable_indices->VVEL) ? q2_velocity[1] :
                                   (n == variable_indices->IE) ? fluid_model.ComputeInternalEnergyFromDensityAndPressure(q2_density, q2_pressure, gamma) : 0.0;
                break;
            case 3:
                phi(i, j, k, n) = (n == variable_indices->DENS) ? q3_density :
                                   (n == variable_indices->UVEL) ? q3_velocity[0] :
                                   (n == variable_indices->VVEL) ? q3_velocity[1] :
                                   (n == variable_indices->IE) ? fluid_model.ComputeInternalEnergyFromDensityAndPressure(q3_density, q3_pressure, gamma) : 0.0;
                break;
            case 4:
                phi(i, j, k, n) = (n == variable_indices->DENS) ? q4_density :
                                   (n == variable_indices->UVEL) ? q4_velocity[0] :
                                   (n == variable_indices->VVEL) ? q4_velocity[1] :
                                   (n == variable_indices->IE) ? fluid_model.ComputeInternalEnergyFromDensityAndPressure(q4_density, q4_density, gamma) : 0.0;
                break;
        }
    }

    // Quadrant-specific parameters
    Set::Scalar q1_density = 1.0;
    Set::Scalar q2_density = 0.125;
    Set::Scalar q3_density = 0.25;
    Set::Scalar q4_density = 0.5;

    Set::Scalar q1_pressure = 1.0;
    Set::Scalar q2_pressure = 0.1;
    Set::Scalar q3_pressure = 0.5;
    Set::Scalar q4_pressure = 0.2;

    Set::Vector q1_velocity = {0.0, 0.0};
    Set::Vector q2_velocity = {0.0, 0.0};
    Set::Vector q3_velocity = {0.0, 0.0};
    Set::Vector q4_velocity = {0.0, 0.0};
};

} // namespace IC

#endif

