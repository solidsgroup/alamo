#ifndef INTEGRATOR_SCIMITARX_H // Include guards
#define INTEGRATOR_SCIMITARX_H // 

// Standard library includes
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <vector>
#include <stdexcept>
#include <memory>

#include "Integrator/Integrator.H"
#include "IO/ParmParse.H"
#include "BC/BC.H"
#include "BC/Nothing.H"
#include "IC/IC.H"
#include "IC/Sphere.H"
#include "IC/Constant.H"
#include "IC/SodShock.H"
#include "IC/Laminate.H"
#include "IC/Expression.H"
#include "Util/ScimitarX_Util.H"
#include "Numeric/TimeStepper.H"
#include "Numeric/FluxHandler.H"
#include "Numeric/WENOReconstruction.H"
#include "Numeric/IntegratorVariableAccessLayer.H"
#include "Numeric/NumericFactory.H"

namespace Numeric {
    template<typename IntegratorDerivedClass> class FluxHandler;
    template<typename IntegratorDerivedClass> class FluxMethod;
    template<typename IntegratorDerivedClass> class FluxReconstruction;
    template<typename IntegratorDerivedClass> class FirstOrderReconstruction;
    template<typename IntegratorDerivedClass, typename VariantTraits> class WENOReconstruction;
    template<typename IntegratorDerivedClass> class LocalLaxFriedrichsMethod;
    template<typename IntegratorDerivedClass> class TimeSteppingScheme;
    template<typename IntegratorDerivedClass> class EulerForwardScheme;
    template<typename IntegratorDerivedClass> class RK3Scheme;
}

namespace Integrator
{
class ScimitarX : virtual public Integrator
{
public:

    ScimitarX() {};
    ScimitarX(IO::ParmParse& pp);
    static void Parse(ScimitarX& value, IO::ParmParse& pp);

    virtual ~ScimitarX()
    {
        delete bc_PVec;
        delete bc_Pressure;

        delete ic_PVec;
        delete ic_Pressure;

    }

    // Enum for Solver Types
    enum class SolverType {
        SolveCompressibleEuler,
        SolveElastoPlastic,
        SolveFiveEquationModel
    };


    // Enum for Euler Field Variables
    enum class EulerFieldVariables {
        DENS,   // Density
        UVEL,   // U-velocity
        VVEL,   // V-velocity
        WVEL,   // W-velocity
        IE      // Specific Internal Energy
    };


    // Maps for SolverTypes string to Enum conversions
    static const std::map<std::string, SolverType> stringToSolverType;

    // Helper function to convert an enum class to int
    template <typename Enum>
    static constexpr int toInt(Enum e) noexcept {
        return static_cast<int>(e);
    }

    // Singleton Class for SolverType Variable Management
    class SolverTypeManager {
    public:
        // Get the Singleton Instance
        static SolverTypeManager& getInstance() {
            static SolverTypeManager instance; // Thread-safe in C++11 and later
            return instance;
        }

        // Get Variables Mapped to Solver Type
        const std::vector<ScimitarX::EulerFieldVariables>& getVariables(ScimitarX::SolverType solverType) const {
            auto it = solverTypeVariables.find(solverType);
            if (it == solverTypeVariables.end()) {
                throw std::runtime_error("Invalid SolverType: No variables mapped.");
            }
            return it->second;
        }

        Util::ScimitarX_Util::getVariableIndex computeAndAssignVariableIndices(ScimitarX::SolverType& solverType) {
            Util::ScimitarX_Util::getVariableIndex result;
            std::set<int> allVariables;

            // Retrieve variables for the given solver type
            const auto& variables = getVariables(solverType);
            for (const auto& var : variables) {
                allVariables.insert(ScimitarX::toInt(var));
            }

            // Assign indices dynamically
            int index = 0;
            for (int var : allVariables) {
                result.variableIndexMap[var] = index++;
            }

            // Cache toInt results for clarity
            const int dens_idx = ScimitarX::toInt(ScimitarX::EulerFieldVariables::DENS);
            const int uvel_idx = ScimitarX::toInt(ScimitarX::EulerFieldVariables::UVEL);
#if AMREX_SPACEDIM >= 2
            const int vvel_idx = ScimitarX::toInt(ScimitarX::EulerFieldVariables::VVEL);
#endif
#if AMREX_SPACEDIM == 3
            const int wvel_idx = ScimitarX::toInt(ScimitarX::EulerFieldVariables::WVEL);
#endif
            const int ie_idx = ScimitarX::toInt(ScimitarX::EulerFieldVariables::IE);

            result.NVAR_MAX = index;
            result.DENS = result.variableIndexMap[dens_idx];
            result.UVEL = result.variableIndexMap[uvel_idx];
#if AMREX_SPACEDIM >= 2
            result.VVEL = result.variableIndexMap[vvel_idx];
#endif
#if AMREX_SPACEDIM == 3
            result.WVEL = result.variableIndexMap[wvel_idx];
#endif
            result.IE = result.variableIndexMap[ie_idx];

            return result;
        }

    private:
        SolverTypeManager() {
            // Initialize the map of solver type to variables
            solverTypeVariables[ScimitarX::SolverType::SolveCompressibleEuler] = {
                ScimitarX::EulerFieldVariables::DENS,
                ScimitarX::EulerFieldVariables::UVEL,
#if AMREX_SPACEDIM >= 2
                ScimitarX::EulerFieldVariables::VVEL,
#endif
#if AMREX_SPACEDIM == 3
                ScimitarX::EulerFieldVariables::WVEL,
#endif
                ScimitarX::EulerFieldVariables::IE
            };
        }

        SolverTypeManager(const SolverTypeManager&) = delete; // Prevent copy
        SolverTypeManager& operator=(const SolverTypeManager&) = delete; // Prevent assignment

        std::map<ScimitarX::SolverType, std::vector<ScimitarX::EulerFieldVariables>> solverTypeVariables;
    };

    static Util::ScimitarX_Util::getVariableIndex variableIndex; // 

    static IO::ParmParse setupPVecBoundaryConditions(IO::ParmParse& pp, const Util::ScimitarX_Util::getVariableIndex& variableIndex);


protected:

// Solver setup and control functions
    void Initialize(int lev) override;
    void TimeStepBegin(Set::Scalar time, int iter) override;
    void TimeStepComplete(Set::Scalar time, int lev) override;
    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override;
    void AdvanceInTimeWithoutStiffTerms(int lev, Set::Scalar time, Set::Scalar dt);
    void TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;
    void Regrid(int lev, Set::Scalar time) override;
    Set::Scalar GetTimeStep();
    void ComputeAndSetNewTimeStep();
    void ComputeConservedVariables(int lev);
    void UpdateSolutions(int lev);
    void ApplyBoundaryConditions(int lev, Set::Scalar time);

protected:
// Template function to specialize GetFluxVectors based on SolverType
template <SolverType solverType>
void ComputeConservedVariables(int lev) {

    //for (amrex::MFIter mfi(*QVec_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    for (amrex::MFIter mfi(*QVec_mf[lev], false); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto const& p_arr = PVec_mf.Patch(lev, mfi);
        auto const& q_arr = QVec_mf.Patch(lev, mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            Set::Scalar rho = p_arr(i, j, k, variableIndex.DENS);
            Set::Scalar uvel = p_arr(i, j, k, variableIndex.UVEL);
            Set::Scalar vvel = p_arr(i, j, k, variableIndex.VVEL);
#if (AMREX_SPACEDIM == 3)
            Set::Scalar wvel = p_arr(i, j, k, variableIndex.WVEL);
#endif
            Set::Scalar internal_energy = p_arr(i, j, k, variableIndex.IE);

            // Total internal energy
            Set::Scalar kinetic_energy = 0.5 * (uvel * uvel + vvel * vvel
#if (AMREX_SPACEDIM == 3)
                                            + wvel * wvel
#endif
                                            );
            Set::Scalar total_internal_energy = internal_energy + kinetic_energy;

            // Set Q vector components (solver-type-specific handling can go here)
            q_arr(i, j, k, variableIndex.DENS) = rho;
            q_arr(i, j, k, variableIndex.UVEL) = rho * uvel;
            q_arr(i, j, k, variableIndex.VVEL) = rho * vvel;
#if (AMREX_SPACEDIM == 3)
            q_arr(i, j, k, variableIndex.WVEL) = rho * wvel;
#endif
            q_arr(i, j, k, variableIndex.IE) = rho * total_internal_energy;

        });
        
    }

}

    // Template function for updating conservative variables to primitive variables
    template <SolverType solverType>
    void UpdateSolutions(int lev) {

        //for (amrex::MFIter mfi(*QVec_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        for (amrex::MFIter mfi(*QVec_mf[lev], false); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.validbox();

            auto const& q_arr = QVec_mf.Patch(lev, mfi);  // Conservative variables array (QVec)
            auto const& p_arr = PVec_mf.Patch(lev, mfi);  // Primitive variables array (PVec)
            auto const& pressure_arr = Pressure_mf.Patch(lev, mfi);  // Pressure array


            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // Extract conservative variables from QVec
                Set::Scalar rho = q_arr(i, j, k, variableIndex.DENS);      // Density


                Set::Scalar uvel = q_arr(i, j, k, variableIndex.UVEL) / rho;  // Velocity in x-direction
                Set::Scalar vvel = q_arr(i, j, k, variableIndex.VVEL) / rho;  // Velocity in y-direction
#if (AMREX_SPACEDIM == 3)
                Set::Scalar wvel = q_arr(i, j, k, variableIndex.WVEL) / rho;  // Velocity in z-direction
#endif


                Set::Scalar total_internal_energy = q_arr(i, j, k, variableIndex.IE) / rho;  // Total specific energy

                // Compute kinetic energy
                Set::Scalar kinetic_energy = 0.5 * (uvel * uvel + vvel * vvel
#if (AMREX_SPACEDIM == 3)
                                                + wvel * wvel
#endif
                                                );

                // Compute internal energy
                Set::Scalar internal_energy = (total_internal_energy - kinetic_energy);

                // Assume ideal gas law: p = (gamma - 1) * rho * internal_energy
                const Set::Scalar gamma = 1.4;  // Ideal gas constant for air (can be adjusted)
                Set::Scalar pressure = (gamma - 1.0) * rho * internal_energy;

                // Update the primitive variables array (PVec)
                p_arr(i, j, k, variableIndex.DENS) = rho;
                p_arr(i, j, k, variableIndex.UVEL) = uvel;
                p_arr(i, j, k, variableIndex.VVEL) = vvel;
#if (AMREX_SPACEDIM == 3)
                p_arr(i, j, k, variableIndex.WVEL) = wvel;
#endif
                // Also update internal energy
                p_arr(i, j, k, variableIndex.IE) = internal_energy;

                pressure_arr(i, j, k) = pressure;  // Store pressure in Pressure_mf for post-processing


            });

        }
    }

protected:
    Set::Field<Set::Scalar> QVec_mf;          // Q Vector holding all the field variables in Conservative Form (current timestep)
    Set::Field<Set::Scalar> QVec_old_mf;      // Q Vector holding all the field variables in Conservative Form (previous timestep)
    Set::Field<Set::Scalar> XFlux_mf;         // X-dir Flux Vector holding all the field variables (current timestep)
    Set::Field<Set::Scalar> YFlux_mf;         // Y-dir Flux Vector holding all the field variables (current timestep)
    Set::Field<Set::Scalar> ZFlux_mf;         // Z-dir Flux Vector holding all the field variables (current timestep)
//    Set::Field<Set::Scalar> Source_mf;         // Source term Vector for all the field variables (current timestep)
    Set::Field<Set::Scalar> PVec_mf;      // Primitive Vector holding all the field variables in primitive form (current timestep)
    Set::Field<Set::Scalar> Pressure_mf;  // Multifab for pressure variable

    int number_of_components;            // Number of components
    int number_of_ghost_cells = 1;           // Number of ghost cells

    //Add FluxHandler member to ScimitarX class declaration
    std::shared_ptr<Numeric::FluxHandler<ScimitarX>> fluxHandler;

    //Add timeStepper member to ScimitarX class declaration
    std::shared_ptr<Numeric::TimeStepper<ScimitarX>> timeStepper;
    
    // Remove duplicate function in ScimitarX.H, keep only one definition:
    std::shared_ptr<Numeric::SolverCapabilities> GetSolverCapabilities() const {
         // Create a CompressibleEulerCapabilities instance
         return std::make_shared<Numeric::CompressibleEuler::CompressibleEulerCapabilities>();
    }

    template <typename IntegratorDerivedClass>
    friend class Numeric::FluxReconstruction;

    template <typename IntegratorDerivedClass>
    friend class Numeric::FluxMethod;    

    template <typename IntegratorDerivedClass>
    friend class Numeric::FluxHandler;    

    template <typename IntegratorDerivedClass>
    friend class Numeric::FirstOrderReconstruction;

    template <typename IntegratorDerivedClass, typename VariantTraits>
    friend class Numeric::WENOReconstruction;

    template <typename IntegratorDerivedClass>
    friend class Numeric::LocalLaxFriedrichsMethod;

    template <typename IntegratorDerivedClass>
    friend class Numeric::TimeSteppingScheme;

    template <typename IntegratorDerivedClass>
    friend class Numeric::EulerForwardScheme;

    template <typename IntegratorDerivedClass>
    friend class Numeric::RK3Scheme; 
 
    friend class Numeric::CompressibleEuler::CompressibleEulerVariableAccessor;

private:

    Set::Scalar refinement_threshold;

    IC::IC<Set::Scalar>* ic_PVec = nullptr;            // Object used to initialize Q Vector
    IC::IC<Set::Scalar>* ic_Pressure = nullptr;       // Object used to initialize Pressure
    BC::BC<Set::Scalar>* bc_PVec = nullptr;            // Object used to provide Boundary Condition for Q Vector
    BC::BC<Set::Scalar>* bc_Pressure = nullptr;            // Object used to provide Boundary Condition for Q Vector
    BC::Nothing bc_nothing;

    SolverType solverType;                //solverType 

  // Solver Configuration
    Numeric::FluxReconstructionType reconstruction_method = 
        Numeric::FluxReconstructionType::FirstOrder;
    Numeric::FluxScheme flux_scheme = 
        Numeric::FluxScheme::LocalLaxFriedrichs;
    Numeric::TimeSteppingSchemeType temporal_scheme = 
        Numeric::TimeSteppingSchemeType::ForwardEuler;
    Numeric::ReconstructionMode variable_space = 
        Numeric::ReconstructionMode::Primitive;

    Numeric::WenoVariant  weno_variant = Numeric::WenoVariant::WENOJS5; 

    std::shared_ptr<Numeric::GenericVariableAccessor> variable_accessor;

    Set::Scalar cflNumber;        // CFL number (read from input)

}; //ScimitarX Class

    inline const std::map<std::string, ScimitarX::SolverType> ScimitarX::stringToSolverType = {
        {"SolveCompressibleEuler", ScimitarX::SolverType::SolveCompressibleEuler},
        {"SolveElastoPlastic", ScimitarX::SolverType::SolveElastoPlastic},
        {"SolveFiveEquationModel", ScimitarX::SolverType::SolveFiveEquationModel}
    };


} // namespace Integrator
#endif
