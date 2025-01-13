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

    std::string constructKey(const std::string& prefix, const std::string& component_name, const std::string& property, const std::string& side) {
        return prefix + component_name + "." + property + "." + side;
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

    // Flux Reconstruction Schemes
    enum class FluxReconstruction {
        FirstOrder,               // Basic first-order reconstruction
        ThirdOrderENO,            // Third-order Essentially Non-Oscillatory (ENO)
        FifthOrderWENO,           // Fifth-order Weighted Essentially Non-Oscillatory (WENO)
        FifthOrderWENOWithChar    // Fifth-order WENO with Characteristic decomposition
    };

    // Flux Schemes
    enum class FluxScheme {
        LocalLaxFriedrichs,   // Local Lax-Friedrichs (Rusanov) flux
        HLLCFlux             // Harten-Lax-van Leer Contact (HLLC) flux
    };

    enum class TimeSteppingScheme {
        ForwardEuler,           // First-order explicit Euler method
        RK3                    // Third-order Runge-Kutta method
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

    // Subclass to manage feature maps
    class FeatureMaps {
    public:
        // Singleton access to ensure one instance
        static FeatureMaps& getInstance() {
            static FeatureMaps instance;
            return instance;
        }

        // Accessor methods for the maps
        const std::map<std::string, FluxReconstruction>& getFluxReconstructionMap() const {
            return fluxReconstructionMap;
        }

        const std::map<std::string, FluxScheme>& getFluxSchemeMap() const {
            return fluxSchemeMap;
        }

        const std::map<std::string, TimeSteppingScheme>& getTimeSteppingSchemeMap() const {
            return timeSteppingSchemeMap;
        }

    private:
        // Private constructor to initialize the maps
        FeatureMaps() {
            fluxReconstructionMap = {
                {"FirstOrder", FluxReconstruction::FirstOrder},
                {"ThirdOrderENO", FluxReconstruction::ThirdOrderENO},
                {"FifthOrderWENO", FluxReconstruction::FifthOrderWENO},
                {"FifthOrderWENOWithChar", FluxReconstruction::FifthOrderWENOWithChar}
            };

            fluxSchemeMap = {
                {"LocalLaxFriedrichs", FluxScheme::LocalLaxFriedrichs},
                {"HLLCFlux", FluxScheme::HLLCFlux}
            };

            timeSteppingSchemeMap = {
                {"ForwardEuler", TimeSteppingScheme::ForwardEuler},
                {"RK3", TimeSteppingScheme::RK3}
            };
        }

        // Disable copy and assignment
        FeatureMaps(const FeatureMaps&) = delete;
        FeatureMaps& operator=(const FeatureMaps&) = delete;

        // The maps for different schemes
        std::map<std::string, FluxReconstruction> fluxReconstructionMap;
        std::map<std::string, FluxScheme> fluxSchemeMap;
        std::map<std::string, TimeSteppingScheme> timeSteppingSchemeMap;
    };

    // Accessors to the singleton
    static const FeatureMaps& getFeatureMaps() {
        return FeatureMaps::getInstance();
    }


protected:

// Solver setup and control functions
    void Initialize(int lev) override;
    void TimeStepBegin(Set::Scalar time, int iter) override;
    void TimeStepComplete(Set::Scalar time, int lev) override;
    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override {
        Util::Message(INFO, "Starting ScimitarX::Advance for level: " + std::to_string(lev));

        // Advance the solution without stiff source terms
        AdvanceInTimeWithoutStiffTerms(lev, time, dt);

        Util::Message(INFO, "Completed ScimitarX::Advance for level: " + std::to_string(lev));
    }
    // Function to advance hyperbolic balance equations without stiff terms
    void AdvanceInTimeWithoutStiffTerms(int lev, Set::Scalar time, Set::Scalar dt) {
        Util::Message(INFO, "Starting AdvanceInTimeWithoutStiffTerms for level: " + std::to_string(lev));

    // Determine the solver type and handle different governing equations
    switch (solverType) {
        case SolverType::SolveCompressibleEuler: {
            Util::Message(INFO, "SolverType: Compressible Euler");
            break;
        }
        case SolverType::SolveElastoPlastic: {
            Util::Message(INFO, "SolverType: Elasto-Plastic");
            break;
        }
        case SolverType::SolveFiveEquationModel: {
            Util::Message(INFO, "SolverType: Five-Equation Model");
            break;
        }
        default:
            Util::Abort(__FILE__, __func__, __LINE__, "Unknown SolverType.");
    }

    // Ensure that the fluxHandler and timeStepper are properly initialized
   /* if (!fluxHandler) {
        Util::Abort(__FILE__, __func__, __LINE__, "FluxHandler is not initialized.");
    }
    if (!timeStepper) {
        Util::Abort(__FILE__, __func__, __LINE__, "TimeStepper is not initialized.");
    }*/

    // Determine the time-stepping scheme
    switch (temporal_scheme) {
        case TimeSteppingScheme::ForwardEuler: {
            Util::Message(INFO, "TimeSteppingScheme: Forward Euler");
           
            Util::ScimitarX_Util::Debug debug; 

            fluxHandler->SetReconstruction(std::make_shared<Numeric::FirstOrderReconstruction<ScimitarX>>());
            fluxHandler->SetFluxMethod(std::make_shared<Numeric::LocalLaxFriedrichsMethod<ScimitarX>>());

            timeStepper->SetTimeSteppingScheme(std::make_shared<Numeric::EulerForwardScheme<ScimitarX>>());

            int numStages = timeStepper->GetNumberOfStages();
            // One-stage loop for Forward Euler
            for (int stage = 0; stage < numStages; ++stage) {
                

                //debug.DebugGhostCellsAndBoundary(*Pressure_mf[lev], geom[lev], lev, "Pressure", "Ghost Cell Debug", true);

                // 1. Compute Conserved Variables
                ComputeConservedVariables<SolverType::SolveCompressibleEuler>(lev);

                // 2. Perform flux reconstruction and compute fluxes in all directions
                fluxHandler->ConstructFluxes(lev, this);

                ApplyBoundaryConditions(lev, time);
                
                // 3. Compute sub-step using the chosen time-stepping scheme
                timeStepper->ComputeSubStep(lev, dt, stage, this);

                // 4. Update solution from conservative to primitive variables
                UpdateSolutions<SolverType::SolveCompressibleEuler>(lev);

                //Util::Message(INFO, "Before Applying Boundary Condition");

                //ApplyBoundaryConditions(lev, time);

                //Util::Message(INFO, "After Applying Boundary Condition");
                //debug.DebugGhostCellsAndBoundary(*Pressure_mf[lev], geom[lev], lev, "QVec", "Ghost Cell Debug", true);

                // 5. Swap the old QVec Fab with new one so that we can use the new one for next substep
                std::swap(*QVec_old_mf[lev], *QVec_mf[lev]);
            }
            break;
        }

        case TimeSteppingScheme::RK3: {
            Util::Message(INFO, "TimeSteppingScheme: RK3");
            fluxHandler->SetReconstruction(std::make_shared<Numeric::FirstOrderReconstruction<ScimitarX>>());
            fluxHandler->SetFluxMethod(std::make_shared<Numeric::LocalLaxFriedrichsMethod<ScimitarX>>());

            timeStepper->SetTimeSteppingScheme(std::make_shared<Numeric::EulerForwardScheme<ScimitarX>>());

            int numStages = timeStepper->GetNumberOfStages();

            for (int stage = 0; stage < numStages; ++stage) {
                ComputeConservedVariables<SolverType::SolveCompressibleEuler>(lev);
                fluxHandler->ConstructFluxes(lev, this);
                timeStepper->ComputeSubStep(lev, dt, stage, this);
                UpdateSolutions<SolverType::SolveCompressibleEuler>(lev);
                ApplyBoundaryConditions(lev, time);
                std::swap(*QVec_old_mf[lev], *QVec_mf[lev]);
            }
            break;
        }

        default:
            Util::Abort(__FILE__, __func__, __LINE__, "Unknown TimeSteppingScheme.");
    }

        Util::Message(INFO, "Completed AdvanceInTimeWithoutStiffTerms for level: " + std::to_string(lev));
    }

    void TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;
    void Regrid(int lev, Set::Scalar time) override;
    Set::Scalar GetTimeStep();
    Set::Scalar ComputeAndSetNewTimeStep();
    void ComputeConservedVariables(int lev);
    void UpdateSolutions(int lev);
    void ApplyBoundaryConditions(int lev, Set::Scalar time);

protected:
// Template function to specialize GetFluxVectors based on SolverType
template <SolverType solverType>
void ComputeConservedVariables(int lev) {
    Util::Message(INFO, "Starting GetFluxVectors for SolverType-specific computation");

       Util::ScimitarX_Util::Debug debug;

        //debug.SetTargetDebugLocation(0.5, 0.2, 0.0, geom[lev]);
        debug.SetTargetDebugLocationIndices(48, -2, 0);

    //for (amrex::MFIter mfi(*QVec_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    for (amrex::MFIter mfi(*QVec_mf[lev], false); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto const& p_arr = PVec_mf.Patch(lev, mfi);
        auto const& q_arr = QVec_mf.Patch(lev, mfi);
        auto const& pressure_arr = Pressure_mf.Patch(lev,mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            Set::Scalar rho = p_arr(i, j, k, variableIndex.DENS);
            Set::Scalar uvel = p_arr(i, j, k, variableIndex.UVEL);
            Set::Scalar vvel = p_arr(i, j, k, variableIndex.VVEL);
#if (AMREX_SPACEDIM == 3)
            Set::Scalar wvel = p_arr(i, j, k, variableIndex.WVEL);
#endif
            Set::Scalar ie = p_arr(i, j, k, variableIndex.IE);

            // Total internal energy
            Set::Scalar kinetic_energy = 0.5 * (uvel * uvel + vvel * vvel
#if (AMREX_SPACEDIM == 3)
                                             + wvel * wvel
#endif
                                             );
            Set::Scalar tot_int_energy = ie + kinetic_energy;

            debug.DebugValuesIfTarget(i, j, k, uvel, "u-velocity", "ComputeConservedVariables", false);
            debug.DebugValuesIfTarget(i, j, k, vvel, "v-velocity", "ComputeConservedVariables", false);

            // Set Q vector components (solver-type-specific handling can go here)
            q_arr(i, j, k, variableIndex.DENS) = rho;
            q_arr(i, j, k, variableIndex.UVEL) = rho * uvel;
            q_arr(i, j, k, variableIndex.VVEL) = rho * vvel;
#if (AMREX_SPACEDIM == 3)
            q_arr(i, j, k, variableIndex.WVEL) = rho * wvel;
#endif
            q_arr(i, j, k, variableIndex.IE) = rho * tot_int_energy;

        });
    
    
    }

    Util::Message(INFO, "Flux Vectors computed successfully for SolverType-specific routine");
}

    // Template function for updating conservative variables to primitive variables
    template <SolverType solverType>
    void UpdateSolutions(int lev) {
        // Access geometry at the specified level
        const Set::Scalar* DX = geom[lev].CellSize();
        Util::ScimitarX_Util::Debug debug;

        //debug.SetTargetDebugLocation(0.5, 0.2, 0.0, geom[lev]);
        debug.SetTargetDebugLocationIndices(48, -2, 0);

        //for (amrex::MFIter mfi(*QVec_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        for (amrex::MFIter mfi(*QVec_mf[lev], false); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.validbox();

            auto const& q_arr = QVec_mf.Patch(lev, mfi);  // Conservative variables array (QVec)
            auto const& p_arr = PVec_mf.Patch(lev, mfi);  // Primitive variables array (PVec)
            auto const& pressure_arr = Pressure_mf.Patch(lev, mfi);  // Pressure array


            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // Extract conservative variables from QVec
                Set::Scalar rho = q_arr(i, j, k, variableIndex.DENS);      // Density

                //debug.CheckNaNAndAbort(i, j, k, rho, "Density", "UpdateSolutions", true);

                Set::Scalar uvel = q_arr(i, j, k, variableIndex.UVEL) / rho;  // Velocity in x-direction
                Set::Scalar vvel = q_arr(i, j, k, variableIndex.VVEL) / rho;  // Velocity in y-direction
#if (AMREX_SPACEDIM == 3)
                Set::Scalar wvel = q_arr(i, j, k, variableIndex.WVEL) / rho;  // Velocity in z-direction
#endif

                debug.DebugValuesIfTarget(i, j, k, uvel, "u-velocity", "UpdateSolutions", false);

                Set::Scalar total_energy = q_arr(i, j, k, variableIndex.IE);  // Total specific energy

                // Compute kinetic energy
                Set::Scalar kinetic_energy = 0.5 * (uvel * uvel + vvel * vvel
#if (AMREX_SPACEDIM == 3)
                                                + wvel * wvel
#endif
                                                );

                // Compute internal energy
                Set::Scalar internal_energy = (total_energy - rho*kinetic_energy)/rho;

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
                pressure_arr(i, j, k) = pressure;  // Store pressure in Pressure_mf for post-processing

                // Also update internal energy explicitly if needed
                p_arr(i, j, k, variableIndex.IE) = internal_energy;

                debug.DebugValuesIfTarget(i, j, k, internal_energy, "internal energy", "GetFluxVector", false);

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

    Numeric::FluxHandler<ScimitarX>* fluxHandler = nullptr;  // Flux handler instance
    Numeric::TimeStepper<ScimitarX>* timeStepper = nullptr;  // Time stepper instance

    template <typename T>
    friend class Numeric::FluxReconstruction;

    template <typename T>
    friend class Numeric::FluxMethod;    

    template <typename T>
    friend class Numeric::FluxHandler;    

    template <typename T>
    friend class Numeric::FirstOrderReconstruction;

    template <typename T>
    friend class Numeric::LocalLaxFriedrichsMethod;

    template <typename T>
    friend class Numeric::TimeSteppingScheme;

    template <typename T>
    friend class Numeric::EulerForwardScheme;

    template <typename T>
    friend class Numeric::RK3Scheme;    

private:

    Set::Scalar refinement_threshold;

    IC::IC* ic_PVec = nullptr;            // Object used to initialize Q Vector
    IC::IC* ic_Pressure = nullptr;       // Object used to initialize Pressure
    BC::BC<Set::Scalar>* bc_PVec = nullptr;            // Object used to provide Boundary Condition for Q Vector
    BC::BC<Set::Scalar>* bc_Pressure = nullptr;            // Object used to provide Boundary Condition for Q Vector
    BC::Nothing bc_nothing;

    SolverType solverType;                //solverType 

    FluxReconstruction reconstruction_method;  // Selected flux reconstruction scheme
    FluxScheme flux_scheme;                    // Selected flux scheme
    TimeSteppingScheme temporal_scheme;         // Selected time-stepping scheme

    Set::Scalar cflNumber;        // CFL number (read from input)

}; //ScimitarX Class

    inline const std::map<std::string, ScimitarX::SolverType> ScimitarX::stringToSolverType = {
        {"SolveCompressibleEuler", ScimitarX::SolverType::SolveCompressibleEuler},
        {"SolveElastoPlastic", ScimitarX::SolverType::SolveElastoPlastic},
        {"SolveFiveEquationModel", ScimitarX::SolverType::SolveFiveEquationModel}
    };


} // namespace Integrator
#endif

