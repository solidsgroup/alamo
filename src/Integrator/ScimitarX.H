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
#include "IC/Shock.H"
#include "IC/Laminate.H"
#include "IC/Expression.H"
#include "Util/ScimitarX_Util.H"
#include "Numeric/NumericTypes.H"
#include "Numeric/NumericFactory.H"
#include "Numeric/SolverCapabilities.H"
#include "Numeric/IntegratorVariableAccessLayer.H"
#include "Numeric/FluxHandler.H"
#include "Numeric/WENOReconstruction.H"
#include "Numeric/TimeStepper.H"

namespace Numeric {
    template<typename IntegratorDerivedClass> class FluxHandler;
    template<typename IntegratorDerivedClass> class FluxMethod;
    template<typename IntegratorDerivedClass> class FluxReconstruction;
    template<typename IntegratorDerivedClass> class FirstOrderReconstruction;
    template<typename IntegratorDerivedClass, typename VariantTraits> class WENOReconstruction;
    template<typename IntegratorDerivedClass> class LocalLaxFriedrichsMethod;
    template<typename IntegratorDerivedClass> class HLLCMethod;
    template<typename IntegratorDerivedClass> class TimeSteppingScheme;
    template<typename IntegratorDerivedClass> class EulerForwardScheme;
    template<typename IntegratorDerivedClass> class RK3Scheme;
}

namespace Integrator
{
class ScimitarX : virtual public Integrator
{
public:
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

    // Enum for ElatoPlastic Field Variables
    enum class ElastoPlasticFieldVariables {
        DENS,         // Density
        UVEL,         // U-velocity
        VVEL,         // V-velocity
        WVEL,         // W-velocity
        IE,           // Specific Internal Energy    
        SXX,          //Sxx scalar
        SXY,          //Sxy scalar
        SYY,          //Syy scalar
        SXZ,          //Sxz scalar
        SYZ,          //Syz scalar
        SZZ,          //Szz scalar
        EPSBAR,       //Effective Plastic Strain
        IE_ELASTIC    //Elastic energy for temperature calculations
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


    // Get Variables Mapped to Solver Type - simplified with type-safe vector
    template<typename EnumType>
    const std::vector<EnumType>& getVariablesForType(SolverType solverType) const {
        switch (solverType) {
            case SolverType::SolveCompressibleEuler:
                return getEulerVariables();
            case SolverType::SolveElastoPlastic:
                return getElastoPlasticVariables();
            case SolverType::SolveFiveEquationModel:
                throw std::runtime_error("Five equation model not yet implemented");
            default:
                throw std::runtime_error("Invalid SolverType: No variables mapped.");
            }
        }

    Numeric::GenericVariableAccessor::VariableIndices computeAndAssignVariableIndices(SolverType& solverType) {
        Numeric::GenericVariableAccessor::VariableIndices result;
        
        switch (solverType) {
            case SolverType::SolveCompressibleEuler:
                computeEulerIndices(result);
                break;
                
            case SolverType::SolveElastoPlastic:
                computeElastoPlasticIndices(result);
                break;
                
            case SolverType::SolveFiveEquationModel:
                throw std::runtime_error("Five equation model not yet implemented");
                
            default:
                throw std::runtime_error("Unsupported solver type for variable indexing");
        }

        return result;
    }

    private:

    SolverTypeManager() {
        // Initialize the vectors of variables for each solver type
        initializeEulerVariables();
        initializeElastoPlasticVariables();
    }

    SolverTypeManager(const SolverTypeManager&) = delete; // Prevent copy
    SolverTypeManager& operator=(const SolverTypeManager&) = delete; // Prevent assignment

    void initializeEulerVariables() {
        eulerVariables = {
            EulerFieldVariables::DENS,
            EulerFieldVariables::UVEL,
#if AMREX_SPACEDIM >= 2
            EulerFieldVariables::VVEL,
#endif
#if AMREX_SPACEDIM == 3
            EulerFieldVariables::WVEL,
#endif
            EulerFieldVariables::IE
        };
    }

    void initializeElastoPlasticVariables() {
        elastoPlasticVariables = {
            ElastoPlasticFieldVariables::DENS,
            ElastoPlasticFieldVariables::UVEL,
#if AMREX_SPACEDIM >= 2
            ElastoPlasticFieldVariables::VVEL,
#endif
#if AMREX_SPACEDIM == 3
            ElastoPlasticFieldVariables::WVEL,
#endif
            ElastoPlasticFieldVariables::IE,
            ElastoPlasticFieldVariables::SXX,
            ElastoPlasticFieldVariables::SYY,
            ElastoPlasticFieldVariables::SZZ,
            ElastoPlasticFieldVariables::SXY,
#if AMREX_SPACEDIM == 3
            ElastoPlasticFieldVariables::SXZ,
            ElastoPlasticFieldVariables::SYZ,
#endif
            ElastoPlasticFieldVariables::EPSBAR,
            ElastoPlasticFieldVariables::IE_ELASTIC
        };
    }

    void computeEulerIndices(Numeric::GenericVariableAccessor::VariableIndices& result) {
        std::set<int> allVariables;
        
        // Add all Euler variables to the set
        for (const auto& var : eulerVariables) {
            allVariables.insert(toInt(var));
        }

        // Assign indices dynamically
        int index = 0;
        for (int var : allVariables) {
            result.variableIndexMap[var] = index++;
        }

        // Set up specific indices
        result.NVAR_MAX = index;
        result.DENS = result.variableIndexMap[toInt(EulerFieldVariables::DENS)];
        result.UVEL = result.variableIndexMap[toInt(EulerFieldVariables::UVEL)];
#if AMREX_SPACEDIM >= 2
        result.VVEL = result.variableIndexMap[toInt(EulerFieldVariables::VVEL)];
#endif
#if AMREX_SPACEDIM == 3
        result.WVEL = result.variableIndexMap[toInt(EulerFieldVariables::WVEL)];
#endif
        result.IE = result.variableIndexMap[toInt(EulerFieldVariables::IE)];
    }


    void computeElastoPlasticIndices(Numeric::GenericVariableAccessor::VariableIndices& result) {
        std::set<int> allVariables;
        
        // Add all ElastoPlastic variables to the set
        for (const auto& var : elastoPlasticVariables) {
            allVariables.insert(toInt(var));
        }

        // Assign indices dynamically
        int index = 0;
        for (int var : allVariables) {
            result.variableIndexMap[var] = index++;
        }

        // Set up specific indices
        result.NVAR_MAX = index;
        result.DENS = result.variableIndexMap[toInt(ElastoPlasticFieldVariables::DENS)];
        result.UVEL = result.variableIndexMap[toInt(ElastoPlasticFieldVariables::UVEL)];
#if AMREX_SPACEDIM >= 2
        result.VVEL = result.variableIndexMap[toInt(ElastoPlasticFieldVariables::VVEL)];
#endif
#if AMREX_SPACEDIM == 3
        result.WVEL = result.variableIndexMap[toInt(ElastoPlasticFieldVariables::WVEL)];
#endif
        result.IE = result.variableIndexMap[toInt(ElastoPlasticFieldVariables::IE)];
        
        // ElastoPlastic specific variables
        result.SXX = result.variableIndexMap[toInt(ElastoPlasticFieldVariables::SXX)];
        result.SYY = result.variableIndexMap[toInt(ElastoPlasticFieldVariables::SYY)];
        result.SZZ = result.variableIndexMap[toInt(ElastoPlasticFieldVariables::SZZ)];
        result.SXY = result.variableIndexMap[toInt(ElastoPlasticFieldVariables::SXY)];
#if AMREX_SPACEDIM == 3
        result.SXZ = result.variableIndexMap[toInt(ElastoPlasticFieldVariables::SXZ)];
        result.SYZ = result.variableIndexMap[toInt(ElastoPlasticFieldVariables::SYZ)];
#endif
        result.EPSBAR = result.variableIndexMap[toInt(ElastoPlasticFieldVariables::EPSBAR)];
        result.IE_ELASTIC = result.variableIndexMap[toInt(ElastoPlasticFieldVariables::IE_ELASTIC)];
    }
    

    const std::vector<EulerFieldVariables>& getEulerVariables() const {
        return eulerVariables;
    }

    const std::vector<ElastoPlasticFieldVariables>& getElastoPlasticVariables() const {
        return elastoPlasticVariables;
    }

    // Member variables
    std::vector<EulerFieldVariables> eulerVariables;
    std::vector<ElastoPlasticFieldVariables> elastoPlasticVariables;

    };

    static Numeric::GenericVariableAccessor::VariableIndices variableIndex; //variableIndex is now a static member of ScimitarX 

    static IO::ParmParse setupPVecBoundaryConditions(IO::ParmParse& pp, const Numeric::GenericVariableAccessor::VariableIndices& variableIndex);

    // Constructors and destructor
    ScimitarX();
    ScimitarX(IO::ParmParse& pp);
    static void Parse(ScimitarX& value, IO::ParmParse& pp);
    virtual ~ScimitarX()
    {
        delete bc_PVec;
        delete bc_Pressure;

        delete ic_PVec;
        delete ic_Pressure;

    }

    void RegisterSolverCapabilities();
protected:
    // Required override methods from Integrator base class
    void Initialize(int lev) override;
    void TimeStepBegin(Set::Scalar time, int iter) override;
    void TimeStepComplete(Set::Scalar time, int lev) override;
    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override;
    void TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;
    void Regrid(int lev, Set::Scalar time) override;

    // ScimitarX specific methods
    void AdvanceInTimeWithoutStiffTerms(int lev, Set::Scalar time, Set::Scalar dt);
    Set::Scalar GetTimeStep();
    void ComputeAndSetNewTimeStep();
    void ApplyBoundaryConditions(int lev, Set::Scalar time);

    // Template methods for solver-type specific operations
    template <SolverType solverType>
    void ComputeConservedVariables(int lev){

    //for (amrex::MFIter mfi(*QVec_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    for (amrex::MFIter mfi(*QVec_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx_with_ghosts = mfi.growntilebox();
        auto const& p_arr = PVec_mf.Patch(lev, mfi);
        auto const& q_arr = QVec_mf.Patch(lev, mfi);

        amrex::ParallelFor(bx_with_ghosts, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
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
        for (amrex::MFIter mfi(*QVec_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
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

    // Method to get solver capabilities (centralized)
    std::shared_ptr<Numeric::SolverCapabilities> GetSolverCapabilities() const;

protected:
    // Field data structures
    Set::Field<Set::Scalar> QVec_mf;          // Conservative variables (current timestep)
    Set::Field<Set::Scalar> QVec_old_mf;      // Conservative variables (previous timestep)
    Set::Field<Set::Scalar> XFlux_mf;         // X-dir Flux Vector
    Set::Field<Set::Scalar> YFlux_mf;         // Y-dir Flux Vector
    Set::Field<Set::Scalar> ZFlux_mf;         // Z-dir Flux Vector
    Set::Field<Set::Scalar> PVec_mf;          // Primitive variables
    Set::Field<Set::Scalar> Pressure_mf;      // Pressure variable

    // Core solver parameters
    int number_of_components;                 // Number of components
    int number_of_ghost_cells;            // Number of ghost cells
    SolverType solverType;                    // Current solver type
    Set::Scalar cflNumber;                    // CFL number (read from input)

    Set::Scalar refinement_threshold;         // Threshold for AMR refinement

    // Refinement threshold values
    Set::Scalar density_refinement_criterion = 0.01;
    Set::Scalar velocity_refinement_criterion = 0.01;
    Set::Scalar pressure_refinement_criterion = 0.01;
    Set::Scalar vorticity_refinement_criterion = 0.01;
    
    // Enable/disable flags for each refinement criterion
    bool enable_density_refinement = false;
    bool enable_velocity_refinement = false;
    bool enable_pressure_refinement = false;
    bool enable_vorticity_refinement = false;
    
    // Numeric handlers (using shared_ptr for better lifetime management)
    std::shared_ptr<Numeric::FluxHandler<ScimitarX>> fluxHandler;
    std::shared_ptr<Numeric::TimeStepper<ScimitarX>> timeStepper;
    std::shared_ptr<Numeric::GenericVariableAccessor> variable_accessor;
    std::shared_ptr<Numeric::SolverCapabilities> solverCapabilities;

    // Initial and boundary conditions
    IC::IC<Set::Scalar>* ic_PVec = nullptr;            // Object used to initialize Primitive Vector
    IC::IC<Set::Scalar>* ic_Pressure = nullptr;        // Object used to initialize Pressure
    BC::BC<Set::Scalar>* bc_PVec = nullptr;            // Boundary Condition for Primitive Vector
    BC::BC<Set::Scalar>* bc_Pressure = nullptr;        // Boundary Condition for Pressure
    BC::Nothing bc_nothing;                            // Default "do-nothing" boundary condition

    // Numeric method configuration 
    Numeric::FluxReconstructionType reconstruction_method;
    Numeric::FluxScheme flux_scheme;
    Numeric::TimeSteppingSchemeType temporal_scheme;
    Numeric::ReconstructionMode variable_space;
    Numeric::WenoVariant weno_variant;

private:
    // Helper method to setup numeric components based on solverType
    void SetupNumericComponents();
    
    // Helper to validate and configure numeric methods
    void ValidateAndSetupNumerics();

    // Friend declarations for Numeric namespace access
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
    friend class Numeric::HLLCMethod;
    
    template <typename IntegratorDerivedClass>
    friend class Numeric::TimeSteppingScheme;

    template <typename IntegratorDerivedClass>
    friend class Numeric::EulerForwardScheme;

    template <typename IntegratorDerivedClass>
    friend class Numeric::RK3Scheme; 
 
    friend class Numeric::CompressibleEuler::CompressibleEulerVariableAccessor;
            
    
}; //ScimitarX Class

// Static member initialization for string to solver type mapping
inline const std::map<std::string, ScimitarX::SolverType> ScimitarX::stringToSolverType = {
    {"SolveCompressibleEuler", ScimitarX::SolverType::SolveCompressibleEuler},
    {"SolveElastoPlastic", ScimitarX::SolverType::SolveElastoPlastic},
    {"SolveFiveEquationModel", ScimitarX::SolverType::SolveFiveEquationModel}
};

} // namespace Integrator
#endif // INTEGRATOR_SCIMITARX_H
