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

    // Maps for string to Enum conversions
    static const std::map<std::string, SolverType> stringToSolverType;

    // Enum for Euler Field Variables
    enum class EulerFieldVariables {
        DENS,   // Density
        UVEL,   // U-velocity
        VVEL,   // V-velocity
        WVEL,   // W-velocity
        IE      // Specific Internal Energy
    };

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
            const int vvel_idx = ScimitarX::toInt(ScimitarX::EulerFieldVariables::VVEL);
            const int wvel_idx = ScimitarX::toInt(ScimitarX::EulerFieldVariables::WVEL);
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

    void Initialize(int lev) override;
    void TimeStepBegin(Set::Scalar a_time, int a_iter) override;
    void TimeStepComplete(Set::Scalar time, int lev) override;
    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override;
    void TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/) override;
    void Regrid(int lev, Set::Scalar time) override;

protected:
    Set::Field<Set::Scalar> QVec_mf;          // Q Vector holding all the field variables in Conservative Form (current timestep)
    Set::Field<Set::Scalar> QVec_old_mf;      // Q Vector holding all the field variables in Conservative Form (previous timestep)
    Set::Field<Set::Scalar> XFlux_mf;         // X-dir Flux Vector holding all the field variables (current timestep)
    Set::Field<Set::Scalar> YFlux_mf;         // Y-dir Flux Vector holding all the field variables (current timestep)
    Set::Field<Set::Scalar> ZFlux_mf;         // Z-dir Flux Vector holding all the field variables (current timestep)
    Set::Field<Set::Scalar> Source_mf;         // Source term Vector for all the field variables (current timestep)
    Set::Field<Set::Scalar> PVec_mf;      // Primitive Vector holding all the field variables in primitive form (current timestep)
    Set::Field<Set::Scalar> Pressure_mf;  // Multifab for pressure variable

private:
    int number_of_components;            // Number of components
    int number_of_ghost_cells = 4;           // Number of ghost cells

    Set::Scalar refinement_threshold;

    IC::IC* ic_PVec = nullptr;            // Object used to initialize Q Vector
    IC::IC* ic_Pressure = nullptr;       // Object used to initialize Pressure
    BC::BC<Set::Scalar>* bc_PVec = nullptr;            // Object used to provide Boundary Condition for Q Vector
    BC::BC<Set::Scalar>* bc_Pressure = nullptr;            // Object used to provide Boundary Condition for Q Vector
    BC::Nothing bc_nothing;

    SolverType solverType;                //solverType 

}; //ScimitarX Class

inline const std::map<std::string, ScimitarX::SolverType> ScimitarX::stringToSolverType = {
    {"SolveCompressibleEuler", ScimitarX::SolverType::SolveCompressibleEuler},
    {"SolveElastoPlastic", ScimitarX::SolverType::SolveElastoPlastic},
    {"SolveFiveEquationModel", ScimitarX::SolverType::SolveFiveEquationModel}
};

} // namespace Integrator
#endif

