#ifndef UTIL_SCIMITARX_SETINDX_H
#define UTIL_SCIMITARX_SETINDX_H

#include <map>
#include <vector>
#include <set>
#include <stdexcept>
#include <iostream>

#include "Util/Util.H"

namespace Util
{

    // Enum for Solver Types
    enum class SolverType {
        SolveCompressibleEuler,
        SolveElastoPlastic,
        SolveFiveEquationModel
    };

    // Maps for string to Enum conversions
    std::map<std::string, SolverType> stringToSolverType{
        {"SolveCompressibleEuler", SolverType::SolveCompressibleEuler},
        {"SolveElastoPlastic", SolverType::SolveElastoPlastic},
        {"SolveFiveEquationModel", SolverType::SolveFiveEquationModel}
    };

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
    constexpr int toInt(Enum e) noexcept {
        return static_cast<int>(e);
    }

    // Singleton Class for SolverType Variable Management
    class ScimitarX_SetIndex {
    public:
        // Get the Singleton Instance
        static ScimitarX_SetIndex& getInstance() {
            static ScimitarX_SetIndex instance; // Thread-safe in C++11 and later
            return instance;
        }

        // Get Variables Mapped to Solver Type
        const std::vector<Util::EulerFieldVariables>& getVariables(Util::SolverType solverType) {
            auto it = solverTypeVariables.find(solverType);
            if (it == solverTypeVariables.end()) {
                throw std::runtime_error("Invalid SolverType: No variables mapped.");
            }
            return it->second;
        }

        // Compute and Assign Indices to Variables
        struct VariableIndicesResult {
            std::map<int, int> variableIndexMap; // Mapping from variable enums to indices
            int NVAR_MAX = 0;                    // Total number of variables
            int DENS = -1;                       // Index for Density
            int UVEL = -1;                       // Index for U-velocity
            int VVEL = -1;                       // Index for V-velocity
            int WVEL = -1;                       // Index for W-velocity
            int IE = -1;                         // Index for Internal Energy
        };

        VariableIndicesResult computeAndAssignVariableIndices(Util::SolverType solverType) {
            VariableIndicesResult result;
            std::set<int> allVariables;

            // Retrieve variables for the given solver type
            const auto& variables = getVariables(solverType);
            for (const auto& var : variables) {
                allVariables.insert(toInt(var));
            }

            // Assign indices dynamically
            int index = 0;
            for (int var : allVariables) {
                result.variableIndexMap[var] = index++;
            }

            result.NVAR_MAX = index;
            result.DENS = result.variableIndexMap[toInt(Util::EulerFieldVariables::DENS)];
            result.UVEL = result.variableIndexMap[toInt(Util::EulerFieldVariables::UVEL)];
#if AMREX_SPACEDIM >= 2
            result.VVEL = result.variableIndexMap[toInt(Util::EulerFieldVariables::VVEL)];
#endif
#if AMREX_SPACEDIM == 3
            result.WVEL = result.variableIndexMap[toInt(Util::EulerFieldVariables::WVEL)];
#endif
            result.IE = result.variableIndexMap[toInt(Util::EulerFieldVariables::IE)];

            return result;
        }

    private:
        ScimitarX_SetIndex() {
            // Initialize the map of solver type to variables
            solverTypeVariables[Util::SolverType::SolveCompressibleEuler] = {
                Util::EulerFieldVariables::DENS,
                Util::EulerFieldVariables::UVEL,
#if AMREX_SPACEDIM >= 2
                Util::EulerFieldVariables::VVEL,
#endif
#if AMREX_SPACEDIM == 3
                Util::EulerFieldVariables::WVEL,
#endif
                Util::EulerFieldVariables::IE
            };
        }

        ScimitarX_SetIndex(const ScimitarX_SetIndex&) = delete; // Prevent copy
        ScimitarX_SetIndex& operator=(const ScimitarX_SetIndex&) = delete; // Prevent assignment

        std::map<Util::SolverType, std::vector<Util::EulerFieldVariables>> solverTypeVariables;
    };

} // namespace Util

#endif // UTIL_SCIMITARX_SETINDX_H

