#include <string>
#include <iostream>
#include <vector>
#include <array>
#include "Util/Util.H"
#include "Set/Set.H"
#include "ScimitarX_Util.H"

namespace Util
{
    void ScimitarX_Util::Debug::DebugComputeInternalEnergyFromDensityAndPressure(int i, int j, int k, int lev,
                                                                            Set::Scalar pressure, Set::Scalar density, Set::Scalar gamma,
                                                                            bool abort_on_checkpoint, const std::string& mode,
                                                                            const std::string& context) {
    std::ostringstream message;
    message << "Debugging ComputeInternalEnergyFromDensityAndPressure (" << context << ") at cell: (" 
            << i << ", " << j << ", " << k << ") at level " << lev << "\n";
    message << "Context: " << context << "\n";
    message << "Density: " << density << "\n";
    message << "Pressure: " << pressure << "\n";
    message << "Gamma: " << gamma << "\n";

    if (mode == "print") {
        Util::Message(INFO, message.str());
    } else if (mode == "check") {
        if (density <= 0.0 || pressure < 0.0 || std::isnan(pressure) || std::isnan(density) || std::isnan(gamma)) {
            message << "Error: Invalid values detected in context: " << context << "\n";
            if (abort_on_checkpoint) {
                Util::Abort(INFO, message.str(), "Checkpoint reached. Aborting.");
            } else {
                Util::Warning(INFO, message.str());
            }
        } else {
            Util::Message(INFO, "Values are valid: " + message.str());
        }
    } else {
        Util::Warning(INFO, "Unknown debug mode specified. Defaulting to 'print'.");
        Util::Message(INFO, message.str());
    }

    if (abort_on_checkpoint) {
        Util::Abort(INFO, message.str(), "Checkpoint reached. Aborting as requested.");
    }
}

    amrex::IntVect Util::ScimitarX_Util::Debug::target_location = amrex::IntVect::TheZeroVector();  // Definition

} // namespace 
