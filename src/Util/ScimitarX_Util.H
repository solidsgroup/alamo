#ifndef UTIL_SCIMITARX_UTIL_H
#define UTIL_SCIMITARX_UTIL_H

#include <string>
#include <iostream>
#include <vector>
#include <array>
#include "Util/Util.H"

namespace Util
{
class ScimitarX_Util
{
public:
    ScimitarX_Util () {};

    struct Constants {
        static constexpr double zeroval = 0.0;
        static constexpr double one = 1.0;
        static constexpr double two = 2.0;
        static constexpr double three = 3.0;
        static constexpr double four = 4.0;
        static constexpr double five = 5.0;
        static constexpr double twothird = 0.6666666666666666667;
        static constexpr double half = 0.5;
        static constexpr double third = twothird * half;
        static constexpr double sixth = twothird * half * half;
        static constexpr double eighth = 0.125;
        static constexpr double sixteenth = 0.0625;
        static constexpr double twelveth = 0.08333333333;
        static constexpr double fourth = 0.25;
        static constexpr double tenth = 0.1;
        static constexpr double threetenth = 0.3;
        static constexpr double threefourth = 0.75;
        static constexpr double PI = 3.14159265358979323846264338327950288;
        static constexpr double TINYVAL = 1e-10;
        static constexpr double HUGEVAL = 1e+10;
    };

    // NOTE: When modifying solver types, ensure the corresponding variables are added here to enable dynamic indexing.
    struct getVariableIndex {
        std::map<int, int> variableIndexMap; // Mapping from variable enums to indices
        int NVAR_MAX = 0;                    // Total number of variables
        int DENS = -1;                       // Index for Density
        int UVEL = -1;                       // Index for U-velocity
        int VVEL = -1;                       // Index for V-velocity
        int WVEL = -1;                       // Index for W-velocity
        int IE = -1;                         // Index for Internal Energy
    };

};
}

#endif // UTIL_SCIMITARX_UTIL_H

