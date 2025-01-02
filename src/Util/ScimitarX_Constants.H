#ifndef UTIL_SCIMITARX_CONSTANTS_H
#define UTIL_SCIMITARX_CONSTANTS_H

#include <string>
#include <iostream>
#include <vector>
#include <array>
#include "Util/Util.H"

namespace Util
{
class ScimitarX_Constants
{
public:
    ScimitarX_Constants () {};

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

    static constexpr double rkCoefficients[4][3] = {
        {one, zeroval, one},
        {half, half, half},
        {threefourth, fourth, fourth},
        {third, twothird, twothird}
    };
};
}


#endif
