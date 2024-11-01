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

    constexpr double zeroval = 0.0;
    constexpr double one = 1.0;
    constexpr double two = 2.0;
    constexpr double three = 3.0;
    constexpr double four = 4.0;
    constexpr double five = 5.0;
    constexpr double twothird = 0.6666666666666666667;
    constexpr double half = 0.5;
    constexpr double third = two3rd * half;
    constexpr double sixth = two3rd * half * half;
    constexpr double eighth = 0.125;
    constexpr double sixteenth = 0.0625;
    constexpr double twelveth = 0.08333333333;
    constexpr double fourth = 0.25;
    constexpr double tenth = 0.1;
    constexpr double threetenth = 0.3;
    constexpr double threefourth = 0.75;
    constexpr double PI = 3.14159265358979323846264338327950288;
    constexpr double TINYVAL = 1e-10;
    constexpr double HUGEVAL = 1e+10;

    constexpr double rkCoefficients[4][3] = {
        {Constants::one, Constants::zeroval, Constants::one},
        {Constants::half, Constants::half, Constants::half},
        {Constants::threefourth, Constants::fourth, Constants::fourth},
        {Constants::third, Constants::twothird, Constants::twothird}
    };
};
}


#endif
