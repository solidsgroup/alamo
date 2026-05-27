#include "BC.H"
#include "AMReX_BC_TYPES.H"

namespace BC
{
namespace BCUtil
{
int ReadString(std::string bcstring)
{
    // From <AMReX_BC_TYPES.H>
    if (bcstring == "BOGUS_BC")         return amrex::BCType::mathematicalBndryTypes::bogus;
    if (bcstring == "INT_DIR")          return amrex::BCType::mathematicalBndryTypes::int_dir;
    if (bcstring == "REFLECT_ODD")      return amrex::BCType::mathematicalBndryTypes::reflect_odd;
    if (bcstring == "REFLECT_EVEN")     return amrex::BCType::mathematicalBndryTypes::reflect_even;
    if (bcstring == "FOEXTRAP")         return amrex::BCType::mathematicalBndryTypes::foextrap;
    if (bcstring == "EXT_DIR")          return amrex::BCType::mathematicalBndryTypes::ext_dir;
    if (bcstring == "HOEXTRAP")         return amrex::BCType::mathematicalBndryTypes::hoextrap;

    // From <AMReX_LO_BCTYPES.H>
    if (bcstring == "interior" )        return (int)amrex::LinOpBCType::interior;
    if (bcstring == "Dirichlet" ||
        bcstring == "dirichlet")        return (int)amrex::LinOpBCType::Dirichlet;
    if (bcstring == "Neumann" ||
        bcstring == "neumann")          return (int)amrex::LinOpBCType::Neumann;
    if (bcstring == "reflect_odd")      return (int)amrex::LinOpBCType::reflect_odd;
    if (bcstring == "Marshak")          return (int)amrex::LinOpBCType::Marshak;
    if (bcstring == "SanchezPomraning") return (int)amrex::LinOpBCType::SanchezPomraning;
    if (bcstring == "inflow")           return (int)amrex::LinOpBCType::inflow;
    if (bcstring == "Periodic" ||
        bcstring == "periodic")         return (int)amrex::LinOpBCType::Periodic;
    return 0;
}
bool IsPeriodic(int bctype)
{
    ///\todo We need to clean up these operators
    if (bctype == (int)amrex::BCType::mathematicalBndryTypes::int_dir) return true;
    if (bctype == (int)amrex::LinOpBCType::interior) return true;
    if (bctype == (int)amrex::LinOpBCType::Periodic) return true;
    else return false;
}
bool IsNeumann(int bctype)
{
    //if (bctype == INT_DIR) return true;
    //if (bctype == Interior) return true;
    if (bctype == (int)amrex::LinOpBCType::Neumann) return true;
    else return false;
}
bool IsDirichlet(int bctype)
{
    if (bctype == (int)amrex::BCType::mathematicalBndryTypes::ext_dir) return true;
    if (bctype == (int)amrex::LinOpBCType::Dirichlet) return true;
    else return false;
}
bool IsReflectEven(int bctype)
{
    if (bctype ==  (int)amrex::BCType::mathematicalBndryTypes::reflect_even) return true;
    else return false;
}
bool IsReflectOdd(int bctype)
{
    if (bctype == (int)amrex::BCType::mathematicalBndryTypes::reflect_odd) return true;
    else return false;
}
}
}
