#include "BC.H"

namespace BC
{
namespace BCUtil
{
int ReadString(std::string bcstring)
{
	// From <AMReX_BC_TYPES.H>
   if (bcstring == "BOGUS_BC")         return BOGUS_BC;
   if (bcstring == "INT_DIR")          return INT_DIR;
   if (bcstring == "REFLECT_ODD")      return REFLECT_ODD;
   if (bcstring == "INT_DIR")          return INT_DIR;
   if (bcstring == "REFLECT_EVEN")     return REFLECT_EVEN;
   if (bcstring == "FOEXTRAP")         return FOEXTRAP;
   if (bcstring == "EXT_DIR")          return EXT_DIR;
   if (bcstring == "HOEXTRAP")         return HOEXTRAP;
   if (bcstring == "Interior")         return Interior;
   if (bcstring == "Inflow")           return Inflow;
   if (bcstring == "Outflow")          return Outflow;
   if (bcstring == "Symmetry")         return Symmetry;
   if (bcstring == "SlipWall")         return SlipWall;
   if (bcstring == "NoSlipWall")       return NoSlipWall;

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
	if (bctype == INT_DIR) return true;
	if (bctype == Interior) return true;
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
	if (bctype == EXT_DIR) return true;
	if (bctype == (int)amrex::LinOpBCType::Dirichlet) return true;
	else return false;
}
bool IsReflectEven(int bctype)
{
	if (bctype == REFLECT_EVEN) return true;
	else return false;
}
bool IsReflectOdd(int bctype)
{
	if (bctype == REFLECT_ODD) return true;
	else return false;
}
}
}
