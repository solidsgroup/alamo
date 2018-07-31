#ifndef OPERATOR_H_
#define OPERATOR_H_

#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>
#include <limits>

#include "BC/BC.H"
#include "Util/Util.H"

#include "Operator2.H"

using namespace amrex;

namespace Operator
{
class Operator
	: public Operator2
{
public:

	/// \brief Base constructor
	Operator () {};
	/// \brief Typical constructor
	Operator (const Vector<Geometry>& a_geom,
		  const Vector<BoxArray>& a_grids,
		  const Vector<DistributionMapping>& a_dmap,
		  BC::BC& a_bc,
		  const LPInfo& a_info);
	/// \brief Destructor (does nothing)
	virtual ~Operator ();

	Operator (const Operator&) = delete;
	Operator (Operator&&) = delete;
	Operator& operator= (const Operator&) = delete;
	Operator& operator= (Operator&&) = delete;

	/// \brief Relay to call define for parent class
	void define (const Vector<Geometry>& a_geom,
		     const Vector<BoxArray>& a_grids,
		     const Vector<DistributionMapping>& a_dmap,
		     BC::BC& a_bc,
		     const LPInfo& a_info = LPInfo());
  
protected:

};
}
#endif
