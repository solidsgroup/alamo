#ifndef MODEL_SOLID_H_
#define MODEL_SOLID_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>

#include "Set/Set.H"

namespace Model
{
namespace Solid
{
class Solid
{
public:
	Solid() {} ;

	// virtual Set::Matrix operator () (Set::Matrix &eps);// = 0;
	// virtual Set::Vector operator () (std::array<Set::Matrix,AMREX_SPACEDIM> &gradeps);// = 0;
	// virtual Solid & operator + (const Solid & rhs);// = 0;
	// virtual Solid & operator - (const Solid & rhs);// = 0;
	// virtual Solid & operator / (const Set::Scalar alpha);// = 0;
};

}
}


#endif


