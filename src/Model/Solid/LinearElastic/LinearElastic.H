#ifndef MODEL_SOLID_LINEARELASTIC_H_
#define MODEL_SOLID_LINEARELASTIC_H_

#include "Model/Solid/Solid.H"
#include "Set/Set.H"

namespace Model
{
namespace Solid
{
namespace LinearElastic
{
///
/// The LinearElastic model provides a structure for computing
/// \f[\mathbb{C}_{ijkl}\varepsilon_{kl}\f]
/// that contains the data stored in \f$\mathbb{C}\f$ in an efficient way.
/// For a general material \f$\mathbb{C}\f$ contains 81 constants.
/// This can be reduced by making assumptions about the material.
///
template<Set::Sym sym>
class LinearElastic
{
public:
	LinearElastic() {};
	virtual ~LinearElastic() {};
	///
	/// Returns
	/// \f[\mathbb{C}_{ijkl}u_{k,l}\f]
	///
	virtual Set::Matrix operator () (Set::Matrix &eps, bool a_homogeneous=true) const = 0;
	///
	/// Returns
	/// \f[\mathbb{C}_{ijkl}u_{k,lj}\f]
	///
	virtual Set::Vector operator () (Set::Matrix3 &gradeps, bool a_homogeneous=true) = 0;
	
protected:
};

//template<Set::Sym sym>
//std::ostream& operator<< (std::ostream& os, const LinearElastic<sym>& b);

}
}
}

#endif


