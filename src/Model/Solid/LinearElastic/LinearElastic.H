#ifndef MODEL_SOLID_ELASTIC_H_
#define MODEL_SOLID_ELASTIC_H_

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
class LinearElastic : public Solid
{
public:
	LinearElastic() {};
	virtual ~LinearElastic() {};
	///
	/// Returns
	/// \f[\mathbb{C}_{ijkl}u_{k,l}\f]
	///
	virtual Set::Matrix operator () (Set::Matrix &eps) const = 0;
	///
	/// Returns
	/// \f[\mathbb{C}_{ijkl}u_{k,lj}\f]
	///
	virtual Set::Vector operator () (std::array<Set::Matrix,AMREX_SPACEDIM> &gradeps) = 0;
};

std::ostream&
operator<< (std::ostream& os,
	    const LinearElastic&    b);

}
}
}

#endif


