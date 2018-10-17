#ifndef MODEL_SOLID_ELASTIC_H_
#define MODEL_SOLID_ELASTIC_H_

#include "Model/Solid/Solid.H"
#include "Set/Set.H"

namespace Model
{
namespace Solid
{
namespace Elastic
{
///
/// The Elastic model provides a structure for computing
/// \f[\mathbb{C}_{ijkl}\varepsilon_{kl}\f]
/// that contains the data stored in \f$\mathbb{C}\f$ in an efficient way.
/// For a general material \f$\mathbb{C}\f$ contains 81 constants.
/// This can be reduced by making assumptions about the material.
///
class Elastic : public Solid
{
public:
	Elastic() {};
	virtual Set::Matrix operator () (Set::Matrix &eps) const = 0;
};

std::ostream&
operator<< (std::ostream& os,
	    const Elastic&    b);

}
}
}

#endif


