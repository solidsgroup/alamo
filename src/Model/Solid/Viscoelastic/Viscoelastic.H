#ifndef MODEL_SOLID_VISCOELASTIC_H_
#define MODEL_SOLID_VISCOELASTIC_H_

#include "Model/Solid/Solid.H"
#include "Set/Set.H"

namespace Model
{
namespace Solid
{
namespace Viscoelastic
{
///
/// The Elastic model provides a structure for computing
/// \f[\mathbb{C}_{ijkl}\varepsilon_{kl}\f]
/// that contains the data stored in \f$\mathbb{C}\f$ in an efficient way.
/// For a general material \f$\mathbb{C}\f$ contains 81 constants.
/// This can be reduced by making assumptions about the material.
///
class Viscoelastic : public Solid
{
public:
	Viscoelastic() {};
	virtual ~Viscoelastic() {};
};

std::ostream&
operator<< (std::ostream& os,
	    const Viscoelastic&    b);

}
}
}

#endif


