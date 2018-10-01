#ifndef MODEL_SOLID_ELASTIC_H_
#define MODEL_SOLID_ELASTIC_H_

#include "Model/Solid/Solid.H"

namespace Model
{
namespace Solid
{
namespace Elastic
{
class Elastic : public Solid
{
public:
	Elastic(Set::Scalar a_mu, Set::Scalar a_lambda): mu(a_mu), lambda(a_lambda) {};

	int operator () (int i, int j, int k, int l);
	
	Set::Scalar C(int i, int j, int k, int l);

	int NComp();

	

private:
	Set::Scalar mu, lambda;
};

}
}
}

#endif


