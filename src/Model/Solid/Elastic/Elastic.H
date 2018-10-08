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
class Elastic : public Solid
{
public:
	Elastic() {};

	void Print()
	{

#if AMREX_SPACEDIM == 2
		std::array<Set::Matrix,3> gradu, eps;
		gradu[0] << 1,0, 0,0; eps[0] = 0.5*(gradu[0] + gradu[0].transpose());
		gradu[1] << 0,0, 0,1; eps[1] = 0.5*(gradu[1] + gradu[1].transpose());
		gradu[2] << 0,1, 0,0; eps[2] = 0.5*(gradu[2] + gradu[2].transpose());
		for (int i = 0; i < 6; i++)
		{
			for (int j = 0; j < 6; j++)
			{
				Set::Scalar comp = (eps[i].transpose() * (*this)(eps[j])).trace();
				std::cout << (fabs(comp)>1E-10 ? comp : 0)  << "\t";
			}
			std::cout << std::endl;
		}
#endif
#if AMREX_SPACEDIM == 3
		std::array<Set::Matrix,6> gradu, eps;
		gradu[0] << 1,0,0, 0,0,0, 0,0,0; eps[0] = 0.5*(gradu[0] + gradu[0].transpose());
		gradu[1] << 0,0,0, 0,1,0, 0,0,0; eps[1] = 0.5*(gradu[1] + gradu[1].transpose());
		gradu[2] << 0,0,0, 0,0,0, 0,0,1; eps[2] = 0.5*(gradu[2] + gradu[2].transpose());
		gradu[3] << 0,0,0, 0,0,1, 0,0,0; eps[3] = 0.5*(gradu[3] + gradu[3].transpose());
		gradu[4] << 0,0,1, 0,0,0, 0,0,0; eps[4] = 0.5*(gradu[4] + gradu[4].transpose());
		gradu[5] << 0,1,0, 0,0,0, 0,0,0; eps[5] = 0.5*(gradu[5] + gradu[5].transpose());
		for (int i = 0; i < 6; i++)
		{
			for (int j = 0; j < 6; j++)
			{
				Set::Scalar comp = (eps[i].transpose() * (*this)(eps[j])).trace();
				std::cout << (fabs(comp)>1E-10 ? comp : 0)  << "\t";
			}
			std::cout << std::endl;
		}
#endif

	}
	virtual Set::Matrix operator () (Set::Matrix &eps) = 0;
	

};

}
}
}

#endif


