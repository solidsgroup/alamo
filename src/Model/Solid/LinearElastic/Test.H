#ifndef MODEL_SOLID_LINEARELASTIC_H
#define MODEL_SOLID_LINEARELASTIC_H
#include "LinearElastic.H"
#include "Util/Util.H"
#include "Set/Set.H"


namespace Model
{
namespace Solid
{
namespace LinearElastic
{
namespace Test
{

template <class T>
bool Consistency(int verbose = 0)
{
	
}

template <class T>
bool MinorSymmetry1(int verbose = 0)
{
	T model;
	for (int i = 0; i < 20; i++)
	{
		model.Randomize();
		Set::Matrix eps1 = Set::Matrix::Random();
		Set::Matrix eps2 = eps1.transpose();
		Set::Matrix sig1 = model(eps1);
		Set::Matrix sig2 = model(eps2);
		Set::Matrix diff = sig1 - sig2;
		Set::Scalar relnorm = diff.norm() / (0.5 * sig1.norm() + 0.5 * sig2.norm());
		if (relnorm > 1E-8)
		{
			if (verbose>0)
			{
				Util::Message(INFO,"Test failed: relnorm = ", relnorm);
				Util::Message(INFO,"eps1 = \n", eps1);
				Util::Message(INFO,"sig1 = \n", sig1);
				Util::Message(INFO,"eps2 = \n", eps2);
				Util::Message(INFO,"sig2 = \n", sig2);
				Util::Message(INFO,"diff = \n", diff);
			}
			return false;
		}
	}
	return true;
}

template <class T>
bool MinorSymmetry2(int verbose = 0)
{
	T model;
	for (int i = 0; i < 20; i++)
	{
		model.Randomize();
		Set::Matrix eps = Set::Matrix::Random();
		Set::Matrix sig = model(eps);
		Set::Matrix diff = sig - sig.transpose();
		Set::Scalar relnorm = diff.norm() / sig.norm();
		if (relnorm > 1E-8)
		{
			if (verbose>0)
			{
				Util::Message(INFO,"Test failed: relnorm = ", relnorm);
				Util::Message(INFO,"eps  = \n", eps);
				Util::Message(INFO,"sig  = \n", sig);
				Util::Message(INFO,"diff = \n", diff);
			}
			return false;
		}
	}
	return true;
}

template <class T>
bool MajorSymmetry(int verbose = 0)
{
	T model;
	for (int i = 0; i < 20; i++)
	{
		model.Randomize();
		Set::Matrix eps1 = Set::Matrix::Random();
		Set::Matrix eps2 = Set::Matrix::Random();

		Set::Scalar energy1 = (eps1.transpose() * model(eps2)).trace();
		Set::Scalar energy2 = (eps2.transpose() * model(eps1)).trace();

		Set::Scalar relnorm = fabs((energy1 - energy2) / (energy1 + energy2));

		if (relnorm > 1E-8)
		{
			return false;
		}
	}
	return true;
}


}
}
}
}



#endif


