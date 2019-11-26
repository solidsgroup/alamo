#ifndef MODEL_SOLID_LINEARELASTIC_TEST_H
#define MODEL_SOLID_LINEARELASTIC_TEST_H
#include "LinearElastic.H"
#include "Util/Util.H"
#include "Set/Set.H"


namespace Model
{
namespace Solid
{
namespace LinearElastic
{
template<class T>
class Test
{
public:
	int Consistency(int verbose = 0)
	{
		T model;
		for (int ctr = 0; ctr < 20; ctr++)
		{
			model.Randomize();
			Set::Matrix3 gradgradu = Set::Matrix3::Random();

			Set::Vector f1 = model(gradgradu);

			Set::Vector f2 = Set::Vector::Zero();
			for (int j = 0; j < AMREX_SPACEDIM; j++)
			{
				Set::Matrix gradu = Set::Matrix::Zero();
			
				for (int k = 0; k < AMREX_SPACEDIM; k++)
					for (int l = 0; l < AMREX_SPACEDIM; l++)
						gradu(k,l) = gradgradu[k](l,j);
				Set::Matrix sig = model(gradu);
				f2 += sig.col(j);
			}

			if ((f1-f2).norm() > 1E-8)
			{
				if (verbose>1)
				{
					Util::Message(INFO,"Test failed");
					Util::Message(INFO,"Implemented = ", f1.transpose());
					Util::Message(INFO,"Calculated  = ", f2.transpose());
				}
				return 1;
			}
		}
		return 0;	
	}

	int MinorSymmetry1(int verbose = 0)
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
				if (verbose>1)
				{
					Util::Message(INFO,"Test failed: relnorm = ", relnorm);
					Util::Message(INFO,"eps1 = \n", eps1);
					Util::Message(INFO,"sig1 = \n", sig1);
					Util::Message(INFO,"eps2 = \n", eps2);
					Util::Message(INFO,"sig2 = \n", sig2);
					Util::Message(INFO,"diff = \n", diff);
				}
				return 1;
			}
		}
		return 0;
	}

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
				if (verbose>1)
				{
					Util::Message(INFO,"Test failed: relnorm = ", relnorm);
					Util::Message(INFO,"eps  = \n", eps);
					Util::Message(INFO,"sig  = \n", sig);
					Util::Message(INFO,"diff = \n", diff);
				}
				return 1;
			}
		}
		return 0;
	}

	int MajorSymmetry(int verbose = 0)
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
				if (verbose>1)
				{
					Util::Message(INFO,"Test failed: relnorm = ", relnorm);
					Util::Message(INFO,"eps1  = \n", eps1);
					Util::Message(INFO,"eps2  = \n", eps2);
					Util::Message(INFO,"energy1 = ", energy1);
					Util::Message(INFO,"energy2 = ", energy2);
				}

				return 1;
			}
		}
		return 0;
	}


};
}
}
}



#endif


