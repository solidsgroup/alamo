#include "Util/Util.H"

#include "Model/Solid/LinearElastic/Test.H"
#include "Model/Solid/LinearElastic/Isotropic.H"
#include "Model/Solid/LinearElastic/Cubic.H"
#include "Model/Solid/LinearElastic/Degradable/Isotropic.H"

#include "Operator/Test.H"
#include "Operator/Elastic.H"

#include "Numeric/Interpolator/Test.H"
#include "Numeric/Interpolator/Linear.H"

int main (int argc, char* argv[])
{
	Util::Initialize(argc, argv);

	int failed = 0;

	{
		Model::Solid::LinearElastic::Test<Model::Solid::LinearElastic::Isotropic> test;
		failed += Util::Test::Message("Model::Solid::LinearElastic::Consistency<Isotropic>",    test.Consistency(0));
		failed += Util::Test::Message("Model::Solid::LinearElastic::MinorSymmetry1<Isotropic>", test.MinorSymmetry1(0));
		failed += Util::Test::Message("Model::Solid::LinearElastic::MinorSymmetry2<Isotropic>", test.MinorSymmetry2(0));
		failed += Util::Test::Message("Model::Solid::LinearElastic::MajorSymmetry<Isotropic>",  test.MajorSymmetry(0));
	}

	{
		Numeric::Interpolator::Test<Numeric::Interpolator::Linear<Set::Scalar> > test;
		failed += Util::Test::Message("Numeric::Interpolator::Match<Numeric::Interpolator::Linear>",test.Match(0));
	}

	{
		Operator::Test<Operator::Elastic<Model::Solid::LinearElastic::Isotropic> > test;
		failed += Util::Test::Message("Operator::RefluxTest", test.RefluxTest(0));
	}

	/// \todo Don't include in test.cc until pass/fail is computed accurately.
	// { 
	// 	Operator::Test<Operator::Elastic<Model::Solid::LinearElastic::Degradable::Isotropic> > test;
	// 	failed += Util::Test::Message("SpatiallyVaryingC Test", test.SpatiallyVaryingCTest(1));
	// }


	Util::Message(INFO,failed," tests failed");

	Util::Finalize();
	return failed;
}
