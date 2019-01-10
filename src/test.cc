#include "Util/Util.H"

#include "Model/Solid/LinearElastic/Test.H"
#include "Model/Solid/LinearElastic/Isotropic.H"
#include "Model/Solid/LinearElastic/Cubic.H"
#include "Model/Solid/LinearElastic/Laplacian.H"
#include "Model/Solid/LinearElastic/Degradable/Isotropic.H"

#include "Test/Operator/Elastic/Analytic.H"
#include "Test/Operator/Elastic/Reflux.H"

#include "Operator/Elastic.H"

#include "Numeric/Interpolator/Test.H"
#include "Numeric/Interpolator/Linear.H"

int main (int argc, char* argv[])
{
	Util::Initialize(argc, argv);


	int failed = 0;

	{
		Model::Solid::LinearElastic::Test<Model::Solid::LinearElastic::Isotropic> test;
		
		Util::Test::Message(          "Model::Solid::LinearElastic<Isotropic>");
		failed += Util::Test::Message("  ├ Consistency",    test.Consistency(0));
		failed += Util::Test::Message("  ├ MinorSymmetry1", test.MinorSymmetry1(0));
		failed += Util::Test::Message("  ├ MinorSymmetry2", test.MinorSymmetry2(0));
		failed += Util::Test::Message("  └ MajorSymmetry",  test.MajorSymmetry(0));
	}

	{
		Model::Solid::LinearElastic::Test<Model::Solid::LinearElastic::Cubic> test;
		Util::Test::Message(          "Model::Solid::LinearElastic<Cubic>");
		failed += Util::Test::Message("  ├ Consistency",    test.Consistency(2));
		failed += Util::Test::Message("  ├ MinorSymmetry1", test.MinorSymmetry1(2));
		failed += Util::Test::Message("  ├ MinorSymmetry2", test.MinorSymmetry2(2));
		failed += Util::Test::Message("  └ MajorSymmetry",  test.MajorSymmetry(2));
	}

	{
		Model::Solid::LinearElastic::Test<Model::Solid::LinearElastic::Laplacian> test;
		Util::Test::Message(          "Model::Solid::LinearElastic<Laplacian>");
		failed += Util::Test::Message("  └ Consistency",    test.Consistency(2));
	}

	{
		Numeric::Interpolator::Test<Numeric::Interpolator::Linear<Set::Scalar> > test;
		Util::Test::Message(          "Numeric::Interpolator<Linear>");
		failed += Util::Test::Message("  └ Match",test.Match(0));
	}

	{
		Test::Operator::Elastic::Reflux<Operator::Elastic<Model::Solid::LinearElastic::Isotropic> > test;
	 	failed += Util::Test::Message("  └ Operator::RefluxTest", test.RefluxTest(0));
	}

	{
		Test::Operator::Elastic::Analytic test;
		test.Define(16,1);
		Util::Test::Message(          "Elastic Operator Trig Test 32x32, 1 level");
		failed += Util::Test::Message("  └ Component 0, period=1",test.TrigTest(0,0,1));

		test.Define(16,2);
		Util::Test::Message(          "Elastic Operator Trig Test 32x32, 2 levels");
		failed += Util::Test::Message("  └ Component 0, period=1",test.TrigTest(0,0,1));
	}


	Util::Message(INFO,failed," tests failed");

	Util::Finalize();
	return failed;
}
