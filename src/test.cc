#include <stdlib.h>

#include "Util/Util.H"

#include "Model/Solid/LinearElastic/Test.H"
#include "Model/Solid/LinearElastic/Isotropic.H"
#include "Model/Solid/LinearElastic/Cubic.H"
#include "Model/Solid/LinearElastic/Laplacian.H"
#include "Model/Solid/LinearElastic/Degradable/Isotropic.H"
#include "Model/Solid/LinearElastic/CrystalPlastic.H"

#include "Test/Numeric/Stencil.H"
#include "Test/Operator/Elastic.H"

#include "Operator/Elastic.H"

#include "Numeric/Interpolator/Test.H"
#include "Numeric/Interpolator/Linear.H"

int main (int argc, char* argv[])
{
	Util::Initialize(argc, argv);

	int failed = 0;

	Util::Test::Message("Model::Solid::LinearElastic<Cubic>");
	{
		int subfailed = 0;
		Model::Solid::LinearElastic::Test<Model::Solid::LinearElastic::Cubic> test;
		subfailed += Util::Test::SubMessage("Consistency",    test.Consistency(2));
		subfailed += Util::Test::SubMessage("MinorSymmetry1", test.MinorSymmetry1(2));
		subfailed += Util::Test::SubMessage("MinorSymmetry2", test.MinorSymmetry2(2));
		subfailed += Util::Test::SubMessage("MajorSymmetry",  test.MajorSymmetry(2));
		failed += Util::Test::SubFinalMessage(subfailed);
	}

	Util::Test::Message("Model::Solid::LinearElastic<Degradable::Isotropic>");
	{
		int subfailed = 0;
		Model::Solid::LinearElastic::Test<Model::Solid::LinearElastic::Degradable::Isotropic> test;
		subfailed += Util::Test::SubMessage("Consistency",    test.Consistency(2));
		subfailed += Util::Test::SubMessage("MinorSymmetry1", test.MinorSymmetry1(2));
		subfailed += Util::Test::SubMessage("MinorSymmetry2", test.MinorSymmetry2(2));
		subfailed += Util::Test::SubMessage("MajorSymmetry",  test.MajorSymmetry(2));
		failed += Util::Test::SubFinalMessage(subfailed);
	}

	Util::Test::Message("Model::Solid::LinearElastic<Laplacian>");
	{
		int subfailed = 0;
		Model::Solid::LinearElastic::Test<Model::Solid::LinearElastic::Laplacian> test;
		subfailed += Util::Test::SubMessage("Consistency",    test.Consistency(2));
		failed += Util::Test::SubFinalMessage(subfailed);
	}

	Util::Test::Message("Numeric::Interpolator<Linear>");
	{
		int subfailed = 0;
		Numeric::Interpolator::Test<Numeric::Interpolator::Linear<Set::Scalar> > test;
		subfailed += Util::Test::SubMessage("Match",test.Match(0));
		failed += Util::Test::SubFinalMessage(subfailed);
	}

	Util::Test::Message("Numeric::Stencil test");
	{
		int subfailed = 0;
		Test::Numeric::Stencil test;
		test.Define(32);
		subfailed += Util::Test::SubMessage("1-0-0",test.Derivative<1,0,0>(0));
		subfailed += Util::Test::SubMessage("0-1-0",test.Derivative<0,1,0>(0));
		subfailed += Util::Test::SubMessage("2-0-0",test.Derivative<2,0,0>(0));
		subfailed += Util::Test::SubMessage("0-2-0",test.Derivative<0,2,0>(0));
		subfailed += Util::Test::SubMessage("1-1-0",test.Derivative<1,1,0>(0));
		subfailed += Util::Test::SubMessage("4-0-0",test.Derivative<4,0,0>(0));
		subfailed += Util::Test::SubMessage("0-4-0",test.Derivative<0,4,0>(0));
		subfailed += Util::Test::SubMessage("3-1-0",test.Derivative<3,1,0>(0));
		subfailed += Util::Test::SubMessage("1-3-0",test.Derivative<1,3,0>(0));
		subfailed += Util::Test::SubMessage("2-2-0",test.Derivative<2,2,0>(0));
#if AMREX_SPACEDIM>2
		subfailed += Util::Test::SubMessage("0-0-1",test.Derivative<0,0,1>(0));
		subfailed += Util::Test::SubMessage("0-0-2",test.Derivative<0,0,2>(0));
		subfailed += Util::Test::SubMessage("1-0-1",test.Derivative<1,0,1>(0));
		subfailed += Util::Test::SubMessage("0-1-1",test.Derivative<0,1,1>(0));
#endif
		failed += Util::Test::SubFinalMessage(subfailed);
	}

	Util::Test::Message("Elastic Operator Trig Test 32x32");
	{
		int subfailed = 0;
		Test::Operator::Elastic test;
		test.Define(32,1);
		subfailed += Util::Test::SubMessage("1 level,  Component 0, period=1",test.TrigTest(0,0,1));
		test.Define(32,2);
		subfailed += Util::Test::SubMessage("2 levels, Reflux test",          test.RefluxTest(0));
		subfailed += Util::Test::SubMessage("2 levels, Component 0, period=1",test.TrigTest(0,0,1));
		test.Define(32,3);
		subfailed += Util::Test::SubMessage("3 levels, Reflux test",          test.RefluxTest(0));
		subfailed += Util::Test::SubMessage("3 levels, Component 0, period=1",test.TrigTest(0,0,1));
		failed += Util::Test::SubFinalMessage(subfailed);

	}

	Util::Test::Message("Elastic Operator Uniaxial Test 32x32");
	{
		int subfailed = 0;
		Test::Operator::Elastic test;
	        test.Define(32,1);
		subfailed += Util::Test::SubMessage("1 level, Component 0, period=1",test.UniaxialTest(0,0));
		test.Define(32,2);
		subfailed += Util::Test::SubMessage("2 levels, Component 0, period=1",test.UniaxialTest(0,0));
		test.Define(32,3);
		subfailed += Util::Test::SubMessage("3 levels, Component 0, period=1",test.UniaxialTest(0,0));
		failed += Util::Test::SubFinalMessage(subfailed);
	}
	

	Util::Message(INFO,failed," tests failed");

	Util::Finalize();
	return failed;
}
