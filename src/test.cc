#include <stdlib.h>

#include "Util/Util.H"

#include "Model/Solid/LinearElastic/Test.H"
#include "Model/Solid/LinearElastic/Degradable/Isotropic.H"
#include "Model/Solid/LinearElastic/CrystalPlastic.H"

#include "Test/Numeric/Stencil.H"
#include "Test/Operator/Elastic.H"
#include "Test/Set/Matrix4.H"

#include "Operator/Elastic.H"

#include "Numeric/Interpolator/Test.H"
#include "Numeric/Interpolator/Linear.H"

#include "Model/Solid/Elastic/Elastic.H"
#include "Model/Solid/Elastic/NeoHookean.H"
#include "Model/Solid/Linear/Isotropic.H"
#include "Model/Solid/Linear/Laplacian.H"

int main (int argc, char* argv[])
{
	Util::Initialize(argc, argv);

	int failed = 0;

	Util::Test::Message("Model::Solid::Linear::Isotropic");
	{
		int subfailed = 0;
		subfailed += Util::Test::SubMessage("DerivativeTest1", Model::Solid::Solid<Set::Sym::Isotropic>::DerivativeTest1<Model::Solid::Linear::Isotropic>(true));
		subfailed += Util::Test::SubMessage("DerivativeTest2", Model::Solid::Solid<Set::Sym::Isotropic>::DerivativeTest2<Model::Solid::Linear::Isotropic>(true));
		failed += Util::Test::SubFinalMessage(subfailed);
	}
	Util::Test::Message("Model::Solid::Linear::Isotropic");
	{
		int subfailed = 0;
		failed += Util::Test::SubFinalMessage(subfailed);
	}

	Util::Test::Message("Set::Matrix4");
	{
		int subfailed = 0;
		Test::Set::Matrix4<2,Set::Sym::Full> test_2d_full;
		subfailed += Util::Test::SubMessage("2D - Full", test_2d_full.SymmetryTest(0));
		Test::Set::Matrix4<3,Set::Sym::Full> test_3d_full;
		subfailed += Util::Test::SubMessage("3D - Full", test_3d_full.SymmetryTest(0));
		Test::Set::Matrix4<3,Set::Sym::MajorMinor> test_3d_majorminor;
		subfailed += Util::Test::SubMessage("3D - MajorMinor", test_3d_majorminor.SymmetryTest(0));
	}

	Util::Test::Message("Model::Solid::Linear::Laplacian");
	{
		int subfailed = 0;
		Model::Solid::LinearElastic::Test<Model::Solid::Linear::Laplacian> test;
		subfailed += Util::Test::SubMessage("Consistency",    test.Consistency(2));
		subfailed += Util::Test::SubMessage("MajorSymmetry",  test.MajorSymmetry(2));
		subfailed += Util::Test::SubMessage("DerivativeTest1", Model::Solid::Solid<Set::Sym::Diagonal>::DerivativeTest1<Model::Solid::Linear::Laplacian>(true));
		subfailed += Util::Test::SubMessage("DerivativeTest2", Model::Solid::Solid<Set::Sym::Diagonal>::DerivativeTest2<Model::Solid::Linear::Laplacian>(true));
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
		// first order
		subfailed += Util::Test::SubMessage("1-0-0",test.Derivative<1,0,0>(0));
		subfailed += Util::Test::SubMessage("0-1-0",test.Derivative<0,1,0>(0));
		// second order
		subfailed += Util::Test::SubMessage("2-0-0",test.Derivative<2,0,0>(0));
		subfailed += Util::Test::SubMessage("0-2-0",test.Derivative<0,2,0>(0));
		subfailed += Util::Test::SubMessage("0-0-1",test.Derivative<0,2,0>(0));
		subfailed += Util::Test::SubMessage("1-1-0",test.Derivative<1,1,0>(0));
		// fourth order
		subfailed += Util::Test::SubMessage("3-1-0",test.Derivative<3,1,0>(0));
		subfailed += Util::Test::SubMessage("1-3-0",test.Derivative<1,3,0>(0));
		subfailed += Util::Test::SubMessage("2-2-0",test.Derivative<2,2,0>(0));
		subfailed += Util::Test::SubMessage("4-0-0",test.Derivative<4,0,0>(0));
		subfailed += Util::Test::SubMessage("0-4-0",test.Derivative<0,4,0>(0));
#if AMREX_SPACEDIM>2
		// first order
		subfailed += Util::Test::SubMessage("0-0-1",test.Derivative<0,0,1>(0));
		// second order
		subfailed += Util::Test::SubMessage("0-0-2",test.Derivative<0,0,2>(0));
		subfailed += Util::Test::SubMessage("1-0-1",test.Derivative<1,0,1>(0));
		subfailed += Util::Test::SubMessage("0-1-1",test.Derivative<0,1,1>(0));
		// fourth order
		subfailed += Util::Test::SubMessage("0-0-4",test.Derivative<0,0,4>(0));
		subfailed += Util::Test::SubMessage("0-1-3",test.Derivative<0,1,3>(0));
		subfailed += Util::Test::SubMessage("0-3-1",test.Derivative<0,3,1>(0));
		subfailed += Util::Test::SubMessage("3-0-1",test.Derivative<3,0,1>(0));
		subfailed += Util::Test::SubMessage("1-0-3",test.Derivative<1,0,3>(0));
		subfailed += Util::Test::SubMessage("0-2-2",test.Derivative<0,2,2>(0));
		subfailed += Util::Test::SubMessage("2-0-2",test.Derivative<2,0,2>(0));
		subfailed += Util::Test::SubMessage("2-1-1",test.Derivative<2,1,1>(0));
		subfailed += Util::Test::SubMessage("1-2-1",test.Derivative<1,2,1>(0));
		subfailed += Util::Test::SubMessage("1-1-2",test.Derivative<1,1,2>(0));
#endif
		failed += Util::Test::SubFinalMessage(subfailed);
	}

	Util::Test::Message("Elastic Operator Trig Test 32^n");
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

	Util::Test::Message("Elastic Operator Uniaxial Test 32^n");
	{
		int subfailed = 0;
		Test::Operator::Elastic test;
	    test.Define(32,1);
		subfailed += Util::Test::SubMessage("1 level,  Component 0",test.UniaxialTest(0,0));
		test.Define(32,2);
		subfailed += Util::Test::SubMessage("2 levels, Component 0",test.UniaxialTest(0,0));
		test.Define(32,3);
		subfailed += Util::Test::SubMessage("3 levels, Component 0",test.UniaxialTest(0,0));
		test.Define(32,2,AMREX_SPACEDIM,test.Grid::YZ);
		subfailed += Util::Test::SubMessage("2 non-centered levels, Component 0",test.UniaxialTest(0,0));
		failed += Util::Test::SubFinalMessage(subfailed);
	}
	

	Util::Message(INFO,failed," tests failed");

	Util::Finalize();
	return failed;
}
