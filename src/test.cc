#include "Util/Util.H"

#include "Model/Solid/LinearElastic/Test.H"
#include "Model/Solid/LinearElastic/Isotropic.H"
#include "Model/Solid/LinearElastic/Cubic.H"

#include "Operator/Test.H"
#include "Operator/Elastic.H"

int main (int argc, char* argv[])
{
	Util::Initialize(argc, argv);

	int failed = 0;

	{
		Model::Solid::LinearElastic::Test<Model::Solid::LinearElastic::Isotropic> test;
		failed += Util::Test::Message("Consistency<Isotropic>",    test.Consistency(2));
		failed += Util::Test::Message("MinorSymmetry1<Isotropic>", test.MinorSymmetry1(2));
		failed += Util::Test::Message("MinorSymmetry2<Isotropic>", test.MinorSymmetry2(2));
		failed += Util::Test::Message("MajorSymmetry<Isotropic>",  test.MajorSymmetry(2));
	}

	{
		Operator::Test<Operator::Elastic<Model::Solid::LinearElastic::Isotropic> > test;
		failed += Util::Test::Message("RefluxTest", test.RefluxTest(2));
	}


	Util::Message(INFO,failed," tests failed");

	Util::Finalize();
	return 1;
}
