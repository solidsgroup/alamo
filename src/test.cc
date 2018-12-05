#include "Util/Util.H"

#include "Model/Solid/LinearElastic/Test.H"
#include "Model/Solid/LinearElastic/Isotropic.H"
#include "Model/Solid/LinearElastic/Cubic.H"

#include "Operator/Test.H"

int main (int argc, char* argv[])
{
	Util::Initialize(argc, argv);

	int failed = 0;

	{
		using namespace Model::Solid::LinearElastic;
		failed += Util::Test::Message("Consistency<Isotropic>",    Test::Consistency<Isotropic>(2));
		failed += Util::Test::Message("MinorSymmetry1<Isotropic>", Test::MinorSymmetry1<Isotropic>(2));
		failed += Util::Test::Message("MinorSymmetry2<Isotropic>", Test::MinorSymmetry2<Isotropic>(2));
		failed += Util::Test::Message("MajorSymmetry<Isotropic>",  Test::MajorSymmetry<Isotropic>(2));
	}

	Util::Message(INFO,failed," tests failed");

	Util::Finalize();
	return 1;
}
