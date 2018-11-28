#include "Util/Util.H"
#include "Model/Solid/LinearElastic/Test.H"
#include "Model/Solid/LinearElastic/Isotropic.H"

int main (int argc, char* argv[])
{
	Util::Initialize(argc, argv);

	int failed = 0;

	{
		using namespace Model::Solid::LinearElastic;

		failed += Util::Test::Message("MinorSymmetry1<Isotropic>", Test::MinorSymmetry1<Isotropic>(2));
		failed += Util::Test::Message("MinorSymmetry2<Isotropic>", Test::MinorSymmetry2<Isotropic>(2));
		failed += Util::Test::Message("MajorSymmetry<Isotropic>",  Test::MajorSymmetry<Isotropic>(2));


		Isotropic model;
		model.Randomize();

		std::array<Set::Matrix,AMREX_SPACEDIM> gradgradu;
		for (int i = 0; i < AMREX_SPACEDIM; i++) gradgradu[i] = Set::Matrix::Random();

		Set::Vector f = model(gradgradu);

		// f_i = C_{ijkl} u_{k,lj} = C_{ijkl} Delta u_{k,l}

		// sig_{ij} = C_{ijkl} u_{k,l}
		

	}

	Util::Message(INFO,failed," tests failed");

	Util::Finalize();
	return 1;
}
