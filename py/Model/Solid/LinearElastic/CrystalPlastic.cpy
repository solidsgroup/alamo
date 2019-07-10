#include "Model/Solid/LinearElastic/CrystalPlastic.H"

void exportModelCrystalPlastic()
{
	bp::object moduleModel(bp::handle<>(bp::borrowed(PyImport_AddModule("alamo.Model"))));
	bp::scope().attr("Model") = moduleModel;
	bp::scope scopeModel = moduleModel;
	
	bp::object moduleModelSolid(bp::handle<>(bp::borrowed(PyImport_AddModule("alamo.Model.Solid"))));
	bp::scope().attr("Solid") = moduleModelSolid;
	bp::scope scopeModelSolid = moduleModelSolid;

	bp::object moduleModelSolidCrystalPlastic(bp::handle<>(bp::borrowed(PyImport_AddModule("alamo.Model.Solid.CrystalPlastic"))));
	bp::scope().attr("CrystalPlastic") = moduleModelSolidCrystalPlastic; 
	bp::scope scopeModelSolidCrystalPlastic = moduleModelSolidCrystalPlastic;
	
	using namespace Model::Solid::CrystalPlastic;

	bp::scope the_scope = bp::class_<CrystalPlastic,boost::noncopyable>("CrystalPlastic")
		.def("UpdateSigma", &CrystalPlastic::UpdateSigma)
		.def("AdvanceEsp", &CrystalPlastic::AdvanceEsp)
		.def("Update", &CrystalPlastic::update)
		.def("Setdt", &CrystalPlastic::Setdt)
		.def("Randomize", &CrystalPlastic::Randomize)
		
		;
}
