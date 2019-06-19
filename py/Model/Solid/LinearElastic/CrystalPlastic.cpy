#include "Model/Solid/LinearElastic/CrystalPlastic.H"

void exportModelSolidLinearElastic()
{
	bp::object moduleModel(bp::handle<>(bp::borrowed(PyImport_AddModule("alamo.Model"))));
	bp::scope().attr("Model") = moduleModel;
	bp::scope scopeModel = moduleModel;
	
	bp::object moduleModelSolid(bp::handle<>(bp::borrowed(PyImport_AddModule("alamo.Model.Solid"))));
	bp::scope().attr("Solid") = moduleModelSolid;
	bp::scope scopeModelSolid = moduleModelSolid;

	bp::object moduleModelSolidLinearElastic(bp::handle<>(bp::borrowed(PyImport_AddModule("alamo.Model.Solid.LinearElastic"))));
	bp::scope().attr("LinearElastic") = moduleModelSolidLinearElastic;
	bp::scope scopeModelSolidLinearElastic = moduleModelSolidLinearElastic;
	
	using namespace Model::Solid::CrystalPlastic;

	//Set::Matrix          (Isotropic::*DW_1) (const Set::Matrix)                         = &SH::W;
	//Set::Scalar          (SH::*W_2) (const Set::Scalar,const Set::Scalar) const = &SH::W;
	//Set::Scalar          (SH::*W_3) (const Set::Vector) const                   = &SH::W;
	//std::array<double,2> (SH::*DW_2)(const Set::Scalar,const Set::Scalar) const = &SH::DW;
	//Set::Vector          (SH::*DW_3)(const Set::Vector) const                   = &SH::DW;
	//Set::Scalar          (SH::*DW_4)(const Set::Vector,const Set::Vector) const = &SH::DW;
	//Set::Scalar          (SH::*DDW)(const Set::Vector,const Set::Vector) const  = &SH::DDW;
	bp::scope the_scope = bp::class_<Isotropic,boost::noncopyable>("Isotropic")
		.def("Define",&Isotropic::Define)
		.def("Randomize",&Isotropic::Randomize)
		.def("W",&Isotropic::W)
		.def("DW",&Isotropic::DW)
		//.def("W",W_1)
		//.def("W",W_2)
		//.def("W",W_3)
		//.def("DW",DW_2)
		//.def("DW",DW_3)
		//.def("DW",DW_4)
		//.def("DDW",DDW)
		;
}
