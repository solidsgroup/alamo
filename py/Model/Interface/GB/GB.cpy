#include "Model/Interface/GB/SH.H"

void exportModelInterfaceGB()
{
	bp::object moduleModel(bp::handle<>(bp::borrowed(PyImport_AddModule("alamo.Model"))));
	bp::scope().attr("Model") = moduleModel;
	bp::scope scopeModel = moduleModel;
	bp::object moduleModelInterface(bp::handle<>(bp::borrowed(PyImport_AddModule("alamo.Model.Interface"))));
	bp::scope().attr("Interface") = moduleModelInterface;
	bp::scope scopeModelInterface = moduleModelInterface;
	bp::object moduleModelInterfaceGB(bp::handle<>(bp::borrowed(PyImport_AddModule("alamo.Model.Interface.GB"))));
	bp::scope().attr("GB") = moduleModelInterfaceGB;
	bp::scope scopeModelInterfaceGB = moduleModelInterfaceGB;

	using namespace Model::Interface::GB;

	Set::Scalar (SH::*W_1) (Set::Scalar)                   = &SH::W;
	Set::Scalar (SH::*W_2) (Set::Scalar,Set::Scalar) const = &SH::W;
	Set::Scalar (SH::*W_3) (Set::Vector) const             = &SH::W;
	Set::Vector (SH::*DW_3)(Set::Vector) const             = &SH::DW;
	bp::scope the_scope = bp::class_<SH,boost::noncopyable>("SH")
		.def("Define",&SH::Define)
		.def("W",W_1)
		.def("W",W_2)
		.def("W",W_3)
		.def("DW",DW_3)
		;
}
