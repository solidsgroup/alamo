#include "Model/Interface/GB/SH.H"
#include "Model/Interface/GB/Sin.H"
#include "Model/Interface/GB/AbsSin.H"

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

	{
		Set::Scalar          (SH::*W_1) (const Set::Scalar)                         = &SH::W;
		Set::Scalar          (SH::*W_2) (const Set::Scalar,const Set::Scalar) const = &SH::W;
		Set::Scalar          (SH::*W_3) (const Set::Vector) const                   = &SH::W;
		std::array<double,2> (SH::*DW_2)(const Set::Scalar,const Set::Scalar) const = &SH::DW;
		Set::Vector          (SH::*DW_3)(const Set::Vector) const                   = &SH::DW;
		Set::Scalar          (SH::*DW_4)(const Set::Vector,const Set::Vector) const = &SH::DW;
		Set::Scalar          (SH::*DDW)(const Set::Vector,const Set::Vector) const  = &SH::DDW;
		bp::scope the_scope = bp::class_<SH,boost::noncopyable>("SH")
			.def("Define",&SH::Define)
			.def("W",W_1)
			.def("W",W_2)
			.def("W",W_3)
			.def("DW",DW_2)
			.def("DW",DW_3)
			.def("DW",DW_4)
			.def("DDW",DDW)
			;
	}
	{
		bp::scope the_scope = bp::class_<Sin,boost::noncopyable>("Sin")
			.def("Define",&Sin::Define)
			.def("W",&Sin::W)
			.def("DW",&Sin::DW)
			.def("DDW",&Sin::DDW)
			;
	}
	{
		bp::scope the_scope = bp::class_<AbsSin,boost::noncopyable>("AbsSin")
			.def("Define",&AbsSin::Define)
			.def("W",&AbsSin::W)
			.def("DW",&AbsSin::DW)
			.def("DDW",&AbsSin::DDW)
			;
	}

}
