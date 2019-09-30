#include "FiniteElement/Q4.H"

void exportElementQ4()
{
	bp::object moduleFiniteElement(bp::handle<>(bp::borrowed(PyImport_AddModule("alamo.FiniteElement"))));
	bp::scope().attr("FiniteElement") = moduleFiniteElement;
	bp::scope scopeFiniteElement = moduleFiniteElement;

	{
		using namespace FiniteElement;
		bp::scope the_scope = bp::class_<Q4,boost::noncopyable>("Q4")
		.def("Define",&Q4::Define)
		.def("Phi",   &Q4::Phi)
		.def("DPhi",  &Q4::DPhi)
		.def("Qpt",   &Q4::Qpt)
		.def("Qwt",   &Q4::Qwt)
		;
	}
};
