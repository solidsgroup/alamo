#ifndef PY_UTIL_H
#define PY_UTIL_H
void exportUtil()
{
	bp::object utilModule(bp::handle<>(bp::borrowed(PyImport_AddModule("alamo.Util"))));
	bp::scope().attr("Util") = utilModule;
	bp::scope util_scope = utilModule;

	void (*Initialize1)() = &Util::Initialize;
	bp::def("Initialize",Initialize1);
	bp::def("Finalize",&Util::Finalize);
}
#endif