void exportTestOperatorElastic()
{
	bp::object testModule(bp::handle<>(bp::borrowed(PyImport_AddModule("alamo.Test"))));
	bp::scope().attr("Test") = testModule;
	bp::scope test_scope = testModule;

	bp::object testOperatorModule(bp::handle<>(bp::borrowed(PyImport_AddModule("alamo.Test.Operator"))));
	bp::scope().attr("Operator") = testOperatorModule;
	bp::scope test_operator_scope = testOperatorModule;

	using namespace Test::Operator;

	void (Elastic::*Define_1)(amrex::IntVect,int,int,Elastic::Grid) = &Elastic::Define;
	bp::scope the_scope = bp::class_<Elastic,boost::noncopyable>("Elastic")
		.def("Define",Define_1)
		.def("RefluxTest",&Elastic::RefluxTest)
		.def("TrigTest",&Elastic::TrigTest)
		.def("UniaxialTest",&Elastic::UniaxialTest)
		.def("setMaxCoarseningLevel",&Elastic::setMaxCoarseningLevel)
		.def("setFixedIter",&Elastic::setFixedIter)
		.def("setMaxIter",&Elastic::setMaxIter)
		.def("setMaxFmgIter",&Elastic::setMaxFmgIter)
		.def("setBottomMaxIter",&Elastic::setBottomMaxIter)
		.def("setBounds",&Elastic::setBounds)
		.def("setAgglomeration",&Elastic::setAgglomeration)
		.def("setConsolidation",&Elastic::setConsolidation)
		.def("setTolRel",&Elastic::setTolRel)
		.def("setTolAbs",&Elastic::setTolAbs)
		;

	bp::enum_<Elastic::Grid>("Grid")
		.value("X",   Elastic::Grid::X)
		.value("Y",   Elastic::Grid::Y)
		.value("Z",   Elastic::Grid::Z)
		.value("YZ",  Elastic::Grid::YZ)
		.value("ZX",  Elastic::Grid::ZX)
		.value("XY",  Elastic::Grid::XY)
		.value("XYZ", Elastic::Grid::XYZ);

}
