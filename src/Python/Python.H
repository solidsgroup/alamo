#include <boost/python.hpp>
using namespace boost::python;
#include "Test/Operator/Elastic.H"
#include "Util/Util.H"
BOOST_PYTHON_MODULE(alamo)
{
	void (*Initialize1)() = &Util::Initialize;
	def("Initialize",Initialize1);

	using namespace Test::Operator;
	void (Elastic::*Define1)(int,int) = &Elastic::Define;
	void (Elastic::*setBounds1)(double,double,double) = &Elastic::setBounds;
	class_<Elastic,boost::noncopyable>("Elastic")
		.def("Define",Define1)
		.def("RefluxTest",&Elastic::RefluxTest)
		.def("TrigTest",&Elastic::TrigTest)
		.def("mytest",&Elastic::mytest)
		.def("setMaxCoarseningLevel",&Elastic::setMaxCoarseningLevel)
		.def("setFixedIter",&Elastic::setFixedIter)
		.def("setMaxIter",&Elastic::setMaxIter)
		.def("setMaxFmgIter",&Elastic::setMaxFmgIter)
		.def("setBottomMaxIter",&Elastic::setBottomMaxIter)
		.def("setBounds",setBounds1)
		;

	//def("invite", invite); // Even better, invite() can also be made a member of module!!!
}
