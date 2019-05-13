#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
using namespace boost::python;
#include "Test/Operator/Elastic.H"
#include "Util/Util.H"

// converter from vector to python list
template<class T>
struct VecToList
{
	static PyObject* convert(const std::vector<T>& vec)
	{
		boost::python::list* l = new boost::python::list();
		for(size_t i = 0; i < vec.size(); i++) {
			l->append(vec[i]);
		}
		return l->ptr();
	}
};

BOOST_PYTHON_MODULE(alamo)
{

	// converters from vectors to python lists
	boost::python::to_python_converter<std::vector<double, std::allocator<double> >, VecToList<double> >();
	boost::python::to_python_converter<std::vector<std::string, std::allocator<std::string> >, VecToList<std::string> >();


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
		.def("data",&Elastic::data)
		;
}
