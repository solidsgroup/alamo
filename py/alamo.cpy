#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/arg_from_python.hpp>
namespace bp = boost::python;
#include "Util/Util.H"
#include "Test/Operator/Elastic.H"
//#include "Test/Operator.cpy"


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

template<class T,int n>
struct ArrToList
{
	static PyObject* convert(const std::array<T,n>& vec)
	{
		boost::python::list* l = new boost::python::list();
		for(size_t i = 0; i < vec.size(); i++) {
			l->append(vec[i]);
		}
		return l->ptr();
	}
};

struct SetVectorToList
{
	static PyObject* convert(const Set::Vector& vec)
	{
		boost::python::list* l = new boost::python::list();
		for(size_t i = 0; i < vec.size(); i++) {
			l->append(vec(i));
		}
		return l->ptr();
	}
};

template<typename containedType>
struct ListToVec{
	ListToVec(){ bp::converter::registry::push_back(&convertible,&construct,bp::type_id<std::vector<containedType> >()); }
	static void* convertible(PyObject* obj_ptr){
		if(!PySequence_Check(obj_ptr) || !PyObject_HasAttrString(obj_ptr,"__len__")) return 0;
		return obj_ptr;
	}
	static void construct(PyObject* obj_ptr, bp::converter::rvalue_from_python_stage1_data* data){
		void* storage=((bp::converter::rvalue_from_python_storage<std::vector<containedType> >*)(data))->storage.bytes;
		 new (storage) std::vector<containedType>();
		 std::vector<containedType>* v=(std::vector<containedType>*)(storage);
		 int l=PySequence_Size(obj_ptr); if(l<0) abort();
		 v->reserve(l); for(int i=0; i<l; i++) { v->push_back(bp::extract<containedType>(PySequence_GetItem(obj_ptr,i))); }
		 data->convertible=storage;
	}
};

template<typename containedType,int dim>
struct ListToArr{
	ListToArr(){ bp::converter::registry::push_back(&convertible,&construct,bp::type_id<std::array<containedType,dim> >()); }
	static void* convertible(PyObject* obj_ptr){
		if(!PySequence_Check(obj_ptr) || !PyObject_HasAttrString(obj_ptr,"__len__")) return 0;
		return obj_ptr;
	}
	static void construct(PyObject* obj_ptr, bp::converter::rvalue_from_python_stage1_data* data){
		void* storage=((bp::converter::rvalue_from_python_storage<std::array<containedType,dim> >*)(data))->storage.bytes;
		new (storage) std::array<containedType,dim>();
		std::array<containedType,dim>* v=(std::array<containedType,dim>*)(storage);
		int l=PySequence_Size(obj_ptr); if(l<0) abort();
		/*v->reserve(l);*/ for(int i=0; i<l; i++) { (*v)[i] = (bp::extract<containedType>(PySequence_GetItem(obj_ptr,i))); }
		data->convertible=storage;
	}
};

struct ListToIntVec{
	ListToIntVec(){ bp::converter::registry::push_back(&convertible,&construct,bp::type_id<amrex::IntVect >()); }
	static void* convertible(PyObject* obj_ptr){
		if(!PySequence_Check(obj_ptr) || !PyObject_HasAttrString(obj_ptr,"__len__")) return 0;
		return obj_ptr;
	}
	static void construct(PyObject* obj_ptr, bp::converter::rvalue_from_python_stage1_data* data){
		void* storage=((bp::converter::rvalue_from_python_storage<amrex::IntVect>*)(data))->storage.bytes;
		new (storage) amrex::IntVect();
		amrex::IntVect* v=(amrex::IntVect*)(storage);
		int l=PySequence_Size(obj_ptr); if(l<0) abort();
		for(int i=0; i<l; i++) {(*v)[i] = (bp::extract<int>(PySequence_GetItem(obj_ptr,i))); }
		data->convertible=storage;
	}
};

struct ListToSetVector{
	ListToSetVector(){ bp::converter::registry::push_back(&convertible,&construct,bp::type_id<Set::Vector>()); }
	static void* convertible(PyObject* obj_ptr){
		if(!PySequence_Check(obj_ptr) || !PyObject_HasAttrString(obj_ptr,"__len__")) return 0;
		return obj_ptr;
	}
	static void construct(PyObject* obj_ptr, bp::converter::rvalue_from_python_stage1_data* data){
		void* storage=((bp::converter::rvalue_from_python_storage<Set::Vector>*)(data))->storage.bytes;
		new (storage) Set::Vector();
		Set::Vector* v=(Set::Vector*)(storage);
		int l=PySequence_Size(obj_ptr); if(l<0) abort();
		for(int i=0; i<l; i++) {(*v)[i] = bp::extract<double>(PySequence_GetItem(obj_ptr,i)); }
		data->convertible=storage;
	}
};

#include "Util.cpy"
#include "Model/Solid/LinearElastic/LinearElastic.cpy"
#include "Model/Solid/LinearElastic/CrystalPlastic.cpy"
#include "Model/Interface/GB/GB.cpy"
#include "Test/Operator/Elastic.cpy"

BOOST_PYTHON_MODULE(alamo)
{
	Py_Initialize();
	np::initialize();

	// converters from vectors to python lists
	boost::python::to_python_converter<std::vector<int, std::allocator<int> >, VecToList<int> >();
	boost::python::to_python_converter<std::vector<double, std::allocator<double> >, VecToList<double> >();
	boost::python::to_python_converter<std::array<double, 2>, ArrToList<double,2> >();
	boost::python::to_python_converter<std::vector<std::string, std::allocator<std::string> >, VecToList<std::string> >();
	boost::python::to_python_converter<Set::Vector, SetVectorToList >();

	ListToIntVec();
	ListToSetVector();
	ListToArr<int,AMREX_SPACEDIM>();
	ListToArr<Set::Scalar,AMREX_SPACEDIM>();
	ListToVec<int>();
	ListToVec<Set::Scalar>();
	ListToVec<std::string>();
	//NdarrayToSetMatrix();

	pygen::convertMatrix<Set::Matrix>(false);

	bp::object package = bp::scope();
	package.attr("__path__") = "alamo";

	exportUtil();
	exportModelSolidLinearElastic();
	exportModelCrystalPlastic();
	exportModelInterfaceGB();
	exportTestOperatorElastic();
}
