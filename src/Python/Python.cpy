#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/arg_from_python.hpp>

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

template<typename containedType>
struct ListToVec{
	ListToVec(){ converter::registry::push_back(&convertible,&construct,type_id<std::vector<containedType> >()); }
	static void* convertible(PyObject* obj_ptr){
		if(!PySequence_Check(obj_ptr) || !PyObject_HasAttrString(obj_ptr,"__len__")) return 0;
		return obj_ptr;
	}
	static void construct(PyObject* obj_ptr, converter::rvalue_from_python_stage1_data* data){
		 void* storage=((converter::rvalue_from_python_storage<std::vector<containedType> >*)(data))->storage.bytes;
		 new (storage) std::vector<containedType>();
		 std::vector<containedType>* v=(std::vector<containedType>*)(storage);
		 int l=PySequence_Size(obj_ptr); if(l<0) abort();
		 v->reserve(l); for(int i=0; i<l; i++) { v->push_back(extract<containedType>(PySequence_GetItem(obj_ptr,i))); }
		 data->convertible=storage;
	}
};

template<typename containedType,int dim>
struct ListToArr{
	ListToArr(){ converter::registry::push_back(&convertible,&construct,type_id<std::array<containedType,dim> >()); }
	static void* convertible(PyObject* obj_ptr){
		if(!PySequence_Check(obj_ptr) || !PyObject_HasAttrString(obj_ptr,"__len__")) return 0;
		return obj_ptr;
	}
	static void construct(PyObject* obj_ptr, converter::rvalue_from_python_stage1_data* data){
		void* storage=((converter::rvalue_from_python_storage<std::array<containedType,dim> >*)(data))->storage.bytes;
		new (storage) std::array<containedType,dim>();
		std::array<containedType,dim>* v=(std::array<containedType,dim>*)(storage);
		int l=PySequence_Size(obj_ptr); if(l<0) abort();
		/*v->reserve(l);*/ for(int i=0; i<l; i++) { (*v)[i] = (extract<containedType>(PySequence_GetItem(obj_ptr,i))); }
		data->convertible=storage;
	}
};

struct ListToIntVec{
	ListToIntVec(){ converter::registry::push_back(&convertible,&construct,type_id<amrex::IntVect >()); }
	static void* convertible(PyObject* obj_ptr){
		if(!PySequence_Check(obj_ptr) || !PyObject_HasAttrString(obj_ptr,"__len__")) return 0;
		return obj_ptr;
	}
	static void construct(PyObject* obj_ptr, converter::rvalue_from_python_stage1_data* data){
		void* storage=((converter::rvalue_from_python_storage<amrex::IntVect>*)(data))->storage.bytes;
		new (storage) amrex::IntVect();
		amrex::IntVect* v=(amrex::IntVect*)(storage);
		int l=PySequence_Size(obj_ptr); if(l<0) abort();
		for(int i=0; i<l; i++) {(*v)[i] = (extract<int>(PySequence_GetItem(obj_ptr,i))); }
		data->convertible=storage;
	}
};


BOOST_PYTHON_MODULE(Operator)
{

	// converters from vectors to python lists
	boost::python::to_python_converter<std::vector<int, std::allocator<int> >, VecToList<int> >();
	boost::python::to_python_converter<std::vector<double, std::allocator<double> >, VecToList<double> >();
	boost::python::to_python_converter<std::vector<std::string, std::allocator<std::string> >, VecToList<std::string> >();

	ListToIntVec();
	ListToArr<int,AMREX_SPACEDIM>();
	ListToArr<Set::Scalar,AMREX_SPACEDIM>();
	ListToVec<int>();
	ListToVec<Set::Scalar>();
	ListToVec<std::string>();

	void (*Initialize1)() = &Util::Initialize;
	def("Initialize",Initialize1);

	using namespace Test::Operator;

	void (Elastic::*Define_1)(int,int) = &Elastic::Define;
	void (Elastic::*Define_2)(amrex::IntVect,int,int,Elastic::Grid) = &Elastic::Define;
	boost::python::scope the_scope = class_<Elastic,boost::noncopyable>("Elastic")
		.def("Define",Define_1)
		.def("Define",Define_2)
		.def("RefluxTest",&Elastic::RefluxTest)
		.def("TrigTest",&Elastic::TrigTest)
		.def("mytest",&Elastic::mytest)
		.def("setMaxCoarseningLevel",&Elastic::setMaxCoarseningLevel)
		.def("setFixedIter",&Elastic::setFixedIter)
		.def("setMaxIter",&Elastic::setMaxIter)
		.def("setMaxFmgIter",&Elastic::setMaxFmgIter)
		.def("setBottomMaxIter",&Elastic::setBottomMaxIter)
		.def("setBounds",&Elastic::setBounds)
		;

	enum_<Elastic::Grid>("Grid")
		.value("X",   Elastic::Grid::X)
		.value("Y",   Elastic::Grid::Y)
		.value("Z",   Elastic::Grid::Z)
		.value("YZ",  Elastic::Grid::YZ)
		.value("ZX",  Elastic::Grid::ZX)
		.value("XY",  Elastic::Grid::XY)
		.value("XYZ", Elastic::Grid::XYZ);

}
