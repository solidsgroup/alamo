// #include <AMReX_ParallelDescriptor.H>
#include "Linear.H"

// template<class T>
// Numeric::Interpolator::Linear::Linear<T>::Linear(const amrex::Vector<T> _data_points, const amrex::Vector<Set::Scalar> _interval_points):
// 	data_points(_data_points), interval_points(_interval_points){}

// template<class T>
// T
// Numeric::Interpolator::Linear::Linear<T>::operator() (Set::Scalar point)
// {
// 	if(data_points.size() != interval_points.size())
// 	{
// 		// some sort of an error or warning here
// 	}

// 	T data;
	
// 	for (int i = 0; i<interval_points.size(); i++)
// 	{
// 		if(point < interval_points[0])
// 		{
// 			data = data_points[0];
// 			break;
// 		}
// 		if (i < interval_points.size()-1 &&  interval_points[i] < point && point < interval_points[i+1])
// 		{
// 			data = data_points[i] +
// 				(point - interval_points[i]) * (data_points[i+1] - data_points[i]) /
// 				(interval_points[i+1] - interval_points[i]);
// 			break;
// 		}
// 		else
// 		{
// 			data = data_points[interval_points.size()-1];
// 			break;
// 		}
// 	}
// 	return data;
// }
