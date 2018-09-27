#ifndef NUMERIC_INTERPOLATOR_LINEAR_LINEAR_H_
#define NUMERIC_INTERPOLATOR_LINEAR_LINEAR_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include "Set/Set.H"
#include "Numeric/Interpolator/Interpolator.H"

namespace Numeric
{
namespace Interpolator
{
namespace Linear
{

template <class T>
class Linear : public Interpolator<T>
{
public:
 	Linear(const amrex::Vector<T> _data_points, const amrex::Vector<Set::Scalar> _interval_points) :
		data_points(_data_points), interval_points(_interval_points) {};

 	T operator() (Set::Scalar point)
	{
		if(data_points.size() != interval_points.size())
		{
			///\todo some sort of an error or warning here
		}

		T data;
	
		for (int i = 0; i<interval_points.size(); i++)
		{
			if(point < interval_points[0])
			{
				data = data_points[0];
				break;
			}
			if (i < interval_points.size()-1 &&  interval_points[i] < point && point < interval_points[i+1])
			{
				data = data_points[i] +
					(point - interval_points[i]) * (data_points[i+1] - data_points[i]) /
					(interval_points[i+1] - interval_points[i]);
				break;
			}
			else
			{
				data = data_points[interval_points.size()-1];
				break;
			}
		}
		return data;
	}

protected:
 	amrex::Vector<T> data_points;
 	amrex::Vector<Set::Scalar> interval_points;
};
}
}
}

#endif
