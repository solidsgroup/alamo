#ifndef NUMERIC_INTERPOLATOR_LINEAR_H_
#define NUMERIC_INTERPOLATOR_LINEAR_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include "Set/Set.H"
#include "Util/Util.H"
#include "Numeric/Interpolator/Interpolator.H"

namespace Numeric
{
namespace Interpolator
{
template <class T>
class Linear : public Interpolator<T>
{
public:
 	Linear(const amrex::Vector<T> _data_points, const amrex::Vector<Set::Scalar> _interval_points) :
		data_points(_data_points), interval_points(_interval_points)
	{
		if(data_points.size() != interval_points.size())
			Util::Abort(INFO,"Data points and interval points have different sizes");
	};

 	T operator() (Set::Scalar point)
	{

		// Do this if point is below interpolation range
		if(point <= interval_points[0])
		{
			return data_points[0];
		}
		// Do this if point is above interpolation range
		if(point >= interval_points[interval_points.size()-1])
		{
			return  data_points[interval_points.size()-1];
		}

		// Do this if point is within interpolation range
		for (int i = 0; i<interval_points.size(); i++)
		{
			if (interval_points[i] < point && point < interval_points[i+1])
			{
				return data_points[i] +
						 (point - interval_points[i]) * (data_points[i+1] - data_points[i]) /
						 (interval_points[i+1] - interval_points[i]);

			}
		}
	}


protected:
 	amrex::Vector<T> data_points;
 	amrex::Vector<Set::Scalar> interval_points;
};
}
}

#endif
