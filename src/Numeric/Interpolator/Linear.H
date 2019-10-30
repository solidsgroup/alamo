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
 	Linear() {};
 	~Linear() {};
 	Linear(const std::vector<T> _data_points, const std::vector<Set::Scalar> _interval_points)
	{
		define(_data_points,_interval_points);
	};
 	void define(const std::vector<T> _data_points, const std::vector<Set::Scalar> _interval_points) 
	{
		data_points = _data_points;
		interval_points = _interval_points;
		if(data_points.size() != interval_points.size())
			Util::Abort(INFO,"Data points and interval points have different sizes");
	};

 	T operator() (const Set::Scalar point) const
	{

		// Do this if point is below interpolation range
		if(point < interval_points[0])
		{
			return data_points[0];
		}
		// Do this if point is within interpolation range
		for (unsigned int i = 0; i<interval_points.size(); i++)
		{
			if (interval_points[i] <= point && point <= interval_points[i+1])
			{
				return data_points[i] +
						 (point - interval_points[i]) * (data_points[i+1] - data_points[i]) /
						 (interval_points[i+1] - interval_points[i]);

			}
		}

		// Otherwise point must be outside interpolation range, so
		// return the value for the highest point
		return  data_points[interval_points.size()-1];
	}


protected:
 	std::vector<T> data_points;
 	std::vector<Set::Scalar> interval_points;
};
}
}

#endif
