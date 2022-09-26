#ifndef NUMERIC_INTERPOLATOR_LINEAR_H_
#define NUMERIC_INTERPOLATOR_LINEAR_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include "Set/Set.H"
#include "Util/Util.H"
#include "Numeric/Interpolator/Interpolator.H"
#include "IO/ParmParse.H"

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
    Linear(const Set::Scalar &a_val)
    {
        interval_points.push_back(0.0);
        data_points.push_back(a_val);
    }
    Linear(const std::string a_str)
    {
        define(a_str);
    }
    Linear(const std::vector<T> _data_points, const std::vector<Set::Scalar> _interval_points)
    {
        define(_data_points,_interval_points);
    };
    
    void define(const std::string a_str)
    {
        std::string str = a_str;
        
        if (Util::String::Contains(str,"(") || Util::String::Contains(str,":") || Util::String::Contains(str,")") || Util::String::Contains(str,","))
        {
            if (!Util::String::Contains(str,"(")) Util::Abort(INFO,"Mismatched parentheses while trying to parse ",str);
            if (!Util::String::Contains(str,")")) Util::Abort(INFO,"Mismatched parentheses while trying to parse ",str);
            if (!Util::String::Contains(str,":")) Util::Abort(INFO,"Missing : in interpolator while trying to parse ",str);
            Util::String::ReplaceAll(str,"(","");
            Util::String::ReplaceAll(str,")","");

            std::vector<std::string> splitstr = Util::String::Split(str,':');
            std::string str_time = splitstr[0];
            std::string str_val  = splitstr[1];
            
            std::vector<std::string> str_time_arr = Util::String::Split(str_time,',');
            std::vector<std::string> str_val_arr  = Util::String::Split(str_val,',');

            if (str_time_arr.size() != str_val_arr.size()) Util::Abort(INFO,"Mismatched number of time values vs array values while trying to parse ", str);

            for (unsigned int i = 0; i < str_time_arr.size(); i++)
            {
                interval_points.push_back(std::stod(str_time_arr[i]));
                if (i > 0) if (interval_points[i] < interval_points[i-1]) Util::Abort(INFO,"Time series must monotonically increase - error while trying to parse ", str);
                data_points.push_back(std::stod(str_val_arr[i]));
            }
        }
        else
        {
            interval_points.push_back(0.0);
            data_points.push_back(std::stod(a_str));
        }
    }

    void define(const std::vector<T> _data_points, const std::vector<Set::Scalar> _interval_points) 
    {
        data_points = _data_points;
        interval_points = _interval_points;
        if(data_points.size() != interval_points.size())
            Util::Abort(INFO,"Data points and interval points have different sizes");
    };

    T operator() (const Set::Scalar point) const
    {
        // If there is only one point stored, return it.
        if (data_points.size() == 1) return data_points[0];

        // Do this if point is below interpolation range
        if(point < interval_points[0])
        {
            return data_points[0];
        }

        Set::Scalar interval = interval_points[1]-interval_points[0];
        int start = 0;
        if (uniform)
        {
            start = (int)( (point - interval_points[0]) / interval);
            start -= 1;
            start = std::max(start,0);
            start = std::min(start,(int)(interval_points.size()-2));
        }

        // Do this if point is within interpolation range
        for (unsigned int i = start; i<interval_points.size(); i++)
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

    static void Parse(Linear<T> & value, IO::ParmParse & pp)
    {
        std::string in;
        pp.query("str",in);
        value.define(in);
    }

protected:
    std::vector<T> data_points;
    std::vector<Set::Scalar> interval_points;
    bool uniform = true;

public:
    template<class U> friend Linear<U> operator * (const Set::Scalar alpha, const Linear<U> &b);
};

template<class T>
inline Linear<T> operator * (const Set::Scalar alpha, const Linear<T> &b)
{
    Linear<T> ret;
    for (unsigned int i = 0; i < b.data_points.size(); i++) 
        ret.data_points.push_back(alpha * b.data_points[i]);
    ret.interval_points = b.interval_points;
    ret.uniform = b.uniform;
    return ret;
}

}
}

#endif
