//
// This is a general-purpose routine for specifying time-interpolated
// quantities in an input file.
// Interpolators obey the following syntax:
// 
// .. code-block:: bash
// 
//    (t1,t2,..,tn;v1,v2,...,vn)
//
// where :code:`tn` are times and :code:`vn` are values at those times.
//
// :bdg-danger-line:`Note: do not include whitespace in an interpolator, this will mess with the parser.`
//
// Interpolators can usually be used in place of regular numbers, but only
// when supported by the calling parse function.
// 
// For instance, the `BC::Operator::Elastic::TensionTest` allows the user to specify
// fixed displacement; for instance, to specify a value of 0.1:
// 
// .. code-block:: bash
// 
//    bc.tenstion_test.disp = 0.1
// 
// However, it may be desirable to change the value over time.
// In this case, one may specify an interpolator string instead:
//
// .. code-block:: bash
// 
//    bc.tenstion_test.disp = (1.0,1.5,2.0;0.0,0.1,0.0)
// 
// This will cause the displacement value to be zero for :code:`t<1.0`; 
// ramp to :code:`0.1` as :code:`1.0<t<2.0`; 
// ramp back to :code:`0.0` as :code:`2.0<t<3.0`; 
// and then remain at the constant value of :code:`3.0` for :code:`t>3.0`.
//

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
        pp.query("str",in); // Interpolator string used when Parsed from queryclass.
        value.define(in);
    }


    static Linear<T> Read(std::string filename, int derivative = 0, bool copyInputFile = true);

protected:
    std::vector<T> data_points;
    std::vector<Set::Scalar> interval_points;
    bool uniform = true;

public:
    template<class U> friend Linear<U> operator * (const Set::Scalar alpha, const Linear<U> &b);
};

template<class T>
AMREX_FORCE_INLINE
Linear<T> operator * (const Set::Scalar alpha, const Linear<T> &b)
{
    Linear<T> ret;
    for (unsigned int i = 0; i < b.data_points.size(); i++) 
        ret.data_points.push_back(alpha * b.data_points[i]);
    ret.interval_points = b.interval_points;
    ret.uniform = b.uniform;
    return ret;
}

template<>
inline Linear<Set::Scalar> Linear<Set::Scalar>::Read(std::string filename, int derivative, bool copyInputFile)
{
    std::ifstream input;
    input.open(filename);
    std::string line;
    std::vector<Set::Scalar> theta, w;
    // parse input 
    while (std::getline(input, line))
    {        std::vector<std::string> dat = Util::String::Split(line);
        theta.push_back(std::stof(dat[0]));
        w.push_back(std::stof(dat[1]));
    }

    if (copyInputFile)  Util::CopyFileToOutputDir(filename);

    if (derivative == 0)
    {
        return Linear<Set::Scalar>(w,theta);
    }
    else if (derivative == 1)
    {
        std::vector<Set::Scalar> thetasmall, dw;
        for (unsigned int i = 1; i < theta.size() - 1; i++)
        {
            thetasmall.push_back(theta[i]);
            dw.push_back((w[i + 1] - w[i - 1]) / (theta[i + 1] - theta[i - 1]));
        }
        return Linear<Set::Scalar>(dw,thetasmall);
    }
    else if (derivative == 2)
    {
        std::vector<Set::Scalar> thetasmall, ddw;
        for (unsigned int i = 1; i < theta.size() - 1; i++)
        {
            thetasmall.push_back(theta[i]);
            ddw.push_back((w[i + 1] - 2.0 * w[i] + w[i - 1]) / ((theta[i + 1] - theta[i]) * (theta[i] - theta[i - 1])));
        }
        return Linear<Set::Scalar>(ddw,thetasmall);
    }
    else Util::Abort(INFO,"derivative must be 0,1,2 but got ",derivative);
    exit(-1);
}


}
}

#endif
