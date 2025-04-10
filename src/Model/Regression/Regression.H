#ifndef MODEL_REGRESSION_REGRESSION_H
#define MODEL_REGRESSION_REGRESSION_H

#include "Set/Base.H"
#include "Model/Model.H"

namespace Model
{
namespace Regression
{

class Regression
{
public:
    Regression() {} ;
    virtual ~Regression(){};

    virtual void
    setPressure(Set::Scalar const pressure)
    {};

    virtual Set::Scalar
    rate(Set::Scalar const ap_frac) = 0;

    const bool thermal_dependent = false;
    const bool pressure_dependent = false;
    const bool heterogeneous = false;
};
}
}

#endif
