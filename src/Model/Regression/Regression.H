#ifndef MODEL_REGRESSION_REGRESSION_H
#define MODEL_REGRESSION_REGRESSION_H

#include "Set/Base.H"
#include "IO/ParmParse.H"


namespace Model
{
namespace Regression
{

class Regression
{
public:
    virtual Set::Scalar operator() (Set::Scalar phi, Set::Scalar T)  = 0;

    bool needs_normal = false;

    Set::Scalar P = NAN;
    virtual void set_pressure(const Set::Scalar P) 
    {
        this->P = P;
    }
};






}
}




#endif
