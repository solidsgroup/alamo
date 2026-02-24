#ifndef TEST_MODEL_INTERFACE_GB_H
#define TEST_MODEL_INTERFACE_GB_H

#include "Set/Set.H"
#include "Model/Interface/GB/GB.H"

namespace Test
{
namespace Model
{
namespace Interface
{
namespace GB
{
template <class T>
class GB
{
public:
    GB() {};
    bool DerivativeTest1(int verbose)
    {
        int failed = 0;
        amrex::Real small = 1E-12;
        amrex::Real tolerance = 1E-3;
   
        T model;
        model.Randomize();

        for (int i = 0; i<20; i++)
        {
            amrex::Real theta = 2.0*Set::Constant::Pi*((amrex::Real)rand()/(amrex::Real)RAND_MAX);

            amrex::Real numerical_DW = (model.W(theta+small) - model.W(theta-small))/(2.0*small);
            amrex::Real exact_DW     = model.DW(theta);
            if (fabs(numerical_DW-exact_DW)/(numerical_DW+exact_DW) > tolerance)
                failed += 1;
            if (verbose)
            {
                Util::Message(INFO,"Theta:         " , theta);
                Util::Message(INFO,"DW Exact:      " , exact_DW);
                Util::Message(INFO,"DW Numerical:  " , numerical_DW);
            }

        }
        return failed;
    };
    bool DerivativeTest2(int verbose)
    {
        int failed = 0;
        amrex::Real small = 1E-8;
        amrex::Real tolerance = 1E-6;
   
        T model;
        model.Randomize();

        for (int i = 0; i<20; i++)
        {
            amrex::Real theta = 2.0*Set::Constant::Pi*((amrex::Real)rand()/(amrex::Real)RAND_MAX);

            amrex::Real numerical_DDW = (model.DW(theta+small) - model.DW(theta-small))/(2.0*small);
            amrex::Real exact_DDW = model.DDW(theta);
            if (fabs((numerical_DDW-exact_DDW)/(numerical_DDW+exact_DDW)) > tolerance)
                failed += 1;

            if (verbose)
            {
                Util::Message(INFO,"Theta:         " , theta);
                Util::Message(INFO,"DDW Exact:     " , exact_DDW);
                Util::Message(INFO,"DDW Numerical: " , numerical_DDW);
            }

        }
        return failed;
    };
};
}
}
}
}

#endif
