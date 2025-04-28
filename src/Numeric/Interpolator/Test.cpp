#include "Test.H"
#include "Linear.H"
#include "Util/Util.H"
#include "Set/Set.H"
#include "Util/Color.H"

namespace Numeric
{
namespace Interpolator
{

template<>
int Test<Linear<Set::Scalar> >::Match(int verbose)
{
    // Test function:
    //
    // y(x) = { -1           x  < -1
        //        { 2*x + 1      -1 < x  < 0
    //        { 1-x          0  < x  < 1
    //        { 0            1  < x

    std::vector<Set::Scalar> xs;
    std::vector<Set::Scalar> ys;
    xs.push_back(-1.0); ys.push_back(-1.0);
    xs.push_back(0.0);  ys.push_back(1.0);
    xs.push_back(1.0);  ys.push_back(0.0);

    Linear<Set::Scalar> interp(ys,xs);

    Set::Scalar normsq = 0.0;
    const Set::Scalar dx = 0.01;
    for (Set::Scalar x = -2.0; x < -1.0; x+=dx)
    {
        Set::Scalar exact = -1.0;
        normsq += pow((interp(x) - exact)/dx,2.0);
        if (verbose>0) Util::Message(INFO,(normsq>1E-8 ? Color::FG::Red : Color::Reset), "x = ", x , "\texact = ", exact, "\tinterp = ", interp(x), " normsq = ", normsq,Color::Reset);
    }
    for (Set::Scalar x = -1.0; x < -0; x+=dx)
    {
        Set::Scalar exact = 2*x + 1;
        normsq += pow((interp(x) - exact)/dx,2.0);
        if (verbose>0) Util::Message(INFO,(normsq>1E-8 ? Color::FG::Red : Color::Reset), "x = ", x , "\texact = ", exact, "\tinterp = ", interp(x), " normsq = ", normsq,Color::Reset);
    }
    for (Set::Scalar x = 0; x < 1.0; x+=dx)
    {
        Set::Scalar exact = 1.0 - x;
        normsq += pow((interp(x) - exact)/dx,2.0);
        if (verbose>0) Util::Message(INFO,(normsq>1E-8 ? Color::FG::Red : Color::Reset), "x = ", x , "\texact = ", exact, "\tinterp = ", interp(x), " normsq = ", normsq,Color::Reset);
    }
    for (Set::Scalar x = 1.0; x < 2.0; x+=dx)
    {
        Set::Scalar exact = 0.0;
        normsq += pow((interp(x) - exact)/dx,2.0);
        if (verbose>0) Util::Message(INFO,(normsq>1E-8 ? Color::FG::Red : Color::Reset), "x = ", x , "\texact = ", exact, "\tinterp = ", interp(x), " normsq = ", normsq,Color::Reset);
    }

    if (normsq > 1E-8) return 1;
    else return 0;
}
}
}
