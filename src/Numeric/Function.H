#ifndef NUMERIC_FUNCTION_H
#define NUMERIC_FUNCTION_H


namespace Numeric
{
namespace Function
{

template<int d> class Polynomial;

template<>
class Polynomial<3>
{
public:
    Polynomial(Set::Scalar a0, Set::Scalar a1, Set::Scalar a2, Set::Scalar a3)
    {
        a[0] = a0; a[1] = a1; a[2] = a2; a[3] = a3;
    }
    AMREX_FORCE_INLINE
        Set::Scalar operator() (const Set::Scalar x) const
    {
        return a[0] + x * (a[1] + x * (a[2] + x * a[3]));
    }
private:
    std::array<Set::Scalar, 4> a;
};

template<>
class Polynomial<4>
{
public:
    Polynomial(Set::Scalar a0, Set::Scalar a1, Set::Scalar a2, Set::Scalar a3, Set::Scalar a4)
    {
        a[0] = a0; a[1] = a1; a[2] = a2; a[3] = a3; a[4] = a4;
    }
    AMREX_FORCE_INLINE
        Set::Scalar operator() (const Set::Scalar x) const
    {
        return a[0] + x * (a[1] + x * (a[2] + x * (a[3] + x * a[4])));
    }
    Polynomial<3> D() const
    {
        return Polynomial<3>(a[1], 2. * a[2], 3 * a[3], 4 * a[4]);
    }
private:
    std::array<Set::Scalar, 5> a;
};

}
}






#endif
