#ifndef IC_CONSTANT_H_
#define IC_CONSTANT_H_

#include "IC/IC.H"
#include "IO/ParmParse.H"
#include "Util/Util.H"

/// \class Constant
/// \brief 
///
namespace IC
{
class Constant : public IC
{
public:
    Constant (amrex::Vector<amrex::Geometry> &_geom) :IC::IC(_geom) {}
    Constant (amrex::Vector<amrex::Geometry> &_geom, IO::ParmParse &pp) : IC::IC(_geom) 
    {pp.queryclass(*this);}
    Constant (amrex::Vector<amrex::Geometry> &_geom, IO::ParmParse &pp, std::string name) : IC::IC(_geom) 
    {pp.queryclass(name,*this);}

    Constant(amrex::Vector<amrex::Geometry> &_geom, std::vector<amrex::Real> a_value)
    : IC(_geom), m_value(a_value)
    {}
    ~Constant() { };
    void Add(const int &lev, Set::Field<Set::Scalar> &a_field)
    {
        Util::Assert(INFO,TEST((m_value.size() == 1 || (int)m_value.size() == (int)a_field[lev]->nComp())));
        for (amrex::MFIter mfi(*a_field[lev],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.growntilebox();
            amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
            for (int m = 0; m < a_field[lev]->nComp(); m++)
                amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {
                field(i,j,k,m) += m_value.size() == 1 ? m_value[0] : m_value[m];
            });  
        }
    }
    static void Parse(Constant & value, IO::ParmParse & pp)
    {
        // Default: set equal to one everywhere
        value.m_value.clear();
        value.m_value.push_back(1.0);
        pp.queryarr("value",value.m_value); // Array of constant values. The number of values should equal either 1 or N where N is the number of fab components
    }
private:
    std::vector<Set::Scalar> m_value;
};
}
#endif
