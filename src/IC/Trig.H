#ifndef IC_TRIG_H_
#define IC_TRIG_H_

#include <complex>

#include "AMReX_Vector.H"
#include "IC/IC.H"
#include "Util/Util.H"
#include "Set/Set.H"

namespace IC
{
/// \brief Initialize using a trigonometric series
class Trig : public IC
{
public:
    Trig (amrex::Vector<amrex::Geometry> &_geom,
        Set::Scalar _alpha = 1.0,
        AMREX_D_DECL(std::complex<int> _n1 = 0,
                std::complex<int> _n2 = 0,
                std::complex<int> _n3 = 0),
        int _dim = AMREX_SPACEDIM) :
        IC(_geom)
    {
        Define(_alpha,AMREX_D_DECL(_n1,_n2,_n3),_dim);
    }


    void Define(Set::Scalar _alpha = 1.0,
            AMREX_D_DECL(std::complex<int> _n1 = 0,
                std::complex<int> _n2 = 0,
                std::complex<int> _n3 = 0),
            int _dim = AMREX_SPACEDIM)
    {
        alpha = _alpha;
        AMREX_D_TERM(n1 = _n1;, n2 = _n2;, n3 = _n3;);
        dim = _dim;
        AMREX_D_DECL(phi1 = std::atan2(n1.imag(),n1.real()),
                phi2 = std::atan2(n2.imag(),n2.real()),
                phi3 = std::atan2(n3.imag(),n3.real()));
    }


    void Add(const int &lev, Set::Field<Set::Scalar> &field)
    {
        bool cellcentered = (field[0]->boxArray().ixType() == amrex::IndexType(amrex::IntVect::TheCellVector()));

        const amrex::Real AMREX_D_DECL(L1 = geom[lev].ProbHi()[0] - geom[lev].ProbLo()[0],
                            L2 = geom[lev].ProbHi()[1] - geom[lev].ProbLo()[1],
                            L3 = geom[lev].ProbHi()[2] - geom[lev].ProbLo()[2]);
                   
        for (amrex::MFIter mfi(*field[lev],true); mfi.isValid(); ++mfi)
        {
            const amrex::Box& box = mfi.validbox();

            amrex::BaseFab<amrex::Real> &field_box = (*field[lev])[mfi];

            AMREX_D_TERM(for (int i = box.loVect()[0] - field[lev]->nGrow(); i<=box.hiVect()[0] + field[lev]->nGrow(); i++),
                    for (int j = box.loVect()[1] - field[lev]->nGrow(); j<=box.hiVect()[1] + field[lev]->nGrow(); j++),
                    for (int k = box.loVect()[2] - field[lev]->nGrow(); k<=box.hiVect()[2] + field[lev]->nGrow(); k++))
            {
                Set::Scalar AMREX_D_DECL(x1,x2,x3);
                if (cellcentered)
                {
                    AMREX_D_TERM(x1 = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];,
                            x2 = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];,
                            x3 = geom[lev].ProbLo()[2] + ((amrex::Real)(k) + 0.5) * geom[lev].CellSize()[2];);
                }
                else
                {
                    AMREX_D_TERM(x1 = geom[lev].ProbLo()[0] + ((amrex::Real)(i)) * geom[lev].CellSize()[0];,
                            x2 = geom[lev].ProbLo()[1] + ((amrex::Real)(j)) * geom[lev].CellSize()[1];,
                            x3 = geom[lev].ProbLo()[2] + ((amrex::Real)(k)) * geom[lev].CellSize()[2];);
                }

                Set::Scalar trigfn = 1.0;
#if AMREX_SPACEDIM > 0
                if (dim > 0)
                    trigfn *= (fabs(std::cos(phi1))*std::cos(n1.real()*Set::Constant::Pi*x1 / L1) +
                            fabs(std::sin(phi1))*std::sin(n1.imag()*Set::Constant::Pi*x1 / L1));
#endif
#if AMREX_SPACEDIM > 1
                if (dim > 1)
                    trigfn *= (fabs(std::cos(phi2))*std::cos(n2.real()*Set::Constant::Pi*x2 / L2) +
                            fabs(std::sin(phi2))*std::sin(n2.imag()*Set::Constant::Pi*x2 / L2));
#endif
#if AMREX_SPACEDIM > 2
                if (dim > 2)
                    trigfn *= (fabs(std::cos(phi3))*std::cos(n3.real()*Set::Constant::Pi*x3 / L3) +
                            fabs(std::sin(phi3))*std::sin(n3.imag()*Set::Constant::Pi*x3 / L3));
#endif
                        
                field_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),comp) += alpha * trigfn;

            }
        }
    }
    using IC::Add;

private:
    int dim = AMREX_SPACEDIM;
    Set::Scalar alpha;
    std::complex<int> AMREX_D_DECL(n1, n2, n3);
    Set::Scalar AMREX_D_DECL(phi1=0.0,phi2=0.0,phi3=0.0);
};
}
#endif
