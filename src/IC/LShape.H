#ifndef IC_LSHAPE_H_
#define IC_LSHAPE_H_

#include "IC/IC.H"
#include "Util/Util.H"
#include "IO/ParmParse.H"

/// \class L-shape
/// \brief Initialize an L-shape inclusion
namespace IC
{
class LShape : public IC
{
public:
    enum Mollifier {Dirac, Gaussian};

    LShape (amrex::Vector<amrex::Geometry> &_geom) : IC(_geom) {}
  
    void Add(const int &lev,Set::Field<Set::Scalar> &a_field)
    {
        for (amrex::MFIter mfi(*a_field[lev],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            bx.grow(a_field[lev]->nGrow());
            amrex::IndexType type = a_field[lev]->ixType();

            amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
            amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {

                    Set::Vector x = Set::Position(i,j,k,geom[lev],type);

                    Set::Scalar min_value = field(i,j,k);

                    Set::Scalar value = 1- 0.5*(1+std::erf(x[0]/eps[0])) * 0.5*(1+std::erf(x[1]/eps[1]));
                    value = in_value + value*(out_value - in_value);
                    min_value = value < min_value ? value : min_value;

                    field(i,j,k) = min_value;
                    if (field(i,j,k) < std::min(in_value,out_value)) field(i,j,k) = std::min(in_value,out_value);
                    if (field(i,j,k) > std::max(in_value,out_value)) field(i,j,k) = std::max(in_value,out_value);
                });
        }
        a_field[lev]->FillBoundary();
    };
  
private:
    int dim = AMREX_SPACEDIM;
    Set::Scalar in_value = 0.0, out_value = 1.0;
    amrex::Vector<Set::Scalar> eps;
    Mollifier moll = Mollifier::Dirac;

public:
    // This IC initializes a single-component fab using the formula
    // 
    // :math:`\phi = \begin{cases} \alpha_{in}  & (\mathbf{x}-\mathbf{x}_0)^T\mathbf{A}(\mathbf{x}-\mathbf{x}_0) \le r /
    //                          \\ \alpha_{out} & \text{else} \end{cases}`
    //
    // The boundary can be mollified using an error function with parameter :math:`\varepsilon`
    static void Parse(LShape & value, IO::ParmParse & pp)
    {
       amrex::Vector<Set::Scalar> eps;
        pp.queryarr("eps", eps);

        if(eps.size() < 1) value.eps.push_back(1.e-5);
        for (int i =0; i < AMREX_SPACEDIM; i++)
        {
            if (eps[i] <= 0.0) eps[i] = 1.e-5;
            value.eps.push_back(eps[i]);
        }
        //if(value.eps <= 0.) value.eps = 1.e-5;

        pp.query("in_value", value.in_value);
        pp.query("out_value", value.out_value);

        std::string mollifier;
        pp.query("mollifier",mollifier);
        if(mollifier == "Dirac" || mollifier == "dirac") 
            value.moll = Mollifier::Dirac;
        else 
            value.moll = Mollifier::Gaussian;
    }
};
}
#endif
