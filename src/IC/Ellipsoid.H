// This IC initializes a single-component fab using the formula
// 
// .. math::
//
//    \phi = \begin{cases} 
//                \alpha_{in}  & (\mathbf{x}-\mathbf{x}_0)^T\mathbf{A}(\mathbf{x}-\mathbf{x}_0) \le r  \\
//                \alpha_{out} & \text{else} 
//           \end{cases}
//
// The boundary can be mollified using an error function with parameter :math:`\varepsilon`
//
// .. WARNING::
//
//    This IC is redundant with :ref:`IC::Ellipse` and :ref:`IC::Expression` and will probably be consolidated
//

#ifndef IC_ELLIPSOID_H_
#define IC_ELLIPSOID_H_

#include "IC/IC.H"
#include "Util/Util.H"
#include "IO/ParmParse.H"

/// \class Ellipsoid
/// \brief Initialize an ellipsoidal inclusion
namespace IC
{
class Ellipsoid : public IC
{
public:
    enum Mollifier {Dirac, Gaussian};

    Ellipsoid (amrex::Vector<amrex::Geometry> &_geom) : IC(_geom) {}
  
    void Add(const int &lev,Set::Field<Set::Scalar> &a_field, Set::Scalar)
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
                    for (int i=0 ; i < center.size(); i++)
                    {
                        Set::Scalar value = 0.;
                        Set::Scalar norm = (A[i]*(x-center[i])).lpNorm<2>();
                        value = 0.5 + 0.5*std::erf(((x-center[i]).transpose() * A[i] * (x-center[i]) - 1.0) / eps[i] / norm);
                        value = in_value + value*(out_value - in_value);
                        min_value = value < min_value ? value : min_value;
                    }
                    field(i,j,k) = min_value;
                    if (field(i,j,k) < std::min(in_value,out_value)) field(i,j,k) = std::min(in_value,out_value);
                    if (field(i,j,k) > std::max(in_value,out_value)) field(i,j,k) = std::max(in_value,out_value);
                });
        }
        a_field[lev]->FillBoundary();
    };
  
private:
    Set::Scalar in_value = 0.0, out_value = 1.0;
    amrex::Vector<Set::Vector> center; 
    amrex::Vector<Set::Scalar> eps;
    amrex::Vector<Set::Matrix> A;
    Mollifier moll = Mollifier::Dirac;

public:
    static void Parse(Ellipsoid & value, IO::ParmParse & pp)
    {
        amrex::Vector<Set::Scalar> center;
        if(pp.contains("center")) pp_queryarr("center",center); // Center of the ellipse :math:`\mathbf{x}_0`

        if(center.size() < AMREX_SPACEDIM) value.center.push_back(Set::Vector::Zero());
        else
        {
            for (int i = 0; i<center.size(); i+=AMREX_SPACEDIM)
                value.center.push_back(Set::Vector(AMREX_D_DECL(center[i],center[i+1],center[i+2])));
        }
        
        amrex::Vector<Set::Scalar> A;
        if (pp.contains("A"))
        {
            pp_queryarr("A", A); // Matrix defining elipse radii and orientation  
            if(A.size() < AMREX_SPACEDIM*AMREX_SPACEDIM) value.A.push_back(Set::Matrix::Identity());
            else
            {
                int i=0, j=0;
                Set::Matrix temp = Set::Matrix::Zero();
                for (i=0; i<A.size(); i+=AMREX_SPACEDIM)
                {
                    for (j=0; j<AMREX_SPACEDIM; j++)
                        temp(i,j) = A[i+j];
                    
                    if((i+j) % (AMREX_SPACEDIM*AMREX_SPACEDIM) == 0) 
                    {
                        value.A.push_back(temp);
                        temp = Set::Matrix::Zero();
                    }
                }
            }
        }
        if (pp.contains("radius"))
        { 
            amrex::Vector<Set::Scalar> a_radius; 
            pp_queryarr("radius",a_radius); // "Vector of radii (use instead of A)"

            if (a_radius.size() < AMREX_SPACEDIM) value.A.push_back(Set::Matrix::Identity());
            else
            {
                for (int i = 0; i<a_radius.size(); i+=AMREX_SPACEDIM)
                {
                    Set::Matrix temp = Set::Matrix::Zero();
                    AMREX_D_TERM( temp(0,0) = 1./(a_radius[i]*a_radius[i]) ;, temp(1,1) = 1./(a_radius[i+1]*a_radius[i+1]) ;, temp(2,2) = 1./(a_radius[i+2]*a_radius[i+2]) ; );
                    value.A.push_back(temp);
                }
            }
        }

        amrex::Vector<Set::Scalar> eps;
        pp_queryarr("eps", eps); // Mollifying value for erf

        if(eps.size() < 1) value.eps.push_back(1.e-5);
        for (int i =0; i < value.A.size(); i++)
        {
            if (eps[i] <= 0.0) eps[i] = 1.e-5;
            value.eps.push_back(eps[i]);
        }
        //if(value.eps <= 0.) value.eps = 1.e-5;

        pp_query_default("in_value",  value.in_value,  0.0); // Value of field inside ellipse
        pp_query_default("out_value", value.out_value, 1.0); // Value of field outside ellipse

        std::string mollifier;
        pp_query("mollifier",mollifier); // Type of mollifier to use (options: dirac, [gaussian])
        if(mollifier == "Dirac" || mollifier == "dirac") 
            value.moll = Mollifier::Dirac;
        else 
            value.moll = Mollifier::Gaussian;
    }
};
}
#endif
