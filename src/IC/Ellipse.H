// If :code:`number_of_inclusions` is specified, then multiple ellipses are specified.
// In this case, each parameter must have number_of_inclusion*M values, where M is the
// number of values specified for the single ellipse case.
//

#ifndef IC_ELLIPSE_H_
#define IC_ELLIPSE_H_

#include "Set/Set.H"
#include "IC/IC.H"
#include "IO/ParmParse.H"

namespace IC
{
class Ellipse : public IC<Set::Scalar>
{
public:
    static constexpr const char* name = "ellipse";

    enum Mollifier {Dirac, Gaussian};

    Ellipse (amrex::Vector<amrex::Geometry> &_geom) : IC(_geom) {}
    Ellipse (amrex::Vector<amrex::Geometry> &_geom, IO::ParmParse &pp, std::string name) : Ellipse(_geom)
    {pp_queryclass(name,*this);}
    
    void Add(const int &lev, Set::Field<Set::Scalar> &a_field, Set::Scalar)
    {
        Set::Vector DX(geom[lev].CellSize());
        amrex::IndexType type = a_field[lev]->ixType();
        int ncomp = a_field[lev]->nComp();
            
        for (amrex::MFIter mfi(*a_field[lev],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx;
            if (type == amrex::IndexType::TheNodeType()) bx = mfi.grownnodaltilebox();
            if (type == amrex::IndexType::TheCellType()) bx = mfi.growntilebox();
            
            amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
            amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {

                Set::Vector x;
                // NODE
                if (type == amrex::IndexType::TheNodeType())
                {
                    AMREX_D_TERM(x(0) = geom[lev].ProbLo()[0] + ((amrex::Real)(i)) * geom[lev].CellSize()[0];,
                                x(1) = geom[lev].ProbLo()[1] + ((amrex::Real)(j)) * geom[lev].CellSize()[1];,
                                x(2) = geom[lev].ProbLo()[2] + ((amrex::Real)(k)) * geom[lev].CellSize()[2];);
                }
                else if (type == amrex::IndexType::TheCellType())
                {
                    AMREX_D_TERM(x(0) = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];,
                                x(1) = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];,
                                x(2) = geom[lev].ProbLo()[2] + ((amrex::Real)(k) + 0.5) * geom[lev].CellSize()[2];);
                }

                if(number_of_inclusions == 0)
                {   
                    Set::Scalar norm = (A[0]*(x-center[0])).lpNorm<2>();
                    field(i,j,k,0) = 0.5 - 0.5*std::erf(((x-center[0]).transpose() * A[0] * (x-center[0]) - 1.0) / eps[0] / norm);
                    if (invert) field(i,j,k,0) = 1.0 - field(i,j,k,0);
                    if (ncomp > 1) field(i,j,k,1) = 1.0 - field(i,j,k,0);
                }
                else
                {
                    // Base matrix is 0, rest of the inclusions are numbered from 1.
                    Set::Scalar value = 0.0;
                    for (int m = 0; m < number_of_inclusions; m++)
                    {
                        Set::Scalar norm = (A[m]*(x-center[m])).lpNorm<2>();
                        value += 0.5 - 0.5*std::erf(((x-center[m]).transpose() * A[m] * (x-center[m]) - 1.0) / eps[m] / norm);
                        // if (field(i,j,k,m+1) < 0.) field(i,j,k,m+1) = 0.;
                        // if (field(i,j,k,m+1) > 1.) field(i,j,k,m+1) = 1.;
                        // value += field(i,j,k,m+1);
                    }
                    field (i,j,k,invert) = 1.0 - value;
                    field (i,j,k,1-invert) = value;
                    if (field(i,j,k,invert) < 0.) field(i,j,k,invert) = 0.;
                    if (field(i,j,k,invert) > 1.) field(i,j,k,invert) = 1.;
                    if (field(i,j,k,1-invert) < 0.) field(i,j,k,1-invert) = 0.;
                    if (field(i,j,k,1-invert) > 1.) field(i,j,k,1-invert) = 1.;
                }
                
            });
        }
        a_field[lev]->FillBoundary();
    }
    
private:
    int number_of_inclusions = 1;
    amrex::Vector<Set::Vector> center;
    amrex::Vector<Set::Matrix> A;
    amrex::Vector<Set::Scalar> eps;
    int invert = 0;

public:
    static void Parse(Ellipse & value, IO::ParmParse & pp)
    {
        amrex::Vector<Set::Scalar> x0;
        if(!(pp.contains("number_of_inclusions")))
        {
            value.center.resize(0);
            value.A.resize(0);
            value.eps.resize(0);
            
            value.number_of_inclusions = 0;
            pp_queryarr("x0",x0); // Coorinates of ellipse center
            value.center.push_back(Set::Vector(AMREX_D_DECL(x0[0],x0[1],x0[2])));

            Set::Scalar _eps;
            pp_query_default("eps",_eps,0.0); // Diffuse boundary thickness
            value.eps.push_back(_eps);

            Set::Matrix _A = Set::Matrix::Zero();
            if (pp.contains("A"))
            {
                pp_queryarr("A",_A); // DxD square matrix defining an ellipse. 
                value.A.push_back(_A);
            }
            else if (pp.contains("a"))
            {
                Set::Matrix _A = Set::Matrix::Zero();
                Set::Vector a = Set::Vector::Ones();
                pp_queryarr("a",a); // If :code:`A` is not defined, then assume a sphere with radius :code:`a`
                for (int d = 0; d < AMREX_SPACEDIM; d++) _A(d,d) = 1./a(d)/a(d);
                value.A.push_back(_A);
            }
        }
        else
        {
            pp_query("number_of_inclusions",value.number_of_inclusions); // Number of ellipses
            if(value.number_of_inclusions < 1) Util::Abort(INFO, "number of inclusions have to be at least 1");
            
            value.center.resize(0);
            value.A.resize(0);
            value.eps.resize(0);

            if (pp.contains("center") && pp.contains("x0")) Util::Abort(INFO,"Cannot specify both center (depricated) and x0");
            pp_queryarr("center", x0); // center of the ellipse
            pp_queryarr("x0", x0); // center of the ellipse
            
            if(x0.size() != value.number_of_inclusions*AMREX_SPACEDIM){
                Util::Message(INFO, value.number_of_inclusions*AMREX_SPACEDIM);
                Util::Message(INFO, x0.size());
                Util::Abort(INFO, "Need centers for all the inclusions");}
            for (int i = 0; i < x0.size(); i+= AMREX_SPACEDIM)
                value.center.push_back(Set::Vector(AMREX_D_DECL(x0[i],x0[i+1],x0[i+2])));

            if (pp.contains("A"))
            {
                amrex::Vector<Set::Scalar> _A;
                pp_queryarr("A", _A); // either a vector containing ellipse radii, or a matrix defining the ellipse
                if(_A.size() != value.number_of_inclusions*AMREX_SPACEDIM*AMREX_SPACEDIM && _A.size() != AMREX_SPACEDIM*AMREX_SPACEDIM)
                    Util::Abort(INFO, "Invalid value of A for ellipse initialization");
                if(_A.size() ==  AMREX_SPACEDIM*AMREX_SPACEDIM)
                {
                    Set::Matrix _A1 = Set::Matrix::Zero();
                    pp_queryarr("A",_A1); // Same
                    for (int i = 0; i < value.number_of_inclusions; i++)
                        value.A.push_back(_A1);
                } 
                else
                {
                    Set::Matrix _A1 = Set::Matrix::Zero();
                    for (int i = 0; i < value.number_of_inclusions; i+= AMREX_SPACEDIM*AMREX_SPACEDIM)
                    {
                        AMREX_D_PICK(   _A1(0,0) = _A[i];
                                        ,
                                        _A1(0,0) = _A[i]; _A1(0,1) = _A[i+1];
                                        _A1(1,0) = _A[i+2]; _A1(1,1) = _A[i+3];
                                        ,
                                        _A1(0,0) = _A[i+0]; _A1(0,1) = _A[i+1]; _A1(0,2) = _A[i+2];
                                        _A1(1,0) = _A[i+3]; _A1(1,1) = _A[i+4]; _A1(1,2) = _A[i+5];
                                        _A1(2,0) = _A[i+6]; _A1(2,1) = _A[i+7]; _A1(2,2) = _A[i+8];
                                    );
                        value.A.push_back(_A1);
                    }
                }
            }
            else if (pp.contains("radius"))
            {
                amrex::Vector<Set::Scalar> _r;
                pp_queryarr("radius", _r); // Array of radii [depricated]

                if(_r.size() != value.number_of_inclusions*AMREX_SPACEDIM && _r.size() != AMREX_SPACEDIM)
                    Util::Abort(INFO, "Invalid value of radius for ellipse initialization");
                if(_r.size() ==  AMREX_SPACEDIM)
                {
                    Set::Matrix _A1 = Set::Matrix::Zero();
                    for (int i = 0; i< AMREX_SPACEDIM; i++) _A1(i,i) = 1.0/(_r[i]*_r[i]);
                    for (int i = 0; i < value.number_of_inclusions; i++)
                        value.A.push_back(_A1);
                } 
                else
                {
                    for (int i = 0; i < value.number_of_inclusions; i++)
                    {
                        Set::Matrix _A1 = Set::Matrix::Zero();
                        for (int d = 0; d <  AMREX_SPACEDIM; d++) _A1(d,d) = 1.0/_r[i+d]/_r[i+d];
                        value.A.push_back(_A1);
                    }
                }
            }

            amrex::Vector<Set::Scalar> _eps;
            pp_queryarr("eps",_eps); // Regularization for smooth boundary
            if(_eps.size() != 1 && _eps.size() != value.number_of_inclusions)
                Util::Abort(INFO, "Incorrect eps. Check the number of values specified");
            
            if (_eps.size() == 1)
                for (int i = 0; i < value.number_of_inclusions; i++) value.eps.push_back(_eps[0]);
            
            else
                for (int i = 0; i < value.number_of_inclusions; i++) value.eps.push_back(_eps[i]);
        }
        // Flip the inside and the outside
        pp_query("invert",value.invert);
    }
};
}
#endif
