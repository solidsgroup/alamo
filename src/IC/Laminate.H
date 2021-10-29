#ifndef IC_LAMINATE_H_
#define IC_LAMINATE_H_

#include "IC/IC.H"
#include "Util/Util.H"
#include "IO/ParmParse.H"

/// \class Laminate
/// \brief Initialize Laminates in a matrix
namespace IC
{
class Laminate : public IC
{
public:
    enum Mollifier {Dirac, Gaussian};

    Laminate (amrex::Vector<amrex::Geometry> &_geom) : IC(_geom) {}
    
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

                    // The base matrix is always field(i,j,k,0)
                    // Inclusions start with field(i,j,k,1)

                    Set::Scalar value = 0.0;
                    for (int m = 0; m < number_of_inclusions; m++)
                    {
                        Set::Scalar t = std::abs((x-center[m]).transpose()*orientation[m]) - (thickness[m]/2.0);
                        field (i,j,k,m+1) = 0.5 - 0.5*std::erf(t / eps[m]);
                        if (field(i,j,k,m+1) < 0.) field(i,j,k,m+1) = 0.;
                        if (field(i,j,k,m+1) > 1.) field(i,j,k,m+1) = 1.;
                        value += field(i,j,k,m+1);
                    }
                    field (i,j,k,0) = 1.0 - value;
                    if (field(i,j,k,0) < 0.) field(i,j,k,0) = 0.;
                    if (field(i,j,k,0) > 1.) field(i,j,k,0) = 1.;
                });
        }
        a_field[lev]->FillBoundary();
    };
  
private:
    int number_of_inclusions = 1;
    amrex::Vector<Set::Vector> center; 
    amrex::Vector<Set::Vector> orientation;
    amrex::Vector<Set::Vector> normal;
    amrex::Vector<Set::Scalar> eps;
    amrex::Vector<Set::Scalar> thickness;
    Mollifier moll = Mollifier::Gaussian;

public:
    static void Parse(Laminate & value, IO::ParmParse & pp)
    {
        // How many laminates (MUST be greater than or equal to 1). Default = 1
        pp.query("number_of_inclusions", value.number_of_inclusions);
        
        if (value.number_of_inclusions < 1)
            Util::Abort(INFO, "Number of inclusions must be at least 1. Aborting.");
        
        amrex::Vector<Set::Scalar> a_center;
        // (x,y,[z]) values for the center point of the laminate
        if(pp.contains("center")) pp.queryarr("center",a_center);

        if(a_center.size() != value.number_of_inclusions*AMREX_SPACEDIM) value.center.push_back(Set::Vector::Zero());
        else
        {
            for (int i = 0; i<a_center.size(); i+=AMREX_SPACEDIM)
                value.center.push_back(Set::Vector(AMREX_D_DECL(a_center[i],a_center[i+1],a_center[i+2])));
        }
        
        amrex::Vector<Set::Scalar> a_thickness; 
        // thickness of the laminate
        if (pp.contains("thickness")) pp.queryarr("thickness", a_thickness);

        if(a_thickness.size() != value.number_of_inclusions && a_thickness.size() != 1)
            Util::Abort(INFO, "Thickness of each inclusion must be specified");
        
        if (a_thickness.size() == 1) 
        {
            if(a_thickness[0] <= 0.0) Util::Abort(INFO, "Invalid value of inclusion thickness");
            for (int i = 0; i<value.number_of_inclusions; i++) value.thickness.push_back(a_thickness[0]);
        }
        
        else
        {
            for (int i = 0; i<value.number_of_inclusions; i++)
            {
                if(a_thickness[i] <= 0.0) Util::Abort(INFO, "Invalid value of inclusion ", i+1, " thickness");
                value.thickness.push_back(a_thickness[i]);
            }
        }
        
        amrex::Vector<Set::Scalar> a_orientation;
        // Vector normal to the interface of the laminate
        pp.queryarr("orientation",a_orientation);

        if(a_orientation.size() != value.number_of_inclusions*AMREX_SPACEDIM && a_orientation.size() != AMREX_SPACEDIM) 
            Util::Abort(INFO, "Orientation of each inclusion must be specified");

        if(a_orientation.size() == AMREX_SPACEDIM)
            for (int i = 0; i < value.number_of_inclusions; i++)
                value.orientation.push_back(Set::Vector(AMREX_D_DECL(a_orientation[0], a_orientation[1], a_orientation[2])));
        
        else
            for (int i=0; i < a_orientation.size(); i+=AMREX_SPACEDIM)
                value.orientation.push_back(Set::Vector(AMREX_D_DECL(a_orientation[i], a_orientation[i+1], a_orientation[i+2])));
        
        for (int i =0; i < value.orientation.size(); i++)
            if(value.orientation[i].lpNorm<2>() <= 0.) value.orientation[i] = Set::Vector::Random();

        for (int i=0; i<value.orientation.size(); i++)
            value.orientation[i] = value.orientation[i] / value.orientation[i].lpNorm<2>();

        value.normal.resize(value.number_of_inclusions);
        for (int i = 0; i<value.orientation.size(); i++)
        {
            value.normal[i] = Set::Vector::Zero();
            if(value.orientation[i](0) != 0.)
            {
                AMREX_D_TERM(value.normal[i](0) = 1.;, value.normal[i](1) = 1.;, value.normal[i](2) = 1.;);
                value.normal[i](0) = -(AMREX_D_TERM(0.,+value.orientation[i](1),+value.orientation[i](2)))/value.orientation[i](0);
                value.normal[i] = value.normal[i]/value.normal[i].lpNorm<2>();
            }
            else if(value.orientation[i](1) != 0.)
            {
                AMREX_D_TERM(value.normal[i](0) = 1.;, value.normal[i](1) = 1.;, value.normal[i](2) = 1.;);
                value.normal[i](1) = -(AMREX_D_TERM(value.orientation[i](0), + 0.0, + value.orientation[i](2)))/value.orientation[i](1);
                value.normal[i] = value.normal[i]/value.normal[i].lpNorm<2>();
            }
        }

        amrex::Vector<Set::Scalar> a_eps;
        // Diffuse thickness
        pp.queryarr("eps", a_eps);

        if(a_eps.size() < 1) 
            for (int i = 0; i< value.number_of_inclusions; i++)
                value.eps.push_back(1.e-5);
        if(a_eps.size() == 1)
        {
            if (a_eps[0] < 0.0) 
            {
                Util::Warning(INFO, "Invalid value of laminate.eps. Resetting to 1e-5");
                a_eps[0] = 1.e-5;
            }
            for (int i = 0; i< value.number_of_inclusions; i++)
                value.eps.push_back(a_eps[0]);
        }
        else
        {
            for (int i = 0; i< value.number_of_inclusions; i++)
                value.eps.push_back(a_eps[i]);
        }
        
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
