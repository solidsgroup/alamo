#ifndef IC_BMP_H
#define IC_BMP_H
#include <cmath>

#include "IC/IC.H"
#include "Util/Util.H"
#include "Util/BMP.H"
#include "Set/Set.H"
#include "IO/ParmParse.H"

namespace IC
{
class BMP : public IC
{
public:
    //enum Type {XYZ, XY, YZ, XZ};
    enum Fit {Stretch,FitWidth,FitHeight};
    enum Channel {R=0, G=1, B=2};

    BMP (amrex::Vector<amrex::Geometry> &_geom) : IC(_geom)  {}

    BMP (amrex::Vector<amrex::Geometry> &_geom, IO::ParmParse &pp, std::string name) : BMP(_geom)
    {
        pp.queryclass(name,*this);
    }

    void Define(std::string bmpfilename)
    {
        bmp.Define(bmpfilename);//"Interface_Blur2.bmp");
    }
  
    void Add(const int &lev, Set::Field<Set::Scalar> &a_field)
    {        
        
        Set::Vector DX(geom[lev].CellSize());
        amrex::Box domain = geom[lev].Domain();

        //Set::Scalar width  = geom[lev].ProbHi()[0] - geom[lev].ProbLo()[0];
        //Set::Scalar height = geom[lev].ProbHi()[1] - geom[lev].ProbLo()[1];

        Set::Scalar img_width = (Set::Scalar)(bmp.nx-1);
        Set::Scalar img_height = (Set::Scalar)(bmp.ny-1);
        Set::Scalar img_dx = 1.0;
        Set::Scalar img_dy = 1.0;

        amrex::IndexType type = a_field[lev]->ixType();
        domain.convert(type);

        Set::Vector domlo(AMREX_D_DECL(geom[lev].ProbLo()[0],geom[lev].ProbLo()[1],0.0));
        Set::Vector domhi(AMREX_D_DECL(geom[lev].ProbHi()[0],geom[lev].ProbHi()[1],0.0));


        for (amrex::MFIter mfi(*a_field[lev],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx;
            if (type == amrex::IndexType::TheNodeType()) bx = mfi.grownnodaltilebox();
            if (type == amrex::IndexType::TheCellType()) bx = mfi.growntilebox();
            bx = bx & domain;


            amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
            amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) 
            {
                Set::Vector x = Set::Vector::Zero();
                // NODE
                if (type == amrex::IndexType::TheNodeType())
                {
                    x(0) = domlo(0) + ((amrex::Real)(i)) * geom[lev].CellSize()[0];
                    x(1) = domlo(1) + ((amrex::Real)(j)) * geom[lev].CellSize()[1];
                }
                else if (type == amrex::IndexType::TheCellType())
                {
                    x(0) = domlo(0) + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];
                    x(1) = domlo(1) + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];
                }

                //Set::Scalar x, y;
                Set::Vector ximg;

                if (fit == Fit::Stretch)
                {
                    ximg(0) = (x(0) - domlo(0)) / (domhi(0) - domlo(0)); 
                    ximg(1) = (x(1) - domlo(1)) / (domhi(1) - domlo(1)); 
                }
                else if (fit == Fit::FitWidth)
                {
                    Set::Scalar aspect_ratio = img_width / img_height;
                    ximg(0) = (x(0) - domlo(0)) / (domhi(0) - domlo(0)); 
                    ximg(1) = (x(1) - domlo(1)) / (domhi(1) - domlo(1)); 
                    ximg(1) -= 0.5 - 0.5 / aspect_ratio;
                    ximg(1) *= aspect_ratio;
                }
                else if (fit == Fit::FitHeight)
                {
                    Set::Scalar aspect_ratio = img_height / img_width;
                    ximg(0) = (x(0) - domlo(0)) / (domhi(0) - domlo(0)); 
                    ximg(1) = (x(1) - domlo(1)) / (domhi(1) - domlo(1)); 
                    ximg(0) -= 0.5 - 0.5 / aspect_ratio;
                    ximg(0) *= aspect_ratio;
                }

                ximg(0) = std::min(ximg(0),1.0); ximg(1) = std::min(ximg(1),1.0);
                ximg(0) = std::max(ximg(0),0.0); ximg(1) = std::max(ximg(1),0.0);
                
                ximg(0) *= img_width;
                ximg(1) *= img_height;

                int I = (int)(ximg(0));
                int J = (int)(ximg(1));

                Set::Scalar x1 = I*img_dx, x2 = (I+1)*img_dx;
                Set::Scalar y1 = J*img_dy, y2 = (J+1)*img_dy;

                if (I < bmp.nx-1 && J < bmp.ny - 1)
                {
                    Set::Scalar fQ11 = ((Set::Scalar)(bmp(I,  J  )[channel]) - min) / (max - min);
                    Set::Scalar fQ12 = ((Set::Scalar)(bmp(I,  J+1)[channel]) - min) / (max - min);
                    Set::Scalar fQ21 = ((Set::Scalar)(bmp(I+1,J  )[channel]) - min) / (max - min);
                    Set::Scalar fQ22 = ((Set::Scalar)(bmp(I+1,J+1)[channel]) - min) / (max - min);

                    field(i,j,k) =
                        (fQ11*(x2-ximg(0))*(y2-ximg(1)) + fQ21*(ximg(0)-x1)*(y2-ximg(1)) + fQ12*(x2-ximg(0))*(ximg(1)-y1) + fQ22*(ximg(0)-x1)*(ximg(1)-y1)) / (img_dx * img_dy);
                }
                else if (I == bmp.nx-1 && J <  bmp.ny - 1)
                {
                    Set::Scalar fQ11 = ((Set::Scalar)(bmp(I  ,J  )[channel]) - min) / (max - min);
                    Set::Scalar fQ12 = ((Set::Scalar)(bmp(I  ,J+1)[channel]) - min) / (max - min);
                    field(i,j,k) = fQ11 + (fQ12-fQ11) * (ximg(1) - y1);
                }
                else if (I <  bmp.nx-1 && J == bmp.ny - 1)
                {
                    Set::Scalar fQ11 = ((Set::Scalar)(bmp(I  ,J  )[channel]) - min) / (max - min);
                    Set::Scalar fQ21 = ((Set::Scalar)(bmp(I+1,J  )[channel]) - min) / (max - min);
                    field(i,j,k) = fQ11 + (fQ21-fQ11) * (ximg(0) - x1);
                }
                else if (I ==  bmp.nx-1 && J == bmp.ny - 1)
                {
                    Set::Scalar fQ11 = ((Set::Scalar)(bmp(I  ,J  )[channel]) - min) / (max - min);
                    field(i,j,k) = fQ11;
                }
                else
                {
                    field(i,j,k) = 0.0;
                }

                if (field.nComp() > 1) field(i,j,k,1) = 1.0 - field(i,j,k,0);

            });
        }
        a_field[lev]->FillBoundary();
    };
  
private:
    Util::BMP bmp;
    Fit fit = Fit::Stretch;
    Channel channel = Channel::G;
    Set::Scalar min=0.0, max=255.0;

public:
    static void Parse(BMP & value, IO::ParmParse & pp)
    {
        std::string filename;
        // BMP filename.
        // Note that in GIMP, you must select "do not write color space information"
        // and "24 bit R8 G8 B8" when exporting the BMP file.
        pp.query("filename",filename);
        value.bmp.Define(filename);

        std::string fit = "stretch";
        pp.query("fit",fit);
        if (fit=="stretch")        value.fit = Fit::Stretch;
        else if (fit=="fitheight") value.fit = Fit::FitHeight;
        else if (fit=="fitwidth")  value.fit = Fit::FitWidth;
        else Util::Abort(INFO,"Invalid value for bmp fit - should be stretch/fitheight/fitwidth but received '",fit,"'");

        std::string channel = "g";
        pp.query("channel",channel);
        if      (channel=="r" || channel=="R") value.channel = Channel::R;
        else if (channel=="g" || channel=="G") value.channel = Channel::G;
        else if (channel=="b" || channel=="B") value.channel = Channel::B;
        else Util::Abort(INFO,"Invalid value for bmp channel - should be r/g/b but received '",channel,"'");

        value.min = (Set::Scalar) value.bmp.min()[value.channel];
        value.max = (Set::Scalar) value.bmp.max()[value.channel];
        pp.query("min",value.min);
        pp.query("max",value.max);
}    
};
}
#endif
