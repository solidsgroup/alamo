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


        for (amrex::MFIter mfi(*a_field[lev],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            bx.grow(a_field[lev]->nGrow());
            amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
            amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {

                Set::Scalar x, y;

                if (fit == Fit::Stretch)
                {
                    x = (Set::Scalar)i / (Set::Scalar)(domain.hiVect()[0]);
                    y = (Set::Scalar)j / (Set::Scalar)(domain.hiVect()[1]);
                }
                else if (fit == Fit::FitWidth)
                {
                    Set::Scalar aspect_ratio = img_width / img_height;
                    x = (Set::Scalar)i / (Set::Scalar)(domain.hiVect()[0]);
                    y = (Set::Scalar)j / (Set::Scalar)(domain.hiVect()[1]);
                    y -= 0.5 - 0.5 / aspect_ratio;
                    y *= aspect_ratio;
                }
                else if (fit == Fit::FitHeight)
                {
                    Set::Scalar aspect_ratio = img_height / img_width;
                    x = (Set::Scalar)i / (Set::Scalar)(domain.hiVect()[0]);
                    y = (Set::Scalar)j / (Set::Scalar)(domain.hiVect()[1]);
                    x -= 0.5 - 0.5 / aspect_ratio;
                    x *= aspect_ratio;
                }

                x = std::min(x,1.0); y = std::min(y,1.0);
                x = std::max(x,0.0); y = std::max(y,0.0);
                
                Set::Scalar img_x = (img_width)  * x;
                Set::Scalar img_y = (img_height) * y;

                int I = (int)(img_x);
                int J = (int)(img_y);

                Set::Scalar x1 = I*img_dx, x2 = (I+1)*img_dx;
                Set::Scalar y1 = J*img_dy, y2 = (J+1)*img_dy;

                if (I < bmp.nx-1 && J < bmp.ny - 1)
                {
                    Set::Scalar fQ11 = ((Set::Scalar)(bmp(I,  J  )[channel]) - min) / (max - min);
                    Set::Scalar fQ12 = ((Set::Scalar)(bmp(I,  J+1)[channel]) - min) / (max - min);
                    Set::Scalar fQ21 = ((Set::Scalar)(bmp(I+1,J  )[channel]) - min) / (max - min);
                    Set::Scalar fQ22 = ((Set::Scalar)(bmp(I+1,J+1)[channel]) - min) / (max - min);

                    field(i,j,k) =
                        (fQ11*(x2-img_x)*(y2-img_y) + fQ21*(img_x-x1)*(y2-img_y) + fQ12*(x2-img_x)*(img_y-y1) + fQ22*(img_x-x1)*(img_y-y1)) / (img_dx * img_dy);
                }
                else if (I == bmp.nx-1 && J <  bmp.ny - 1)
                {
                    Set::Scalar fQ11 = ((Set::Scalar)(bmp(I  ,J  )[channel]) - min) / (max - min);
                    Set::Scalar fQ12 = ((Set::Scalar)(bmp(I  ,J+1)[channel]) - min) / (max - min);
                    field(i,j,k) = fQ11 + (fQ12-fQ11) * (img_y - y1);
                }
                else if (I <  bmp.nx-1 && J == bmp.ny - 1)
                {
                    Set::Scalar fQ11 = ((Set::Scalar)(bmp(I  ,J  )[channel]) - min) / (max - min);
                    Set::Scalar fQ21 = ((Set::Scalar)(bmp(I+1,J  )[channel]) - min) / (max - min);
                    field(i,j,k) = fQ11 + (fQ21-fQ11) * (img_x - x1);
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
