//
// Initialize a field using a bitmap image. (2D only)
//
// Note that in GIMP, you must select "do not write color space information"
// and "24 bit R8 G8 B8" when exporting the BMP file.
//

#ifndef IC_Quaternions_H
#define IC_Quaternions_H

#include <cmath>

#include "IC/IC.H"
#include "Util/Util.H"
#include "Util/BMP4.H"  // Updated to include BMP4.H instead of BMP.H
#include "Set/Set.H"
#include "IO/ParmParse.H"

namespace IC
{
class Quaternions : public IC  // Updated class name from EulerAngles to Quaternions
{
public:
    enum Fit {Stretch, FitWidth, FitHeight, Coord};
    enum Channel {R=0, G=1, B=2, A=3};  // Channels for RGBA or quaternion components

    Quaternions(amrex::Vector<amrex::Geometry> &_geom) : IC(_geom)  {}

    Quaternions(amrex::Vector<amrex::Geometry> &_geom, IO::ParmParse &pp, std::string name) : Quaternions(_geom)
    {
        pp_queryclass(name, *this);
    }

    void Define(std::string bmpfilename)
    {
        bmp.Define(bmpfilename);  // Calls the updated BMP4 Define function
    }

    void Add(const int &lev, Set::Field<Set::Scalar> &a_field, Set::Scalar)
    {        
        Util::Assert(INFO, TEST(a_field[0]->nComp() == 4));  // Ensure the field has 4 components for quaternions
        Set::Vector DX(geom[lev].CellSize());
        amrex::Box domain = geom[lev].Domain();

        Set::Scalar img_width = static_cast<Set::Scalar>(bmp.nx - 1);
        Set::Scalar img_height = static_cast<Set::Scalar>(bmp.ny - 1);
        Set::Scalar img_dx = 1.0;
        Set::Scalar img_dy = 1.0;

        amrex::IndexType type = a_field[lev]->ixType();
        domain.convert(type);

        Set::Vector domlo(AMREX_D_DECL(geom[lev].ProbLo()[0], geom[lev].ProbLo()[1], 0.0));
        Set::Vector domhi(AMREX_D_DECL(geom[lev].ProbHi()[0], geom[lev].ProbHi()[1], 0.0));

        for (amrex::MFIter mfi(*a_field[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx;
            if (type == amrex::IndexType::TheNodeType()) bx = mfi.grownnodaltilebox();
            if (type == amrex::IndexType::TheCellType()) bx = mfi.growntilebox();
            bx = bx & domain;

            amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Vector x = Set::Vector::Zero();
                // NODE
                if (type == amrex::IndexType::TheNodeType())
                {
                    x(0) = domlo(0) + static_cast<amrex::Real>(i) * geom[lev].CellSize()[0];
                    x(1) = domlo(1) + static_cast<amrex::Real>(j) * geom[lev].CellSize()[1];
                }
                else if (type == amrex::IndexType::TheCellType())
                {
                    x(0) = domlo(0) + (static_cast<amrex::Real>(i) + 0.5) * geom[lev].CellSize()[0];
                    x(1) = domlo(1) + (static_cast<amrex::Real>(j) + 0.5) * geom[lev].CellSize()[1];
                }

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
                else if (fit == Fit::Coord)
                {
                    ximg(0) = (x(0) - coord_lo(0)) / (coord_hi(0) - coord_lo(0));
                    ximg(1) = (x(1) - coord_lo(1)) / (coord_hi(1) - coord_lo(1));
                }

                ximg(0) = std::min(ximg(0), 1.0); ximg(1) = std::min(ximg(1), 1.0);
                ximg(0) = std::max(ximg(0), 0.0); ximg(1) = std::max(ximg(1), 0.0);

                ximg(0) *= img_width;
                ximg(1) *= img_height;

                int I = static_cast<int>(ximg(0));
                int J = static_cast<int>(ximg(1));
                
                Set::Scalar x1 = I*img_dx, x2 = (I+1)*img_dx;
                Set::Scalar y1 = J*img_dy, y2 = (J+1)*img_dy;

                // Handle interpolation within the image bounds
                if (I < bmp.nx - 1 && J < bmp.ny - 1)
                {
                    // Bilinear interpolation on the four corners of the cell
                    Set::Scalar fQ11R = (static_cast<Set::Scalar>(bmp(I,  J  )[channel1]) - min[0]) / (max[0] - min[0]);
                    Set::Scalar fQ12R = (static_cast<Set::Scalar>(bmp(I,  J+1)[channel1]) - min[0]) / (max[0] - min[0]);
                    Set::Scalar fQ21R = (static_cast<Set::Scalar>(bmp(I+1,J  )[channel1]) - min[0]) / (max[0] - min[0]);
                    Set::Scalar fQ22R = (static_cast<Set::Scalar>(bmp(I+1,J+1)[channel1]) - min[0]) / (max[0] - min[0]);

                    Set::Scalar fQ11G = (static_cast<Set::Scalar>(bmp(I,  J  )[channel2]) - min[1]) / (max[1] - min[1]);
                    Set::Scalar fQ12G = (static_cast<Set::Scalar>(bmp(I,  J+1)[channel2]) - min[1]) / (max[1] - min[1]);
                    Set::Scalar fQ21G = (static_cast<Set::Scalar>(bmp(I+1,J  )[channel2]) - min[1]) / (max[1] - min[1]);
                    Set::Scalar fQ22G = (static_cast<Set::Scalar>(bmp(I+1,J+1)[channel2]) - min[1]) / (max[1] - min[1]);

                    Set::Scalar fQ11B = (static_cast<Set::Scalar>(bmp(I,  J  )[channel3]) - min[2]) / (max[2] - min[2]);
                    Set::Scalar fQ12B = (static_cast<Set::Scalar>(bmp(I,  J+1)[channel3]) - min[2]) / (max[2] - min[2]);
                    Set::Scalar fQ21B = (static_cast<Set::Scalar>(bmp(I+1,J  )[channel3]) - min[2]) / (max[2] - min[2]);
                    Set::Scalar fQ22B = (static_cast<Set::Scalar>(bmp(I+1,J+1)[channel3]) - min[2]) / (max[2] - min[2]);

                    Set::Scalar fQ11A = (static_cast<Set::Scalar>(bmp(I,  J  )[channel4]) - min[3]) / (max[3] - min[3]);
                    Set::Scalar fQ12A = (static_cast<Set::Scalar>(bmp(I,  J+1)[channel4]) - min[3]) / (max[3] - min[3]);
                    Set::Scalar fQ21A = (static_cast<Set::Scalar>(bmp(I+1,J  )[channel4]) - min[3]) / (max[3] - min[3]);
                    Set::Scalar fQ22A = (static_cast<Set::Scalar>(bmp(I+1,J+1)[channel4]) - min[3]) / (max[3] - min[3]);

                    // Interpolated quaternion components (normalized)
                    Set::Scalar quatR_norm = (fQ11R*(x2-ximg(0))*(y2-ximg(1)) + fQ21R*(ximg(0)-x1)*(y2-ximg(1)) 
                                               + fQ12R*(x2-ximg(0))*(ximg(1)-y1) + fQ22R*(ximg(0)-x1)*(ximg(1)-y1)) / (img_dx * img_dy);
                    Set::Scalar quatG_norm = (fQ11G*(x2-ximg(0))*(y2-ximg(1)) + fQ21G*(ximg(0)-x1)*(y2-ximg(1)) 
                                               + fQ12G*(x2-ximg(0))*(ximg(1)-y1) + fQ22G*(ximg(0)-x1)*(ximg(1)-y1)) / (img_dx * img_dy);
                    Set::Scalar quatB_norm = (fQ11B*(x2-ximg(0))*(y2-ximg(1)) + fQ21B*(ximg(0)-x1)*(y2-ximg(1)) 
                                               + fQ12B*(x2-ximg(0))*(ximg(1)-y1) + fQ22B*(ximg(0)-x1)*(ximg(1)-y1)) / (img_dx * img_dy);
                    Set::Scalar quatA_norm = (fQ11A*(x2-ximg(0))*(y2-ximg(1)) + fQ21A*(ximg(0)-x1)*(y2-ximg(1)) 
                                               + fQ12A*(x2-ximg(0))*(ximg(1)-y1) + fQ22A*(ximg(0)-x1)*(ximg(1)-y1)) / (img_dx * img_dy);

                    // Assign quaternion components to the field
                    field(i,j,k,0) = quatR_norm;  // w
                    field(i,j,k,1) = quatG_norm;  // x
                    field(i,j,k,2) = quatB_norm;  // y
                    field(i,j,k,3) = quatA_norm;  // z
                }
                // Handle edge cases when on the right edge but not the bottom
                else if (I == bmp.nx-1 && J < bmp.ny - 1)
                {
                    Set::Scalar fQ11R = (static_cast<Set::Scalar>(bmp(I,  J  )[channel1]) - min[0]) / (max[0] - min[0]);
                    Set::Scalar fQ12R = (static_cast<Set::Scalar>(bmp(I,  J+1)[channel1]) - min[0]) / (max[0] - min[0]);
                    Set::Scalar quatR_norm = fQ11R + (fQ12R-fQ11R) * (ximg(1) - y1);

                    Set::Scalar fQ11G = (static_cast<Set::Scalar>(bmp(I,  J  )[channel2]) - min[1]) / (max[1] - min[1]);
                    Set::Scalar fQ12G = (static_cast<Set::Scalar>(bmp(I,  J+1)[channel2]) - min[1]) / (max[1] - min[1]);
                    Set::Scalar quatG_norm = fQ11G + (fQ12G-fQ11G) * (ximg(1) - y1);

                    Set::Scalar fQ11B = (static_cast<Set::Scalar>(bmp(I,  J  )[channel3]) - min[2]) / (max[2] - min[2]);
                    Set::Scalar fQ12B = (static_cast<Set::Scalar>(bmp(I,  J+1)[channel3]) - min[2]) / (max[2] - min[2]);
                    Set::Scalar quatB_norm = fQ11B + (fQ12B-fQ11B) * (ximg(1) - y1);

                    Set::Scalar fQ11A = (static_cast<Set::Scalar>(bmp(I,  J  )[channel4]) - min[3]) / (max[3] - min[3]);
                    Set::Scalar fQ12A = (static_cast<Set::Scalar>(bmp(I,  J+1)[channel4]) - min[3]) / (max[3] - min[3]);
                    Set::Scalar quatA_norm = fQ11A + (fQ12A-fQ11A) * (ximg(1) - y1);

                    // Assign quaternion components to the field
                    field(i,j,k,0) = quatR_norm;  // w
                    field(i,j,k,1) = quatG_norm;  // x
                    field(i,j,k,2) = quatB_norm;  // y
                    field(i,j,k,3) = quatA_norm;  // z
                }
                // Handle edge cases when on the bottom edge but not the right
                else if (I <  bmp.nx-1 && J == bmp.ny - 1)
                {
                    Set::Scalar fQ11R = (static_cast<Set::Scalar>(bmp(I,  J  )[channel1]) - min[0]) / (max[0] - min[0]);
                    Set::Scalar fQ21R = (static_cast<Set::Scalar>(bmp(I+1,J  )[channel1]) - min[0]) / (max[0] - min[0]);
                    Set::Scalar quatR_norm = fQ11R + (fQ21R-fQ11R) * (ximg(0) - x1);

                    Set::Scalar fQ11G = (static_cast<Set::Scalar>(bmp(I,  J  )[channel2]) - min[1]) / (max[1] - min[1]);
                    Set::Scalar fQ21G = (static_cast<Set::Scalar>(bmp(I+1,J  )[channel2]) - min[1]) / (max[1] - min[1]);
                    Set::Scalar quatG_norm = fQ11G + (fQ21G-fQ11G) * (ximg(0) - x1);

                    Set::Scalar fQ11B = (static_cast<Set::Scalar>(bmp(I,  J  )[channel3]) - min[2]) / (max[2] - min[2]);
                    Set::Scalar fQ21B = (static_cast<Set::Scalar>(bmp(I+1,J  )[channel3]) - min[2]) / (max[2] - min[2]);
                    Set::Scalar quatB_norm = fQ11B + (fQ21B-fQ11B) * (ximg(0) - x1);

                    Set::Scalar fQ11A = (static_cast<Set::Scalar>(bmp(I,  J  )[channel4]) - min[3]) / (max[3] - min[3]);
                    Set::Scalar fQ21A = (static_cast<Set::Scalar>(bmp(I+1,J  )[channel4]) - min[3]) / (max[3] - min[3]);
                    Set::Scalar quatA_norm = fQ11A + (fQ21A-fQ11A) * (ximg(0) - x1);

                    // Assign quaternion components to the field
                    field(i,j,k,0) = quatR_norm;  // w
                    field(i,j,k,1) = quatG_norm;  // x
                    field(i,j,k,2) = quatB_norm;  // y
                    field(i,j,k,3) = quatA_norm;  // z
                }
                // Handle the corner case (bottom-right)
                else if (I ==  bmp.nx-1 && J == bmp.ny - 1)
                {
                    Set::Scalar quatR_norm = (static_cast<Set::Scalar>(bmp(I, J)[channel1]) - min[0]) / (max[0] - min[0]);
                    
                    
                    Set::Scalar quatG_norm = (static_cast<Set::Scalar>(bmp(I, J)[channel2]) - min[1]) / (max[1] - min[1]);

                    
                    Set::Scalar quatB_norm = (static_cast<Set::Scalar>(bmp(I, J)[channel3]) - min[2]) / (max[2] - min[2]);

                    
                    Set::Scalar quatA_norm = (static_cast<Set::Scalar>(bmp(I, J)[channel4]) - min[3]) / (max[3] - min[3]);


                    // Assign quaternion components to the field
                    field(i,j,k,0) = quatR_norm;  // w
                    field(i,j,k,1) = quatG_norm;  // x
                    field(i,j,k,2) = quatB_norm;  // y
                    field(i,j,k,3) = quatA_norm;  // z
                }
                else
                {
                    // Handle cases where we are outside the bounds or invalid
                    field(i,j,k,0) = 0.0;  // Identity quaternion w
                    field(i,j,k,1) = 0.0;  // Identity quaternion x
                    field(i,j,k,2) = 0.0;  // Identity quaternion y
                    field(i,j,k,3) = 0.0;  // Identity quaternion z
                }
                
                //if (field.nComp() > 1) field(i,j,k,1) = 1.0 - field(i,j,k,0);
            });
        }
        a_field[lev]->FillBoundary();
    };

private:
    Util::BMP4 bmp;
    Fit fit = Fit:: FitWidth;
    Channel channel1 = Channel::R;
    Channel channel2 = Channel::G;
    Channel channel3 = Channel::B;
    Channel channel4 = Channel::A;  // Added channel for the 4th component
    Set::Scalar min[4] = {0.0, 0.0, 0.0, 0.0};  // Adjusted for 4 channels
    Set::Scalar max[4] = {255.0, 255.0, 255.0, 255.0};  // Adjusted for 4 channels
    Set::Vector coord_lo = Set::Vector::Zero();
    Set::Vector coord_hi = Set::Vector::Zero();

public:
    static void Parse(Quaternions & value, IO::ParmParse & pp)
    {
        std::string filename;
        pp_query_file("filename", filename);
        value.bmp.Define(filename);

        std::string fit;
        pp_query_validate("fit", fit, { "stretch", "fitheight", "fitwidth" });
        if (fit == "stretch")        value.fit = Fit::Stretch;
        else if (fit == "fitheight") value.fit = Fit::FitHeight;
        else if (fit == "fitwidth")  value.fit = Fit::FitWidth;
        else if (fit == "coord")
        {
            value.fit = Fit::Coord;
            pp_queryarr("coord.lo", value.coord_lo); // Location of lower-left corner in the domain
            pp_queryarr("coord.hi", value.coord_hi); // Location of upper-right corner in the domain
        }
        else Util::Abort(INFO, "Invalid value for bmp fit - should be stretch/fitheight/fitwidth but received '", fit, "'");

        // Use RGBA channels as quaternion components
        value.channel1 = Channel::R;
        value.channel2 = Channel::G;
        value.channel3 = Channel::B;
        value.channel4 = Channel::A;  // Added for 4th channel handling

        value.min[0] = static_cast<Set::Scalar>(value.bmp.min()[value.channel1]);
        value.max[0] = static_cast<Set::Scalar>(value.bmp.max()[value.channel1]);

        value.min[1] = static_cast<Set::Scalar>(value.bmp.min()[value.channel2]);
        value.max[1] = static_cast<Set::Scalar>(value.bmp.max()[value.channel2]);

        value.min[2] = static_cast<Set::Scalar>(value.bmp.min()[value.channel3]);
        value.max[2] = static_cast<Set::Scalar>(value.bmp.max()[value.channel3]);

        value.min[3] = static_cast<Set::Scalar>(value.bmp.min()[value.channel4]);  // Min for the 4th channel
        value.max[3] = static_cast<Set::Scalar>(value.bmp.max()[value.channel4]);  // Max for the 4th channel
        
        pp_query("minR", value.min[0]);
        pp_query("maxR", value.max[0]);
        pp_query("minG", value.min[1]);
        pp_query("maxG", value.max[1]);
        pp_query("minB", value.min[2]);
        pp_query("maxB", value.max[2]);
        pp_query("minA", value.min[3]);  // Query for the 4th channel
        pp_query("maxA", value.max[3]);  // Query for the 4th channel
    }    
};
}

#endif

