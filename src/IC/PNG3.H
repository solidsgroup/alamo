// 
// Initialize a field using a PNG image. (2D only)
//

#ifndef IC_PNG3_H
#define IC_PNG3_H
#include <cmath>

#ifndef ALAMO_NOPNG
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <png.h>
#endif

#include "IC/IC.H"
#include "Util/Util.H"
#include "Util/BMP.H"
#include "Set/Set.H"
#include "IO/ParmParse.H"

namespace IC
{
class PNG3 : public IC
{
public:
    //enum Type {XYZ, XY, YZ, XZ};
    enum Fit { Stretch, FitWidth, FitHeight, Coord };
    enum Channel { R = 0, G = 1, B = 2};

    PNG3(amrex::Vector<amrex::Geometry>& _geom) : IC(_geom) {}

    PNG3(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name) : PNG3(_geom)
    {
        pp_queryclass(name, *this);
    }

    void Define(std::string filename)
    {
#ifndef ALAMO_NOPNG

        FILE* fp = std::fopen(filename.c_str(), "rb");

        if (fp == NULL) Util::Abort(INFO, "Cannot find file", filename);

        png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
        if (!png) Util::Abort(INFO);

        png_infop info = png_create_info_struct(png);
        if (!info) Util::Abort(INFO);

        if (setjmp(png_jmpbuf(png))) Util::Abort(INFO);


        png_init_io(png, fp);

        png_read_info(png, info);

        png_width = png_get_image_width(png, info);
        png_height = png_get_image_height(png, info);
        color_type = png_get_color_type(png, info);
        bit_depth = png_get_bit_depth(png, info);

        if (bit_depth == 16)
            png_set_strip_16(png);

        if (color_type == PNG_COLOR_TYPE_PALETTE)
            png_set_palette_to_rgb(png);

        if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
            png_set_expand_gray_1_2_4_to_8(png);

        if (png_get_valid(png, info, PNG_INFO_tRNS))
            png_set_tRNS_to_alpha(png);

        if (color_type == PNG_COLOR_TYPE_RGB ||
            color_type == PNG_COLOR_TYPE_GRAY ||
            color_type == PNG_COLOR_TYPE_PALETTE)
            png_set_filler(png, 0xFF, PNG_FILLER_AFTER);

        if (color_type == PNG_COLOR_TYPE_GRAY ||
            color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
            png_set_gray_to_rgb(png);

        png_read_update_info(png, info);

        if (row_pointers) Util::Abort(INFO);

        row_pointers = (png_bytep*)malloc(sizeof(png_bytep) * png_height);
        for (int y = 0; y < png_height; y++) {
            row_pointers[y] = (png_byte*)malloc(png_get_rowbytes(png, info));
        }

        png_read_image(png, row_pointers);

        fclose(fp);

        png_destroy_read_struct(&png, &info, NULL);
#else
        Util::Abort(INFO,"PNG is disabled");
#endif
    }

    void Add(const int& lev, Set::Field<Set::Scalar>& a_field, Set::Scalar)
    {

#ifndef ALAMO_NOPNG
        Set::Vector DX(geom[lev].CellSize());
        amrex::Box domain = geom[lev].Domain();

        if (!row_pointers) Util::Abort(INFO, "Running IC without initialization...");

        Set::Scalar img_width = (Set::Scalar)(png_width - 1);
        Set::Scalar img_height = (Set::Scalar)(png_height - 1);
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
                    x(0) = domlo(0) + ((amrex::Real)(i)) * geom[lev].CellSize()[0];
                    x(1) = domlo(1) + ((amrex::Real)(j)) * geom[lev].CellSize()[1];
                }
                else if (type == amrex::IndexType::TheCellType())
                {
                    x(0) = domlo(0) + ((amrex::Real)(i)+0.5) * geom[lev].CellSize()[0];
                    x(1) = domlo(1) + ((amrex::Real)(j)+0.5) * geom[lev].CellSize()[1];
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
                else if (fit == Fit::Coord)
                {
                    ximg(0) = (x(0) - coord_lo(0)) / (coord_hi(0) - coord_lo(0));
                    ximg(1) = (x(1) - coord_lo(1)) / (coord_hi(1) - coord_lo(1));
                }

                ximg(0) = std::min(ximg(0), 1.0); ximg(1) = std::min(ximg(1), 1.0);
                ximg(0) = std::max(ximg(0), 0.0); ximg(1) = std::max(ximg(1), 0.0);

                ximg(0) *= img_width;
                ximg(1) *= img_height;

                int I = (int)(ximg(0));
                int J = (int)(ximg(1));

                Set::Scalar x1 = I * img_dx, x2 = (I + 1) * img_dx;
                Set::Scalar y1 = J * img_dy, y2 = (J + 1) * img_dy;


                if (I > 0 && I < png_width - 1 &&
                    J>0 && J < png_height - 1)
                {
                    png_bytep px_sw = &(row_pointers[J][I * 4]);
                    png_bytep px_se = &(row_pointers[J][(I + 1) * 4]);
                    png_bytep px_nw = &(row_pointers[J + 1][I * 4]);
                    png_bytep px_ne = &(row_pointers[J + 1][(I + 1) * 4]);


                    Set::Scalar fQ11R = ((Set::Scalar)(px_sw[channel1]) - min[0]) / (max[0] - min[0]);
                    Set::Scalar fQ12R = ((Set::Scalar)(px_nw[channel1]) - min[0]) / (max[0] - min[0]);
                    Set::Scalar fQ21R = ((Set::Scalar)(px_se[channel1]) - min[0]) / (max[0] - min[0]);
                    Set::Scalar fQ22R = ((Set::Scalar)(px_ne[channel1]) - min[0]) / (max[0] - min[0]);
                    
                    Set::Scalar fQ11G = ((Set::Scalar)(px_sw[channel2]) - min[1]) / (max[1] - min[1]);
                    Set::Scalar fQ12G = ((Set::Scalar)(px_nw[channel2]) - min[1]) / (max[1] - min[1]);
                    Set::Scalar fQ21G = ((Set::Scalar)(px_se[channel2]) - min[1]) / (max[1] - min[1]);
                    Set::Scalar fQ22G = ((Set::Scalar)(px_ne[channel2]) - min[1]) / (max[1] - min[1]);
                    
                    Set::Scalar fQ11B = ((Set::Scalar)(px_sw[channel3]) - min[2]) / (max[2] - min[2]);
                    Set::Scalar fQ12B = ((Set::Scalar)(px_nw[channel3]) - min[2]) / (max[2] - min[2]);
                    Set::Scalar fQ21B = ((Set::Scalar)(px_se[channel3]) - min[2]) / (max[2] - min[2]);
                    Set::Scalar fQ22B = ((Set::Scalar)(px_ne[channel3]) - min[2]) / (max[2] - min[2]);

                    

                    field(i, j, k,0) =
                        (fQ11R * (x2 - ximg(0)) * (y2 - ximg(1)) +
                            fQ21R * (ximg(0) - x1) * (y2 - ximg(1)) +
                            fQ12R * (x2 - ximg(0)) * (ximg(1) - y1) +
                            fQ22R * (ximg(0) - x1) * (ximg(1) - y1)) / (img_dx * img_dy);
                            
                            
                    field(i, j, k,1) =
                        (fQ11G * (x2 - ximg(0)) * (y2 - ximg(1)) +
                            fQ21G * (ximg(0) - x1) * (y2 - ximg(1)) +
                            fQ12G * (x2 - ximg(0)) * (ximg(1) - y1) +
                            fQ22G * (ximg(0) - x1) * (ximg(1) - y1)) / (img_dx * img_dy);
                            
                            
                    field(i, j, k,2) =
                        (fQ11B * (x2 - ximg(0)) * (y2 - ximg(1)) +
                            fQ21B * (ximg(0) - x1) * (y2 - ximg(1)) +
                            fQ12B * (x2 - ximg(0)) * (ximg(1) - y1) +
                            fQ22B * (ximg(0) - x1) * (ximg(1) - y1)) / (img_dx * img_dy);
                            



                }
                else if ((I == 0 || I == png_width - 1) && J < png_height - 1)
                {
                    png_bytep px_sw = &(row_pointers[J][I * 4]);
                    png_bytep px_nw = &(row_pointers[J + 1][I * 4]);

                    Set::Scalar fQ11R = ((Set::Scalar)(px_sw[channel1]) - min[0]) / (max[0] - min[0]);
                    Set::Scalar fQ12R = ((Set::Scalar)(px_nw[channel1]) - min[0]) / (max[0] - min[0]);
                    
                    
                    Set::Scalar fQ11G = ((Set::Scalar)(px_sw[channel2]) - min[1]) / (max[1] - min[1]);
                    Set::Scalar fQ12G = ((Set::Scalar)(px_nw[channel2]) - min[1]) / (max[1] - min[1]);
                    
                    
                    Set::Scalar fQ11B = ((Set::Scalar)(px_sw[channel3]) - min[2]) / (max[2] - min[2]);
                    Set::Scalar fQ12B = ((Set::Scalar)(px_nw[channel3]) - min[2]) / (max[2] - min[2]);
                    

                    
                    
                    field(i, j, k, 0) = fQ11R + (fQ12R - fQ11R) * (ximg(1) - y1);
                    field(i, j, k, 1) = fQ11G + (fQ12G - fQ11G) * (ximg(1) - y1);
                    field(i, j, k, 2) = fQ11B + (fQ12B - fQ11B) * (ximg(1) - y1);

                    
                }
                else if (I < png_width - 1 && (J == 0 || J == png_height - 1))
                {
                    png_bytep px_sw = &(row_pointers[J][I * 4]);
                    png_bytep px_se = &(row_pointers[J][(I + 1) * 4]);

                    Set::Scalar fQ11R = ((Set::Scalar)(px_sw[channel1]) - min[0]) / (max[0] - min[0]);
                    Set::Scalar fQ21R = ((Set::Scalar)(px_se[channel1]) - min[0]) / (max[0] - min[0]);
                    
                    Set::Scalar fQ11G = ((Set::Scalar)(px_sw[channel2]) - min[1]) / (max[1] - min[1]);
                    Set::Scalar fQ21G = ((Set::Scalar)(px_se[channel2]) - min[1]) / (max[1] - min[1]);
                    
                    Set::Scalar fQ11B = ((Set::Scalar)(px_sw[channel3]) - min[2]) / (max[2] - min[2]);
                    Set::Scalar fQ21B = ((Set::Scalar)(px_se[channel3]) - min[2]) / (max[2] - min[2]);

                    
                                                                                
                    field(i, j, k, 0) = fQ11R + (fQ21R - fQ11R) * (ximg(0) - x1);
                    field(i, j, k, 1) = fQ11G + (fQ21G - fQ11G) * (ximg(0) - x1);
                    field(i, j, k, 2) = fQ11B + (fQ21B - fQ11B) * (ximg(0) - x1);
                                     
                }
                else if (I == png_width - 1 && J == png_height - 1)
                {
                    png_bytep px_sw = &(row_pointers[J][I * 4]);

                    Set::Scalar fQ11R = ((Set::Scalar)(px_sw[channel1]) - min[0]) / (max[0] - min[0]);
                    
                    Set::Scalar fQ11G = ((Set::Scalar)(px_sw[channel2]) - min[1]) / (max[1] - min[1]);
                    
                    Set::Scalar fQ11B = ((Set::Scalar)(px_sw[channel3]) - min[2]) / (max[2] - min[2]);
                    
       
                    field(i, j, k, 0) = fQ11R;
                    
                    field(i, j, k, 1) = fQ11G;
                    
                    field(i, j, k, 2) = fQ11B;
                    
            
                    
                }
                else
                {
                    field(i, j, k, 0) = 0.0;
                    field(i, j, k, 1) = 0.0;
                    field(i, j, k, 2) = 0.0;
   
                }

                //if (field.nComp() > 1) field(i, j, k, 1) = 1.0 - field(i, j, k, 0);

            });
        }
        a_field[lev]->FillBoundary();
#else
        Util::Abort(INFO,"PNG is disabled");
#endif
        
    };

private:
#ifndef ALAMO_NOPNG
    int png_width, png_height;
    png_byte color_type;
    png_byte bit_depth;
    png_bytep* row_pointers = NULL;
    Set::Vector coord_lo = Set::Vector::Zero();
    Set::Vector coord_hi = Set::Vector::Zero();

    //Util::BMP bmp;
    Fit fit = Fit::Stretch;
    Channel channel1 = Channel::R;
    Channel channel2 = Channel::G;
    Channel channel3 = Channel::B;
  
    
    Set::Scalar min[4] = {0.0, 0.0, 0.0};
    Set::Scalar max[4] = {255.0, 255.0, 255.0};
    //Set::Scalar min = NAN, max = NAN;
#endif

public:
    static void Parse(PNG3& value, IO::ParmParse& pp)
    {
#ifndef ALAMO_NOPNG
        std::string filename;
        pp_query_file("filename", filename); // BMP filename.
        value.Define(filename);
        //value.bmp.Define(filename);

        std::string fit = "stretch";
        pp_query("fit", fit); // How to fit. {options: stretch, fitheight, fitwidth}
        if (fit == "stretch")        value.fit = Fit::Stretch;
        else if (fit == "fitheight") value.fit = Fit::FitHeight;
        else if (fit == "fitwidth")  value.fit = Fit::FitWidth;
        else if (fit == "coord")
        {
            value.fit = Fit::Coord;
            pp_queryarr("coord.lo", value.coord_lo); // Lower-left coordinates of image in domain
            pp_queryarr("coord.hi", value.coord_hi); // Upper-right coordinates of image in domain
        }
        else Util::Abort(INFO, "Invalid value for png fit - should be stretch/fitheight/fitwidth/coord but received '", fit, "'");
        
        value.channel1 = Channel::R;
        value.channel2 = Channel::G;
        value.channel3 = Channel::B;

        
        value.min[0] = 0.0;
        value.max[0] = 255.0;

        value.min[1] = 0.0;
        value.max[1] = 255.0;

        value.min[2] = 0.0;
        value.max[2] = 255.0;;
        

        
        pp_query("minR",value.min[0], 0.0 );
        pp_query("maxR",value.max[0], 255.0);
        pp_query("minG",value.min[1], 0.0 );
        pp_query("maxG",value.max[1], 255.0);
        pp_query("minB",value.min[2], 0.0 );
        pp_query("maxB",value.max[2], 255.0);        

#else
        Util::Abort(INFO,"PNG is disabled");
#endif
    }
};
}
#endif
