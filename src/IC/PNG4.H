// 
// Initialize a field using a PNG image. (2D only)
//

#ifndef IC_PNG4_H
#define IC_PNG4_H
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
class PNG4 : public IC
{
public:
    enum Fit { Stretch, FitWidth, FitHeight, Coord };
    enum Channel { R = 0, G = 1, B = 2, A = 3 };

    PNG4(amrex::Vector<amrex::Geometry>& _geom) : IC(_geom) {}

    PNG4(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name) : PNG4(_geom)
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
        
        // Flip the row pointers to invert the image vertically (newly added since image is getting inverted)
        for (int y = 0; y < png_height / 2; y++) {
            png_bytep temp = row_pointers[y];
            row_pointers[y] = row_pointers[png_height - 1 - y];
            row_pointers[png_height - 1 - y] = temp;
        }
        
        // Initialize min and max for each channel with extreme values
        for (int c = 0; c < 4; ++c) {
            min[c] = 255.0;
            max[c] = 0.0;
        }

        // First pass: Calculate actual min/max values for each channel (R, G, B, A)
        for (int y = 0; y < png_height; y++) {
            png_bytep row = row_pointers[y];
            for (int x = 0; x < png_width; x++) {
                png_bytep px = &(row[x * 4]);  // Assuming 4 bytes per pixel (RGBA)

                // Update min/max for each channel
                min[0] = std::min(min[0], (Set::Scalar)px[Channel::R]);
                max[0] = std::max(max[0], (Set::Scalar)px[Channel::R]);
                min[1] = std::min(min[1], (Set::Scalar)px[Channel::G]);
                max[1] = std::max(max[1], (Set::Scalar)px[Channel::G]);
                min[2] = std::min(min[2], (Set::Scalar)px[Channel::B]);
                max[2] = std::max(max[2], (Set::Scalar)px[Channel::B]);
                min[3] = std::min(min[3], (Set::Scalar)px[Channel::A]);
                max[3] = std::max(max[3], (Set::Scalar)px[Channel::A]);
            }
        }  
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

                if (I >= 0 && I < png_width && J >= 0 && J < png_height)
                {
                    png_bytep px = &(row_pointers[J][I * 4]);

                    // Assign quaternion components
                    field(i, j, k, 0) = ((Set::Scalar)(px[channel1]) - min[0]) / (max[0] - min[0]);
                    field(i, j, k, 1) = ((Set::Scalar)(px[channel2]) - min[1]) / (max[1] - min[1]);
                    field(i, j, k, 2) = ((Set::Scalar)(px[channel3]) - min[2]) / (max[2] - min[2]);
                    field(i, j, k, 3) = ((Set::Scalar)(px[channel4]) - min[3]) / (max[3] - min[3]);

                    // Normalize the quaternion to maintain unit length
                    Set::Scalar q0 = field(i, j, k, 0);
                    Set::Scalar q1 = field(i, j, k, 1);
                    Set::Scalar q2 = field(i, j, k, 2);
                    Set::Scalar q3 = field(i, j, k, 3);

                    Set::Scalar norm = std::sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);

                    if (norm > 0) {
                        field(i, j, k, 0) = q0 / norm;
                        field(i, j, k, 1) = q1 / norm;
                        field(i, j, k, 2) = q2 / norm;
                        field(i, j, k, 3) = q3 / norm;
                    }
                }
                else
                {
                    field(i, j, k, 0) = 0.0;
                    field(i, j, k, 1) = 0.0;
                    field(i, j, k, 2) = 0.0;
                    field(i, j, k, 3) = 0.0;
                }
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

    Fit fit = Fit::Stretch;
    Channel channel1 = Channel::R;
    Channel channel2 = Channel::G;
    Channel channel3 = Channel::B;
    Channel channel4 = Channel::A;
    
    Set::Scalar min[4] = {0.0, 0.0, 0.0, 0.0};
    Set::Scalar max[4] = {255.0, 255.0, 255.0, 255.0};
#endif

public:
    static void Parse(PNG4& value, IO::ParmParse& pp)
    {
#ifndef ALAMO_NOPNG
        std::string filename;
        pp_query_file("filename", filename);
        value.Define(filename);

        std::string fit = "stretch";
        pp_query("fit", fit);
        if (fit == "stretch") value.fit = Fit::Stretch;
        else if (fit == "fitheight") value.fit = Fit::FitHeight;
        else if (fit == "fitwidth") value.fit = Fit::FitWidth;
        else if (fit == "coord")
        {
            value.fit = Fit::Coord;
            pp_queryarr("coord.lo", value.coord_lo);
            pp_queryarr("coord.hi", value.coord_hi);
        }
        else Util::Abort(INFO, "Invalid value for png fit");

        value.channel1 = Channel::R;
        value.channel2 = Channel::G;
        value.channel3 = Channel::B;
        value.channel4 = Channel::A;

        pp_query("minR", value.min[0]);
        pp_query("maxR", value.max[0]);
        pp_query("minG", value.min[1]);
        pp_query("maxG", value.max[1]);
        pp_query("minB", value.min[2]);
        pp_query("maxB", value.max[2]);
        pp_query("minA", value.min[3]);
        pp_query("maxA", value.max[3]);

#else
        Util::Abort(INFO,"PNG is disabled");
#endif
    }
};
}
#endif

