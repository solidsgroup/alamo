
#include <string.h>

#include "Set/Set.H"
#include "IO/ParmParse.H"

#ifndef ALAMO_NOPNG
#include <png.h>
#endif

namespace Util
{
class PNG
{
public:
    PNG() {} 
    PNG(IO::ParmParse &pp,std::string name) 
    {pp.queryclass(name, *this);}
    ~PNG() {} 

    //enum Type {XYZ, XY, YZ, XZ};
    enum Fit { Stretch, FitWidth, FitHeight, Coord };
    enum Channel { R = 0, G = 1, B = 2, A = 3 };

    static void Parse(PNG &value, IO::ParmParse &pp)
    {
#ifndef ALAMO_NOPNG

        std::string filename;
        pp_query_file("filename", filename); // BMP filename.

        std::string fit = "stretch";
        // how to position the image
        pp_query_validate("fit", fit, {"stretch","fitheight","fitwidth","coord"}); 
        if (fit == "stretch")        value.fit = Fit::Stretch;
        else if (fit == "fitheight") value.fit = Fit::FitHeight;
        else if (fit == "fitwidth")  value.fit = Fit::FitWidth;
        else if (fit == "coord")
        {
            value.fit = Fit::Coord;
            pp_queryarr("coord.lo", value.coord_lo); // Lower-left coordinates of image in domain
            pp_queryarr("coord.hi", value.coord_hi); // Upper-right coordinates of image in domain
        }
        else Util::Abort(INFO);

        // Desired minimum value to scale pixels by
        pp_query_default("min", value.min, 0.0  );
        // Desired maximum value to scale pixels by
        pp_query_default("max", value.max, 255.0);

        //
        //
        // Now, read in the PNG
        //
        //


        FILE* fp = std::fopen(filename.c_str(), "rb");

        if (fp == NULL) Util::Abort(INFO, "Cannot find file", filename);

        png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
        if (!png) Util::Abort(INFO);

        png_infop info = png_create_info_struct(png);
        if (!info) Util::Abort(INFO);

        if (setjmp(png_jmpbuf(png))) Util::Abort(INFO);


        png_init_io(png, fp);

        png_read_info(png, info);

        value.png_width = png_get_image_width(png, info);
        value.png_height = png_get_image_height(png, info);
        value.color_type = png_get_color_type(png, info);
        value.bit_depth = png_get_bit_depth(png, info);

        if (value.bit_depth == 16)
            png_set_strip_16(png);

        if (value.color_type == PNG_COLOR_TYPE_PALETTE)
            png_set_palette_to_rgb(png);

        if (value.color_type == PNG_COLOR_TYPE_GRAY && value.bit_depth < 8)
            png_set_expand_gray_1_2_4_to_8(png);

        if (png_get_valid(png, info, PNG_INFO_tRNS))
            png_set_tRNS_to_alpha(png);

        if (value.color_type == PNG_COLOR_TYPE_RGB ||
            value.color_type == PNG_COLOR_TYPE_GRAY ||
            value.color_type == PNG_COLOR_TYPE_PALETTE)
            png_set_filler(png, 0xFF, PNG_FILLER_AFTER);

        if (value.color_type == PNG_COLOR_TYPE_GRAY ||
            value.color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
            png_set_gray_to_rgb(png);

        png_read_update_info(png, info);

        if (value.row_pointers) Util::Abort(INFO);

        value.row_pointers = (png_bytep*)malloc(sizeof(png_bytep) * value.png_height);
        for (int y = 0; y < value.png_height; y++) {
            value.row_pointers[y] = (png_byte*)malloc(png_get_rowbytes(png, info));
        }

        png_read_image(png, value.row_pointers);

        fclose(fp);

        png_destroy_read_struct(&png, &info, NULL);
#else
        Util::Abort(INFO,"PNG is disabled");
#endif

    }
    void setDomain(Set::Vector &_lo, Set::Vector &_hi)
    {
#ifndef ALAMO_NOPNG
        domain_lo = _lo;
        domain_hi = _hi;
        domain_defined = false;
#endif
    }

    [[nodiscard]]
    std::array<int, 4> operator() (int I, int J)
    {
#ifndef ALAMO_NOPNG
        std::array<int, 4> ret;
        png_bytep px = &row_pointers[J][I * 4];
        ret = { px[0], px[1], px[2], px[3] };
        return ret;
#endif
    }

    [[nodiscard]]
    std::array<Set::Scalar, 4> operator()(Set::Vector &x)
    {
#ifndef ALAMO_NOPNG
        Util::Assert(INFO, TEST(row_pointers), "Running IC without initialization...");

        std::array<Set::Scalar, 4> ret;

        Set::Scalar img_width = (Set::Scalar)(png_width - 1);
        Set::Scalar img_height = (Set::Scalar)(png_height - 1);
        Set::Scalar img_dx = 1.0;
        Set::Scalar img_dy = 1.0;
        
        //Set::Scalar x, y;
        Set::Vector ximg;

        if (fit == Fit::Stretch)
        {
            ximg(0) = (x(0) - domain_lo(0)) / (domain_hi(0) - domain_lo(0));
            ximg(1) = (x(1) - domain_lo(1)) / (domain_hi(1) - domain_lo(1));
        }
        else if (fit == Fit::FitWidth)
        {
            Set::Scalar aspect_ratio = img_width / img_height;
            ximg(0) = (x(0) - domain_lo(0)) / (domain_hi(0) - domain_lo(0));
            ximg(1) = (x(1) - domain_lo(1)) / (domain_hi(1) - domain_lo(1));
            ximg(1) -= 0.5 - 0.5 / aspect_ratio;
            ximg(1) *= aspect_ratio;
        }
        else if (fit == Fit::FitHeight)
        {
            Set::Scalar aspect_ratio = img_height / img_width;
            ximg(0) = (x(0) - domain_lo(0)) / (domain_hi(0) - domain_lo(0));
            ximg(1) = (x(1) - domain_lo(1)) / (domain_hi(1) - domain_lo(1));
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


        for (int channel = 0; channel < 4; channel++)
        {

            if (I > 0 && I < png_width - 1 &&
                J>0 && J < png_height - 1)
            {
                png_bytep px_sw = &(row_pointers[J][I * 4]);
                png_bytep px_se = &(row_pointers[J][(I + 1) * 4]);
                png_bytep px_nw = &(row_pointers[J + 1][I * 4]);
                png_bytep px_ne = &(row_pointers[J + 1][(I + 1) * 4]);


                Set::Scalar fQ11 = ((Set::Scalar)(px_sw[channel]) - min) / (max - min);
                Set::Scalar fQ12 = ((Set::Scalar)(px_nw[channel]) - min) / (max - min);
                Set::Scalar fQ21 = ((Set::Scalar)(px_se[channel]) - min) / (max - min);
                Set::Scalar fQ22 = ((Set::Scalar)(px_ne[channel]) - min) / (max - min);

                ret[channel] = (
                    fQ11 * (x2 - ximg(0)) * (y2 - ximg(1)) +
                    fQ21 * (ximg(0) - x1) * (y2 - ximg(1)) +
                    fQ12 * (x2 - ximg(0)) * (ximg(1) - y1) +
                    fQ22 * (ximg(0) - x1) * (ximg(1) - y1)) / (img_dx * img_dy);

            }
            else if ((I == 0 || I == png_width - 1) && J < png_height - 1)
            {
                png_bytep px_sw = &(row_pointers[J][I * 4]);
                png_bytep px_nw = &(row_pointers[J + 1][I * 4]);

                Set::Scalar fQ11 = ((Set::Scalar)(px_sw[channel]) - min) / (max - min);
                Set::Scalar fQ12 = ((Set::Scalar)(px_nw[channel]) - min) / (max - min);
                ret[channel] = fQ11 + (fQ12 - fQ11) * (ximg(1) - y1);
            }
            else if (I < png_width - 1 && (J == 0 || J == png_height - 1))
            {
                png_bytep px_sw = &(row_pointers[J][I * 4]);
                png_bytep px_se = &(row_pointers[J][(I + 1) * 4]);

                Set::Scalar fQ11 = ((Set::Scalar)(px_sw[channel]) - min) / (max - min);
                Set::Scalar fQ21 = ((Set::Scalar)(px_se[channel]) - min) / (max - min);
                ret[channel] = fQ11 + (fQ21 - fQ11) * (ximg(0) - x1);
            }
            else if (I == png_width - 1 && J == png_height - 1)
            {
                png_bytep px_sw = &(row_pointers[J][I * 4]);

                Set::Scalar fQ11 = ((Set::Scalar)(px_sw[channel]) - min) / (max - min);
                ret[channel] = fQ11;
            }
            else
            {
                ret[channel] = 0.0;
            }
        }
        return ret;
#endif
    }

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
    Set::Scalar min = NAN, max = NAN;

    bool domain_defined = false;
    Set::Vector domain_lo = Set::Vector::Zero();
    Set::Vector domain_hi = Set::Vector::Zero();
#endif

};

}
