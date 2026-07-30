
#include <cstring>

#include <AMReX_GpuContainers.H>

#include "Set/Set.H"
#include "IO/ParmParse.H"
#include "Util/Util.H"
#include "Unit/Unit.H"

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
    ~PNG() = default;

    //enum Type {XYZ, XY, YZ, XZ};
    enum Fit { Stretch, FitWidth, FitHeight, Coord };
    enum Channel { R = 0, G = 1, B = 2, A = 3 };

#ifndef ALAMO_NOPNG
    struct DeviceView
    {
        const unsigned char* pixels = nullptr;
        int width = 0;
        int height = 0;
        Set::Vector coord_lo = Set::Vector::Zero();
        Set::Vector coord_hi = Set::Vector::Zero();
        Fit fit = Fit::Stretch;
        Set::Scalar min = NAN;
        Set::Scalar max = NAN;
        int rotation = 0;
        bool flip_x = false;
        bool flip_y = false;
        Set::Vector domain_lo = Set::Vector::Zero();
        Set::Vector domain_hi = Set::Vector::Zero();

        AMREX_GPU_HOST_DEVICE
        unsigned char Pixel(int i, int j, int channel) const noexcept
        {
            return pixels[4 * (j * width + i) + channel];
        }

        AMREX_GPU_HOST_DEVICE
        std::array<Set::Scalar, 4> operator()(
            const Set::Vector& x) const noexcept
        {
            std::array<Set::Scalar, 4> ret{};
            const Set::Scalar img_width =
                static_cast<Set::Scalar>(width - 1);
            const Set::Scalar img_height =
                static_cast<Set::Scalar>(height - 1);
            const bool sideways = rotation == 90 || rotation == 270;
            const Set::Scalar displayed_width =
                sideways ? img_height : img_width;
            const Set::Scalar displayed_height =
                sideways ? img_width : img_height;

            Set::Vector ximg;
            if (fit == Fit::Stretch)
            {
                ximg(0) =
                    (x(0) - domain_lo(0)) / (domain_hi(0) - domain_lo(0));
                ximg(1) =
                    (x(1) - domain_lo(1)) / (domain_hi(1) - domain_lo(1));
            }
            else if (fit == Fit::FitWidth)
            {
                const Set::Scalar aspect_ratio =
                    displayed_width / displayed_height;
                ximg(0) =
                    (x(0) - domain_lo(0)) / (domain_hi(0) - domain_lo(0));
                ximg(1) =
                    (x(1) - domain_lo(1)) / (domain_hi(1) - domain_lo(1));
                ximg(1) -= 0.5 - 0.5 / aspect_ratio;
                ximg(1) *= aspect_ratio;
            }
            else if (fit == Fit::FitHeight)
            {
                const Set::Scalar aspect_ratio =
                    displayed_height / displayed_width;
                ximg(0) =
                    (x(0) - domain_lo(0)) / (domain_hi(0) - domain_lo(0));
                ximg(1) =
                    (x(1) - domain_lo(1)) / (domain_hi(1) - domain_lo(1));
                ximg(0) -= 0.5 - 0.5 / aspect_ratio;
                ximg(0) *= aspect_ratio;
            }
            else
            {
                ximg(0) = (x(0) - coord_lo(0)) / (coord_hi(0) - coord_lo(0));
                ximg(1) = (x(1) - coord_lo(1)) / (coord_hi(1) - coord_lo(1));
            }

            if (flip_x) ximg(0) = 1.0 - ximg(0);
            if (flip_y) ximg(1) = 1.0 - ximg(1);

            const Set::Vector unrotated = ximg;
            if (rotation == 90)
            {
                ximg(0) = unrotated(1);
                ximg(1) = 1.0 - unrotated(0);
            }
            else if (rotation == 180)
            {
                ximg(0) = 1.0 - unrotated(0);
                ximg(1) = 1.0 - unrotated(1);
            }
            else if (rotation == 270)
            {
                ximg(0) = 1.0 - unrotated(1);
                ximg(1) = unrotated(0);
            }

            ximg(0) = ximg(0) < 0.0 ? 0.0 : (ximg(0) > 1.0 ? 1.0 : ximg(0));
            ximg(1) = ximg(1) < 0.0 ? 0.0 : (ximg(1) > 1.0 ? 1.0 : ximg(1));
            ximg(0) *= img_width;
            ximg(1) *= img_height;

            const int I = static_cast<int>(ximg(0));
            const int J = static_cast<int>(ximg(1));
            const Set::Scalar x1 = static_cast<Set::Scalar>(I);
            const Set::Scalar x2 = static_cast<Set::Scalar>(I + 1);
            const Set::Scalar y1 = static_cast<Set::Scalar>(J);
            const Set::Scalar y2 = static_cast<Set::Scalar>(J + 1);

            for (int channel = 0; channel < 4; ++channel)
            {
                if (I > 0 && I < width - 1 && J > 0 && J < height - 1)
                {
                    const Set::Scalar fQ11 =
                        (static_cast<Set::Scalar>(Pixel(I, J, channel)) - min) /
                        (max - min);
                    const Set::Scalar fQ12 =
                        (static_cast<Set::Scalar>(Pixel(I, J + 1, channel)) -
                         min) /
                        (max - min);
                    const Set::Scalar fQ21 =
                        (static_cast<Set::Scalar>(Pixel(I + 1, J, channel)) -
                         min) /
                        (max - min);
                    const Set::Scalar fQ22 =
                        (static_cast<Set::Scalar>(
                             Pixel(I + 1, J + 1, channel)) -
                         min) /
                        (max - min);
                    ret[channel] =
                        fQ11 * (x2 - ximg(0)) * (y2 - ximg(1)) +
                        fQ21 * (ximg(0) - x1) * (y2 - ximg(1)) +
                        fQ12 * (x2 - ximg(0)) * (ximg(1) - y1) +
                        fQ22 * (ximg(0) - x1) * (ximg(1) - y1);
                }
                else if ((I == 0 || I == width - 1) && J < height - 1)
                {
                    const Set::Scalar fQ11 =
                        (static_cast<Set::Scalar>(Pixel(I, J, channel)) - min) /
                        (max - min);
                    const Set::Scalar fQ12 =
                        (static_cast<Set::Scalar>(Pixel(I, J + 1, channel)) -
                         min) /
                        (max - min);
                    ret[channel] = fQ11 + (fQ12 - fQ11) * (ximg(1) - y1);
                }
                else if (I < width - 1 && (J == 0 || J == height - 1))
                {
                    const Set::Scalar fQ11 =
                        (static_cast<Set::Scalar>(Pixel(I, J, channel)) - min) /
                        (max - min);
                    const Set::Scalar fQ21 =
                        (static_cast<Set::Scalar>(Pixel(I + 1, J, channel)) -
                         min) /
                        (max - min);
                    ret[channel] = fQ11 + (fQ21 - fQ11) * (ximg(0) - x1);
                }
                else
                {
                    ret[channel] =
                        (static_cast<Set::Scalar>(Pixel(I, J, channel)) - min) /
                        (max - min);
                }
            }
            return ret;
        }
    };
#endif

    static void Parse(PNG &value, IO::ParmParse &pp)
    {
#ifndef ALAMO_NOPNG

        std::string filename;
        pp.query_file("filename", filename); // BMP filename.

        // how to position the image
        pp.query_switch("fit", {
                {"stretch", [&]() {
                    value.fit = Fit::Stretch;
                }},
                {"fitheight", [&]() {
                    value.fit = Fit::FitHeight;
                }},
                {"fitwidth", [&]() {
                    value.fit = Fit::FitWidth;
                }},
                {"coord", [&]() {
                    value.fit = Fit::Coord;
                    // Lower-left coordinates of image in domain
                    pp.queryarr_required("coord.lo", value.coord_lo,Unit::Length());
                    // Upper-right coordinates of image in domain
                    pp.queryarr_required("coord.hi", value.coord_hi,Unit::Length());
                    for (int d = 0; d < 2; ++d)
                        if (!(value.coord_hi(d) > value.coord_lo(d)))
                            Util::Exception(
                                INFO,
                                "PNG coord.hi must be greater than coord.lo "
                                "in each direction; use flip.x or flip.y "
                                "to reverse an image");
                }}
            });

        // Desired minimum value to scale pixels by
        pp.query_default("min", value.min, 0.0  );
        // Desired maximum value to scale pixels by
        pp.query_default("max", value.max, 255.0);

        // Clockwise rotation applied while sampling the image.
        pp.query_default("rotation", value.rotation, 0);
        if (value.rotation != 0 && value.rotation != 90 &&
            value.rotation != 180 && value.rotation != 270)
            Util::Exception(
                INFO, "PNG rotation must be one of 0, 90, 180, or 270");
        pp.query_default("flip.x", value.flip_x, false);
        pp.query_default("flip.y", value.flip_y, false);

        if (IO::ParmParse::InTraversalMode()) return;

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

        png_bytep* row_pointers =
            (png_bytep*)malloc(sizeof(png_bytep) * value.png_height);
        for (int y = 0; y < value.png_height; y++) {
            row_pointers[y] = (png_byte*)malloc(png_get_rowbytes(png, info));
        }

        png_read_image(png, row_pointers);

        value.pixels.resize(
            static_cast<std::size_t>(value.png_width) * value.png_height * 4);
        for (int y = 0; y < value.png_height; ++y)
        {
            std::memcpy(
                value.pixels.data() +
                    static_cast<std::size_t>(y) * value.png_width * 4,
                row_pointers[y],
                static_cast<std::size_t>(value.png_width) * 4);
            free(row_pointers[y]);
        }
        free(row_pointers);

        fclose(fp);

        png_destroy_read_struct(&png, &info, NULL);
#else
        Util::IgnoreUnused(value,pp);
        Util::Abort(INFO,"PNG is disabled");
#endif

    }
    void setDomain(Set::Vector &_lo, Set::Vector &_hi)
    {
#ifndef ALAMO_NOPNG
        domain_lo = _lo;
        domain_hi = _hi;
        domain_defined = false;
#else
        Util::IgnoreUnused(_lo,_hi);
#endif
    }

    [[nodiscard]]
    std::array<int, 4> operator() (int I, int J)
    {
#ifndef ALAMO_NOPNG
        const unsigned char* px =
            pixels.data() + 4 * (J * png_width + I);
        return {px[0], px[1], px[2], px[3]};
#else
        Util::IgnoreUnused(I,J);
#endif
    }

    [[nodiscard]]
    std::array<Set::Scalar, 4> operator()(Set::Vector &x)
    {
#ifndef ALAMO_NOPNG
        Util::Assert(
            INFO, TEST(!pixels.empty()), "Running IC without initialization...");
        return GetDeviceView()(x);
#else
        Util::IgnoreUnused(x);
        return {NAN,NAN,NAN,NAN};
#endif
    }

#ifndef ALAMO_NOPNG
    DeviceView GetDeviceView() const noexcept
    {
        DeviceView view;
        view.pixels = pixels.data();
        view.width = png_width;
        view.height = png_height;
        view.coord_lo = coord_lo;
        view.coord_hi = coord_hi;
        view.fit = fit;
        view.min = min;
        view.max = max;
        view.rotation = rotation;
        view.flip_x = flip_x;
        view.flip_y = flip_y;
        view.domain_lo = domain_lo;
        view.domain_hi = domain_hi;
        return view;
    }
#endif

private:

#ifndef ALAMO_NOPNG
    int png_width = 0;
    int png_height = 0;
    png_byte color_type = 0;
    png_byte bit_depth = 0;
    amrex::Gpu::ManagedDeviceVector<unsigned char> pixels;
    Set::Vector coord_lo = Set::Vector::Zero();
    Set::Vector coord_hi = Set::Vector::Zero();

    //Util::BMP bmp;
    Fit fit = Fit::Stretch;
    Set::Scalar min = NAN, max = NAN;
    int rotation = 0;
    bool flip_x = false;
    bool flip_y = false;

    bool domain_defined = false;
    Set::Vector domain_lo = Set::Vector::Zero();
    Set::Vector domain_hi = Set::Vector::Zero();
#endif

};

}
