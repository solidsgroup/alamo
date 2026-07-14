//
// Initialize a field using a PNG image. (2D only)
// Nearest-neighbor sampler ONLY.
//
#ifndef IC_PNG4_H
#define IC_PNG4_H
#include <cmath>
#include <algorithm>

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
class PNG4 : public IC<Set::Scalar>
{
public:
    enum Fit { Stretch, FitWidth, FitHeight, Coord };
    enum Channel { R = 0, G = 1, B = 2, A = 3 };

    PNG4(amrex::Vector<amrex::Geometry>& _geom) : IC(_geom) {}
    PNG4(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name) : PNG4(_geom)
    {
        pp_queryclass(name, *this);
    }

    ~PNG4() {
#ifndef ALAMO_NOPNG
        if (row_pointers) {
            for (int y = 0; y < png_height; ++y) free(row_pointers[y]);
            free(row_pointers);
            row_pointers = nullptr;
        }
#endif
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

        png_width  = png_get_image_width(png, info);
        png_height = png_get_image_height(png, info);
        color_type = png_get_color_type(png, info);
        bit_depth  = png_get_bit_depth(png, info);

        if (bit_depth == 16) png_set_strip_16(png);
        if (color_type == PNG_COLOR_TYPE_PALETTE) png_set_palette_to_rgb(png);
        if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8) png_set_expand_gray_1_2_4_to_8(png);
        if (png_get_valid(png, info, PNG_INFO_tRNS)) png_set_tRNS_to_alpha(png);

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

        // Flip vertically (to match expected orientation)
        for (int y = 0; y < png_height / 2; y++) {
            png_bytep temp = row_pointers[y];
            row_pointers[y] = row_pointers[png_height - 1 - y];
            row_pointers[png_height - 1 - y] = temp;
        }

        // Compute per-channel min/max across the image
        for (int c = 0; c < 4; ++c) { min[c] = 255.0; max[c] = 0.0; }
        for (int y = 0; y < png_height; y++) {
            png_bytep row = row_pointers[y];
            for (int x = 0; x < png_width; x++) {
                png_bytep px = &(row[x * 4]);
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
        amrex::Box domain = geom[lev].Domain();
        if (!row_pointers) Util::Abort(INFO, "Running IC without initialization...");

        // Image coordinate helpers
        const Set::Scalar img_w = (Set::Scalar)(png_width  - 1);
        const Set::Scalar img_h = (Set::Scalar)(png_height - 1);

        amrex::IndexType type = a_field[lev]->ixType();
        domain.convert(type);

        // Precompute geometry scalars (safer for GPU lambdas)
        const Set::Scalar domlo_x = geom[lev].ProbLo()[0];
        const Set::Scalar domlo_y = geom[lev].ProbLo()[1];
        const Set::Scalar domhi_x = geom[lev].ProbHi()[0];
        const Set::Scalar domhi_y = geom[lev].ProbHi()[1];
        const Set::Scalar dx      = geom[lev].CellSize()[0];
        const Set::Scalar dy      = geom[lev].CellSize()[1];

        for (amrex::MFIter mfi(*a_field[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx;
            if (type == amrex::IndexType::TheNodeType()) bx = mfi.grownnodaltilebox();
            if (type == amrex::IndexType::TheCellType()) bx = mfi.growntilebox();
            bx = bx & domain;

            amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);

            // Capture everything needed by value (AMReX GPU lambdas require PODs)
            const int png_w = png_width;
            const int png_h = png_height;
            png_bytep* rows = row_pointers;

            const Fit fit_local = fit;
            const Channel ch1 = channel1, ch2 = channel2, ch3 = channel3, ch4 = channel4;

            const Set::Scalar minR = min[0], maxR = max[0];
            const Set::Scalar minG = min[1], maxG = max[1];
            const Set::Scalar minB = min[2], maxB = max[2];
            const Set::Scalar minA = min[3], maxA = max[3];

            const Set::Scalar coord_lo_x = coord_lo(0);
            const Set::Scalar coord_lo_y = coord_lo(1);
            const Set::Scalar coord_hi_x = coord_hi(0);
            const Set::Scalar coord_hi_y = coord_hi(1);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                // Physical to normalized domain coords
                Set::Scalar x = 0.0, y = 0.0;
                if (type == amrex::IndexType::TheNodeType()) {
                    x = domlo_x + ((amrex::Real)(i)) * dx;
                    y = domlo_y + ((amrex::Real)(j)) * dy;
                } else { // Cell-centered
                    x = domlo_x + ((amrex::Real)(i)+0.5) * dx;
                    y = domlo_y + ((amrex::Real)(j)+0.5) * dy;
                }

                Set::Scalar u = 0.0, v = 0.0; // normalized [0,1] before scaling to pixel space
                if (fit_local == Fit::Stretch) {
                    u = (x - domlo_x) / (domhi_x - domlo_x);
                    v = (y - domlo_y) / (domhi_y - domlo_y);
                } else if (fit_local == Fit::FitWidth) {
                    Set::Scalar aspect = img_w / img_h;
                    u = (x - domlo_x) / (domhi_x - domlo_x);
                    v = (y - domlo_y) / (domhi_y - domlo_y);
                    v -= 0.5 - 0.5 / aspect;
                    v *= aspect;
                } else if (fit_local == Fit::FitHeight) {
                    Set::Scalar aspect = img_h / img_w;
                    u = (x - domlo_x) / (domhi_x - domlo_x);
                    v = (y - domlo_y) / (domhi_y - domlo_y);
                    u -= 0.5 - 0.5 / aspect;
                    u *= aspect;
                } else { // Coord
                    u = (x - coord_lo_x) / (coord_hi_x - coord_lo_x);
                    v = (y - coord_lo_y) / (coord_hi_y - coord_lo_y);
                }

                // Clamp to [0,1], then scale to pixel space
                u = amrex::min<Set::Scalar>(amrex::max<Set::Scalar>(u, 0.0), 1.0) * img_w;
                v = amrex::min<Set::Scalar>(amrex::max<Set::Scalar>(v, 0.0), 1.0) * img_h;

                // --- Nearest neighbor sampling ---
                int Ii = static_cast<int>(u + 0.5);
                int Jj = static_cast<int>(v + 0.5);
                Ii = Ii < 0 ? 0 : (Ii > png_w - 1 ? png_w - 1 : Ii);
                Jj = Jj < 0 ? 0 : (Jj > png_h - 1 ? png_h - 1 : Jj);

                png_bytep px = &(rows[Jj][Ii * 4]);

                // Safe denominators
                const Set::Scalar dR = (maxR > minR) ? (maxR - minR) : 1.0;
                const Set::Scalar dG = (maxG > minG) ? (maxG - minG) : 1.0;
                const Set::Scalar dB = (maxB > minB) ? (maxB - minB) : 1.0;
                const Set::Scalar dA = (maxA > minA) ? (maxA - minA) : 1.0;

                Set::Scalar r01 = ((Set::Scalar)(px[ch1]) - minR) / dR;
                Set::Scalar g01 = ((Set::Scalar)(px[ch2]) - minG) / dG;
                Set::Scalar b01 = ((Set::Scalar)(px[ch3]) - minB) / dB;
                Set::Scalar a01 = ((Set::Scalar)(px[ch4]) - minA) / dA;

                // Map [0,1] -> [-1,1]
                field(i, j, k, 0) = 2 * r01 - 1;
                field(i, j, k, 1) = 2 * g01 - 1;
                field(i, j, k, 2) = 2 * b01 - 1;
                field(i, j, k, 3) = 2 * a01 - 1;
            });
        }

        a_field[lev]->FillBoundary();
#else
        Util::Abort(INFO,"PNG is disabled");
#endif
    };

private:
#ifndef ALAMO_NOPNG
    int png_width = 0, png_height = 0;
    png_byte color_type = 0;
    png_byte bit_depth  = 0;
    png_bytep* row_pointers = NULL;

    Set::Vector coord_lo = Set::Vector::Zero();
    Set::Vector coord_hi = Set::Vector::Zero();

    Fit     fit      = Fit::Stretch;
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
        if      (fit == "stretch")   value.fit = Fit::Stretch;
        else if (fit == "fitheight") value.fit = Fit::FitHeight;
        else if (fit == "fitwidth")  value.fit = Fit::FitWidth;
        else if (fit == "coord") {
            value.fit = Fit::Coord;
            pp_queryarr("coord.lo", value.coord_lo);
            pp_queryarr("coord.hi", value.coord_hi);
        } else {
            Util::Abort(INFO, "Invalid value for png fit - should be stretch/fitheight/fitwidth/coord but received '", fit, "'");
        }

        // Channel mapping (defaults RGBA)
        value.channel1 = Channel::R;
        value.channel2 = Channel::G;
        value.channel3 = Channel::B;
        value.channel4 = Channel::A;

        // Optional explicit min/max overrides (keeps your behavior)
        pp_query("minR", value.min[0]); pp_query("maxR", value.max[0]);
        pp_query("minG", value.min[1]); pp_query("maxG", value.max[1]);
        pp_query("minB", value.min[2]); pp_query("maxB", value.max[2]);
        pp_query("minA", value.min[3]); pp_query("maxA", value.max[3]);

        // NOTE: no 'interp' knob anymore; always nearest.
#else
        Util::Abort(INFO,"PNG is disabled");
#endif
    }
};
}
#endif