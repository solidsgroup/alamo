#ifndef IC_CAHNHILLARDMIXTURE_H_
#define IC_CAHNHILLARDMIXTURE_H_

#include "Set/Set.H"
#include "IC/IC.H"
#include "IO/ParmParse.H"

namespace IC
{
class CahnHillardMixture : public IC
{
public:
    CahnHillardMixture(amrex::Vector<amrex::Geometry>& _geom) : IC(_geom) {}
    CahnHillardMixture(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name) : CahnHillardMixture(_geom)
    {
        pp.queryclass(name, *this);
    }

    void Add(const int& lev, Set::Field<Set::Scalar>& a_field, Set::Field<Set::Scalar>& b_field, Set::Field<Set::Scalar>& c_field, Set::Scalar)
    {
        Set::Vector DX(geom[lev].CellSize());
        amrex::IndexType type = a_field[lev]->ixType();
        int ncomp = a_field[lev]->nComp();

        for (amrex::MFIter mfi(*a_field[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx;
            if (type == amrex::IndexType::TheNodeType()) bx = mfi.grownnodaltilebox();
            if (type == amrex::IndexType::TheCellType()) bx = mfi.growntilebox();

            amrex::Array4<Set::Scalar> const& field_a = a_field[lev]->array(mfi);
            amrex::Array4<Set::Scalar> const& field_b = b_field[lev]->array(mfi);
            amrex::Array4<Set::Scalar> const& field_c = c_field[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

                if (mode == "perturbation")
                {
                    if (j % period < period / 2.0)
                    {
                        field_a(i, j, k) = ratio_a;
                        field_b(i, j, k) = ratio_b;
                        field_c(i, j, k) = 1.0 - ratio_a - ratio_b;
                    }
                    else
                    {
                        field_a(i, j, k) = ratio_b;
                        field_b(i, j, k) = ratio_a;
                        field_c(i, j, k) = 1.0 - ratio_a - ratio_b;
                    }

                    if (j > base)
                    {
                        field_a(i, j, k) = ratio_c;
                        field_b(i, j, k) = ratio_c;
                        field_c(i, j, k) = 1.0 - ratio_c - ratio_c;
                    }
                }
                else if (mode == "seed")
                {
                    Set::Scalar x_pos;
                    Set::Scalar y_pos;
                    Set::Scalar rtmp;
                    Set::Scalar phi;
                    int N_frac = (Ny_ic + 1) / 10;

                    if (i < (Nx + 1) / 4)
                    {
                        x_pos = i * dx;
                        y_pos = j * dx;
                        rtmp = sqrt(y_pos * y_pos + x_pos * x_pos);
                        phi = 0.5 * (1.0 + tanh((rtmp - radius) / xi))
                    }
                    else if (i > (Nx + 1) / 4 && i <= 3 * (Nx + 1) / 4)
                    {
                        x_pos = (i - 2 * radius) * dx;
                        y_pos = j * dx;
                        rtmp = sqrt(y_pos * y_pos + x_pos * x_pos);
                        phi = 0.5 * (1.0 + tanh((rtmp - radius) / xi))
                    }
                    else if (i > 3 * (Nx + 1) / 4 && i <= Nx)
                    {
                        x_pos = (i - 4 * radius) * dx;
                        y_pos = j * dx;
                        rtmp = sqrt(y_pos * y_pos + x_pos * x_pos);
                        phi = 0.5 * (1.0 + tanh((rtmp - radius) / xi))
                    }

                    for (m = 0; m <= ((Ny_ic + 1) / N_frac); m++)
                    {
                        if (phi <= 0.5)
                        {
                            if (((m * N_frac) % 2 == 0) && (j >= (m * N_frac)) && (j < ((m + 1) * N_frac)))
                            {
                                field_a(i, j, k) = ratio_a;
                                field_b(i, j, k) = ratio_b;
                                field_c(i, j, k) = 1.0 - ratio_a - ratio_b;
                            }

                            if (((m * N_frac) % 2 == 1) && (j >= (m * N_frac)) && (j < ((m + 1) * N_frac)))
                            {
                                field_a(i, j, k) = ratio_b;
                                field_b(i, j, k) = ratio_a;
                                field_c(i, j, k) = 1.0 - ratio_a - ratio_b;
                            }
                        }
                        else
                        {
                            field_a(i, j, k) = ratio_c;
                            field_b(i, j, k) = ratio_c;
                            field_c(i, j, k) = 1.0 - ratio_c - ratio_c;
                        }
                    }
                }
                else Util::Abort(INFO, "Invalid Initialization mode ", mode);
            });
        }
        a_field[lev]->FillBoundary();
        b_field[lev]->FillBoundary();
        c_field[lev]->FillBoundary();
    }

private:
    std::string mode = "perturbation"; //seed
    int period;
    int amp;
    int base;
    Set::Scalar ratio_a = 0.96;
    Set::Scalar ratio_b = 0.02;
    Set::Scalar ratio_c = 0.02;
    int Nx = 199;
    int Ny = 44;
    Set::Scalar dx = 1.0;
    Set::Scalar xi = 1.0;
    int Ny_ic = 49;
    int lambda = (Nx + 1) / 2;

public:
    static void Parse(ChanHillardMixture& value, IO::ParmParse& pp)
    {
        if (pp.contains("ic_mode"))
        {
            pp.query("ic_mode", mode);
        }
        pp.query_default("amp", amp, 15);
        pp.query_default("base", base, 20);
        pp.query_default("period", period, 10);
    }

};

}
#endif
