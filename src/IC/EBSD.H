#ifndef IC_EBSD_H
#define IC_EBSD_H

#include "IC/IC.H"
#include "IO/ParmParse.H"
#include "Util/PNG.H"
#include "Model/Solid/Linear/Cubic.H"

namespace IC
{
template<class MODEL>
class EBSD : public IC<MODEL>
{
public:
    static constexpr const char* name = "ebsd";

    EBSD(amrex::Vector<amrex::Geometry> &_geom) : IC<MODEL>(_geom) {}
    EBSD(amrex::Vector<amrex::Geometry> &_geom, IO::ParmParse &pp, std::string name) : EBSD(_geom)
    {
        pp_queryclass(name,*this);
    }

    virtual void
    Add(const int &lev, Set::Field<MODEL> &a_field, Set::Scalar) override
    {
        amrex::Vector<amrex::Geometry> &geom = IC<MODEL>::geom;

        Set::Vector DX(geom[lev].CellSize());
        amrex::Box domain = geom[lev].Domain();

        amrex::IndexType type = a_field[lev]->ixType();
        domain.convert(type);

        Set::Vector domlo(AMREX_D_DECL(geom[lev].ProbLo()[0], geom[lev].ProbLo()[1], 0.0));
        Set::Vector domhi(AMREX_D_DECL(geom[lev].ProbHi()[0], geom[lev].ProbHi()[1], 0.0));

        png.setDomain(domlo, domhi);

        for (amrex::MFIter mfi(*a_field[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx;
            if (type == amrex::IndexType::TheNodeType()) bx = mfi.grownnodaltilebox();
            if (type == amrex::IndexType::TheCellType()) bx = mfi.growntilebox();
            bx = bx & domain;

            amrex::Array4<MODEL> const& field = a_field[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Vector x = Set::Position(i,j,k,geom[lev],type);
                std::array<Set::Scalar, 4> val = png(x);

                Eigen::Quaternion q(val[0], val[1], val[2], val[3]);
                
                field(i,j,k) = base;
                //field(i,j,k).Rotate(q.toRotationMatrix());
            });
        }
        a_field[lev]->FillBoundary();
    }


    static void Parse(EBSD &value, IO::ParmParse &pp)
    {
        // png image storing EBSD data
        pp.queryclass<Util::PNG>(value.png);
        // read in the parameters for the base model
        pp.queryclass<MODEL>(value.base);
    }

private:
    Util::PNG png;
    MODEL base;
};

}



#endif
