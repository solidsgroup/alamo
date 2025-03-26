//
// This is the most commonly used standard boundary condition implementation.
// The name "Constant" refers to the invariance of the BC value or character along each face
// of the simulation domain; however, you may cause a change in the value with time by using a
// :ref:`Numeric::Interpolator::Linear` string.
// 
// The BC on each face is specified with a :code:`type` string that specifies the nature of 
// the BC, and a :code:`val` that specifies the value.
// By default, all types are Dirichlet and all values are 0.0.
// The complete list of BC types are below:
//
// `Dirichlet <https://en.wikipedia.org/wiki/Dirichlet_boundary_condition>`_  boundary condition
// types can be specified with :code:`Dirichlet`, :code:`dirichlet`, :code:`EXT_DIR`.
// `Neumann <https://en.wikipedia.org/wiki/Neumann_boundary_condition>`_ boundary condition
// types are specified with :code:`Neumann`, :code:`neumann`.
// `Periodic <https://en.wikipedia.org/wiki/Periodic_boundary_conditions>`_ boundary conditions
// can be specified with :code:`Periodic`, :code:`periodic`, :code:`INT_DIR`.
// **Important**: Ensure that your geometry is specified to be periodic using :code:`geometry.is_periodic`.
// For instance, if periodic in x, set to :code:`1 0 0`
//
// The BC values can be specified either as a number (e.g. :code:`1.0`) or using a linear interpolator
// string (e.g. :code:`(0,1:2.3,-0.1)`).
//
// The number of values and types must be either 0 (for defaults), 1 (to set the same for all field components)
// or N (where N=number of field components).
//

#ifndef BC_CONSTANT_H_
#define BC_CONSTANT_H_

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCRec.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Array.H>

#include "Set/Set.H"
#include "BC/BC.H"
#include "Numeric/Interpolator/Linear.H"

namespace BC
{
class Constant
    : public BC<Set::Scalar>
{

#if AMREX_SPACEDIM==2
    enum Face {
        XLO, YLO, XHI, YHI,
        INT
    };
#elif AMREX_SPACEDIM==3
    enum Face {
        XLO, YLO, ZLO, XHI, YHI, ZHI, // 6
        INT
    };
#endif

public:
    static constexpr const char* name = "constant";

    //Constant (amrex::Vector<amrex::Geometry> &_geom) : geom(_geom) {};
    Constant(int a_ncomp) : m_ncomp(a_ncomp) {};
    Constant(int a_ncomp, IO::ParmParse& pp, std::string name) : m_ncomp(a_ncomp)
    {
        pp_queryclass(name, *this);
    };
    Constant(int ncomp, amrex::Vector<std::string> bc_hi_str,
        amrex::Vector<std::string> bc_lo_str,
        AMREX_D_DECL(amrex::Vector<amrex::Real> _bc_lo_1,
            amrex::Vector<amrex::Real> _bc_lo_2,
            amrex::Vector<amrex::Real> _bc_lo_3),
        AMREX_D_DECL(amrex::Vector<amrex::Real> _bc_hi_1,
            amrex::Vector<amrex::Real> _bc_hi_2,
            amrex::Vector<amrex::Real> _bc_hi_3));

    static Constant ZeroNeumann(int ncomp = 1)
    {
        Constant ret = Constant(ncomp);

        for (int d = 0; d < m_nfaces; d++)
        {
            ret.m_bc_type[d].clear();
            for (int n = 0; n < ncomp; n++)
            {
                ret.m_bc_type[d].push_back((int)amrex::LinOpBCType::Neumann);
                ret.m_bc_val[d].push_back (0.0);;
            }
        }

        return ret;
    }

    static Constant ZeroDirichlet(int ncomp = 1)
    {
        Constant ret = Constant(ncomp);

        for (int d = 0; d < m_nfaces; d++)
            for (int n = 0; n < ncomp; n++)
            {
                ret.m_bc_type[n][d] = (int)amrex::LinOpBCType::Dirichlet;
                ret.m_bc_val [n][d]  = 0.0;
            }

        return ret;
    }


    virtual ~Constant() {};

    virtual void FillBoundary(amrex::BaseFab<Set::Scalar>& in, const amrex::Box& box,
        int ngrow, int dcomp, int ncomp, amrex::Real time,
        Orientation face = Orientation::All,
        const amrex::Mask* mask = nullptr) override;

    using BC::FillBoundary;

    amrex::BCRec GetBCRec() override;
    virtual amrex::Array<int, AMREX_SPACEDIM> IsPeriodic() override;
    virtual amrex::Periodicity Periodicity() const override;
    virtual amrex::Periodicity Periodicity(const amrex::Box& b) override;



    template<class T>
    const amrex::Array<amrex::Array<T, AMREX_SPACEDIM>, 2> GetBCTypes()
    {
        return { {{AMREX_D_DECL((T)m_bc_type[Face::XLO][0],(T)m_bc_type[Face::YLO][0],(T)m_bc_type[Face::ZLO][0])},
                {AMREX_D_DECL((T)m_bc_type[Face::XLO][0],(T)m_bc_type[Face::YLO][0],(T)m_bc_type[Face::ZLO][0])}} };
    }


private:
#if AMREX_SPACEDIM==2
    static const int m_nfaces = 4;
#elif AMREX_SPACEDIM==3
    static const int m_nfaces = 6;
#endif

    unsigned int m_ncomp = 0;

    //int bc_lo[BL_SPACEDIM];
    //int bc_hi[BL_SPACEDIM];
    //amrex::Vector<amrex::Real> AMREX_D_DECL(bc_lo_1, bc_lo_2, bc_lo_3);
    //amrex::Vector<amrex::Real> AMREX_D_DECL(bc_hi_1, bc_hi_2, bc_hi_3);

    std::array<std::vector<int>, m_nfaces> m_bc_type;
    std::array<std::vector<Numeric::Interpolator::Linear<Set::Scalar>>, m_nfaces> m_bc_val;

public:
    static void Parse(Constant& value, IO::ParmParse& pp)
    {
        std::map<std::string, int> bcmap;
        bcmap["BOGUS_BC"] = amrex::BCType::mathematicalBndryTypes::bogus;
        bcmap["INT_DIR"] = amrex::BCType::mathematicalBndryTypes::int_dir;
        bcmap["REFLECT_ODD"] = amrex::BCType::mathematicalBndryTypes::reflect_odd;
        bcmap["INT_DIR"] = amrex::BCType::mathematicalBndryTypes::int_dir;
        bcmap["REFLECT_EVEN"] = amrex::BCType::mathematicalBndryTypes::reflect_even;
        bcmap["FOEXTRAP"] = amrex::BCType::mathematicalBndryTypes::foextrap;
        bcmap["EXT_DIR"] = amrex::BCType::mathematicalBndryTypes::ext_dir;
        bcmap["HOEXTRAP"] = amrex::BCType::mathematicalBndryTypes::hoextrap;
        bcmap["Interior"] = amrex::BCType::mathematicalBndryTypes::int_dir;
        bcmap["Inflow"] = amrex::BCType::mathematicalBndryTypes::ext_dir;
        bcmap["Outflow"] = amrex::BCType::mathematicalBndryTypes::foextrap;
        bcmap["Symmetry"] = amrex::BCType::mathematicalBndryTypes::reflect_even;
        bcmap["SlipWall"] = amrex::BCType::mathematicalBndryTypes::ext_dir;
        bcmap["NoSlipWall"] = amrex::BCType::mathematicalBndryTypes::ext_dir;
        // From <AMReX_LO_BCTYPES.H>
        bcmap["interior"] = (int)amrex::LinOpBCType::interior;
        bcmap["Dirichlet"] = (int)amrex::LinOpBCType::Dirichlet;
        bcmap["dirichlet"] = (int)amrex::LinOpBCType::Dirichlet;
        bcmap["Neumann"] = (int)amrex::LinOpBCType::Neumann;
        bcmap["NEUMANN"] = (int)amrex::LinOpBCType::Neumann;
        bcmap["neumann"] = (int)amrex::LinOpBCType::Neumann;
        bcmap["reflect_odd"] = (int)amrex::LinOpBCType::reflect_odd;
        bcmap["Marshak"] = (int)amrex::LinOpBCType::Marshak;
        bcmap["SanchezPomraning"] = (int)amrex::LinOpBCType::SanchezPomraning;
        bcmap["inflow"] = (int)amrex::LinOpBCType::inflow;
        bcmap["Periodic"] = (int)amrex::LinOpBCType::Periodic;
        bcmap["periodic"] = (int)amrex::LinOpBCType::Periodic;



        value.m_bc_type[Face::XLO].clear(); value.m_bc_val[Face::XLO].clear();
        value.m_bc_type[Face::XHI].clear(); value.m_bc_val[Face::XHI].clear();
        value.m_bc_type[Face::YLO].clear(); value.m_bc_val[Face::YLO].clear();
        value.m_bc_type[Face::YHI].clear(); value.m_bc_val[Face::YHI].clear();
#if AMREX_SPACEDIM == 3
        value.m_bc_type[Face::ZLO].clear(); value.m_bc_val[Face::ZLO].clear();
        value.m_bc_type[Face::ZHI].clear(); value.m_bc_val[Face::ZHI].clear();
#endif

        // TYPES

        std::vector<std::string> str;
        pp_queryarr_default("type.xlo", str,"dirichlet");  // BC type on the lower x edge (2d) face (3d)
        for (unsigned int i = 0; i < str.size(); i++) if (!bcmap.count(str[i])) Util::Exception(INFO, "Invalid BC: ", str[i]);
        if (str.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_type[Face::XLO].push_back(bcmap[str[i]]);
        else if (str.size() == 1) value.m_bc_type[Face::XLO].resize(value.m_ncomp, bcmap[str[0]]);
        else Util::Exception(INFO, "Incorrect number of ", pp.prefix(), " BC type args: expected ", value.m_ncomp, " or 1 but got ", str.size());
        pp_queryarr_default("type.xhi", str,"dirichlet");  // BC type on the upper x edge (2d) face (3d)
        for (unsigned int i = 0; i < str.size(); i++) if (!bcmap.count(str[i])) Util::Exception(INFO, "Invalid BC: ", str[i]);
        if (str.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_type[Face::XHI].push_back(bcmap[str[i]]);
        else if (str.size() == 1) value.m_bc_type[Face::XHI].resize(value.m_ncomp, bcmap[str[0]]);
        else Util::Exception(INFO, "Incorrect number of ", pp.prefix(), " BC type args: expected ", value.m_ncomp, " or 1 but got ", str.size());
        pp_queryarr_default("type.ylo", str,"dirichlet");  // BC type on the lower y edge (2d) face (3d)
        for (unsigned int i = 0; i < str.size(); i++) if (!bcmap.count(str[i])) Util::Exception(INFO, "Invalid BC: ", str[i]);
        if (str.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_type[Face::YLO].push_back(bcmap[str[i]]);
        else if (str.size() == 1) value.m_bc_type[Face::YLO].resize(value.m_ncomp, bcmap[str[0]]);
        else Util::Exception(INFO, "Incorrect number of ", pp.prefix(), " BC type args: expected ", value.m_ncomp, " or 1 but got ", str.size());
        pp_queryarr_default("type.yhi", str,"dirichlet");  // BC type on the upper y edge (2d) face (3d)
        for (unsigned int i = 0; i < str.size(); i++) if (!bcmap.count(str[i])) Util::Exception(INFO, "Invalid BC: ", str[i]);
        if (str.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_type[Face::YHI].push_back(bcmap[str[i]]);
        else if (str.size() == 1) value.m_bc_type[Face::YHI].resize(value.m_ncomp, bcmap[str[0]]);
        else Util::Exception(INFO, "Incorrect number of ", pp.prefix(), " BC type args: expected ", value.m_ncomp, " or 1 but got ", str.size());
        pp_queryarr_default("type.zlo", str,"dirichlet");  // BC type on the lower z face (processed but ignored in 2d to prevent unused input errors)
#if AMREX_SPACEDIM==3
        for (unsigned int i = 0; i < str.size(); i++) if (!bcmap.count(str[i])) Util::Exception(INFO, "Invalid BC: ", str[i]);
        if (str.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_type[Face::ZLO].push_back(bcmap[str[i]]);
        else if (str.size() == 1) value.m_bc_type[Face::ZLO].resize(value.m_ncomp, bcmap[str[0]]);
        else Util::Exception(INFO, "Incorrect number of ", pp.prefix(), " BC type args: expected ", value.m_ncomp, " or 1 but got ", str.size());
#endif
        pp_queryarr_default("type.zhi", str,"dirichlet");  // BC type on the upper z face (processed but ignored in 2d to prevent unused input errors)
#if AMREX_SPACEDIM==3
        for (unsigned int i = 0; i < str.size(); i++) if (!bcmap.count(str[i])) Util::Exception(INFO, "Invalid BC: ", str[i]);
        if (str.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_type[Face::ZHI].push_back(bcmap[str[i]]);
        else if (str.size() == 1) value.m_bc_type[Face::ZHI].resize(value.m_ncomp, bcmap[str[0]]);
        else Util::Exception(INFO, "Incorrect number of ", pp.prefix(), " BC type args: expected ", value.m_ncomp, " or 1 but got ", str.size());
#endif

        // VALS
        std::vector<std::string> val;
        value.m_bc_val[Face::XLO].clear();
        pp_queryarr_default("val.xlo", val,"0.0");  // BC value on the lower x edge (2d) face (3d)
        if (val.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_val[Face::XLO].push_back(Numeric::Interpolator::Linear<Set::Scalar>(val[i]));
        else if (val.size() == 1) value.m_bc_val[Face::XLO].resize(value.m_ncomp, Numeric::Interpolator::Linear<Set::Scalar>(val[0]));
        else if (val.size() == 0) value.m_bc_val[Face::XLO].resize(value.m_ncomp, 0.0);
        else Util::Exception(INFO, "Incorrect number of ", pp.prefix(), " BC value args: expected ", value.m_ncomp, " or 0 or 1 but got ", val.size());
        value.m_bc_val[Face::XHI].clear();
        pp_queryarr_default("val.xhi", val,"0.0"); // BC value on the upper x edge (2d) face (3d)
        if (val.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_val[Face::XHI].push_back(Numeric::Interpolator::Linear<Set::Scalar>(val[i]));
        else if (val.size() == 1) value.m_bc_val[Face::XHI].resize(value.m_ncomp, Numeric::Interpolator::Linear<Set::Scalar>(val[0]));
        else if (val.size() == 0) value.m_bc_val[Face::XHI].resize(value.m_ncomp, 0.0);
        else Util::Exception(INFO, "Incorrect number of ", pp.prefix(), " BC value args: expected ", value.m_ncomp, " or 0 or 1 but got ", val.size());
        value.m_bc_val[Face::YLO].clear();
        pp_queryarr_default("val.ylo", val,"0.0");  // BC value on the lower y edge (2d) face (3d)
        if (val.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_val[Face::YLO].push_back(Numeric::Interpolator::Linear<Set::Scalar>(val[i]));
        else if (val.size() == 1) value.m_bc_val[Face::YLO].resize(value.m_ncomp, Numeric::Interpolator::Linear<Set::Scalar>(val[0]));
        else if (val.size() == 0) value.m_bc_val[Face::YLO].resize(value.m_ncomp, 0.0);
        else Util::Exception(INFO, "Incorrect number of ", pp.prefix(), " BC value args: expected ", value.m_ncomp, " or 0 or 1 but got ", val.size());
        value.m_bc_val[Face::YHI].clear();
        pp_queryarr_default("val.yhi", val,"0.0"); // BC value on the upper y edge (2d) face (3d)
        if (val.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_val[Face::YHI].push_back(Numeric::Interpolator::Linear<Set::Scalar>(val[i]));
        else if (val.size() == 1) value.m_bc_val[Face::YHI].resize(value.m_ncomp, Numeric::Interpolator::Linear<Set::Scalar>(val[0]));
        else if (val.size() == 0) value.m_bc_val[Face::YHI].resize(value.m_ncomp, 0.0);
        else Util::Exception(INFO, "Incorrect number of ", pp.prefix(), " BC value args: expected ", value.m_ncomp, " or 0 or 1 but got ", val.size());
        pp_queryarr_default("val.zlo", val,"0.0");  // BC value on the lower z face (processed but ignored in 2d to prevent unused input errors)
#if AMREX_SPACEDIM==3
        value.m_bc_val[Face::ZLO].clear();
        if (val.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_val[Face::ZLO].push_back(Numeric::Interpolator::Linear<Set::Scalar>(val[i]));
        else if (val.size() == 1) value.m_bc_val[Face::ZLO].resize(value.m_ncomp, Numeric::Interpolator::Linear<Set::Scalar>(val[0]));
        else if (val.size() == 0) value.m_bc_val[Face::ZLO].resize(value.m_ncomp, 0.0);
        else Util::Exception(INFO, "Incorrect number of ", pp.prefix(), " BC value args: expected ", value.m_ncomp, " or 0 or 1 but got ", val.size());
#endif
        pp_queryarr_default("val.zhi", val,"0.0"); // BC value on the upper z face (processed but ignored in 2d to prevent unused input errors)
#if AMREX_SPACEDIM==3
        value.m_bc_val[Face::ZHI].clear();
        if (val.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_val[Face::ZHI].push_back(Numeric::Interpolator::Linear<Set::Scalar>(val[i]));
        else if (val.size() == 1) value.m_bc_val[Face::ZHI].resize(value.m_ncomp, Numeric::Interpolator::Linear<Set::Scalar>(val[0]));
        else if (val.size() == 0) value.m_bc_val[Face::ZHI].resize(value.m_ncomp, 0.0);
        else Util::Exception(INFO, "Incorrect number of ", pp.prefix(), " BC value args: expected ", value.m_ncomp, " or 0 or 1 but got ", val.size());
#endif
    }

};
}
#endif
