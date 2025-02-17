//
//

#ifndef BC_EXPRESSION_H_
#define BC_EXPRESSION_H_

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCRec.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Array.H>
#include <AMReX_Parser.H>

#include "Set/Set.H"
#include "BC/BC.H"
#include "Numeric/Interpolator/Linear.H"

namespace BC
{
class Expression
    : public BC<Set::Scalar>
{
public:
    static constexpr const char* name = "expression";

    //static constexpr const char* const Elastic::strings[];
    #if AMREX_SPACEDIM==2
    static const constexpr char * const facestr[] = {
        "xlo","ylo","xhi","yhi"
    };
    #elif AMREX_SPACEDIM==3
    static const constexpr char * const facestr[] = {
        "xlo","ylo","zlo","xhi","yhi","zhi",
    };
    #endif

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
    //Constant (amrex::Vector<amrex::Geometry> &_geom) : geom(_geom) {};
    Expression(int a_ncomp) : m_ncomp(a_ncomp) {};
    Expression(int a_ncomp, IO::ParmParse& pp, std::string name) : m_ncomp(a_ncomp)
    {
        pp_queryclass(name, *this);
    };

    virtual ~Expression() {};

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

    std::array<std::vector<int>, m_nfaces> m_bc_type;
    std::array<std::vector<amrex::Parser>, m_nfaces> m_bc_func_parser; 
    std::array<std::vector<amrex::ParserExecutor<4>>, m_nfaces> m_bc_func; 

    //std::array<std::vector<Numeric::Interpolator::Linear<Set::Scalar>>, m_nfaces> m_bc_val;

public:
    static void Parse(Expression& value, IO::ParmParse& pp)
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



        // TYPES

        std::vector<std::string> str;
        pp_queryarr("type.xlo", str);  // BC type on the lower x edge (2d) face (3d)
        for (unsigned int i = 0; i < str.size(); i++) if (!bcmap.count(str[i])) Util::Abort(INFO, "Invalid BC: ", str[i]);
        if (str.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_type[Face::XLO].push_back(bcmap[str[i]]);
        else if (str.size() == 1) value.m_bc_type[Face::XLO].resize(value.m_ncomp, bcmap[str[0]]);
        else Util::Abort(INFO, "Incorrect number of ", pp.prefix(), " BC type args: expected ", value.m_ncomp, " or 1 but got ", str.size());
        pp_queryarr("type.xhi", str);  // BC type on the upper x edge (2d) face (3d)
        for (unsigned int i = 0; i < str.size(); i++) if (!bcmap.count(str[i])) Util::Abort(INFO, "Invalid BC: ", str[i]);
        if (str.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_type[Face::XHI].push_back(bcmap[str[i]]);
        else if (str.size() == 1) value.m_bc_type[Face::XHI].resize(value.m_ncomp, bcmap[str[0]]);
        else Util::Abort(INFO, "Incorrect number of ", pp.prefix(), " BC type args: expected ", value.m_ncomp, " or 1 but got ", str.size());
        pp_queryarr("type.ylo", str);  // BC type on the lower y edge (2d) face (3d)
        for (unsigned int i = 0; i < str.size(); i++) if (!bcmap.count(str[i])) Util::Abort(INFO, "Invalid BC: ", str[i]);
        if (str.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_type[Face::YLO].push_back(bcmap[str[i]]);
        else if (str.size() == 1) value.m_bc_type[Face::YLO].resize(value.m_ncomp, bcmap[str[0]]);
        else Util::Abort(INFO, "Incorrect number of ", pp.prefix(), " BC type args: expected ", value.m_ncomp, " or 1 but got ", str.size());
        pp_queryarr("type.yhi", str);  // BC type on the upper y edge (2d) face (3d)
        for (unsigned int i = 0; i < str.size(); i++) if (!bcmap.count(str[i])) Util::Abort(INFO, "Invalid BC: ", str[i]);
        if (str.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_type[Face::YHI].push_back(bcmap[str[i]]);
        else if (str.size() == 1) value.m_bc_type[Face::YHI].resize(value.m_ncomp, bcmap[str[0]]);
        else Util::Abort(INFO, "Incorrect number of ", pp.prefix(), " BC type args: expected ", value.m_ncomp, " or 1 but got ", str.size());
        pp_queryarr("type.zlo", str);  // BC type on the lower z face (processed but ignored in 2d to prevent unused input errors)
#if AMREX_SPACEDIM==3
        for (unsigned int i = 0; i < str.size(); i++) if (!bcmap.count(str[i])) Util::Abort(INFO, "Invalid BC: ", str[i]);
        if (str.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_type[Face::ZLO].push_back(bcmap[str[i]]);
        else if (str.size() == 1) value.m_bc_type[Face::ZLO].resize(value.m_ncomp, bcmap[str[0]]);
        else Util::Abort(INFO, "Incorrect number of ", pp.prefix(), " BC type args: expected ", value.m_ncomp, " or 1 but got ", str.size());
#endif
        pp_queryarr("type.zhi", str);  // BC type on the upper z face (processed but ignored in 2d to prevent unused input errors)
#if AMREX_SPACEDIM==3
        for (unsigned int i = 0; i < str.size(); i++) if (!bcmap.count(str[i])) Util::Abort(INFO, "Invalid BC: ", str[i]);
        if (str.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_type[Face::ZHI].push_back(bcmap[str[i]]);
        else if (str.size() == 1) value.m_bc_type[Face::ZHI].resize(value.m_ncomp, bcmap[str[0]]);
        else Util::Abort(INFO, "Incorrect number of ", pp.prefix(), " BC type args: expected ", value.m_ncomp, " or 1 but got ", str.size());
#endif

        // VALS
        std::vector<std::string> val;
        for (int face = 0; face != Face::INT; face++)
        {
            std::string querystr = std::string("val.") + std::string(facestr[face]);
            if (pp.contains(querystr.c_str())) pp_queryarr(querystr.c_str(),val);
            else val.resize(value.m_ncomp,"0.0");
            if (val.size() != value.m_ncomp) Util::Abort(INFO,"Incorrect number of expressions specified for ",querystr,": expected 1 or ",value.m_ncomp," but received ",val.size());

            value.m_bc_func_parser[face].clear();
            value.m_bc_func[face].clear();
            for (unsigned int i = 0 ; i < value.m_ncomp; i++)
            {
                value.m_bc_func_parser[face].push_back( amrex::Parser(val[i].c_str()) );
                value.m_bc_func_parser[face][i].registerVariables({"x","y","z","t"});
                value.m_bc_func[face].push_back( value.m_bc_func_parser[face][i].compile<4>() );
            }
        }

#if AMREX_SPACEDIM==2
        // We may wish to use an input file that has 3D BC inputs
        // This will prevent the parser from complaining that there are unused inputs.
        std::vector<std::string> ignore_in_2d = {
            "zlo","zhi",
            "zhixlo","zloxlo","zhixhi","zloxhi","ylozlo","ylozhi","yhizlo","yhizhi",
            "xloylozlo","xloylozhi","xloyhizlo","xloyhizhi","xhiylozlo","xhiylozhi","xhiyhizlo","xhiyhizhi"};
        for (unsigned int n = 0; n < ignore_in_2d.size(); n++)
        {
            std::string querystr = std::string("val.") + ignore_in_2d[n];
            pp.ignore(querystr);
            querystr = std::string("type.") + ignore_in_2d[n];
            pp.ignore(querystr);
        }
#endif
        
    }
};
}
#endif
