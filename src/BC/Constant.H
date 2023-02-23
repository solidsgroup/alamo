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
    enum Face{
        XLO, YLO, XHI, YHI, 
        INT
    };
    #elif AMREX_SPACEDIM==3
    enum Face{
        XLO, YLO, ZLO, XHI, YHI, ZHI, // 6
        INT
    };
    #endif

public:
    //Constant (amrex::Vector<amrex::Geometry> &_geom) : geom(_geom) {};
    Constant(int a_ncomp) : m_ncomp(a_ncomp) {} ;
    Constant(int a_ncomp, IO::ParmParse &pp, std::string name) : m_ncomp(a_ncomp) 
    {pp.queryclass(name,*this);} ;
    Constant (int ncomp, amrex::Vector<std::string> bc_hi_str,
        amrex::Vector<std::string> bc_lo_str,
        AMREX_D_DECL(amrex::Vector<amrex::Real> _bc_lo_1,
                    amrex::Vector<amrex::Real> _bc_lo_2,
                    amrex::Vector<amrex::Real> _bc_lo_3),
        AMREX_D_DECL(amrex::Vector<amrex::Real> _bc_hi_1,
                    amrex::Vector<amrex::Real> _bc_hi_2,
                    amrex::Vector<amrex::Real> _bc_hi_3));

    virtual ~Constant() {};

    virtual void FillBoundary (amrex::BaseFab<Set::Scalar> &in, const amrex::Box &box,
                    int ngrow, int dcomp, int ncomp, amrex::Real time,
                    Orientation face = Orientation::All,
                    const amrex::Mask *mask = nullptr) override;

    using BC::FillBoundary;
    
    amrex::BCRec GetBCRec() override;
    virtual amrex::Array<int,AMREX_SPACEDIM> IsPeriodic() override;
    virtual amrex::Periodicity Periodicity () const override;
    virtual amrex::Periodicity Periodicity (const amrex::Box& b) override;



    template<class T>
    const amrex::Array<amrex::Array<T,AMREX_SPACEDIM>,2> GetBCTypes()
    {
        return {{{AMREX_D_DECL((T)m_bc_type[Face::XLO][0],(T)m_bc_type[Face::YLO][0],(T)m_bc_type[Face::ZLO][0])},
                {AMREX_D_DECL((T)m_bc_type[Face::XLO][0],(T)m_bc_type[Face::YLO][0],(T)m_bc_type[Face::ZLO][0])}}};
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

    std::array<std::vector<int>,                                        m_nfaces> m_bc_type; 
    std::array<std::vector<Numeric::Interpolator::Linear<Set::Scalar>>, m_nfaces> m_bc_val; 

public:
    static void Parse(Constant & value, IO::ParmParse & pp)
    {
        std::map<std::string, int> bcmap;
            bcmap["BOGUS_BC"]         = BOGUS_BC;
            bcmap["INT_DIR"]          = INT_DIR;
            bcmap["REFLECT_ODD"]      = REFLECT_ODD;
            bcmap["INT_DIR"]          = INT_DIR;
            bcmap["REFLECT_EVEN"]     = REFLECT_EVEN;
            bcmap["FOEXTRAP"]         = FOEXTRAP;
            bcmap["EXT_DIR"]          = EXT_DIR;
            bcmap["HOEXTRAP"]         = HOEXTRAP;
            bcmap["Interior"]         = Interior;
            bcmap["Inflow"]           = Inflow;
            bcmap["Outflow"]          = Outflow;
            bcmap["Symmetry"]         = Symmetry;
            bcmap["SlipWall"]         = SlipWall;
            bcmap["NoSlipWall"]       = NoSlipWall;
        // From <AMReX_LO_BCTYPES.H>
        bcmap["interior"]         = (int)amrex::LinOpBCType::interior;
        bcmap["Dirichlet"]        = (int)amrex::LinOpBCType::Dirichlet;
        bcmap["dirichlet"]        = (int)amrex::LinOpBCType::Dirichlet;
        bcmap["Neumann"]          = (int)amrex::LinOpBCType::Neumann;
        bcmap["NEUMANN"]          = (int)amrex::LinOpBCType::Neumann;
        bcmap["neumann"]          = (int)amrex::LinOpBCType::Neumann;
        bcmap["reflect_odd"]      = (int)amrex::LinOpBCType::reflect_odd;
        bcmap["Marshak"]          = (int)amrex::LinOpBCType::Marshak;
        bcmap["SanchezPomraning"] = (int)amrex::LinOpBCType::SanchezPomraning;
        bcmap["inflow"]           = (int)amrex::LinOpBCType::inflow;
        bcmap["Periodic"]         = (int)amrex::LinOpBCType::Periodic;
        bcmap["periodic"]         = (int)amrex::LinOpBCType::Periodic;


        
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
        pp.queryarr("type.xlo",str); 
            for (unsigned int i = 0; i < str.size(); i++) if (!bcmap.count(str[i])) Util::Abort(INFO,"Invalid BC: ", str[i]);
            if (str.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_type[Face::XLO].push_back(bcmap[str[i]]); 
            else if (str.size() == 1) value.m_bc_type[Face::XLO].resize(value.m_ncomp,bcmap[str[0]]);
            else Util::Abort(INFO,"Incorrect number of ",pp.prefix()," BC type args: expected ", value.m_ncomp, " or 1 but got ", str.size());
        pp.queryarr("type.xhi",str);
            for (unsigned int i = 0; i < str.size(); i++) if (!bcmap.count(str[i])) Util::Abort(INFO,"Invalid BC: ", str[i]);
            if (str.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_type[Face::XHI].push_back(bcmap[str[i]]); 
            else if (str.size() == 1) value.m_bc_type[Face::XHI].resize(value.m_ncomp,bcmap[str[0]]);
            else Util::Abort(INFO,"Incorrect number of ",pp.prefix()," BC type args: expected ", value.m_ncomp, " or 1 but got ", str.size());
        pp.queryarr("type.ylo",str);
            for (unsigned int i = 0; i < str.size(); i++) if (!bcmap.count(str[i])) Util::Abort(INFO,"Invalid BC: ", str[i]);
            if (str.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_type[Face::YLO].push_back(bcmap[str[i]]); 
            else if (str.size() == 1) value.m_bc_type[Face::YLO].resize(value.m_ncomp,bcmap[str[0]]);
            else Util::Abort(INFO,"Incorrect number of ",pp.prefix()," BC type args: expected ", value.m_ncomp, " or 1 but got ", str.size());
        pp.queryarr("type.yhi",str);
            for (unsigned int i = 0; i < str.size(); i++) if (!bcmap.count(str[i])) Util::Abort(INFO,"Invalid BC: ", str[i]);
            if (str.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_type[Face::YHI].push_back(bcmap[str[i]]); 
            else if (str.size() == 1) value.m_bc_type[Face::YHI].resize(value.m_ncomp,bcmap[str[0]]);
            else Util::Abort(INFO,"Incorrect number of ",pp.prefix()," BC type args: expected ", value.m_ncomp, " or 1 but got ", str.size());
        pp.queryarr("type.zlo",str);
        #if AMREX_SPACEDIM==3
            for (unsigned int i = 0; i < str.size(); i++) if (!bcmap.count(str[i])) Util::Abort(INFO,"Invalid BC: ", str[i]);
            if (str.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_type[Face::ZLO].push_back(bcmap[str[i]]); 
            else if (str.size() == 1) value.m_bc_type[Face::ZLO].resize(value.m_ncomp,bcmap[str[0]]);
            else Util::Abort(INFO,"Incorrect number of ",pp.prefix()," BC type args: expected ", value.m_ncomp, " or 1 but got ", str.size());
        #endif
        pp.queryarr("type.zhi",str);
        #if AMREX_SPACEDIM==3
            for (unsigned int i = 0; i < str.size(); i++) if (!bcmap.count(str[i])) Util::Abort(INFO,"Invalid BC: ", str[i]);
            if (str.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_type[Face::ZHI].push_back(bcmap[str[i]]); 
            else if (str.size() == 1) value.m_bc_type[Face::ZHI].resize(value.m_ncomp,bcmap[str[0]]);
            else Util::Abort(INFO,"Incorrect number of ",pp.prefix()," BC type args: expected ", value.m_ncomp, " or 1 but got ", str.size());
        #endif

        // VALS
        std::vector<std::string> val;
        value.m_bc_val[Face::XLO].clear();
        pp.queryarr("val.xlo",val); 
            if (val.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_val[Face::XLO].push_back(Numeric::Interpolator::Linear<Set::Scalar>(val[i])); 
            else if (val.size() == 1) value.m_bc_val[Face::XLO].resize(value.m_ncomp,Numeric::Interpolator::Linear<Set::Scalar>(val[0]));
            else if (val.size() == 0) value.m_bc_val[Face::XLO].resize(value.m_ncomp,0.0);
            else Util::Abort(INFO,"Incorrect number of ",pp.prefix()," BC value args: expected ", value.m_ncomp, " or 0 or 1 but got ", val.size());        
        value.m_bc_val[Face::XHI].clear();
        pp.queryarr("val.xhi",val); 
            if (val.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_val[Face::XHI].push_back(Numeric::Interpolator::Linear<Set::Scalar>(val[i])); 
            else if (val.size() == 1) value.m_bc_val[Face::XHI].resize(value.m_ncomp,Numeric::Interpolator::Linear<Set::Scalar>(val[0]));
            else if (val.size() == 0) value.m_bc_val[Face::XHI].resize(value.m_ncomp,0.0);
            else Util::Abort(INFO,"Incorrect number of ",pp.prefix()," BC value args: expected ", value.m_ncomp, " or 0 or 1 but got ", val.size());        
        value.m_bc_val[Face::YLO].clear();
        pp.queryarr("val.ylo",val); 
            if (val.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_val[Face::YLO].push_back(Numeric::Interpolator::Linear<Set::Scalar>(val[i])); 
            else if (val.size() == 1) value.m_bc_val[Face::YLO].resize(value.m_ncomp,Numeric::Interpolator::Linear<Set::Scalar>(val[0]));
            else if (val.size() == 0) value.m_bc_val[Face::YLO].resize(value.m_ncomp,0.0);
            else Util::Abort(INFO,"Incorrect number of ",pp.prefix()," BC value args: expected ", value.m_ncomp, " or 0 or 1 but got ", val.size());        
        value.m_bc_val[Face::YHI].clear();
        pp.queryarr("val.yhi",val); 
            if (val.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_val[Face::YHI].push_back(Numeric::Interpolator::Linear<Set::Scalar>(val[i])); 
            else if (val.size() == 1) value.m_bc_val[Face::YHI].resize(value.m_ncomp,Numeric::Interpolator::Linear<Set::Scalar>(val[0]));
            else if (val.size() == 0) value.m_bc_val[Face::YHI].resize(value.m_ncomp,0.0);
            else Util::Abort(INFO,"Incorrect number of ",pp.prefix()," BC value args: expected ", value.m_ncomp, " or 0 or 1 but got ", val.size());        
        pp.queryarr("val.zlo",val); 
        #if AMREX_SPACEDIM==3
            value.m_bc_val[Face::ZLO].clear();
            if (val.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_val[Face::ZLO].push_back(Numeric::Interpolator::Linear<Set::Scalar>(val[i])); 
            else if (val.size() == 1) value.m_bc_val[Face::ZLO].resize(value.m_ncomp,Numeric::Interpolator::Linear<Set::Scalar>(val[0]));
            else if (val.size() == 0) value.m_bc_val[Face::ZLO].resize(value.m_ncomp,0.0);
            else Util::Abort(INFO,"Incorrect number of ",pp.prefix()," BC value args: expected ", value.m_ncomp, " or 0 or 1 but got ", val.size());        
        #endif
        pp.queryarr("val.zhi",val); 
        #if AMREX_SPACEDIM==3
            value.m_bc_val[Face::ZHI].clear();
            if (val.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_val[Face::ZHI].push_back(Numeric::Interpolator::Linear<Set::Scalar>(val[i])); 
            else if (val.size() == 1) value.m_bc_val[Face::ZHI].resize(value.m_ncomp,Numeric::Interpolator::Linear<Set::Scalar>(val[0]));
            else if (val.size() == 0) value.m_bc_val[Face::ZHI].resize(value.m_ncomp,0.0);
            else Util::Abort(INFO,"Incorrect number of ",pp.prefix()," BC value args: expected ", value.m_ncomp, " or 0 or 1 but got ", val.size());        
        #endif
    }

};
}
#endif
