#ifndef BC_INTBC_H_
#define BC_INTBC_H_

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCRec.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Array.H>
#include "Set/Set.H"
#include "BC/BC.H"
#include "BC/Nothing.H"
#include "IO/ParmParse.H" // Including IO for ParmParse

namespace BC {

/// \class NothingInt
/// \brief Specialized "do-nothing" boundary condition class for integer fields (`Set::IntScalar`).
class NothingInt : public BC<Set::IntScalar> {
public:
    NothingInt() = default;
    virtual ~NothingInt() {}

    /// Override FillBoundary for integer fields
    void FillBoundary(amrex::BaseFab<Set::IntScalar>& in,
                      const amrex::Box& box,
                      int ngrow, int dcomp, int ncomp, amrex::Real time,
                      Orientation face = Orientation::All,
                      const amrex::Mask* mask = nullptr) override;

    /// Override GetBCRec
    amrex::BCRec GetBCRec() override {
        return amrex::BCRec();
    }
};

/// \class ConstantInt
/// \brief Constant boundary condition class for integer fields (`Set::IntScalar`).
class ConstantInt : public BC<Set::IntScalar> {
public:
    ConstantInt(int a_ncomp);
    ConstantInt(int a_ncomp, IO::ParmParse& pp, std::string name);
    ConstantInt(int ncomp, amrex::Vector<std::string> bc_hi_str,
                amrex::Vector<std::string> bc_lo_str,
                AMREX_D_DECL(amrex::Vector<int> _bc_lo_1,
                             amrex::Vector<int> _bc_lo_2,
                             amrex::Vector<int> _bc_lo_3),
                AMREX_D_DECL(amrex::Vector<int> _bc_hi_1,
                             amrex::Vector<int> _bc_hi_2,
                             amrex::Vector<int> _bc_hi_3));

    virtual ~ConstantInt() {}

    /// Override FillBoundary for integer fields
    void FillBoundary(amrex::BaseFab<Set::IntScalar>& in,
                      const amrex::Box& box,
                      int ngrow, int dcomp, int ncomp, amrex::Real time,
                      Orientation face = Orientation::All,
                      const amrex::Mask* mask = nullptr) override;

    /// Override GetBCRec
    amrex::BCRec GetBCRec() override;

    /// Override IsPeriodic
    amrex::Array<int, AMREX_SPACEDIM> IsPeriodic() override;

    /// Override Periodicity functions
    amrex::Periodicity Periodicity() const override;
    amrex::Periodicity Periodicity(const amrex::Box& b) override;

private:
#if AMREX_SPACEDIM == 2
    static const int m_nfaces = 4;
#elif AMREX_SPACEDIM == 3
    static const int m_nfaces = 6;
#endif

    /// Enumeration for domain faces
#if AMREX_SPACEDIM==2
    enum Face {
        XLO, XHI, YLO, YHI,
	INT
       };		
#elif AMREX_SPACEDIM == 3
    enum Face {
	XLO, YLO, ZLO, XHI, YHI, ZHI,    
        INT
        };		
#endif


    unsigned int m_ncomp = 0;
    std::array<std::vector<int>, m_nfaces> m_bc_type;
    std::array<std::vector<int>, m_nfaces> m_bc_val;

public:
    /// Parse input from ParmParse
    static void Parse(ConstantInt& value, IO::ParmParse& pp) {
        std::map<std::string, int> bcmap;
        bcmap["Neumann"] = (int)amrex::LinOpBCType::Neumann;
        bcmap["neumann"] = (int)amrex::LinOpBCType::Neumann;

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

        // XLO
        pp.queryarr("type.xlo", str);
        for (unsigned int i = 0; i < str.size(); i++) 
            if (!bcmap.count(str[i])) Util::Exception(INFO, "Invalid BC: ", str[i]);
        if (str.size() == value.m_ncomp) 
            for (unsigned int i = 0; i < value.m_ncomp; i++) 
                value.m_bc_type[Face::XLO].push_back(bcmap[str[i]]);
        else if (str.size() == 1) 
            value.m_bc_type[Face::XLO].resize(value.m_ncomp, bcmap[str[0]]);
        else 
            Util::Exception(INFO, "Incorrect number of ", pp.prefix(), " BC type args: expected ", value.m_ncomp, " or 1 but got ", str.size());

        // Repeat for other faces

        // XHI
        pp.queryarr("type.xhi", str);
        for (unsigned int i = 0; i < str.size(); i++) 
            if (!bcmap.count(str[i])) Util::Exception(INFO, "Invalid BC: ", str[i]);
        if (str.size() == value.m_ncomp) 
            for (unsigned int i = 0; i < value.m_ncomp; i++) 
                value.m_bc_type[Face::XHI].push_back(bcmap[str[i]]);
        else if (str.size() == 1) 
            value.m_bc_type[Face::XHI].resize(value.m_ncomp, bcmap[str[0]]);
        else 
            Util::Exception(INFO, "Incorrect number of ", pp.prefix(), " BC type args: expected ", value.m_ncomp, " or 1 but got ", str.size());

        // YLO
        pp.queryarr("type.ylo", str);
        for (unsigned int i = 0; i < str.size(); i++) 
            if (!bcmap.count(str[i])) Util::Exception(INFO, "Invalid BC: ", str[i]);
        if (str.size() == value.m_ncomp) 
            for (unsigned int i = 0; i < value.m_ncomp; i++) 
                value.m_bc_type[Face::YLO].push_back(bcmap[str[i]]);
        else if (str.size() == 1) 
            value.m_bc_type[Face::YLO].resize(value.m_ncomp, bcmap[str[0]]);
        else 
            Util::Exception(INFO, "Incorrect number of ", pp.prefix(), " BC type args: expected ", value.m_ncomp, " or 1 but got ", str.size());

        // YHI
        pp.queryarr("type.yhi", str);
        for (unsigned int i = 0; i < str.size(); i++) 
            if (!bcmap.count(str[i])) Util::Exception(INFO, "Invalid BC: ", str[i]);
        if (str.size() == value.m_ncomp) 
            for (unsigned int i = 0; i < value.m_ncomp; i++) 
                value.m_bc_type[Face::YHI].push_back(bcmap[str[i]]);
        else if (str.size() == 1) 
            value.m_bc_type[Face::YHI].resize(value.m_ncomp, bcmap[str[0]]);
        else 
            Util::Exception(INFO, "Incorrect number of ", pp.prefix(), " BC type args: expected ", value.m_ncomp, " or 1 but got ", str.size());

#if AMREX_SPACEDIM == 3
        // ZLO
        pp.queryarr("type.zlo", str);
        for (unsigned int i = 0; i < str.size(); i++) 
            if (!bcmap.count(str[i])) Util::Exception(INFO, "Invalid BC: ", str[i]);
        if (str.size() == value.m_ncomp) 
            for (unsigned int i = 0; i < value.m_ncomp; i++) 
                value.m_bc_type[Face::ZLO].push_back(bcmap[str[i]]);
        else if (str.size() == 1) 
            value.m_bc_type[Face::ZLO].resize(value.m_ncomp, bcmap[str[0]]);
        else 
            Util::Exception(INFO, "Incorrect number of ", pp.prefix(), " BC type args: expected ", value.m_ncomp, " or 1 but got ", str.size());

        // ZHI
        pp.queryarr("type.zhi", str);
        for (unsigned int i = 0; i < str.size(); i++) 
            if (!bcmap.count(str[i])) Util::Exception(INFO, "Invalid BC: ", str[i]);
        if (str.size() == value.m_ncomp) 
            for (unsigned int i = 0; i < value.m_ncomp; i++) 
                value.m_bc_type[Face::ZHI].push_back(bcmap[str[i]]);
        else if (str.size() == 1) 
            value.m_bc_type[Face::ZHI].resize(value.m_ncomp, bcmap[str[0]]);
        else 
            Util::Exception(INFO, "Incorrect number of ", pp.prefix(), " BC type args: expected ", value.m_ncomp, " or 1 but got ", str.size());
#endif

        // VALS
	std::vector<int> val;

	// XLO Face
	value.m_bc_val[Face::XLO].clear();
	pp.queryarr("val.xlo", val);  // BC value on the lower x edge (2d) face (3d)
	if (val.size() == value.m_ncomp) {
	    for (unsigned int i = 0; i < value.m_ncomp; i++) {
		value.m_bc_val[Face::XLO].push_back(val[i]);
	    }
	} else if (val.size() == 1) {
	    value.m_bc_val[Face::XLO].resize(value.m_ncomp, val[0]);
	} else if (val.size() == 0) {
	    value.m_bc_val[Face::XLO].resize(value.m_ncomp, 0);
	} else {
	    Util::Exception(INFO, "Incorrect number of ", pp.prefix(), " BC value args: expected ", value.m_ncomp, " or 0 or 1 but got ", val.size());
	}

	// XHI Face
	value.m_bc_val[Face::XHI].clear();
	pp.queryarr("val.xhi", val);  // BC value on the upper x edge (2d) face (3d)
	if (val.size() == value.m_ncomp) {
	    for (unsigned int i = 0; i < value.m_ncomp; i++) {
		value.m_bc_val[Face::XHI].push_back(val[i]);
	    }
	} else if (val.size() == 1) {
	    value.m_bc_val[Face::XHI].resize(value.m_ncomp, val[0]);
	} else if (val.size() == 0) {
	    value.m_bc_val[Face::XHI].resize(value.m_ncomp, 0);
	} else {
	    Util::Exception(INFO, "Incorrect number of ", pp.prefix(), " BC value args: expected ", value.m_ncomp, " or 0 or 1 but got ", val.size());
	}

	// YLO Face
	value.m_bc_val[Face::YLO].clear();
	pp.queryarr("val.ylo", val);  // BC value on the lower y edge (2d) face (3d)
	if (val.size() == value.m_ncomp) {
	    for (unsigned int i = 0; i < value.m_ncomp; i++) {
		value.m_bc_val[Face::YLO].push_back(val[i]);
	    }
	} else if (val.size() == 1) {
	    value.m_bc_val[Face::YLO].resize(value.m_ncomp, val[0]);
	} else if (val.size() == 0) {
	    value.m_bc_val[Face::YLO].resize(value.m_ncomp, 0);
	} else {
	    Util::Exception(INFO, "Incorrect number of ", pp.prefix(), " BC value args: expected ", value.m_ncomp, " or 0 or 1 but got ", val.size());
	}

	// YHI Face
	value.m_bc_val[Face::YHI].clear();
	pp.queryarr("val.yhi", val);  // BC value on the upper y edge (2d) face (3d)
	if (val.size() == value.m_ncomp) {
	    for (unsigned int i = 0; i < value.m_ncomp; i++) {
		value.m_bc_val[Face::YHI].push_back(val[i]);
	    }
	} else if (val.size() == 1) {
	    value.m_bc_val[Face::YHI].resize(value.m_ncomp, val[0]);
	} else if (val.size() == 0) {
	    value.m_bc_val[Face::YHI].resize(value.m_ncomp, 0);
	} else {
	    Util::Exception(INFO, "Incorrect number of ", pp.prefix(), " BC value args: expected ", value.m_ncomp, " or 0 or 1 but got ", val.size());
	}

#if AMREX_SPACEDIM == 3
	// ZLO Face
	value.m_bc_val[Face::ZLO].clear();
	pp.queryarr("val.zlo", val);  // BC value on the lower z edge (3d)
	if (val.size() == value.m_ncomp) {
	    for (unsigned int i = 0; i < value.m_ncomp; i++) {
		value.m_bc_val[Face::ZLO].push_back(val[i]);
	    }
	} else if (val.size() == 1) {
	    value.m_bc_val[Face::ZLO].resize(value.m_ncomp, val[0]);
	} else if (val.size() == 0) {
	    value.m_bc_val[Face::ZLO].resize(value.m_ncomp, 0);
	} else {
	    Util::Exception(INFO, "Incorrect number of ", pp.prefix(), " BC value args: expected ", value.m_ncomp, " or 0 or 1 but got ", val.size());
	}

	// ZHI Face
	value.m_bc_val[Face::ZHI].clear();
	pp.queryarr("val.zhi", val);  // BC value on the upper z edge (3d)
	if (val.size() == value.m_ncomp) {
	    for (unsigned int i = 0; i < value.m_ncomp; i++) {
		value.m_bc_val[Face::ZHI].push_back(val[i]);
	    }
	} else if (val.size() == 1) {
	    value.m_bc_val[Face::ZHI].resize(value.m_ncomp, val[0]);
	} else if (val.size() == 0) {
	    value.m_bc_val[Face::ZHI].resize(value.m_ncomp, 0);
	} else {
	    Util::Exception(INFO, "Incorrect number of ", pp.prefix(), " BC value args: expected ", value.m_ncomp, " or 0 or 1 but got ", val.size());
	}
#endif

	    }
	};

} // namespace BC

#endif // BC_INTBC_H_

