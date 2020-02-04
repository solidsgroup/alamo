#ifndef BC_CONSTANT_H_
#define BC_CONSTANT_H_

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCRec.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Array.H>

#include "BC/BC.H"
#include "Numeric/Interpolator/Linear.H"

namespace BC
{
class Constant
	: public BC
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
	Constant (int ncomp, amrex::Vector<std::string> bc_hi_str,
		  amrex::Vector<std::string> bc_lo_str,
		  AMREX_D_DECL(amrex::Vector<amrex::Real> _bc_lo_1,
			       amrex::Vector<amrex::Real> _bc_lo_2,
			       amrex::Vector<amrex::Real> _bc_lo_3),
		  AMREX_D_DECL(amrex::Vector<amrex::Real> _bc_hi_1,
			       amrex::Vector<amrex::Real> _bc_hi_2,
			       amrex::Vector<amrex::Real> _bc_hi_3));

	virtual ~Constant() {};

	virtual void FillBoundary (amrex::FArrayBox &in, const amrex::Box &box,
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

	unsigned int m_ncomp;

	//int bc_lo[BL_SPACEDIM];
	//int bc_hi[BL_SPACEDIM];
	//amrex::Vector<amrex::Real> AMREX_D_DECL(bc_lo_1, bc_lo_2, bc_lo_3);
	//amrex::Vector<amrex::Real> AMREX_D_DECL(bc_hi_1, bc_hi_2, bc_hi_3);

	std::array<std::vector<int>,                                        m_nfaces> m_bc_type; 
	std::array<std::vector<Numeric::Interpolator::Linear<Set::Scalar>>, m_nfaces> m_bc_val; 

public:
    static void Parse(Constant & value, amrex::ParmParse & pp)
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
		bcmap["neumann"]          = (int)amrex::LinOpBCType::Neumann;
		bcmap["reflect_odd"]      = (int)amrex::LinOpBCType::reflect_odd;
		bcmap["Marshak"]          = (int)amrex::LinOpBCType::Marshak;
		bcmap["SanchezPomraning"] = (int)amrex::LinOpBCType::SanchezPomraning;
		bcmap["inflow"]           = (int)amrex::LinOpBCType::inflow;
		bcmap["Periodic"]         = (int)amrex::LinOpBCType::Periodic;
		bcmap["periodic"]         = (int)amrex::LinOpBCType::Periodic;


		std::vector<std::string> str;
		pp.queryarr("type.xlo",str);
		Util::Assert(INFO,TEST(str.size() > 0));
		value.m_ncomp = str.size();
		
		// TYPES
		
		pp.queryarr("type.xlo",str); 
			if (str.size() == value.m_ncomp) for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XLO][i] = bcmap[str[i]]; 
			else if (str.size() == 1) value.m_bc_type[Face::XLO].resize(value.m_ncomp,bcmap[str[0]]);
			else Util::Abort(INFO,"Incorrect number of BC args");
		pp.queryarr("type.xhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XHI][i] = bcmap[str[i]];
			if (str.size() == value.m_ncomp) for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XHI][i] = bcmap[str[i]]; 
			else if (str.size() == 1) value.m_bc_type[Face::XHI].resize(value.m_ncomp,bcmap[str[0]]);
			else Util::Abort(INFO,"Incorrect number of BC args");
		pp.queryarr("type.ylo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::YLO][i] = bcmap[str[i]];
			if (str.size() == value.m_ncomp) for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::YLO][i] = bcmap[str[i]]; 
			else if (str.size() == 1) value.m_bc_type[Face::YLO].resize(value.m_ncomp,bcmap[str[0]]);
			else Util::Abort(INFO,"Incorrect number of BC args");
		pp.queryarr("type.yhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::YHI][i] = bcmap[str[i]];
			if (str.size() == value.m_ncomp) for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::YHI][i] = bcmap[str[i]]; 
			else if (str.size() == 1) value.m_bc_type[Face::YHI].resize(value.m_ncomp,bcmap[str[0]]);
			else Util::Abort(INFO,"Incorrect number of BC args");
		#if AMREX_SPACEDIM==3
		pp.queryarr("type.zlo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::ZLO][i] = bcmap[str[i]];
			if (str.size() == value.m_ncomp) for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::ZLO][i] = bcmap[str[i]]; 
			else if (str.size() == 1) value.m_bc_type[Face::ZLO].resize(value.m_ncomp,bcmap[str[0]]);
			else Util::Abort(INFO,"Incorrect number of BC args");
		pp.queryarr("type.zhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::ZHI][i] = bcmap[str[i]];
			if (str.size() == value.m_ncomp) for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::ZHI][i] = bcmap[str[i]]; 
			else if (str.size() == 1) value.m_bc_type[Face::ZHI].resize(value.m_ncomp,bcmap[str[0]]);
			else Util::Abort(INFO,"Incorrect number of BC args");
		#endif

		// VALS
		std::vector<std::string> val(0);
		if (pp.contains("val.xlo")) { pp.queryarr("val.xlo",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XLO][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); }
		Util::Assert(INFO,TEST(value.m_bc_val[Face::XLO].size() == value.m_ncomp || value.m_bc_val[Face::XLO].size() == 0));
		if (pp.contains("val.xhi")) { pp.queryarr("val.xhi",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XHI][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); }
		Util::Assert(INFO,TEST(value.m_bc_val[Face::XHI].size() == value.m_ncomp || value.m_bc_val[Face::XHI].size() == 0));
		if (pp.contains("val.ylo")) { pp.queryarr("val.ylo",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::YLO][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); }
		Util::Assert(INFO,TEST(value.m_bc_val[Face::YLO].size() == value.m_ncomp || value.m_bc_val[Face::YLO].size() == 0));
		if (pp.contains("val.yhi")) { pp.queryarr("val.yhi",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::YHI][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); }
		Util::Assert(INFO,TEST(value.m_bc_val[Face::YHI].size() == value.m_ncomp || value.m_bc_val[Face::YHI].size() == 0));
		#if AMREX_SPACEDIM==3
		if (pp.contains("val.zlo")) { pp.queryarr("val.zlo",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::ZLO][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); }
		Util::Assert(INFO,TEST(value.m_bc_val[Face::ZLO].size() == value.m_ncomp || value.m_bc_val[Face::ZLO].size() == 0));
		if (pp.contains("val.zhi")) { pp.queryarr("val.zhi",val); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::ZHI][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); }
		Util::Assert(INFO,TEST(value.m_bc_val[Face::ZHI].size() == value.m_ncomp || value.m_bc_val[Face::ZHI].size() == 0));
		#endif
    }

};
}
#endif
