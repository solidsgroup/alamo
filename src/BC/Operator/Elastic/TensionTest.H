#ifndef BC_OPERATOR_ELASTIC_TENSIONTEST_H
#define BC_OPERATOR_ELASTIC_TENSIONTEST_H

// #include "Operator/Elastic.H"
#include "Constant.H"

namespace BC
{
namespace Operator
{
namespace Elastic
{
class TensionTest : public Constant
{
public:

    TensionTest(): Constant()
    {};

    TensionTest(IO::ParmParse &pp, std::string name) : TensionTest()
    { pp.queryclass(name,*this); }


    ~TensionTest() {};

    using Constant::Init;
    using Constant::operator();
    using Constant::set;

    enum LoadingMode {TRAC,DISP};

public:
    static void Parse(TensionTest & value, IO::ParmParse & pp)
    {
        std::string type = "uniaxial_stress_clamp";
        std::string disp = "0.0", trac = "0.0";
        LoadingMode mode;

        pp.query("type",type); // Tension test type. (Only option right now is :code:`uniaxial_stress_clamp`)
        
        if (pp.contains("trac") && pp.contains("disp")) 
            Util::Abort(INFO,"Cannot specify both trac and disp");
        else if (pp.contains("disp"))
        {
            pp.query("disp",disp); // Displacement of the test
            mode = LoadingMode::DISP;
        }
        else if (pp.contains("trac"))
        {
            pp.query("trac",trac);
            mode = LoadingMode::TRAC;
        }
        else
            Util::Abort(INFO,"You must specify either trac or disp");

        if (type == "uniaxial_stress_clamp")
        {
            // Add a X displacement in the X direction
            value.m_bc_val[Face::XHI][0]         = Numeric::Interpolator::Linear<Set::Scalar>(disp);
            value.m_bc_val[Face::XHI_YLO][0]     = Numeric::Interpolator::Linear<Set::Scalar>(disp);
            value.m_bc_val[Face::XHI_YHI][0]     = Numeric::Interpolator::Linear<Set::Scalar>(disp);
            #if AMREX_SPACEDIM==3
            value.m_bc_val[Face::ZLO_XHI][0]     = Numeric::Interpolator::Linear<Set::Scalar>(disp);
            value.m_bc_val[Face::ZHI_XHI][0]     = Numeric::Interpolator::Linear<Set::Scalar>(disp);
            value.m_bc_val[Face::XHI_YLO_ZLO][0] = Numeric::Interpolator::Linear<Set::Scalar>(disp);
            value.m_bc_val[Face::XHI_YLO_ZHI][0] = Numeric::Interpolator::Linear<Set::Scalar>(disp);
            value.m_bc_val[Face::XHI_YHI_ZLO][0] = Numeric::Interpolator::Linear<Set::Scalar>(disp);
            value.m_bc_val[Face::XHI_YHI_ZHI][0] = Numeric::Interpolator::Linear<Set::Scalar>(disp);
            #endif

            // Set the Y faces, Z faces, and YZ edges to be traction-free
            for (int i = 0 ; i < AMREX_SPACEDIM ; i++)
            {
                value.m_bc_type[Face::YLO][i]     = Type::Traction;                
                value.m_bc_type[Face::YHI][i]     = Type::Traction;                
                #if AMREX_SPACEDIM==3
                value.m_bc_type[Face::ZLO][i]     = Type::Traction;                
                value.m_bc_type[Face::ZHI][i]     = Type::Traction;                
                value.m_bc_type[Face::YLO_ZLO][i] = Type::Traction;                
                value.m_bc_type[Face::YHI_ZLO][i] = Type::Traction;                
                value.m_bc_type[Face::YLO_ZHI][i] = Type::Traction;                
                value.m_bc_type[Face::YHI_ZHI][i] = Type::Traction;                
                #endif
            }
        }
        else if (type == "uniaxial_stress")
        {
            // Go ahead and set everything to traction first.
            for (int face = 0; face < m_nfaces; face++)
                for (int i = 0; i < AMREX_SPACEDIM; i++)
                    value.m_bc_type[face][i] = Type::Traction;
                    

            // Fix the X displacement in the XMIN
            value.m_bc_type[Face::XLO][0]         = Type::Displacement;
            value.m_bc_type[Face::XLO_YLO][0]     = Type::Displacement;
            value.m_bc_type[Face::XLO_YHI][0]     = Type::Displacement;
            #if AMREX_SPACEDIM == 3
            value.m_bc_type[Face::ZLO_XLO][0]     = Type::Displacement;
            value.m_bc_type[Face::ZHI_XLO][0]     = Type::Displacement;
            value.m_bc_type[Face::XLO_YLO_ZLO][0] = Type::Displacement;
            value.m_bc_type[Face::XLO_YLO_ZHI][0] = Type::Displacement;
            value.m_bc_type[Face::XLO_YHI_ZLO][0] = Type::Displacement;
            value.m_bc_type[Face::XLO_YHI_ZHI][0] = Type::Displacement;
            #endif

            // Add a X displacement in the X direction (same as above)
            if (mode == TensionTest::LoadingMode::DISP)
            {
                value.m_bc_type[Face::XHI][0]         = Type::Displacement;
                value.m_bc_val [Face::XHI][0]         = Numeric::Interpolator::Linear<Set::Scalar>(disp);
                value.m_bc_type[Face::XHI_YLO][0]     = Type::Displacement;
                value.m_bc_val [Face::XHI_YLO][0]     = Numeric::Interpolator::Linear<Set::Scalar>(disp);
                value.m_bc_type[Face::XHI_YHI][0]     = Type::Displacement;
                value.m_bc_val [Face::XHI_YHI][0]     = Numeric::Interpolator::Linear<Set::Scalar>(disp);
            }
            else if (mode == TensionTest::LoadingMode::TRAC)
            {
                value.m_bc_type[Face::XHI][0]         = Type::Traction;
                value.m_bc_val [Face::XHI][0]         = Numeric::Interpolator::Linear<Set::Scalar>(trac);
                value.m_bc_type[Face::XHI_YLO][0]     = Type::Traction;
                value.m_bc_val [Face::XHI_YLO][0]     = Numeric::Interpolator::Linear<Set::Scalar>(trac);
                value.m_bc_type[Face::XHI_YHI][0]     = Type::Traction;
                value.m_bc_val [Face::XHI_YHI][0]     = Numeric::Interpolator::Linear<Set::Scalar>(trac);
            }
            #if AMREX_SPACEDIM == 3
            if (mode == TensionTest::LoadingMode::DISP)
            {
                value.m_bc_type[Face::ZLO_XHI][0]     = Type::Displacement;
                value.m_bc_val [Face::ZLO_XHI][0]     = Numeric::Interpolator::Linear<Set::Scalar>(disp);
                value.m_bc_type[Face::ZHI_XHI][0]     = Type::Displacement;
                value.m_bc_val [Face::ZHI_XHI][0]     = Numeric::Interpolator::Linear<Set::Scalar>(disp);
                value.m_bc_type[Face::XHI_YLO_ZLO][0] = Type::Displacement;
                value.m_bc_val [Face::XHI_YLO_ZLO][0] = Numeric::Interpolator::Linear<Set::Scalar>(disp);
                value.m_bc_type[Face::XHI_YLO_ZHI][0] = Type::Displacement;
                value.m_bc_val [Face::XHI_YLO_ZHI][0] = Numeric::Interpolator::Linear<Set::Scalar>(disp);
                value.m_bc_type[Face::XHI_YHI_ZLO][0] = Type::Displacement;
                value.m_bc_val [Face::XHI_YHI_ZLO][0] = Numeric::Interpolator::Linear<Set::Scalar>(disp);
                value.m_bc_type[Face::XHI_YHI_ZHI][0] = Type::Displacement;
                value.m_bc_val [Face::XHI_YHI_ZHI][0] = Numeric::Interpolator::Linear<Set::Scalar>(disp);
            }
            if (mode == TensionTest::LoadingMode::TRAC)
            {
                value.m_bc_type[Face::ZLO_XHI][0]     = Type::Traction;
                value.m_bc_val [Face::ZLO_XHI][0]     = Numeric::Interpolator::Linear<Set::Scalar>(trac);
                value.m_bc_type[Face::ZHI_XHI][0]     = Type::Traction;
                value.m_bc_val [Face::ZHI_XHI][0]     = Numeric::Interpolator::Linear<Set::Scalar>(trac);
                value.m_bc_type[Face::XHI_YLO_ZLO][0] = Type::Traction;
                value.m_bc_val [Face::XHI_YLO_ZLO][0] = Numeric::Interpolator::Linear<Set::Scalar>(trac);
                value.m_bc_type[Face::XHI_YLO_ZHI][0] = Type::Traction;
                value.m_bc_val [Face::XHI_YLO_ZHI][0] = Numeric::Interpolator::Linear<Set::Scalar>(trac);
                value.m_bc_type[Face::XHI_YHI_ZLO][0] = Type::Traction;
                value.m_bc_val [Face::XHI_YHI_ZLO][0] = Numeric::Interpolator::Linear<Set::Scalar>(trac);
                value.m_bc_type[Face::XHI_YHI_ZHI][0] = Type::Traction;
                value.m_bc_val [Face::XHI_YHI_ZHI][0] = Numeric::Interpolator::Linear<Set::Scalar>(trac);
            }
            #endif

            // Fix Ymin displacement
            value.m_bc_type[Face::YLO][1] = Type::Displacement;
            value.m_bc_type[Face::XLO_YLO][1] = Type::Displacement;
            value.m_bc_type[Face::XHI_YLO][1] = Type::Displacement;
            #if AMREX_SPACEDIM == 3
            value.m_bc_type[Face::YLO_ZLO][1] = Type::Displacement;
            value.m_bc_type[Face::YLO_ZHI][1] = Type::Displacement;
            value.m_bc_type[Face::XLO_YLO_ZLO][1] = Type::Displacement;
            value.m_bc_type[Face::XLO_YLO_ZHI][1] = Type::Displacement;
            value.m_bc_type[Face::XHI_YLO_ZLO][1] = Type::Displacement;
            value.m_bc_type[Face::XHI_YLO_ZHI][1] = Type::Displacement;
            #endif

            #if AMREX_SPACEDIM == 3
            // Fix Zmin displacement
            value.m_bc_type[Face::ZLO][2] = Type::Displacement;
            value.m_bc_type[Face::YLO_ZLO][2] = Type::Displacement;
            value.m_bc_type[Face::YHI_ZLO][2] = Type::Displacement;
            value.m_bc_type[Face::ZLO_XLO][2] = Type::Displacement;
            value.m_bc_type[Face::ZLO_XHI][2] = Type::Displacement;
            value.m_bc_type[Face::XLO_YLO_ZLO][2] = Type::Displacement;
            value.m_bc_type[Face::XLO_YHI_ZLO][2] = Type::Displacement;
            value.m_bc_type[Face::XHI_YLO_ZLO][2] = Type::Displacement;
            value.m_bc_type[Face::XHI_YHI_ZLO][2] = Type::Displacement;
            #endif
            
        }
        else if (type == "uniaxial_strain")
        {
            // Go ahead and set everything to traction first.
            for (int face = 0; face < m_nfaces; face++)
                for (int i = 0; i < AMREX_SPACEDIM; i++)
                    value.m_bc_type[face][i] = Type::Traction;
                    

            // Fix the X displacement in the XMIN
            value.m_bc_type[Face::XLO][0]         = Type::Displacement;
            value.m_bc_type[Face::XLO_YLO][0]     = Type::Displacement;
            value.m_bc_type[Face::XLO_YHI][0]     = Type::Displacement;
            #if AMREX_SPACEDIM == 3
            value.m_bc_type[Face::ZLO_XLO][0]     = Type::Displacement;
            value.m_bc_type[Face::ZHI_XLO][0]     = Type::Displacement;
            value.m_bc_type[Face::XLO_YLO_ZLO][0] = Type::Displacement;
            value.m_bc_type[Face::XLO_YLO_ZHI][0] = Type::Displacement;
            value.m_bc_type[Face::XLO_YHI_ZLO][0] = Type::Displacement;
            value.m_bc_type[Face::XLO_YHI_ZHI][0] = Type::Displacement;
            #endif

            // Add a X displacement in the X direction (same as above)
            value.m_bc_type[Face::XHI][0]         = Type::Displacement;
            value.m_bc_val [Face::XHI][0]         = Numeric::Interpolator::Linear<Set::Scalar>(disp);
            value.m_bc_type[Face::XHI_YLO][0]     = Type::Displacement;
            value.m_bc_val [Face::XHI_YLO][0]     = Numeric::Interpolator::Linear<Set::Scalar>(disp);
            value.m_bc_type[Face::XHI_YHI][0]     = Type::Displacement;
            value.m_bc_val [Face::XHI_YHI][0]     = Numeric::Interpolator::Linear<Set::Scalar>(disp);
            #if AMREX_SPACEDIM == 3
            value.m_bc_type[Face::ZLO_XHI][0]     = Type::Displacement;
            value.m_bc_val [Face::ZLO_XHI][0]     = Numeric::Interpolator::Linear<Set::Scalar>(disp);
            value.m_bc_type[Face::ZHI_XHI][0]     = Type::Displacement;
            value.m_bc_val [Face::ZHI_XHI][0]     = Numeric::Interpolator::Linear<Set::Scalar>(disp);
            value.m_bc_type[Face::XHI_YLO_ZLO][0] = Type::Displacement;
            value.m_bc_val [Face::XHI_YLO_ZLO][0] = Numeric::Interpolator::Linear<Set::Scalar>(disp);
            value.m_bc_type[Face::XHI_YLO_ZHI][0] = Type::Displacement;
            value.m_bc_val [Face::XHI_YLO_ZHI][0] = Numeric::Interpolator::Linear<Set::Scalar>(disp);
            value.m_bc_type[Face::XHI_YHI_ZLO][0] = Type::Displacement;
            value.m_bc_val [Face::XHI_YHI_ZLO][0] = Numeric::Interpolator::Linear<Set::Scalar>(disp);
            value.m_bc_type[Face::XHI_YHI_ZHI][0] = Type::Displacement;
            value.m_bc_val [Face::XHI_YHI_ZHI][0] = Numeric::Interpolator::Linear<Set::Scalar>(disp);
            #endif

            // Fix Ymin displacement
            value.m_bc_type[Face::YLO][1] = Type::Displacement;
            value.m_bc_type[Face::XLO_YLO][1] = Type::Displacement;
            value.m_bc_type[Face::XHI_YLO][1] = Type::Displacement;
            value.m_bc_type[Face::YHI][1] = Type::Displacement;
            value.m_bc_type[Face::XLO_YHI][1] = Type::Displacement;
            value.m_bc_type[Face::XHI_YHI][1] = Type::Displacement;
            #if AMREX_SPACEDIM == 3
            value.m_bc_type[Face::YLO_ZLO][1] = Type::Displacement;
            value.m_bc_type[Face::YLO_ZHI][1] = Type::Displacement;
            value.m_bc_type[Face::YHI_ZLO][1] = Type::Displacement;
            value.m_bc_type[Face::YHI_ZHI][1] = Type::Displacement;
            value.m_bc_type[Face::XLO_YLO_ZLO][1] = Type::Displacement;
            value.m_bc_type[Face::XLO_YLO_ZHI][1] = Type::Displacement;
            value.m_bc_type[Face::XHI_YLO_ZLO][1] = Type::Displacement;
            value.m_bc_type[Face::XHI_YLO_ZHI][1] = Type::Displacement;
            value.m_bc_type[Face::XLO_YHI_ZLO][1] = Type::Displacement;
            value.m_bc_type[Face::XLO_YHI_ZHI][1] = Type::Displacement;
            value.m_bc_type[Face::XHI_YHI_ZLO][1] = Type::Displacement;
            value.m_bc_type[Face::XHI_YHI_ZHI][1] = Type::Displacement;
            #endif

            #if AMREX_SPACEDIM == 3
            value.m_bc_type[Face::ZLO][2] = Type::Displacement;
            value.m_bc_type[Face::ZHI][2] = Type::Displacement;
            value.m_bc_type[Face::YLO_ZLO][2] = Type::Displacement;
            value.m_bc_type[Face::YHI_ZLO][2] = Type::Displacement;
            value.m_bc_type[Face::ZLO_XLO][2] = Type::Displacement;
            value.m_bc_type[Face::ZLO_XHI][2] = Type::Displacement;
            value.m_bc_type[Face::YLO_ZHI][2] = Type::Displacement;
            value.m_bc_type[Face::YHI_ZHI][2] = Type::Displacement;
            value.m_bc_type[Face::ZHI_XLO][2] = Type::Displacement;
            value.m_bc_type[Face::ZHI_XHI][2] = Type::Displacement;
            value.m_bc_type[Face::XLO_YLO_ZLO][2] = Type::Displacement;
            value.m_bc_type[Face::XLO_YHI_ZLO][2] = Type::Displacement;
            value.m_bc_type[Face::XHI_YLO_ZLO][2] = Type::Displacement;
            value.m_bc_type[Face::XHI_YHI_ZLO][2] = Type::Displacement;
            value.m_bc_type[Face::XLO_YLO_ZHI][2] = Type::Displacement;
            value.m_bc_type[Face::XLO_YHI_ZHI][2] = Type::Displacement;
            value.m_bc_type[Face::XHI_YLO_ZHI][2] = Type::Displacement;
            value.m_bc_type[Face::XHI_YHI_ZHI][2] = Type::Displacement;
            #endif
            
        }
        else
        {
            Util::Abort(INFO,"Unknown BC type ", type);
        }
    }
};
}
}
}
#endif
