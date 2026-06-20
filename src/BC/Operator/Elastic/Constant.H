//
// This BC defines boundary conditions that are constant with respect to space.
// However they may change in time using the :ref:`Numeric::Interpolator::Linear` method.
// 
// In 2D, BCs are prescribed for each edge (4) and corner (4) on the boundary, for a total of 8 possible regions.
// In 3D, BCs are prescribed for each face (6), each edge (12), and each corner (8), for a total of 26 possible regions.
// Care must be taken to define BCs for edges/corners consistently with the faces/edges, or you will get poor 
// convergence and inaccurate behavior.
// (See :ref:`BC::Operator::Elastic::TensionTest` for a reduced BC with streamlined load options.)
//
// To define a boundary condition, you must define both a "type" and a "value" in each direction.
// The type can be "displacement", "disp", "neumann", "traction", "trac".
// The value is the corresponding value.
// All BCs are, by default, dirichlet (displacement) with a value of zero.
//

#ifndef BC_OPERATOR_ELASTIC_CONSTANT_H
#define BC_OPERATOR_ELASTIC_CONSTANT_H

// #include "Operator/Elastic.H"
#include "AMReX_Loop.H"
#include "IO/ParmParse.H"
#include "Numeric/Interpolator/Linear.H"
#include "BC/Operator/Elastic/Elastic.H"

#ifdef ALAMO_GPU
#define ALAMO_ELASTIC_BC_FOR amrex::ParallelFor
#define ALAMO_ELASTIC_BC_DEVICE AMREX_GPU_DEVICE
#define ALAMO_ELASTIC_BC_TYPE(face, dir) bc_type[face][dir]
#define ALAMO_ELASTIC_BC_VAL(face, dir) bc_val[face][dir]
#define ALAMO_ELASTIC_BC_PERIODIC(dir) periodic[dir]
#else
#define ALAMO_ELASTIC_BC_FOR amrex::LoopConcurrentOnCpu
#define ALAMO_ELASTIC_BC_DEVICE
#define ALAMO_ELASTIC_BC_TYPE(face, dir) m_bc_type[face][dir]
#define ALAMO_ELASTIC_BC_VAL(face, dir) m_bc_val[face][dir](m_time)
#define ALAMO_ELASTIC_BC_PERIODIC(dir) a_geom.isPeriodic(dir)
#endif

namespace BC
{
namespace Operator
{
namespace Elastic
{
class Constant : public Elastic
{
public:
    static constexpr const char* name = "constant";

    Constant()
    {
        // By default, all boundary conditions are displacement
        for (int face = 0; face < m_nfaces; face++)
            for (int direction = 0; direction < AMREX_SPACEDIM; direction++)
            {
                m_bc_type[face][direction] = Type::Displacement;
                m_bc_val [face][direction] = 0.0;
            }
    };

    Constant(IO::ParmParse &pp, std::string name): Constant()
    { pp_queryclass(name,*this); }

    ~Constant() {};


    void
    Set(const Face face,
        const Direction direction,
        const Type type,
        const Set::Scalar value)
    {
        m_bc_type[face][direction] = type;
        m_bc_val [face][direction] = value;
    }

    void
    Set(const Face /*face*/,
        const Direction /*direction*/,
        const Type /*type*/,
        const Set::Scalar /*value*/,
        amrex::Vector<amrex::MultiFab *> &/*a_rhs*/,
        const amrex::Vector<amrex::Geometry> &/*a_geom*/)
    {
        Util::Abort(INFO,"This has been depricated");
    }

    void
    Set(const Face face,
        const Direction direction,
        const Type type,
        const Set::Scalar value,
        amrex::Vector<amrex::MultiFab> &a_rhs,
        const amrex::Vector<amrex::Geometry> &a_geom)
    {
        amrex::Vector<amrex::MultiFab *> pa_rhs = amrex::GetVecOfPtrs(a_rhs);
        Set(face,direction,type,value,pa_rhs,a_geom);
    }

    void
    Set(const Face face,
        const Direction direction,
        const Type type,
        const Set::Scalar value,
        amrex::Vector<std::unique_ptr<amrex::MultiFab> > &a_rhs,
        const amrex::Vector<amrex::Geometry> &a_geom)
    {
        amrex::Vector<amrex::MultiFab *> pa_rhs = amrex::GetVecOfPtrs(a_rhs);
        Set(face,direction,type,value,pa_rhs,a_geom);
    }

    //using Elastic::Init;
    
    virtual void
    Init(amrex::FabArray<amrex::BaseFab<Set::Vector>> *a_rhs,
        const amrex::Geometry &a_geom,
        bool a_homogeneous = false) const override
    {
            amrex::Box domain(a_geom.growPeriodicDomain(2));
            domain.convert(amrex::IntVect::TheNodeVector());
            const amrex::Dim3 lo= amrex::lbound(domain), hi = amrex::ubound(domain);
#ifdef ALAMO_GPU
            amrex::GpuArray<amrex::GpuArray<Type, AMREX_SPACEDIM>, m_nfaces> bc_type;
            amrex::GpuArray<amrex::GpuArray<Set::Scalar, AMREX_SPACEDIM>, m_nfaces> bc_val;
            for (int face = 0; face < m_nfaces; ++face)
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
                {
                    bc_type[face][dir] = m_bc_type[face][dir];
                    bc_val[face][dir] = m_bc_val[face][dir](m_time);
                }
            const amrex::GpuArray<int, AMREX_SPACEDIM> periodic = a_geom.isPeriodicArray();
#endif
            for (amrex::MFIter mfi(*a_rhs, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.grownnodaltilebox();
                bx = bx & domain;
                amrex::Array4<Set::Vector> const& rhs       = a_rhs->array(mfi);
                ALAMO_ELASTIC_BC_FOR(bx, [=] ALAMO_ELASTIC_BC_DEVICE (int i, int j, int k) noexcept {
                    
                for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
                {
                    Face face = Face::INT;
                        
                    #if AMREX_SPACEDIM == 2
                        
                    if      (i==lo.x && j==lo.y) face = Face::XLO_YLO;
                    else if (i==lo.x && j==hi.y) face = Face::XLO_YHI;
                    else if (i==hi.x && j==lo.y) face = Face::XHI_YLO;
                    else if (i==hi.x && j==hi.y) face = Face::XHI_YHI;

                    else if (i==lo.x) face = Face::XLO;
                    else if (i==hi.x) face = Face::XHI;
                    else if (j==lo.y) face = Face::YLO;
                    else if (j==hi.y) face = Face::YHI;

                    #elif AMREX_SPACEDIM == 3

                    if      (i==lo.x && j==lo.y && k==lo.z) face = Face::XLO_YLO_ZLO;
                    else if (i==lo.x && j==lo.y && k==hi.z) face = Face::XLO_YLO_ZHI;
                    else if (i==lo.x && j==hi.y && k==lo.z) face = Face::XLO_YHI_ZLO;
                    else if (i==lo.x && j==hi.y && k==hi.z) face = Face::XLO_YHI_ZHI;
                    else if (i==hi.x && j==lo.y && k==lo.z) face = Face::XHI_YLO_ZLO;
                    else if (i==hi.x && j==lo.y && k==hi.z) face = Face::XHI_YLO_ZHI;
                    else if (i==hi.x && j==hi.y && k==lo.z) face = Face::XHI_YHI_ZLO;
                    else if (i==hi.x && j==hi.y && k==hi.z) face = Face::XHI_YHI_ZHI;

                    else if (j==lo.y && k==lo.z) face = Face::YLO_ZLO;
                    else if (j==lo.y && k==hi.z) face = Face::YLO_ZHI;
                    else if (j==hi.y && k==lo.z) face = Face::YHI_ZLO;
                    else if (j==hi.y && k==hi.z) face = Face::YHI_ZHI;
                    else if (k==lo.z && i==lo.x) face = Face::ZLO_XLO;
                    else if (k==lo.z && i==hi.x) face = Face::ZLO_XHI;
                    else if (k==hi.z && i==lo.x) face = Face::ZHI_XLO;
                    else if (k==hi.z && i==hi.x) face = Face::ZHI_XHI;
                    else if (i==lo.x && j==lo.y) face = Face::XLO_YLO;
                    else if (i==lo.x && j==hi.y) face = Face::XLO_YHI;
                    else if (i==hi.x && j==lo.y) face = Face::XHI_YLO;
                    else if (i==hi.x && j==hi.y) face = Face::XHI_YHI;

                    else if (i==lo.x) face = Face::XLO;
                    else if (i==hi.x) face = Face::XHI;
                    else if (j==lo.y) face = Face::YLO;
                    else if (j==hi.y) face = Face::YHI;
                    else if (k==lo.z) face = Face::ZLO;
                    else if (k==hi.z) face = Face::ZHI;

                    #endif

                    if (ALAMO_ELASTIC_BC_PERIODIC(0) && ALAMO_ELASTIC_BC_PERIODIC(1))
                    {
                        if (face == Face::XLO || face == Face::XHI ||
                            face == Face::YLO || face == Face::YHI) continue;
                    }
                    if (ALAMO_ELASTIC_BC_PERIODIC(0))
                    {
                        if (face == Face::XLO || face==Face::XHI) continue;
                    }
                    if (ALAMO_ELASTIC_BC_PERIODIC(1))
                    {
                        if (face == Face::YLO || face==Face::YHI) continue;
                    }
                    #if AMREX_SPACEDIM == 3
                    if (ALAMO_ELASTIC_BC_PERIODIC(2))
                    {
                        if (face == Face::ZLO || face==Face::ZHI) continue;
                    }
                    #endif

                    if (!(face == Face::INT))
                    {
                        if (a_homogeneous && ALAMO_ELASTIC_BC_TYPE(face, dir) == Type::Displacement) 
                            rhs(i,j,k)(dir) = 0.0;
                        else 
                            rhs(i,j,k)(dir) = ALAMO_ELASTIC_BC_VAL(face, dir);
                    }
                }
            });
        }                    
    }

    using Elastic::Init;
    virtual void
    Init(amrex::MultiFab * a_rhs,
        const amrex::Geometry &a_geom,
        bool a_homogeneous = false) const override
    {
            amrex::Box domain(a_geom.Domain());
            domain.convert(amrex::IntVect::TheNodeVector());
            const amrex::Dim3 lo= amrex::lbound(domain), hi = amrex::ubound(domain);
#ifdef ALAMO_GPU
            amrex::GpuArray<amrex::GpuArray<Type, AMREX_SPACEDIM>, m_nfaces> bc_type;
            amrex::GpuArray<amrex::GpuArray<Set::Scalar, AMREX_SPACEDIM>, m_nfaces> bc_val;
            for (int face = 0; face < m_nfaces; ++face)
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
                {
                    bc_type[face][dir] = m_bc_type[face][dir];
                    bc_val[face][dir] = m_bc_val[face][dir](m_time);
                }
            const amrex::GpuArray<int, AMREX_SPACEDIM> periodic = a_geom.isPeriodicArray();
#endif
            for (amrex::MFIter mfi(*a_rhs, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.tilebox();
                bx.grow(2);
                bx = bx & domain;
                amrex::Array4<amrex::Real> const& rhs       = a_rhs->array(mfi);
                ALAMO_ELASTIC_BC_FOR(bx, [=] ALAMO_ELASTIC_BC_DEVICE (int i, int j, int k) noexcept {
                    
                for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
                {
                    Face face = Face::INT;
                        
                    #if AMREX_SPACEDIM == 2
                        
                    if      (i==lo.x && j==lo.y) face = Face::XLO_YLO;
                    else if (i==lo.x && j==hi.y) face = Face::XLO_YHI;
                    else if (i==hi.x && j==lo.y) face = Face::XHI_YLO;
                    else if (i==hi.x && j==hi.y) face = Face::XHI_YHI;

                    else if (i==lo.x) face = Face::XLO;
                    else if (i==hi.x) face = Face::XHI;
                    else if (j==lo.y) face = Face::YLO;
                    else if (j==hi.y) face = Face::YHI;

                    #elif AMREX_SPACEDIM == 3

                    if      (i==lo.x && j==lo.y && k==lo.z) face = Face::XLO_YLO_ZLO;
                    else if (i==lo.x && j==lo.y && k==hi.z) face = Face::XLO_YLO_ZHI;
                    else if (i==lo.x && j==hi.y && k==lo.z) face = Face::XLO_YHI_ZLO;
                    else if (i==lo.x && j==hi.y && k==hi.z) face = Face::XLO_YHI_ZHI;
                    else if (i==hi.x && j==lo.y && k==lo.z) face = Face::XHI_YLO_ZLO;
                    else if (i==hi.x && j==lo.y && k==hi.z) face = Face::XHI_YLO_ZHI;
                    else if (i==hi.x && j==hi.y && k==lo.z) face = Face::XHI_YHI_ZLO;
                    else if (i==hi.x && j==hi.y && k==hi.z) face = Face::XHI_YHI_ZHI;

                    else if (j==lo.y && k==lo.z) face = Face::YLO_ZLO;
                    else if (j==lo.y && k==hi.z) face = Face::YLO_ZHI;
                    else if (j==hi.y && k==lo.z) face = Face::YHI_ZLO;
                    else if (j==hi.y && k==hi.z) face = Face::YHI_ZHI;
                    else if (k==lo.z && i==lo.x) face = Face::ZLO_XLO;
                    else if (k==lo.z && i==hi.x) face = Face::ZLO_XHI;
                    else if (k==hi.z && i==lo.x) face = Face::ZHI_XLO;
                    else if (k==hi.z && i==hi.x) face = Face::ZHI_XHI;
                    else if (i==lo.x && j==lo.y) face = Face::XLO_YLO;
                    else if (i==lo.x && j==hi.y) face = Face::XLO_YHI;
                    else if (i==hi.x && j==lo.y) face = Face::XHI_YLO;
                    else if (i==hi.x && j==hi.y) face = Face::XHI_YHI;

                    else if (i==lo.x) face = Face::XLO;
                    else if (i==hi.x) face = Face::XHI;
                    else if (j==lo.y) face = Face::YLO;
                    else if (j==hi.y) face = Face::YHI;
                    else if (k==lo.z) face = Face::ZLO;
                    else if (k==hi.z) face = Face::ZHI;

                    #endif
                    
                    if (ALAMO_ELASTIC_BC_PERIODIC(0) && ALAMO_ELASTIC_BC_PERIODIC(1))
                    {
                        if (face == Face::XLO || face == Face::XHI ||
                            face == Face::YLO || face == Face::YHI) continue;
                    }
                    if (ALAMO_ELASTIC_BC_PERIODIC(0))
                    {
                        if (face == Face::XLO || face==Face::XHI) continue;
                    }
                    if (ALAMO_ELASTIC_BC_PERIODIC(1))
                    {
                        if (face == Face::YLO || face==Face::YHI) continue;
                    }

                    if (!(face == Face::INT))
                    {
                        if (a_homogeneous && ALAMO_ELASTIC_BC_TYPE(face, dir) == Type::Displacement) 
                            rhs(i,j,k,dir) = 0.0;
                        else 
                            rhs(i,j,k,dir) = ALAMO_ELASTIC_BC_VAL(face, dir);
                    }
                }
            });
        }                    
    }

#ifndef ALAMO_GPU
    AMREX_GPU_HOST_DEVICE
    std::array<Type,AMREX_SPACEDIM> getType (
                const int &i, const int &j, const int &k,
                const amrex::Box &domain) override
    {
        return Elastic::typeOf(m_bc_type, i, j, k, domain);
    }

    AMREX_FORCE_INLINE
    AMREX_GPU_HOST_DEVICE
    Set::Vector operator () (const Set::Vector &u,
                const Set::Matrix &gradu,
                const Set::Matrix &sigma,
                const int &i, const int &j, const int &k,
                const amrex::Box &domain) override
    {
        return Elastic::eval(m_bc_type, u, gradu, sigma, i, j, k, domain);
    }

    AMREX_FORCE_INLINE
    AMREX_GPU_HOST_DEVICE
    Set::Vector set(std::array<Type,AMREX_SPACEDIM> &bc_type,
                    const Set::Vector &u, const Set::Matrix &gradu, const Set::Matrix &sigma, Set::Vector n) const
    {
        return Elastic::set(bc_type, u, gradu, sigma, n);
    }
#endif


protected:
    
    std::array<std::array<Numeric::Interpolator::Linear<Set::Scalar>,AMREX_SPACEDIM>, m_nfaces> m_bc_val; 

public:
    static void Parse(Constant & value, IO::ParmParse & pp)
    {
        std::map<std::string, Type> bcmap;
        bcmap["displacement"] = Type::Displacement;
        bcmap["disp"]         = Type::Displacement;
        bcmap["neumann"]      = Type::Neumann;
        bcmap["traction"]     = Type::Traction;
        bcmap["trac"]         = Type::Traction;
        bcmap["periodic"]     = Type::Periodic;
        
        std::vector<std::string> str;
        
        // TYPES
        
        #if AMREX_SPACEDIM==3
        if (pp.contains("type.xloylozlo")) { 
            pp_queryarr("type.xloylozlo",str); // 3D Corner
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XLO_YLO_ZLO][i] = bcmap[str[i]]; } 
        if (pp.contains("type.xloylozhi")) { 
            pp_queryarr("type.xloylozhi",str); // 3D Corner
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XLO_YLO_ZHI][i] = bcmap[str[i]]; }
        if (pp.contains("type.xloyhizlo")) { 
            pp_queryarr("type.xloyhizlo",str); // 3D Corner
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XLO_YHI_ZLO][i] = bcmap[str[i]]; }
        if (pp.contains("type.xloyhizhi")) { 
            pp_queryarr("type.xloyhizhi",str); // 3D Corner
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XLO_YHI_ZHI][i] = bcmap[str[i]]; }
        if (pp.contains("type.xhiylozlo")) { 
            pp_queryarr("type.xhiylozlo",str); // 3D Corner
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XHI_YLO_ZLO][i] = bcmap[str[i]]; }
        if (pp.contains("type.xhiylozhi")) { 
            pp_queryarr("type.xhiylozhi",str); // 3D Corner
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XHI_YLO_ZHI][i] = bcmap[str[i]]; }
        if (pp.contains("type.xhiyhizlo")) { 
            pp_queryarr("type.xhiyhizlo",str); // 3D Corner
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XHI_YHI_ZLO][i] = bcmap[str[i]]; }
        if (pp.contains("type.xhiyhizhi")) { 
            pp_queryarr("type.xhiyhizhi",str); // 3D Corner
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XHI_YHI_ZHI][i] = bcmap[str[i]]; }
        #else // ignore unused inputs
        pp.ignore("type.xloylozlo"); pp.ignore("type.xloylozhi"); pp.ignore("type.xloyhizlo"); pp.ignore("type.xloyhizhi"); 
        pp.ignore("type.xhiylozlo"); pp.ignore("type.xhiylozhi"); pp.ignore("type.xhiyhizlo"); pp.ignore("type.xhiyhizhi");  
        #endif

        #if AMREX_SPACEDIM==3
        if (pp.contains("type.ylozlo")) { 
            pp_queryarr("type.ylozlo",str);  // 3D Edge
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::YLO_ZLO][i] = bcmap[str[i]]; }
        if (pp.contains("type.ylozhi")) { 
            pp_queryarr("type.ylozhi",str);  // 3D Edge
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::YLO_ZHI][i] = bcmap[str[i]]; }
        if (pp.contains("type.yhizlo")) { 
            pp_queryarr("type.yhizlo",str);  // 3D Edge
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::YHI_ZLO][i] = bcmap[str[i]]; }
        if (pp.contains("type.yhizhi")) { 
            pp_queryarr("type.yhizhi",str);  // 3D Edge
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::YHI_ZHI][i] = bcmap[str[i]]; }
        if (pp.contains("type.zloxlo")) { 
            pp_queryarr("type.zloxlo",str);  // 3D Edge
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::ZLO_XLO][i] = bcmap[str[i]]; }
        if (pp.contains("type.zloxhi")) { 
            pp_queryarr("type.zloxhi",str);  // 3D Edge
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::ZLO_XHI][i] = bcmap[str[i]]; }
        if (pp.contains("type.zhixlo")) { 
            pp_queryarr("type.zhixlo",str);  // 3D Edge
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::ZHI_XLO][i] = bcmap[str[i]]; }
        if (pp.contains("type.zhixhi")) { 
            pp_queryarr("type.zhixhi",str);  // 3D Edge
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::ZHI_XHI][i] = bcmap[str[i]]; }
        #else 
        pp.ignore("type.ylozlo"); pp.ignore("type.ylozhi"); pp.ignore("type.yhizlo"); pp.ignore("type.yhizhi"); 
        pp.ignore("type.zloxlo"); pp.ignore("type.zloxhi"); pp.ignore("type.zhixlo"); pp.ignore("type.zhixhi"); 
        #endif
        
        if (pp.contains("type.xloylo")) { 
            pp_queryarr("type.xloylo",str);  // 3D Edge / 2D Corner
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XLO_YLO][i] = bcmap[str[i]]; }
        if (pp.contains("type.xloyhi")) { 
            pp_queryarr("type.xloyhi",str);  // 3D Edge / 2D Corner
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XLO_YHI][i] = bcmap[str[i]]; }
        if (pp.contains("type.xhiylo")) { 
            pp_queryarr("type.xhiylo",str);  // 3D Edge / 2D Corner
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XHI_YLO][i] = bcmap[str[i]]; }
        if (pp.contains("type.xhiyhi")) { 
            pp_queryarr("type.xhiyhi",str);  // 3D Edge / 2D Corner
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XHI_YHI][i] = bcmap[str[i]]; }

        if (pp.contains("type.xlo")) { 
            pp_queryarr("type.xlo",str); // 3D Face / 2D Edge
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XLO][i] = bcmap[str[i]]; } 
        if (pp.contains("type.xhi")) { 
            pp_queryarr("type.xhi",str); // 3D Face / 2D Edge
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XHI][i] = bcmap[str[i]]; }
        if (pp.contains("type.ylo")) { 
            pp_queryarr("type.ylo",str); // 3D Face / 2D Edge
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::YLO][i] = bcmap[str[i]]; }
        if (pp.contains("type.yhi")) { 
            pp_queryarr("type.yhi",str); // 3D Face / 2D Edge
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::YHI][i] = bcmap[str[i]]; }
        #if AMREX_SPACEDIM==3
        if (pp.contains("type.zlo")) { 
            pp_queryarr("type.zlo",str); // 3D Face
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::ZLO][i] = bcmap[str[i]]; } 
        if (pp.contains("type.zhi")) { 
            pp_queryarr("type.zhi",str); // 3D Face
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::ZHI][i] = bcmap[str[i]]; } 
        #else
        pp.ignore("type.zlo"); pp.ignore("type.zhi");
        #endif

        // VALS
        //std::vector<Set::Scalar> val;
        std::vector<std::string> val;
        
        #if AMREX_SPACEDIM==3
        if (pp.contains("val.xloylozlo")) { 
            pp_queryarr("val.xloylozlo",val); // 3D Corner
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XLO_YLO_ZLO][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); } 
        if (pp.contains("val.xloylozhi")) { 
            pp_queryarr("val.xloylozhi",val); // 3D Corner
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XLO_YLO_ZHI][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); }
        if (pp.contains("val.xloyhizlo")) { 
            pp_queryarr("val.xloyhizlo",val); // 3D Corner
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XLO_YHI_ZLO][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); }
        if (pp.contains("val.xloyhizhi")) { 
            pp_queryarr("val.xloyhizhi",val); // 3D Corner
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XLO_YHI_ZHI][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); }
        if (pp.contains("val.xhiylozlo")) { 
            pp_queryarr("val.xhiylozlo",val); // 3D Corner
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XHI_YLO_ZLO][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); }
        if (pp.contains("val.xhiylozhi")) { 
            pp_queryarr("val.xhiylozhi",val); // 3D Corner
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XHI_YLO_ZHI][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); }
        if (pp.contains("val.xhiyhizlo")) { 
            pp_queryarr("val.xhiyhizlo",val); // 3D Corner
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XHI_YHI_ZLO][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); }
        if (pp.contains("val.xhiyhizhi")) { 
            pp_queryarr("val.xhiyhizhi",val); // 3D Corner
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XHI_YHI_ZHI][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); }
        #else
        pp.ignore("val.xloylozlo"); pp.ignore("val.xloylozhi"); pp.ignore("val.xloyhizlo"); pp.ignore("val.xloyhizhi"); 
        pp.ignore("val.xhiylozlo"); pp.ignore("val.xhiylozhi"); pp.ignore("val.xhiyhizlo"); pp.ignore("val.xhiyhizhi");  
        #endif

        #if AMREX_SPACEDIM==3
        if (pp.contains("val.ylozlo")) { 
            pp_queryarr("val.ylozlo",val); // 3D Edge
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::YLO_ZLO][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); } 
        if (pp.contains("val.ylozhi")) { 
            pp_queryarr("val.ylozhi",val); // 3D Edge
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::YLO_ZHI][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); } 
        if (pp.contains("val.yhizlo")) { 
            pp_queryarr("val.yhizlo",val); // 3D Edge
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::YHI_ZLO][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); } 
        if (pp.contains("val.yhizhi")) { 
            pp_queryarr("val.yhizhi",val); // 3D Edge
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::YHI_ZHI][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); } 
        if (pp.contains("val.zloxlo")) { 
            pp_queryarr("val.zloxlo",val); // 3D Edge
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::ZLO_XLO][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); } 
        if (pp.contains("val.zloxhi")) { 
            pp_queryarr("val.zloxhi",val); // 3D Edge
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::ZLO_XHI][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); } 
        if (pp.contains("val.zhixlo")) { 
            pp_queryarr("val.zhixlo",val); // 3D Edge
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::ZHI_XLO][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); }
        if (pp.contains("val.zhixhi")) { 
            pp_queryarr("val.zhixhi",val); // 3D Edge
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::ZHI_XHI][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); } 
        #else
        pp.ignore("val.ylozlo"); pp.ignore("val.ylozhi"); pp.ignore("val.yhizlo"); pp.ignore("val.yhizhi"); 
        pp.ignore("val.zloxlo"); pp.ignore("val.zloxhi"); pp.ignore("val.zhixlo"); pp.ignore("val.zhixhi"); 
        #endif
        if (pp.contains("val.xloylo")) { 
            pp_queryarr("val.xloylo",val); // 3D Edge / 2D Corner
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XLO_YLO][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); } 
        if (pp.contains("val.xloyhi")) { 
            pp_queryarr("val.xloyhi",val); // 3D Edge / 2D Corner
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XLO_YHI][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); } 
        if (pp.contains("val.xhiylo")) { 
            pp_queryarr("val.xhiylo",val); // 3D Edge / 2D Corner
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XHI_YLO][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); } 
        if (pp.contains("val.xhiyhi")) { 
            pp_queryarr("val.xhiyhi",val); // 3D Edge / 2D Corner
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XHI_YHI][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); } 

        if (pp.contains("val.xlo")) { 
            pp_queryarr("val.xlo",val); // 3D Face / 2D Edge
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XLO][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); } 
        if (pp.contains("val.xhi")) { 
            pp_queryarr("val.xhi",val); // 3D Face / 2D Edge
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::XHI][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); } 
        if (pp.contains("val.ylo")) { 
            pp_queryarr("val.ylo",val); // 3D Face / 2D Edge
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::YLO][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); } 
        if (pp.contains("val.yhi")) { 
            pp_queryarr("val.yhi",val); // 3D Face / 2D Edge
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::YHI][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); } 
        #if AMREX_SPACEDIM==3
        if (pp.contains("val.zlo")) { 
            pp_queryarr("val.zlo",val); // 3D Face
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::ZLO][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); } 
        if (pp.contains("val.zhi")) { 
            pp_queryarr("val.zhi",val); // 3D Face
            for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_val[Face::ZHI][i] = Numeric::Interpolator::Linear<Set::Scalar>(val[i]); } 
        #else
        pp.ignore("val.zlo"); pp.ignore("val.zhi");
        #endif

    }
};
}
}
}
#endif
