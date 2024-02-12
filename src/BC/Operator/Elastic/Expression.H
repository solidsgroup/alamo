#ifndef BC_OPERATOR_ELASTIC_EXPRESSION_H
#define BC_OPERATOR_ELASTIC_EXPRESSION_H

// #include "Operator/Elastic.H"
#include "IO/ParmParse.H"
#include "Numeric/Interpolator/Linear.H"
#include "BC/Operator/Elastic/Elastic.H"

#include "AMReX_Parser.H"

namespace BC
{
namespace Operator
{
namespace Elastic
{
class Expression : public Elastic
{
public:
    //static constexpr const char* const Elastic::strings[];
    #if AMREX_SPACEDIM==2
    static const constexpr char * const facestr[] = {
        "xlo","ylo","xhi","yhi","xloylo","xloyhi","xhiylo","xhiyhi","int"
    };
    #elif AMREX_SPACEDIM==3
    static const constexpr char * const facestr[] = {
        "xlo","ylo","zlo","xhi","yhi","zhi",
        "ylozlo","ylozhi","yhizlo","yhizhi",
        "zloxlo","zloxhi","zhixlo","zhixhi",
        "xloylo","xloyhi","xhiylo","xhiyhi",
        "xloylozlo","xloylozhi","xloyhizlo","xloyhizhi",
        "xhiylozlo","xhiylozhi","xhiyhizlo","xhiyhizhi",
        "int"
    };
    #endif

    Expression()
    {
        // By default, all boundary conditions are displacement
        for (int face = 0; face < m_nfaces; face++)
            for (int direction = 0; direction < AMREX_SPACEDIM; direction++)
            {
                m_bc_type[face][direction] = Type::Displacement;
            }
    };
    ~Expression() {};


    using Elastic::Init;

    virtual void
    Init(amrex::FabArray<amrex::BaseFab<Set::Vector>> *a_rhs,
        const amrex::Geometry &a_geom,
        bool a_homogeneous = false) const override
    {
            amrex::Box domain(a_geom.Domain());
            domain.convert(amrex::IntVect::TheNodeVector());
            const amrex::Dim3 lo= amrex::lbound(domain), hi = amrex::ubound(domain);
            for (amrex::MFIter mfi(*a_rhs, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.tilebox();
                bx.grow(2);
                bx = bx & domain;
                amrex::Array4<Set::Vector> const& rhs       = a_rhs->array(mfi);
                amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    
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

                    amrex::IndexType type = a_rhs->ixType();
                    if (!(face == Face::INT))
                    {
                        if (a_homogeneous && m_bc_type[face][dir] == Type::Displacement) 
                            rhs(i,j,k)(dir) = 0.0;
                        else 
                        {
                            Set::Vector x = Set::Position(i, j, k, a_geom, type);
                            #if AMREX_SPACEDIM==1
                            rhs(i,j,k)(dir) = m_bc_func[face][dir](x(0),0.0,0.0,m_time);
                            #elif AMREX_SPACEDIM==2
                            rhs(i,j,k)(dir) = m_bc_func[face][dir](x(0),x(1),0.0,m_time);
                            #elif AMREX_SPACEDIM==3
                            rhs(i,j,k)(dir) = m_bc_func[face][dir](x(0),x(1),x(2),m_time);
                            #endif
                        }
                    }
                }
            });
        }                    
    }

    
    virtual void
    Init(amrex::MultiFab * a_rhs,
        const amrex::Geometry &a_geom,
        bool a_homogeneous = false) const override
    {
            amrex::Box domain(a_geom.Domain());
            domain.convert(amrex::IntVect::TheNodeVector());
            const amrex::Dim3 lo= amrex::lbound(domain), hi = amrex::ubound(domain);
            for (amrex::MFIter mfi(*a_rhs, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.tilebox();
                bx.grow(2);
                bx = bx & domain;
                amrex::Array4<amrex::Real> const& rhs       = a_rhs->array(mfi);
                amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    
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

                    amrex::IndexType type = a_rhs->ixType();
                    if (!(face == Face::INT))
                    {
                        if (a_homogeneous && m_bc_type[face][dir] == Type::Displacement) 
                            rhs(i,j,k,dir) = 0.0;
                        else 
                        {
                            Set::Vector x = Set::Position(i, j, k, a_geom, type);
                            #if AMREX_SPACEDIM==1
                            rhs(i,j,k,dir) = (m_bc_func[face][dir])(x(0),0.0,0.0,m_time);
                            #elif AMREX_SPACEDIM==2
                            rhs(i,j,k,dir) = (m_bc_func[face][dir])(x(0),x(1),0.0,m_time);
                            #elif AMREX_SPACEDIM==3
                            rhs(i,j,k,dir) = (m_bc_func[face][dir])(x(0),x(1),x(2),m_time);
                            #endif
                        }
                    }
                }
            });
        }                    
    }



    virtual
    std::array<Type,AMREX_SPACEDIM> getType (
                const int &i, const int &j, [[maybe_unused]] const int &k,
                const amrex::Box &domain)
    {
        (void)i, (void)j, (void)k;
        amrex::IntVect m(AMREX_D_DECL(i,j,k));
        const amrex::Dim3 lo = amrex::lbound(domain), hi = amrex::ubound(domain);

        // Corners
        #if AMREX_SPACEDIM == 2
        if (m[0] == lo.x && m[1] == lo.y)                  return m_bc_type[Face::XLO_YLO];
        if (m[0] == lo.x && m[1] == hi.y)                  return m_bc_type[Face::XLO_YHI];
        if (m[0] == hi.x && m[1] == lo.y)                  return m_bc_type[Face::XHI_YLO];
        if (m[0] == hi.x && m[1] == hi.y)                  return m_bc_type[Face::XHI_YHI];
        if (m[0] == lo.x)                                  return m_bc_type[Face::XLO];
        if (m[0] == hi.x)                                  return m_bc_type[Face::XHI];
        if (m[1] == lo.y)                                  return m_bc_type[Face::YLO];
        if (m[1] == hi.y)                                  return m_bc_type[Face::YHI];
        
        #elif AMREX_SPACEDIM == 3
        
        if (m[0] == lo.x &&  m[1] == lo.y && m[2] == lo.z) return m_bc_type[Face::XLO_YLO_ZLO];
        if (m[0] == lo.x &&  m[1] == lo.y && m[2] == hi.z) return m_bc_type[Face::XLO_YLO_ZHI];
        if (m[0] == lo.x &&  m[1] == hi.y && m[2] == lo.z) return m_bc_type[Face::XLO_YHI_ZLO];
        if (m[0] == lo.x &&  m[1] == hi.y && m[2] == hi.z) return m_bc_type[Face::XLO_YHI_ZHI];
        if (m[0] == hi.x &&  m[1] == lo.y && m[2] == lo.z) return m_bc_type[Face::XHI_YLO_ZLO];
        if (m[0] == hi.x &&  m[1] == lo.y && m[2] == hi.z) return m_bc_type[Face::XHI_YLO_ZHI];
        if (m[0] == hi.x &&  m[1] == hi.y && m[2] == lo.z) return m_bc_type[Face::XHI_YHI_ZLO];
        if (m[0] == hi.x &&  m[1] == hi.y && m[2] == hi.z) return m_bc_type[Face::XHI_YHI_ZHI];
        
        if (m[1] == lo.y && m[2] == lo.z)                  return m_bc_type[Face::YLO_ZLO];
        if (m[1] == lo.y && m[2] == hi.z)                  return m_bc_type[Face::YLO_ZHI];
        if (m[1] == hi.y && m[2] == lo.z)                  return m_bc_type[Face::YHI_ZLO];
        if (m[1] == hi.y && m[2] == hi.z)                  return m_bc_type[Face::YHI_ZHI];
        if (m[2] == lo.z && m[0] == lo.x)                  return m_bc_type[Face::ZLO_XLO];
        if (m[2] == lo.z && m[0] == hi.x)                  return m_bc_type[Face::ZLO_XHI];
        if (m[2] == hi.z && m[0] == lo.x)                  return m_bc_type[Face::ZHI_XLO];
        if (m[2] == hi.z && m[0] == hi.x)                  return m_bc_type[Face::ZHI_XHI];
        if (m[0] == lo.x && m[1] == lo.y)                  return m_bc_type[Face::XLO_YLO];
        if (m[0] == lo.x && m[1] == hi.y)                  return m_bc_type[Face::XLO_YHI];
        if (m[0] == hi.x && m[1] == lo.y)                  return m_bc_type[Face::XHI_YLO];
        if (m[0] == hi.x && m[1] == hi.y)                  return m_bc_type[Face::XHI_YHI];
        
        if (m[0] == lo.x)                                  return m_bc_type[Face::XLO];
        if (m[0] == hi.x)                                  return m_bc_type[Face::XHI];
        if (m[1] == lo.y)                                  return m_bc_type[Face::YLO];
        if (m[1] == hi.y)                                  return m_bc_type[Face::YHI];
        if (m[2] == lo.z)                                  return m_bc_type[Face::ZLO];
        if (m[2] == hi.z)                                  return m_bc_type[Face::ZHI];
        
        #endif        
        
        return {AMREX_D_DECL(Type::None,Type::None,Type::None)};
    }
    

    AMREX_FORCE_INLINE
    Set::Vector operator () (const Set::Vector &u,
                const Set::Matrix &gradu,
                const Set::Matrix &sigma,
                const int &i, const int &j, const int &k,
                const amrex::Box &domain) override
    {
        (void)i; (void)j; (void)k; // Suppress "unused variable" warnings
        //Set::Vector f;

        const amrex::Dim3 lo= amrex::lbound(domain), hi = amrex::ubound(domain);
        
        amrex::IntVect m(AMREX_D_DECL(i,j,k));

        // Corners
        #if AMREX_SPACEDIM == 2
        
        if (m[0] == lo.x && m[1] == lo.y)                  return set(m_bc_type[Face::XLO_YLO],     u, gradu, sigma, Set::Vector(-SQRT2INV, -SQRT2INV));
        if (m[0] == lo.x && m[1] == hi.y)                  return set(m_bc_type[Face::XLO_YHI],     u, gradu, sigma, Set::Vector(-SQRT2INV, +SQRT2INV));
        if (m[0] == hi.x && m[1] == lo.y)                  return set(m_bc_type[Face::XHI_YLO],     u, gradu, sigma, Set::Vector(+SQRT2INV, -SQRT2INV));
        if (m[0] == hi.x && m[1] == hi.y)                  return set(m_bc_type[Face::XHI_YHI],     u, gradu, sigma, Set::Vector(+SQRT2INV, +SQRT2INV));
        
        if (m[0] == lo.x)                                  return set(m_bc_type[Face::XLO],         u, gradu, sigma, Set::Vector(-1, 0));
        if (m[0] == hi.x)                                  return set(m_bc_type[Face::XHI],         u, gradu, sigma, Set::Vector(+1, 0));
        if (m[1] == lo.y)                                  return set(m_bc_type[Face::YLO],         u, gradu, sigma, Set::Vector( 0,-1));
        if (m[1] == hi.y)                                  return set(m_bc_type[Face::YHI],         u, gradu, sigma, Set::Vector( 0,+1));
        
        #elif AMREX_SPACEDIM == 3
        
        if (m[0] == lo.x &&  m[1] == lo.y && m[2] == lo.z) return set(m_bc_type[Face::XLO_YLO_ZLO], u, gradu, sigma, Set::Vector(-SQRT3INV,-SQRT3INV,-SQRT3INV));
        if (m[0] == lo.x &&  m[1] == lo.y && m[2] == hi.z) return set(m_bc_type[Face::XLO_YLO_ZHI], u, gradu, sigma, Set::Vector(-SQRT3INV,-SQRT3INV,+SQRT3INV));
        if (m[0] == lo.x &&  m[1] == hi.y && m[2] == lo.z) return set(m_bc_type[Face::XLO_YHI_ZLO], u, gradu, sigma, Set::Vector(-SQRT3INV,+SQRT3INV,-SQRT3INV));
        if (m[0] == lo.x &&  m[1] == hi.y && m[2] == hi.z) return set(m_bc_type[Face::XLO_YHI_ZHI], u, gradu, sigma, Set::Vector(-SQRT3INV,+SQRT3INV,+SQRT3INV));
        if (m[0] == hi.x &&  m[1] == lo.y && m[2] == lo.z) return set(m_bc_type[Face::XHI_YLO_ZLO], u, gradu, sigma, Set::Vector(+SQRT3INV,-SQRT3INV,-SQRT3INV));
        if (m[0] == hi.x &&  m[1] == lo.y && m[2] == hi.z) return set(m_bc_type[Face::XHI_YLO_ZHI], u, gradu, sigma, Set::Vector(+SQRT3INV,-SQRT3INV,+SQRT3INV));
        if (m[0] == hi.x &&  m[1] == hi.y && m[2] == lo.z) return set(m_bc_type[Face::XHI_YHI_ZLO], u, gradu, sigma, Set::Vector(+SQRT3INV,+SQRT3INV,-SQRT3INV));
        if (m[0] == hi.x &&  m[1] == hi.y && m[2] == hi.z) return set(m_bc_type[Face::XHI_YHI_ZHI], u, gradu, sigma, Set::Vector(+SQRT3INV,+SQRT3INV,+SQRT3INV));
        
        if (m[1] == lo.y && m[2] == lo.z)                  return set(m_bc_type[Face::YLO_ZLO],     u, gradu, sigma, Set::Vector(0,        -SQRT2INV,-SQRT2INV));
        if (m[1] == lo.y && m[2] == hi.z)                  return set(m_bc_type[Face::YLO_ZHI],     u, gradu, sigma, Set::Vector(0,        -SQRT2INV,+SQRT2INV));
        if (m[1] == hi.y && m[2] == lo.z)                  return set(m_bc_type[Face::YHI_ZLO],     u, gradu, sigma, Set::Vector(0,        +SQRT2INV,-SQRT2INV));
        if (m[1] == hi.y && m[2] == hi.z)                  return set(m_bc_type[Face::YHI_ZHI],     u, gradu, sigma, Set::Vector(0,        +SQRT2INV,+SQRT2INV));
        if (m[2] == lo.z && m[0] == lo.x)                  return set(m_bc_type[Face::ZLO_XLO],     u, gradu, sigma, Set::Vector(-SQRT2INV, 0,       -SQRT2INV));
        if (m[2] == lo.z && m[0] == hi.x)                  return set(m_bc_type[Face::ZLO_XHI],     u, gradu, sigma, Set::Vector(+SQRT2INV, 0,       -SQRT2INV));
        if (m[2] == hi.z && m[0] == lo.x)                  return set(m_bc_type[Face::ZHI_XLO],     u, gradu, sigma, Set::Vector(-SQRT2INV, 0,       +SQRT2INV));
        if (m[2] == hi.z && m[0] == hi.x)                  return set(m_bc_type[Face::ZHI_XHI],     u, gradu, sigma, Set::Vector(+SQRT2INV, 0,       +SQRT2INV));
        if (m[0] == lo.x && m[1] == lo.y)                  return set(m_bc_type[Face::XLO_YLO],     u, gradu, sigma, Set::Vector(-SQRT2INV, -SQRT2INV,0       ));
        if (m[0] == lo.x && m[1] == hi.y)                  return set(m_bc_type[Face::XLO_YHI],     u, gradu, sigma, Set::Vector(-SQRT2INV, +SQRT2INV,0       ));
        if (m[0] == hi.x && m[1] == lo.y)                  return set(m_bc_type[Face::XHI_YLO],     u, gradu, sigma, Set::Vector(+SQRT2INV, -SQRT2INV,0       ));
        if (m[0] == hi.x && m[1] == hi.y)                  return set(m_bc_type[Face::XHI_YHI],     u, gradu, sigma, Set::Vector(+SQRT2INV, +SQRT2INV,0       ));
        
        if (m[0] == lo.x)                                  return set(m_bc_type[Face::XLO],         u, gradu, sigma, Set::Vector(-1, 0, 0));
        if (m[0] == hi.x)                                  return set(m_bc_type[Face::XHI],         u, gradu, sigma, Set::Vector(+1, 0, 0));
        if (m[1] == lo.y)                                  return set(m_bc_type[Face::YLO],         u, gradu, sigma, Set::Vector( 0,-1, 0));
        if (m[1] == hi.y)                                  return set(m_bc_type[Face::YHI],         u, gradu, sigma, Set::Vector( 0,+1, 0));
        if (m[2] == lo.z)                                  return set(m_bc_type[Face::ZLO],         u, gradu, sigma, Set::Vector( 0, 0,-1));
        if (m[2] == hi.z)                                  return set(m_bc_type[Face::ZHI],         u, gradu, sigma, Set::Vector( 0, 0,+1));
        
        #endif
        
        Util::Abort(INFO,"Boundary condition error");
        return Set::Vector::Zero();
    }

    AMREX_FORCE_INLINE
    Set::Vector set(std::array<Type,AMREX_SPACEDIM> &bc_type, 
                    const Set::Vector &u, const Set::Matrix &gradu, const Set::Matrix &sigma, Set::Vector n) const
    {
        Set::Vector f = Set::Vector::Zero();
        for (int i = 0; i < AMREX_SPACEDIM; i++)
        {
            if      (bc_type[i] == Type::Displacement) 
                f(i) = u(i);
            else if (bc_type[i] == Type::Traction)
                f(i) = (sigma*n)(i);
            else if (bc_type[i] == Type::Neumann)
                f(i) = (gradu*n)(i);
            else if (bc_type[i] == Periodic)
                continue;
        }
        return f;
    }

protected:
    
    #if AMREX_SPACEDIM==2
    static const int m_nfaces = 8;
    #elif AMREX_SPACEDIM==3
    static const int m_nfaces = 26;
    #endif

    std::array<std::array<Type,            AMREX_SPACEDIM>, m_nfaces> m_bc_type; 
    //std::array<std::array<FunctionParserAD,AMREX_SPACEDIM>, m_nfaces> m_bc_func; 
    std::array<std::array<amrex::Parser,AMREX_SPACEDIM>, m_nfaces> m_bc_func_parser; 
    std::array<std::array<amrex::ParserExecutor<4>,AMREX_SPACEDIM>, m_nfaces> m_bc_func; 

public:

    //
    // This is basically the same as :ref:`BC::Operator::Elastic::Constant` except that
    // you can use dynamically compiled expressions in space and time to define values.
    //
    // Usage is exactly the same except that the "val" inputs can depend on x, y, z, and t.
    //
    static void Parse(Expression & value, IO::ParmParse & pp)
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
        if (pp.contains("type.xloylozlo")) { pp.queryarr("type.xloylozlo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XLO_YLO_ZLO][i] = bcmap[str[i]]; } // 3D Corner
        if (pp.contains("type.xloylozhi")) { pp.queryarr("type.xloylozhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XLO_YLO_ZHI][i] = bcmap[str[i]]; } // 3D Corner
        if (pp.contains("type.xloyhizlo")) { pp.queryarr("type.xloyhizlo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XLO_YHI_ZLO][i] = bcmap[str[i]]; } // 3D Corner
        if (pp.contains("type.xloyhizhi")) { pp.queryarr("type.xloyhizhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XLO_YHI_ZHI][i] = bcmap[str[i]]; } // 3D Corner
        if (pp.contains("type.xhiylozlo")) { pp.queryarr("type.xhiylozlo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XHI_YLO_ZLO][i] = bcmap[str[i]]; } // 3D Corner
        if (pp.contains("type.xhiylozhi")) { pp.queryarr("type.xhiylozhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XHI_YLO_ZHI][i] = bcmap[str[i]]; } // 3D Corner
        if (pp.contains("type.xhiyhizlo")) { pp.queryarr("type.xhiyhizlo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XHI_YHI_ZLO][i] = bcmap[str[i]]; } // 3D Corner
        if (pp.contains("type.xhiyhizhi")) { pp.queryarr("type.xhiyhizhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XHI_YHI_ZHI][i] = bcmap[str[i]]; } // 3D Corner
        #endif

        #if AMREX_SPACEDIM==3
        if (pp.contains("type.ylozlo")) { pp.queryarr("type.ylozlo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::YLO_ZLO][i] = bcmap[str[i]]; } // 3D Edge
        if (pp.contains("type.ylozhi")) { pp.queryarr("type.ylozhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::YLO_ZHI][i] = bcmap[str[i]]; } // 3D Edge
        if (pp.contains("type.yhizlo")) { pp.queryarr("type.yhizlo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::YHI_ZLO][i] = bcmap[str[i]]; } // 3D Edge
        if (pp.contains("type.yhizhi")) { pp.queryarr("type.yhizhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::YHI_ZHI][i] = bcmap[str[i]]; } // 3D Edge
        if (pp.contains("type.zloxlo")) { pp.queryarr("type.zloxlo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::ZLO_XLO][i] = bcmap[str[i]]; } // 3D Edge
        if (pp.contains("type.zloxhi")) { pp.queryarr("type.zloxhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::ZLO_XHI][i] = bcmap[str[i]]; } // 3D Edge
        if (pp.contains("type.zhixlo")) { pp.queryarr("type.zhixlo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::ZHI_XLO][i] = bcmap[str[i]]; } // 3D Edge
        if (pp.contains("type.zhixhi")) { pp.queryarr("type.zhixhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::ZHI_XHI][i] = bcmap[str[i]]; } // 3D Edge
        #endif
        if (pp.contains("type.xloylo")) { pp.queryarr("type.xloylo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XLO_YLO][i] = bcmap[str[i]]; } // 3D Edge / 2D Corner
        if (pp.contains("type.xloyhi")) { pp.queryarr("type.xloyhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XLO_YHI][i] = bcmap[str[i]]; } // 3D Edge / 2D Corner
        if (pp.contains("type.xhiylo")) { pp.queryarr("type.xhiylo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XHI_YLO][i] = bcmap[str[i]]; } // 3D Edge / 2D Corner
        if (pp.contains("type.xhiyhi")) { pp.queryarr("type.xhiyhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XHI_YHI][i] = bcmap[str[i]]; } // 3D Edge / 2D Corner

        if (pp.contains("type.xlo")) { pp.queryarr("type.xlo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XLO][i] = bcmap[str[i]]; } // 3D Face  2D Edge
        if (pp.contains("type.xhi")) { pp.queryarr("type.xhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::XHI][i] = bcmap[str[i]]; } // 3D Face  2D Edge
        if (pp.contains("type.ylo")) { pp.queryarr("type.ylo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::YLO][i] = bcmap[str[i]]; } // 3D Face  2D Edge
        if (pp.contains("type.yhi")) { pp.queryarr("type.yhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::YHI][i] = bcmap[str[i]]; } // 3D Face  2D Edge
        #if AMREX_SPACEDIM==3
        if (pp.contains("type.zlo")) { pp.queryarr("type.zlo",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::ZLO][i] = bcmap[str[i]]; } // 3D Face
        if (pp.contains("type.zhi")) { pp.queryarr("type.zhi",str); for (int i = 0; i < AMREX_SPACEDIM; i++) value.m_bc_type[Face::ZHI][i] = bcmap[str[i]]; } // 3D Face
        #endif

        // VALS
        std::vector<std::string> val;

        for (int face = 0; face != Face::INT; face++)
        {
            std::string querystr = std::string("val.") + std::string(facestr[face]);
            if (pp.contains(querystr.c_str())) pp.queryarr(querystr.c_str(),val);
            else val.resize(AMREX_SPACEDIM,"0.0");
            
            for (int i = 0 ; i < AMREX_SPACEDIM; i++)
            {
                value.m_bc_func_parser[face][i] = amrex::Parser(val[i].c_str());
                value.m_bc_func_parser[face][i].registerVariables({"x","y","z","t"});
                value.m_bc_func[face][i] = value.m_bc_func_parser[face][i].compile<4>();
            }
        }        
    }
    Expression(IO::ParmParse &pp, std::string name): Expression()
    {pp.queryclass(name,*this);}
};
}
}
}
#endif
