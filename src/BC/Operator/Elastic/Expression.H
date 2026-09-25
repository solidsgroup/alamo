//
// Use expressions to define space and time varying boundary conditions
// for elastic problems. 
//

#ifndef BC_OPERATOR_ELASTIC_EXPRESSION_H
#define BC_OPERATOR_ELASTIC_EXPRESSION_H

#include "IO/ParmParse.H"
#include "BC/Operator/Elastic/Elastic.H"

#include "AMReX_Parser.H"
#include <cmath>
#include <limits>

namespace BC
{
namespace Operator
{
namespace Elastic
{
class Expression : public Elastic
{
public:
    static constexpr const char* name = "expression";

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

    bool HasContact() const { return m_contact; }
    int ContactIterations() const { return m_contact_iterations; }

    void TagCellsForRefinement(amrex::TagBoxArray& a_tags,
        const amrex::Geometry& a_geom, Set::Scalar a_time) const override
    {
        if (!m_contact) return;

        constexpr int normal = AMREX_SPACEDIM - 1;
        const auto prob_lo = a_geom.ProbLoArray();
        const auto prob_hi = a_geom.ProbHiArray();
        const auto dx = a_geom.CellSizeArray();
        const amrex::Box domain = a_geom.Domain();

        for (amrex::MFIter mfi(a_tags); mfi.isValid(); ++mfi)
        {
            const amrex::Box box = mfi.validbox() & domain;
            auto tags = a_tags.array(mfi);
            amrex::LoopOnCpu(box, [&](int i, int j, int k) {
                const int index[3] = {i, j, k};
                Set::Scalar pos[3] = {0.0, 0.0, 0.0};
                for (int d = 0; d < AMREX_SPACEDIM; ++d)
                    pos[d] = prob_lo[d] + (index[d] + 0.5)*dx[d];

                bool refine = m_refinement_region(pos[0], pos[1], pos[2], a_time);

                // Cell-center tagging alone can miss a boundary node near the
                // edge of the contact region. Tag every top cell having at
                // least one contact-region corner so the finest level covers
                // every node that the active-set update may inspect.
                if (!refine && index[normal] == domain.bigEnd(normal))
                {
                    constexpr int tangential_corners = 1 << normal;
                    for (int corner = 0; corner < tangential_corners && !refine; ++corner)
                    {
                        for (int d = 0; d < normal; ++d)
                            pos[d] = prob_lo[d] + (index[d] + ((corner >> d) & 1))*dx[d];
                        pos[normal] = prob_hi[normal];
                        refine = m_region(pos[0], pos[1], pos[2], a_time);
                    }
                }

                if (refine) tags(i,j,k) = amrex::TagBox::SET;
            });
        }
    }

    // Small-slope contact on the upper face in the last coordinate direction.
    // Freeze this active set throughout each equilibrium solve (including MG).
    bool UpdateContact(const amrex::FabArray<amrex::BaseFab<Set::Vector>>& displacement,
                       const amrex::FabArray<amrex::BaseFab<Set::Matrix>>& stress,
                       const amrex::Geometry& geom, bool report = false)
    {
        if (!m_contact) return false;
        constexpr int normal = AMREX_SPACEDIM-1;
        auto domain = amrex::surroundingNodes(geom.Domain());
        int size = 1;
        for (int d = 0; d < normal; ++d)
        {
            m_contact_length[d] = domain.length(d);
            size *= domain.length(d);
        }
        if (m_active.empty()) m_active.resize(size, 0);
        Util::Assert(INFO, TEST(int(m_active.size()) == size), "Contact mesh cannot change during a run");
        std::vector<int> next(size, 0);
        std::vector<Set::Scalar> pressure_field(size, -std::numeric_limits<Set::Scalar>::infinity());
        Set::Scalar max_penetration = 0.0, max_tension = 0.0;
        for (amrex::MFIter mfi(displacement); mfi.isValid(); ++mfi)
        {
            auto box = mfi.validbox() & domain;
            box.setSmall(normal, domain.bigEnd(normal));
            box &= mfi.validbox();
            auto u = displacement.const_array(mfi);
            auto p = stress.const_array(mfi);
            amrex::LoopOnCpu(box, [&](int i, int j, int k) {
                const int index[3] = {i,j,k};
                int offset = 0, stride = 1;
                for (int d = 0; d < normal; ++d)
                {
                    if (index[d] == domain.smallEnd(d) || index[d] == domain.bigEnd(d)) return;
                    offset += (index[d]-domain.smallEnd(d))*stride;
                    stride *= domain.length(d);
                }
                auto pos = Set::Position(i,j,k,geom,amrex::IndexType::TheNodeType());
                const Set::Scalar x = pos(0), y = pos(1);
                const Set::Scalar z = AMREX_D_PICK(0.0,0.0,pos(2));
                if (!m_region(x,y,z,m_time)) return;
                const auto penetration = u(i,j,k)(normal)-m_obstacle(x,y,z,m_time);
                const auto pressure = -p(i,j,k)(normal,normal);
                pressure_field[offset] = pressure;
                max_penetration = std::max(max_penetration, penetration);
                max_tension = std::max(max_tension, -pressure);
                const auto criterion = pressure + m_contact_stiffness*penetration;
                Util::Assert(INFO, TEST(std::isfinite(criterion)), "Nonfinite contact state");
                next[offset] = criterion > m_contact_tolerance ? 1 :
                    (criterion < -m_contact_tolerance ? 0 : m_active[offset]);
            });
        }
        amrex::ParallelDescriptor::ReduceIntMax(next.data(), size);
        amrex::ParallelDescriptor::ReduceRealMax(pressure_field.data(), size);
        // The finest level must cover the entire possible contact region.
        // This avoids mixing coarse/fine reactions or duplicating contact load.
        for (int n = 0; n < size; ++n)
        {
            Set::Scalar pos[3] = {0,0,0};
            int remainder = n;
            bool edge = false;
            for (int d = 0; d < normal; ++d)
            {
                const int index = remainder % domain.length(d);
                remainder /= domain.length(d);
                edge = edge || index == 0 || index == domain.length(d)-1;
                pos[d] = geom.ProbLo(d) + index*geom.CellSize(d);
            }
            pos[normal] = geom.ProbHi(normal);
            if (!edge && m_region(pos[0],pos[1],pos[2],m_time))
                Util::Assert(INFO, TEST(std::isfinite(pressure_field[n])),
                    "The finest grid must cover the entire contact.region on the upper surface");
        }
        amrex::ParallelDescriptor::ReduceRealMax(max_penetration);
        amrex::ParallelDescriptor::ReduceRealMax(max_tension);
        const bool changed = next != m_active;
        m_active.swap(next);
        Set::Scalar load = 0.0, area = 1.0;
        int active = 0;
        for (int d = 0; d < normal; ++d) area *= geom.CellSize(d);
        for (int n = 0; n < size; ++n)
        {
            active += m_active[n];
            if (std::isfinite(pressure_field[n])) load += pressure_field[n]*area;
        }
        if (report) Util::Message(INFO, "Contact time=", m_time, " stable=", !changed,
            " active_nodes=", active, " load=", load,
            " penetration=", max_penetration, " tension=", max_tension);
        return changed;
    }

    void CheckGeometry(const amrex::Geometry& geom) const
    {
        for (const auto& face : m_type_expression)
            for (bool expression : face)
                if (expression)
                    for (int d = 0; d < AMREX_SPACEDIM; ++d)
                        Util::Assert(INFO, TEST(!geom.isPeriodic(d)),
                            "Spatial elastic type expressions currently require nonperiodic geometry");
    }

    virtual void
    Init(amrex::FabArray<amrex::BaseFab<Set::Vector>> *a_rhs,
        const amrex::Geometry &a_geom,
        bool a_homogeneous = false) const override
    {
            m_prob_lo = a_geom.ProbLoArray();
            m_prob_hi = a_geom.ProbHiArray();
            CheckGeometry(a_geom);
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
                        if (a_homogeneous && Types(face, i, j, k, domain)[dir] == Type::Displacement)
                            rhs(i,j,k)(dir) = 0.0;
                        else 
                        {
                            Set::Vector x = Set::Position(i, j, k, a_geom, type);
                            #if AMREX_SPACEDIM==1
                            rhs(i,j,k)(dir) = Value(face,dir,x(0),0.0,0.0);
                            #elif AMREX_SPACEDIM==2
                            rhs(i,j,k)(dir) = Value(face,dir,x(0),x(1),0.0);
                            #elif AMREX_SPACEDIM==3
                            rhs(i,j,k)(dir) = Value(face,dir,x(0),x(1),x(2));
                            #endif
                            if (Types(face,i,j,k,domain)[dir] == Type::Traction)
                                rhs(i,j,k)(dir) *= m_traction_scale;
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
            m_prob_lo = a_geom.ProbLoArray();
            m_prob_hi = a_geom.ProbHiArray();
            CheckGeometry(a_geom);
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
                        if (a_homogeneous && Types(face, i, j, k, domain)[dir] == Type::Displacement)
                            rhs(i,j,k,dir) = 0.0;
                        else 
                        {
                            Set::Vector x = Set::Position(i, j, k, a_geom, type);
                            #if AMREX_SPACEDIM==1
                            rhs(i,j,k,dir) = Value(face,dir,x(0),0.0,0.0);
                            #elif AMREX_SPACEDIM==2
                            rhs(i,j,k,dir) = Value(face,dir,x(0),x(1),0.0);
                            #elif AMREX_SPACEDIM==3
                            rhs(i,j,k,dir) = Value(face,dir,x(0),x(1),x(2));
                            #endif
                            if (Types(face,i,j,k,domain)[dir] == Type::Traction)
                                rhs(i,j,k,dir) *= m_traction_scale;
                        }
                    }
                }
            });
        }                    
    }



    virtual
    std::array<Type,AMREX_SPACEDIM> getType (
                const int &i, const int &j, [[maybe_unused]] const int &k,
                const amrex::Box &domain) override
    {
        (void)i, (void)j, (void)k;
        amrex::IntVect m(AMREX_D_DECL(i,j,k));
        const amrex::Dim3 lo = amrex::lbound(domain), hi = amrex::ubound(domain);

        // Corners
        #if AMREX_SPACEDIM == 2
        if (m[0] == lo.x && m[1] == lo.y)                  return Types(Face::XLO_YLO, i, j, k, domain);
        if (m[0] == lo.x && m[1] == hi.y)                  return Types(Face::XLO_YHI, i, j, k, domain);
        if (m[0] == hi.x && m[1] == lo.y)                  return Types(Face::XHI_YLO, i, j, k, domain);
        if (m[0] == hi.x && m[1] == hi.y)                  return Types(Face::XHI_YHI, i, j, k, domain);
        if (m[0] == lo.x)                                  return Types(Face::XLO, i, j, k, domain);
        if (m[0] == hi.x)                                  return Types(Face::XHI, i, j, k, domain);
        if (m[1] == lo.y)                                  return Types(Face::YLO, i, j, k, domain);
        if (m[1] == hi.y)                                  return Types(Face::YHI, i, j, k, domain);
        
        #elif AMREX_SPACEDIM == 3
        
        if (m[0] == lo.x &&  m[1] == lo.y && m[2] == lo.z) return Types(Face::XLO_YLO_ZLO, i, j, k, domain);
        if (m[0] == lo.x &&  m[1] == lo.y && m[2] == hi.z) return Types(Face::XLO_YLO_ZHI, i, j, k, domain);
        if (m[0] == lo.x &&  m[1] == hi.y && m[2] == lo.z) return Types(Face::XLO_YHI_ZLO, i, j, k, domain);
        if (m[0] == lo.x &&  m[1] == hi.y && m[2] == hi.z) return Types(Face::XLO_YHI_ZHI, i, j, k, domain);
        if (m[0] == hi.x &&  m[1] == lo.y && m[2] == lo.z) return Types(Face::XHI_YLO_ZLO, i, j, k, domain);
        if (m[0] == hi.x &&  m[1] == lo.y && m[2] == hi.z) return Types(Face::XHI_YLO_ZHI, i, j, k, domain);
        if (m[0] == hi.x &&  m[1] == hi.y && m[2] == lo.z) return Types(Face::XHI_YHI_ZLO, i, j, k, domain);
        if (m[0] == hi.x &&  m[1] == hi.y && m[2] == hi.z) return Types(Face::XHI_YHI_ZHI, i, j, k, domain);
        
        if (m[1] == lo.y && m[2] == lo.z)                  return Types(Face::YLO_ZLO, i, j, k, domain);
        if (m[1] == lo.y && m[2] == hi.z)                  return Types(Face::YLO_ZHI, i, j, k, domain);
        if (m[1] == hi.y && m[2] == lo.z)                  return Types(Face::YHI_ZLO, i, j, k, domain);
        if (m[1] == hi.y && m[2] == hi.z)                  return Types(Face::YHI_ZHI, i, j, k, domain);
        if (m[2] == lo.z && m[0] == lo.x)                  return Types(Face::ZLO_XLO, i, j, k, domain);
        if (m[2] == lo.z && m[0] == hi.x)                  return Types(Face::ZLO_XHI, i, j, k, domain);
        if (m[2] == hi.z && m[0] == lo.x)                  return Types(Face::ZHI_XLO, i, j, k, domain);
        if (m[2] == hi.z && m[0] == hi.x)                  return Types(Face::ZHI_XHI, i, j, k, domain);
        if (m[0] == lo.x && m[1] == lo.y)                  return Types(Face::XLO_YLO, i, j, k, domain);
        if (m[0] == lo.x && m[1] == hi.y)                  return Types(Face::XLO_YHI, i, j, k, domain);
        if (m[0] == hi.x && m[1] == lo.y)                  return Types(Face::XHI_YLO, i, j, k, domain);
        if (m[0] == hi.x && m[1] == hi.y)                  return Types(Face::XHI_YHI, i, j, k, domain);
        
        if (m[0] == lo.x)                                  return Types(Face::XLO, i, j, k, domain);
        if (m[0] == hi.x)                                  return Types(Face::XHI, i, j, k, domain);
        if (m[1] == lo.y)                                  return Types(Face::YLO, i, j, k, domain);
        if (m[1] == hi.y)                                  return Types(Face::YHI, i, j, k, domain);
        if (m[2] == lo.z)                                  return Types(Face::ZLO, i, j, k, domain);
        if (m[2] == hi.z)                                  return Types(Face::ZHI, i, j, k, domain);
        
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
        auto types = getType(i, j, k, domain);

        const amrex::Dim3 lo= amrex::lbound(domain), hi = amrex::ubound(domain);
        
        amrex::IntVect m(AMREX_D_DECL(i,j,k));

        // Corners
        #if AMREX_SPACEDIM == 2
        
        if (m[0] == lo.x && m[1] == lo.y)                  return set(types,     u, gradu, sigma, Set::Vector(-SQRT2INV, -SQRT2INV));
        if (m[0] == lo.x && m[1] == hi.y)                  return set(types,     u, gradu, sigma, Set::Vector(-SQRT2INV, +SQRT2INV));
        if (m[0] == hi.x && m[1] == lo.y)                  return set(types,     u, gradu, sigma, Set::Vector(+SQRT2INV, -SQRT2INV));
        if (m[0] == hi.x && m[1] == hi.y)                  return set(types,     u, gradu, sigma, Set::Vector(+SQRT2INV, +SQRT2INV));
        
        if (m[0] == lo.x)                                  return set(types,         u, gradu, sigma, Set::Vector(-1, 0));
        if (m[0] == hi.x)                                  return set(types,         u, gradu, sigma, Set::Vector(+1, 0));
        if (m[1] == lo.y)                                  return set(types,         u, gradu, sigma, Set::Vector( 0,-1));
        if (m[1] == hi.y)                                  return set(types,         u, gradu, sigma, Set::Vector( 0,+1));
        
        #elif AMREX_SPACEDIM == 3
        
        if (m[0] == lo.x &&  m[1] == lo.y && m[2] == lo.z) return set(types, u, gradu, sigma, Set::Vector(-SQRT3INV,-SQRT3INV,-SQRT3INV));
        if (m[0] == lo.x &&  m[1] == lo.y && m[2] == hi.z) return set(types, u, gradu, sigma, Set::Vector(-SQRT3INV,-SQRT3INV,+SQRT3INV));
        if (m[0] == lo.x &&  m[1] == hi.y && m[2] == lo.z) return set(types, u, gradu, sigma, Set::Vector(-SQRT3INV,+SQRT3INV,-SQRT3INV));
        if (m[0] == lo.x &&  m[1] == hi.y && m[2] == hi.z) return set(types, u, gradu, sigma, Set::Vector(-SQRT3INV,+SQRT3INV,+SQRT3INV));
        if (m[0] == hi.x &&  m[1] == lo.y && m[2] == lo.z) return set(types, u, gradu, sigma, Set::Vector(+SQRT3INV,-SQRT3INV,-SQRT3INV));
        if (m[0] == hi.x &&  m[1] == lo.y && m[2] == hi.z) return set(types, u, gradu, sigma, Set::Vector(+SQRT3INV,-SQRT3INV,+SQRT3INV));
        if (m[0] == hi.x &&  m[1] == hi.y && m[2] == lo.z) return set(types, u, gradu, sigma, Set::Vector(+SQRT3INV,+SQRT3INV,-SQRT3INV));
        if (m[0] == hi.x &&  m[1] == hi.y && m[2] == hi.z) return set(types, u, gradu, sigma, Set::Vector(+SQRT3INV,+SQRT3INV,+SQRT3INV));
        
        if (m[1] == lo.y && m[2] == lo.z)                  return set(types,     u, gradu, sigma, Set::Vector(0,        -SQRT2INV,-SQRT2INV));
        if (m[1] == lo.y && m[2] == hi.z)                  return set(types,     u, gradu, sigma, Set::Vector(0,        -SQRT2INV,+SQRT2INV));
        if (m[1] == hi.y && m[2] == lo.z)                  return set(types,     u, gradu, sigma, Set::Vector(0,        +SQRT2INV,-SQRT2INV));
        if (m[1] == hi.y && m[2] == hi.z)                  return set(types,     u, gradu, sigma, Set::Vector(0,        +SQRT2INV,+SQRT2INV));
        if (m[2] == lo.z && m[0] == lo.x)                  return set(types,     u, gradu, sigma, Set::Vector(-SQRT2INV, 0,       -SQRT2INV));
        if (m[2] == lo.z && m[0] == hi.x)                  return set(types,     u, gradu, sigma, Set::Vector(+SQRT2INV, 0,       -SQRT2INV));
        if (m[2] == hi.z && m[0] == lo.x)                  return set(types,     u, gradu, sigma, Set::Vector(-SQRT2INV, 0,       +SQRT2INV));
        if (m[2] == hi.z && m[0] == hi.x)                  return set(types,     u, gradu, sigma, Set::Vector(+SQRT2INV, 0,       +SQRT2INV));
        if (m[0] == lo.x && m[1] == lo.y)                  return set(types,     u, gradu, sigma, Set::Vector(-SQRT2INV, -SQRT2INV,0       ));
        if (m[0] == lo.x && m[1] == hi.y)                  return set(types,     u, gradu, sigma, Set::Vector(-SQRT2INV, +SQRT2INV,0       ));
        if (m[0] == hi.x && m[1] == lo.y)                  return set(types,     u, gradu, sigma, Set::Vector(+SQRT2INV, -SQRT2INV,0       ));
        if (m[0] == hi.x && m[1] == hi.y)                  return set(types,     u, gradu, sigma, Set::Vector(+SQRT2INV, +SQRT2INV,0       ));
        
        if (m[0] == lo.x)                                  return set(types,         u, gradu, sigma, Set::Vector(-1, 0, 0));
        if (m[0] == hi.x)                                  return set(types,         u, gradu, sigma, Set::Vector(+1, 0, 0));
        if (m[1] == lo.y)                                  return set(types,         u, gradu, sigma, Set::Vector( 0,-1, 0));
        if (m[1] == hi.y)                                  return set(types,         u, gradu, sigma, Set::Vector( 0,+1, 0));
        if (m[2] == lo.z)                                  return set(types,         u, gradu, sigma, Set::Vector( 0, 0,-1));
        if (m[2] == hi.z)                                  return set(types,         u, gradu, sigma, Set::Vector( 0, 0,+1));
        
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
                f(i) = m_traction_scale*(sigma*n)(i);
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

    // Physical bounds are shared by all AMR and multigrid levels. The nodal
    // domain supplied by the operator determines the spacing at each level.
    mutable amrex::GpuArray<Set::Scalar, AMREX_SPACEDIM> m_prob_lo{}, m_prob_hi{};
    std::array<std::array<bool, AMREX_SPACEDIM>, m_nfaces> m_type_expression{};
    std::array<std::array<amrex::Parser, AMREX_SPACEDIM>, m_nfaces> m_type_parser;
    std::array<std::array<amrex::ParserExecutor<6>, AMREX_SPACEDIM>, m_nfaces> m_type_func;

    bool m_contact = false;
    Set::Scalar m_traction_scale = 1.0;
    int m_contact_iterations = 40;
    Set::Scalar m_contact_stiffness = 1.0, m_contact_tolerance = 1e-8;
    amrex::Parser m_obstacle_parser, m_region_parser, m_refinement_region_parser;
    amrex::ParserExecutor<4> m_obstacle, m_region, m_refinement_region;
    std::vector<int> m_active;
    std::array<int,AMREX_SPACEDIM> m_contact_length{};

    int ContactAt(Set::Scalar x, Set::Scalar y, Set::Scalar z) const
    {
        if (!m_contact || m_active.empty()) return 0;
        const Set::Scalar pos[3] = {x,y,z};
        constexpr int normal = AMREX_SPACEDIM-1;
        if (std::abs(pos[normal]-m_prob_hi[normal]) > 1e-10*(m_prob_hi[normal]-m_prob_lo[normal])) return 0;
        int offset = 0, stride = 1;
        for (int d = 0; d < normal; ++d)
        {
            int index = std::lround((pos[d]-m_prob_lo[d])/(m_prob_hi[d]-m_prob_lo[d])*(m_contact_length[d]-1));
            if (index < 0 || index >= m_contact_length[d]) return 0;
            offset += index*stride;
            stride *= m_contact_length[d];
        }
        return m_active[offset];
    }

    Set::Scalar Value(Face face, int dir, Set::Scalar x, Set::Scalar y, Set::Scalar z) const
    {
        return m_bc_func[face][dir](x,y,z,m_time,Set::Scalar(ContactAt(x,y,z)),
            m_contact ? m_obstacle(x,y,z,m_time) : 0.0);
    }

    std::array<Type, AMREX_SPACEDIM> Types(
        Face face, int i, int j, int k, const amrex::Box& domain) const
    {
        auto ret = m_bc_type[face];
        const int index[3] = {i,j,k};
        Set::Scalar x[3] = {0,0,0};
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
            x[d] = m_prob_lo[d] + (m_prob_hi[d]-m_prob_lo[d]) *
                (index[d]-domain.smallEnd(d)) / (domain.length(d)-1);
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
            if (m_type_expression[face][d])
            {
                const auto type = m_type_func[face][d](x[0],x[1],x[2],m_time,
                    Set::Scalar(ContactAt(x[0],x[1],x[2])), m_contact ? m_obstacle(x[0],x[1],x[2],m_time) : 0.0);
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                    type == int(Type::Displacement) || type == int(Type::Traction) || type == int(Type::Neumann),
                    "Elastic type expression must return disp, trac, or neumann");
                ret[d] = static_cast<Type>(static_cast<int>(type));
            }
        return ret;
    }

    std::array<std::array<Type,            AMREX_SPACEDIM>, m_nfaces> m_bc_type; 
    //std::array<std::array<FunctionParserAD,AMREX_SPACEDIM>, m_nfaces> m_bc_func; 
    std::array<std::array<amrex::Parser,AMREX_SPACEDIM>, m_nfaces> m_bc_func_parser; 
    std::array<std::array<amrex::ParserExecutor<6>,AMREX_SPACEDIM>, m_nfaces> m_bc_func;

public:

    //
    // This is basically the same as :ref:`BC::Operator::Elastic::Constant` except that
    // you can use dynamically compiled expressions in space and time to define values.
    //
    // Both "val" and "type" accept expressions in x, y, z, and t.
    // Example: type.yhi = trac "if(abs(x)<0.25,disp,trac)"
    //
    static void Parse(Expression & value, IO::ParmParse & pp)
    {
        // Scale BOTH sides of traction equations, without changing physical loads.
        pp.query_default("traction_scale", value.m_traction_scale, 1.0);
        Util::Assert(INFO, TEST(value.m_traction_scale > 0 && std::isfinite(value.m_traction_scale)));
        pp.query_default("contact.enabled", value.m_contact, false);
        if (value.m_contact)
        {
#ifdef AMREX_USE_GPU
            Util::Abort(INFO, "Expression contact currently supports CPU builds only");
#endif
            std::string obstacle, region, refinement_region;
            pp.query_required("contact.obstacle", obstacle);
            pp.query_default("contact.region", region, std::string("1"));
            pp.query_default("contact.refinement_region", refinement_region, std::string("0"));
            pp.query_required("contact.stiffness", value.m_contact_stiffness);
            pp.query_default("contact.tolerance", value.m_contact_tolerance, 1e-8);
            pp.query_default("contact.max_iter", value.m_contact_iterations, 40);
            Util::Assert(INFO, TEST(value.m_contact_stiffness > 0 && value.m_contact_tolerance >= 0 && value.m_contact_iterations > 0));
            value.m_obstacle_parser = amrex::Parser(obstacle);
            value.m_region_parser = amrex::Parser(region);
            value.m_refinement_region_parser = amrex::Parser(refinement_region);
            value.m_obstacle_parser.registerVariables({"x","y","z","t"});
            value.m_region_parser.registerVariables({"x","y","z","t"});
            value.m_refinement_region_parser.registerVariables({"x","y","z","t"});
            value.m_obstacle = value.m_obstacle_parser.compile<4>();
            value.m_region = value.m_region_parser.compile<4>();
            value.m_refinement_region = value.m_refinement_region_parser.compile<4>();
        }
        std::map<std::string, Type> bcmap;
        bcmap["displacement"] = Type::Displacement;
        bcmap["disp"]         = Type::Displacement;
        bcmap["neumann"]      = Type::Neumann;
        bcmap["traction"]     = Type::Traction;
        bcmap["trac"]         = Type::Traction;
        bcmap["periodic"]     = Type::Periodic;
        
        std::vector<std::string> type, val;
        
        //
        // Lambda to parse the type and val strings, and initialize parser.
        // 

        auto set =  [&] (std::vector<std::string> type, std::vector<std::string> val, Face face)  
        {
            if (type.size() == 1)
                type.resize(AMREX_SPACEDIM,type[0]);
            else if (type.size() < AMREX_SPACEDIM)
                Util::Abort(INFO,"Invalid number of types specified for ", facestr[face]);

            if (val.size() == 1)
                val.resize(AMREX_SPACEDIM,val[0]);
            else if (val.size() < AMREX_SPACEDIM)
                Util::Abort(INFO,"Invalid number of vals specified for ", facestr[face]);

            for (int i = 0 ; i < AMREX_SPACEDIM; i++)
            {
                // set enums to appropriate type
                auto known = bcmap.find(type[i]);
                value.m_type_expression[face][i] = known == bcmap.end();
                if (known != bcmap.end()) value.m_bc_type[face][i] = known->second;
                else
                {
                    auto& parser = value.m_type_parser[face][i];
                    parser = amrex::Parser(type[i]);
                    parser.registerVariables({"x","y","z","t","contact","obstacle"});
                    for (const auto& entry : bcmap) parser.setConstant(entry.first, entry.second);
                    value.m_type_func[face][i] = parser.compile<6>();
                }
                
                // initialize the parser
                value.m_bc_func_parser[face][i] = amrex::Parser(val[i].c_str());
                value.m_bc_func_parser[face][i].registerVariables({"x","y","z","t","contact","obstacle"});
                value.m_bc_func[face][i] = value.m_bc_func_parser[face][i].compile<6>();
            }
        };

        pp.queryarr_default("type.xlo",type, {"disp"}); // 3D Face / 2D Edge type
        pp.queryarr_default("val.xlo", val , {"0.0"}); // 3D Face / 2D Edge value
        set(type,val,Face::XLO);
        type.clear(); val.clear();

        pp.queryarr_default("type.xhi",type, {"disp"}); // 3D Face / 2D Edge  type
        pp.queryarr_default("val.xhi" ,val , {"0.0"}); // 3D Face / 2D Edge value
        set(type,val,Face::XHI);
        type.clear(); val.clear();

        pp.queryarr_default("type.ylo",type, {"disp"}); // 3D Face / 2D Edge type
        pp.queryarr_default("val.ylo" ,val , {"0.0"}); // 3D Face / 2D Edge value
        set(type,val,Face::YLO);
        type.clear(); val.clear();

        pp.queryarr_default("type.yhi",type, {"disp"}); // 3D Face / 2D Edge type
        pp.queryarr_default("val.yhi", val , {"0.0"}); // 3D Face / 2D Edge value
        set(type,val,Face::YHI);
        type.clear(); val.clear();

#if AMREX_SPACEDIM==3
        pp.queryarr_default("type.zlo",type, {"disp"}); // 3D Face type
        pp.queryarr_default("val.zlo", val , {"0.0"}); // 3D Face value
        set(type,val,Face::ZLO);
        type.clear(); val.clear();

        pp.queryarr_default("type.zhi",type, {"disp"}); // 3D Face type
        pp.queryarr_default("val.zhi", val , {"0.0"}); // 3D Face value
        set(type,val,Face::ZHI);
        type.clear(); val.clear();
#endif


        pp.queryarr_default("type.xloylo",type, {"disp"}); // 3D Edge / 2D Corner type
        pp.queryarr_default("val.xloylo", val , {"0.0"}); // 3D Edge / 2D Corner value
        set(type,val,Face::XLO_YLO);
        type.clear(); val.clear();

        pp.queryarr_default("type.xloyhi",type, {"disp"}); // 3D Edge / 2D Corner type
        pp.queryarr_default("val.xloyhi", val , {"0.0"}); // 3D Edge / 2D Corner value
        set(type,val,Face::XLO_YHI);
        type.clear(); val.clear();

        pp.queryarr_default("type.xhiylo",type, {"disp"}); // 3D Edge / 2D Corner type
        pp.queryarr_default("val.xhiylo", val , {"0.0"}); // 3D Edge / 2D Corner value
        set(type,val,Face::XHI_YLO);
        type.clear(); val.clear();

        pp.queryarr_default("type.xhiyhi",type, {"disp"}); // 3D Edge / 2D Corner type
        pp.queryarr_default("val.xhiyhi", val , {"0.0"}); // 3D Edge / 2D Corner value
        set(type,val,Face::XHI_YHI);
        type.clear(); val.clear();

#if AMREX_SPACEDIM==3
        pp.queryarr_default("type.xloylozlo",type,{"disp"}); // 3D Corner type
        pp.queryarr_default("val.xloylozlo",val,{"0.0"}); // 3D Corner value
        set(type,val,Face::XLO_YLO_ZLO);
        type.clear(); val.clear();

        pp.queryarr_default("type.xloylozhi",type, {"disp"});  // 3D Corner type
        pp.queryarr_default("val.xloylozhi", val , {"0.0"});  // 3D Corner value
        set(type,val,Face::XLO_YLO_ZHI);
        type.clear(); val.clear();

        pp.queryarr_default("type.xloyhizlo",type, {"disp"}); // 3D Corner type
        pp.queryarr_default("val.xloyhizlo", val , {"0.0"}); // 3D Corner value
        set(type,val,Face::XLO_YHI_ZLO);
        type.clear(); val.clear();

        pp.queryarr_default("type.xloyhizhi",type, {"disp"}); // 3D Corner type
        pp.queryarr_default("val.xloyhizhi", val , {"0.0"}); // 3D Corner value
        set(type,val,Face::XLO_YHI_ZHI);
        type.clear(); val.clear();

        pp.queryarr_default("type.xhiylozlo",type, {"disp"}); // 3D Corner type
        pp.queryarr_default("val.xhiylozlo", val , {"0.0"}); // 3D Corner value
        set(type,val,Face::XHI_YLO_ZLO);
        type.clear(); val.clear();

        pp.queryarr_default("type.xhiylozhi",type, {"disp"}); // 3D Corner type
        pp.queryarr_default("val.xhiylozhi", val , {"0.0"}); // 3D Corner value
        set(type,val,Face::XHI_YLO_ZHI);
        type.clear(); val.clear();

        pp.queryarr_default("type.xhiyhizlo",type, {"disp"}); // 3D Corner type
        pp.queryarr_default("val.xhiyhizlo", val , {"0.0"}); // 3D Corner value
        set(type,val,Face::XHI_YHI_ZLO);
        type.clear(); val.clear();

        pp.queryarr_default("type.xhiyhizhi",type, {"disp"}); // 3D Corner type
        pp.queryarr_default("val.xhiyhizhi", val , {"0.0"}); // 3D Corner value
        set(type,val,Face::XHI_YHI_ZHI);
        type.clear(); val.clear();

        pp.queryarr_default("type.ylozlo",type, {"disp"}); // 3D Edge type
        pp.queryarr_default("val.ylozlo", val , {"0.0"}); // 3D Edge value
        set(type,val,Face::YLO_ZLO);
        type.clear(); val.clear();

        pp.queryarr_default("type.ylozhi",type, {"disp"}); // 3D Edge type
        pp.queryarr_default("val.ylozhi", val , {"0.0"}); // 3D Edge value
        set(type,val,Face::YLO_ZHI);
        type.clear(); val.clear();

        pp.queryarr_default("type.yhizlo",type, {"disp"}); // 3D Edge type
        pp.queryarr_default("val.yhizlo", val , {"0.0"}); // 3D Edge value
        set(type,val,Face::YHI_ZLO);
        type.clear(); val.clear();

        pp.queryarr_default("type.yhizhi",type, {"disp"}); // 3D Edge type
        pp.queryarr_default("val.yhizhi", val , {"0.0"}); // 3D Edge value
        set(type,val,Face::YHI_ZHI);
        type.clear(); val.clear();

        pp.queryarr_default("type.zloxlo",type, {"disp"}); // 3D Edge type
        pp.queryarr_default("val.zloxlo", val , {"0.0"}); // 3D Edge value
        set(type,val,Face::ZLO_XLO);
        type.clear(); val.clear();

        pp.queryarr_default("type.zloxhi",type, {"disp"}); // 3D Edge type
        pp.queryarr_default("val.zloxhi", val , {"0.0"}); // 3D Edge value
        set(type,val,Face::ZLO_XHI);
        type.clear(); val.clear();

        pp.queryarr_default("type.zhixlo",type, {"disp"}); // 3D Edge type
        pp.queryarr_default("val.zhixlo", val , {"0.0"}); // 3D Edge value
        set(type,val,Face::ZHI_XLO);
        type.clear(); val.clear();

        pp.queryarr_default("type.zhixhi",type, {"disp"}); // 3D Edge type
        pp.queryarr_default("val.zhixhi", val , {"0.0"}); // 3D Edge value
        set(type,val,Face::ZHI_XHI);
        type.clear(); val.clear();
#endif



  



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
    Expression(IO::ParmParse &pp, std::string name): Expression()
    {pp.queryclass(name,*this);}
};
}
}
}
#endif
