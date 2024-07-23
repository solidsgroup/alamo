#ifndef MODEL_SOLID_COMPOSITE_H
#define MODEL_SOLID_COMPOSITE_H


#include "Solid.H"
#include "IO/ParmParse.H"

namespace Model
{
namespace Solid
{
template<class SOLID, int NUM>
class Composite: public Solid<SOLID::sym>
{
public:
    static constexpr Set::Sym sym = SOLID::sym;
    static const KinematicVariable kinvar = SOLID::kinvar;
    int order = 1;

    virtual Set::Scalar W(const Set::Matrix& F) const override
    {
        Set::Scalar etasum = 0;
        Set::Scalar wsum = 0;
        for (int i = 0; i < NUM; i++)
        {
            if (order == 1)
            {
                etasum += m_eta[i];
                wsum += m_eta[i] * m_models[i].W(F);
            }
            else if (order == 2)
            {
                etasum += m_eta[i] * m_eta[i];
                wsum += m_eta[i] * m_eta[i] * m_models[i].W(F);
            }
            else if (order == 4)
            {
                etasum += m_eta[i] * m_eta[i] * m_eta[i] * m_eta[i];
                wsum += m_eta[i] * m_eta[i] * m_eta[i] * m_eta[i] * m_models[i].W(F);
            }
        }
        return wsum / etasum;
    };
    virtual Set::Matrix DW(const Set::Matrix& F) const override
    {
        Set::Scalar etasum = 0;
        Set::Matrix dwsum = Set::Matrix::Zero();
        for (int i = 0; i < NUM; i++)
        {
            if (order == 1)
            {
                etasum += m_eta[i];
                dwsum += m_eta[i] * m_models[i].DW(F);
            }
            else if (order == 4)
            {
                etasum += m_eta[i] * m_eta[i] * m_eta[i] * m_eta[i];
                dwsum += m_eta[i] * m_eta[i] * m_eta[i] * m_eta[i] * m_models[i].DW(F);
            }
        }
        return dwsum / etasum;
    };
    virtual Set::Matrix4<AMREX_SPACEDIM, sym> DDW(const Set::Matrix& F) const override
    {
        Set::Scalar etasum = 0;
        Set::Matrix4<AMREX_SPACEDIM, sym> ddwsum = Set::Matrix4<AMREX_SPACEDIM, sym>::Zero();
        for (int i = 0; i < NUM; i++)
        {
            if (order == 1)
            {
                etasum += m_eta[i];
                ddwsum += m_models[i].DDW(F) * m_eta[i];
            }
            else if (order == 4)
            {
                etasum += m_eta[i] * m_eta[i] * m_eta[i] * m_eta[i];
                ddwsum += m_models[i].DDW(F) * m_eta[i] * m_eta[i] * m_eta[i] * m_eta[i];
            }
        }
        return ddwsum / etasum;
    };

    ///
    /// \f[]
    Set::Scalar DWDeta(const Set::Matrix& F, int n) const
    {
        Set::Scalar etasum = 0;
        if (order == 1)
        {
            Util::Abort(INFO, "This implementation does not work");
            for (int i = 0; i < NUM; i++) etasum += m_eta[i];
            return (1. / etasum) * (m_models[n].W(F) - W(F));
        }
        else if (order == 4)
        {
            for (int i = 0; i < NUM; i++) etasum += m_eta[i] * m_eta[i] * m_eta[i] * m_eta[i];
            return (4. * m_eta[n] * m_eta[n] * m_eta[n] / etasum) * (m_models[n].W(F) - W(F));
        }
    }


    std::array<SOLID, NUM> m_models;
    std::array<Set::Scalar, NUM> m_eta;

    static Composite Zero()
    {
        Composite ret;
        for (int i = 0; i < NUM; i++)
        {
            ret.m_models[i] = SOLID::Zero();
            ret.m_eta[i] = 0.0;
        }
        return ret;
    }

    AMREX_FORCE_INLINE
        void operator+=(const Composite& rhs)
    {
        for (int i = 0; i < NUM; i++)
        {
            m_models[i] += rhs.m_models[i];
            m_eta[i] += rhs.m_eta[i];
        }
    }
    AMREX_FORCE_INLINE
        void operator-=(const Composite& rhs)
    {
        for (int i = 0; i < NUM; i++)
        {
            m_models[i] = m_models[i] - rhs.m_models[i];
            m_eta[i] -= rhs.m_eta[i];
        }
    }
    AMREX_FORCE_INLINE
        Composite operator*(const Set::Scalar alpha) const
    {
        Composite ret = Zero();
        for (int i = 0; i < NUM; i++)
        {
            ret.m_models[i] = m_models[i] * alpha;
            ret.m_eta[i] = m_eta[i] * alpha;
        }
        return ret;
    }
    // Composite operator * (const Set::Scalar alpha, const Composite b) {};
    // Composite operator + (const Composite a, const Composite b) {};
    // Composite operator - (const Composite a, const Composite b) {};
public:
    static void Parse(Composite<SOLID, NUM>& value, IO::ParmParse& pp)
    {
        for (int i = 0; i < NUM; i++)
            pp.queryclass("model" + std::to_string(i + 1), value.m_models[i]);
        pp.query("order", value.order); // order to use when mixing
    }
};

template<class SOLID, int NUM>
AMREX_FORCE_INLINE
Composite<SOLID, NUM> operator * (const Set::Scalar alpha, const Composite<SOLID, NUM> b)
{
    return b * alpha;
}

template<class SOLID, int NUM>
AMREX_FORCE_INLINE
Composite<SOLID, NUM> operator + (const Composite<SOLID, NUM> a, const Composite<SOLID, NUM> b)
{
    Composite<SOLID, NUM> ret = a;
    ret += b;
    return ret;
}

template<class SOLID, int NUM>
AMREX_FORCE_INLINE
Composite<SOLID, NUM> operator - (const Composite<SOLID, NUM> a, const Composite<SOLID, NUM> b)
{
    Composite<SOLID, NUM> ret = a;
    ret -= b;
    return ret;
}

}
}








#endif
