#ifndef MODEL_SOLID_FINITE_CRYSTALPLASTIC_H_
#define MODEL_SOLID_FINITE_CRYSTALPLASTIC_H_


#include "IO/ParmParse.H"
#include "Model/Solid/Solid.H"
#include "PseudoLinearCubic.H"

namespace Model
{
namespace Solid
{
namespace Finite
{
class CrystalPlastic : public PseudoLinearCubic
{
public:
    CrystalPlastic() {};
    //CrystalPlastic(Solid<Set::Sym::Major> base) : Solid<Set::Sym::Major>(base) {};
    virtual ~CrystalPlastic() {};

    Set::Scalar W(const Set::Matrix & F) const override
    {
        return PseudoLinearCubic::W(F * Set::reduce(Fp).inverse());
    }
    Set::Matrix DW(const Set::Matrix & F) const override
    {
        Set::Matrix Fpinv = Set::reduce(Fp).inverse();
        return PseudoLinearCubic::DW(F * Fpinv) * Fpinv.transpose();
    }
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Major> DDW(const Set::Matrix & F) const override
    {
        Set::Matrix Fpinv = Set::reduce(Fp).inverse();
        auto DDW = PseudoLinearCubic::DDW(F * Fpinv);
        auto ret = Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Major>::Zero();
        for (int i = 0; i < AMREX_SPACEDIM; i++)
            for (int J = 0; J < AMREX_SPACEDIM; J++)
                for (int k = 0; k < AMREX_SPACEDIM; k++)
                    for (int L = 0; L < AMREX_SPACEDIM; L++)
                    {
                        ret(i,J,k,L) = 0.0;
                        for (int Q = 0; Q < AMREX_SPACEDIM; Q++)
                            for (int R = 0; R < AMREX_SPACEDIM; R++)
                                ret(i,J,k,L) += DDW(i,Q,k,R)*Fpinv(J,Q)*Fpinv(L,R);
                    }
        return ret;
    }
    
    virtual void Advance(Set::Scalar dt, Set::Matrix F, Set::Matrix P, Set::Scalar time) override
    {
        if (time  < tstart) return;

        q.normalize();
        Set::Matrix3d R = q.toRotationMatrix();
        std::array<std::pair<Set::Vector3d,Set::Vector3d>,12> ss = slipSystems(R);



        Set::Matrix3d L = Set::Matrix3d::Zero();

        for (int n = 0; n < 12; n++)
        {
            Set::Vector3d A = ss[n].first;
            Set::Vector3d N = ss[n].second;

            Set::Scalar tau = A.dot(Set::expand(P)*N);

            Set::Scalar sign = tau>0 ? 1.0 : -1.0;

            Set::Scalar gammadot
                = gammadot0 * std::pow(fabs(tau/tau_crss),m_rate_inv) * sign;
            
            L += gammadot * A * N.transpose();

            gamma[n] += gammadot * dt;

        }

        Fp += L * Fp * dt;

    }


public:
    //Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor> CC;
    //static constexpr KinematicVariable kinvar = KinematicVariable::F;
    //Set::Scalar C11=1.68, C12=1.21, C44=0.75;
    //mutable Set::Quaternion q = Set::Quaternion(1.0, 0.0, 0.0, 0.0);

    // Plastic state variables
    Set::Matrix3d Fp = Set::Matrix3d::Identity();
    Eigen::Matrix<Set::Scalar,12,1> gamma = Eigen::Matrix<Set::Scalar,12,1>::Zero();

    // Parameters
    Set::Scalar m_rate_inv = 0.5;
    Set::Scalar tau_crss = 1.0;
    Set::Scalar gammadot0 = 1.0;

    // Control variables
    Set::Scalar tstart = 0.0;

private:
    static std::array<std::pair<Set::Vector3d,Set::Vector3d>,12> slipSystems(Set::Matrix3d R)
    {
        const std::array<Set::Vector3d,12> n0 =
            {
             Set::Vector3d::Constant( 1, 1, 1) / sqrt(3.0),
             Set::Vector3d::Constant( 1, 1, 1) / sqrt(3.0),
             Set::Vector3d::Constant( 1, 1, 1) / sqrt(3.0),
             Set::Vector3d::Constant(-1, 1,-1) / sqrt(3.0),
             Set::Vector3d::Constant(-1, 1,-1) / sqrt(3.0),
             Set::Vector3d::Constant(-1, 1,-1) / sqrt(3.0),
             Set::Vector3d::Constant( 1,-1,-1) / sqrt(3.0),
             Set::Vector3d::Constant( 1,-1,-1) / sqrt(3.0),
             Set::Vector3d::Constant( 1,-1,-1) / sqrt(3.0),
             Set::Vector3d::Constant(-1,-1, 1) / sqrt(3.0),
             Set::Vector3d::Constant(-1,-1, 1) / sqrt(3.0),
             Set::Vector3d::Constant(-1,-1, 1) / sqrt(3.0)
            };
        const std::array<Set::Vector3d,12> a0 =
            {
             Set::Vector3d::Constant( 1, 0,-1) / sqrt(2.0),
             Set::Vector3d::Constant( 0,-1, 1) / sqrt(2.0),
             Set::Vector3d::Constant( 1,-1, 0) / sqrt(2.0),
             Set::Vector3d::Constant( 1, 0,-1) / sqrt(2.0),
             Set::Vector3d::Constant( 1, 1, 0) / sqrt(2.0),
             Set::Vector3d::Constant( 0, 1, 1) / sqrt(2.0),
             Set::Vector3d::Constant( 1, 1, 0) / sqrt(2.0),
             Set::Vector3d::Constant( 0,-1, 1) / sqrt(2.0),
             Set::Vector3d::Constant( 1, 0, 1) / sqrt(2.0),
             Set::Vector3d::Constant( 0, 1, 1) / sqrt(2.0),
             Set::Vector3d::Constant( 1, 0, 1) / sqrt(2.0),
             Set::Vector3d::Constant( 1,-1, 0) / sqrt(2.0)
            };


        static std::array<std::pair<Set::Vector3d,Set::Vector3d>,12> ret;
        for (int n = 0 ; n < 12 ; n ++)
        {
            ret[n].first = R*a0[n];
            ret[n].second = R*n0[n];
        }
        return ret;
    }

public:
    static CrystalPlastic Zero()
    {
        CrystalPlastic ret;
        ret.C11 = 0.0;
        ret.C12 = 0.0;
        ret.C44 = 0.0;
        ret.q = Set::Quaternion(0.0,0.0,0.0,0.0);
        ret.Fp = Set::Matrix3d::Zero();
        ret.gamma = Eigen::Matrix<Set::Scalar,12,1>::Zero();
        ret.m_rate_inv = 0.0;
        ret.tau_crss = 0.0;
        ret.gammadot0 = 0.0;
        return ret;
    }
    static CrystalPlastic Random()
    {
        return Random(Util::Random(), Util::Random(), Util::Random());
    }
    static CrystalPlastic Random(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44)
    {
        CrystalPlastic ret = CrystalPlastic::Zero();
        ret.C11 = C11;
        ret.C12 = C12;
        ret.C44 = C44;
        ret.q = Set::Quaternion::UnitRandom();
        return ret;
    }
    static void Parse(CrystalPlastic & value, IO::ParmParse & pp)
    {
        PseudoLinearCubic::Parse(value,pp);
        // TODO Add inputs for the CP values

        pp.query("tstart",value.tstart);    // Time to activate plastic slip

    }

#define OP_CLASS CrystalPlastic
#define OP_VARS X(C11) X(C12) X(C44) X(q) X(Fp) X(gamma) X(m_rate_inv) X(tau_crss) X(gammadot0) X(tstart)
#include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"





}
}
}


template<>
inline int Set::Field<Model::Solid::Finite::CrystalPlastic>::NComp() const 
{
    return 12;
}

template<>
inline std::string Set::Field<Model::Solid::Finite::CrystalPlastic>::Name(int i) const 
{
    if (i==0) return name + "gamma1";
    if (i==1) return name + "gamma2";
    if (i==2) return name + "gamma3";
    if (i==3) return name + "gamma4";
    if (i==4) return name + "gamma5";
    if (i==5) return name + "gamma6";
    if (i==6) return name + "gamma7";
    if (i==7) return name + "gamma8";
    if (i==8) return name + "gamma9";
    if (i==9) return name + "gamma10";
    if (i==10) return name + "gamma11";
    if (i==11) return name + "gamma12";
    return name;
}

template<>
inline void Set::Field<Model::Solid::Finite::CrystalPlastic>::Copy(int a_lev, amrex::MultiFab &a_dst, int a_dstcomp, int a_nghost) const
{
    for (amrex::MFIter mfi(a_dst, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox(amrex::IntVect(a_nghost));
        if (bx.ok())
        {
            amrex::Array4<const Model::Solid::Finite::CrystalPlastic> const & src = ((*this)[a_lev])->array(mfi);
            amrex::Array4<Set::Scalar> const & dst = a_dst.array(mfi);
            for (int n = 0; n < 12; n++)
            {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    dst(i,j,k,a_dstcomp + n) = src(i,j,k).gamma(n);
                });
            }
        }
    }    
}



#endif
