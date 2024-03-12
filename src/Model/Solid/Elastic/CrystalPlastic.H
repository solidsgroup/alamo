#ifndef MODEL_SOLID_ELASTIC_CRYSTALPLASTIC_H_
#define MODEL_SOLID_ELASTIC_CRYSTALPLASTIC_H_


#include "IO/ParmParse.H"
#include "Model/Solid/Solid.H"
#include "PseudoLinearCubic.H"

namespace Model
{
namespace Solid
{
namespace Elastic
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
    
    virtual void Advance(Set::Scalar dt, Set::Matrix F, Set::Matrix P) override
    {
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
    Set::Matrix3d Fp = Set::Matrix3d::Zero();
    Eigen::Matrix<Set::Scalar,12,1> gamma = Eigen::Matrix<Set::Scalar,12,1>::Zero();

    // Parameters
    Set::Scalar m_rate_inv = 0.5;
    Set::Scalar tau_crss = 1.0;
    Set::Scalar gammadot0 = 1.0;
    

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
    }

#define OP_CLASS CrystalPlastic
#define OP_VARS X(C11) X(C12) X(C44) X(q) X(Fp) X(gamma) X(m_rate_inv) X(tau_crss) X(gammadot0)
#include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"


}
}
}

#endif
