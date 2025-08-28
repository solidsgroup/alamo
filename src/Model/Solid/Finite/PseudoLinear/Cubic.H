#ifndef MODEL_SOLID_FINITE_PSEUDOLINEAR_CUBIC_H_
#define MODEL_SOLID_FINITE_PSEUDOLINEAR_CUBIC_H_

#include "IO/ParmParse.H"
#include "Model/Solid/Solid.H"
#include "Set/Base.H"

namespace Model
{
namespace Solid
{
namespace Finite
{
namespace PseudoLinear
{
class Cubic : public Solid<Set::Sym::Major>
{
public:
    Cubic() {};
    //Cubic(Solid<Set::Sym::Major> base) : Solid<Set::Sym::Major>(base) {};
    virtual ~Cubic() {};

    Set::Scalar W(const Set::Matrix & F) const override
    {
        q.normalize();
        Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor> CC
            = Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor>::Cubic(C11,C12,C44,q);
        Set::Matrix E = 0.5 * (F.transpose() * F - Set::Matrix::Identity());
        return 0.5 * (E * (CC*E).transpose()).trace();
    }
    Set::Matrix DW(const Set::Matrix & F) const override
    {
        q.normalize();
        Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor> CC
            = Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor>::Cubic(C11,C12,C44,q);
        Set::Matrix C = F.transpose() * F;
        Set::Matrix dw = Set::Matrix::Zero();
        for (int i = 0; i < AMREX_SPACEDIM; i++)
        for (int J = 0; J < AMREX_SPACEDIM; J++)
        for (int A = 0; A < AMREX_SPACEDIM; A++)
        for (int R = 0; R < AMREX_SPACEDIM; R++)
        for (int S = 0; S < AMREX_SPACEDIM; S++)
        {
            dw(i,J) += 0.5 * F(i,A)*CC(J,A,R,S) * (C(R,S) - ((R==S) ? 1.0 : 0.0));
        }
        return dw;      
    }
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Major> DDW(const Set::Matrix & F) const override
    {
        q.normalize();
        Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor> CC
            = Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor>::Cubic(C11,C12,C44,q);
        Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Major> ddw;
        Set::Matrix C = F.transpose() * F;
        Set::Matrix I = Set::Matrix::Identity();
        
        for (int i = 0; i < AMREX_SPACEDIM; i++)
        for (int J = 0; J < AMREX_SPACEDIM; J++)
        for (int k = 0; k < AMREX_SPACEDIM; k++)
        for (int L = 0; L < AMREX_SPACEDIM; L++)
        {
            ddw(i,J,k,L) = 0.0;
            for (int A = 0; A < AMREX_SPACEDIM; A++)
                for (int B = 0; B < AMREX_SPACEDIM; B++)
                {
                    ddw(i,J,k,L) += 0.5*I(i,k)*CC(J,L,A,B) * (C(A,B) - I(A,B));
                    ddw(i,J,k,L) += F(i,A)*F(k,B)*CC(A,J,B,L);
                    if (std::isnan(ddw(i,J,k,L))) Util::Abort(INFO);
                }
        }

        return ddw;   
    }
    virtual void Print(std::ostream &out) const override 
    {
        //out << "CC = " << CC;
        out << "C11=" << C11 << " C12=" << C12 << " C44=" << C44 <<
            " q=" << q.w() << " + " << q.x() << "i + " << q.y() << "j + " << q.z() << "k";
    }
    
public:
    //Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor> CC;
    static constexpr KinematicVariable kinvar = KinematicVariable::F;
    Set::Scalar C11=1.68, C12=1.21, C44=0.75;
    mutable Set::Quaternion q = Set::Quaternion(1.0, 0.0, 0.0, 0.0);

public:
    static Cubic Zero()
    {
        Cubic ret;
        ret.C11 = 0.0;
        ret.C12 = 0.0;
        ret.C44 = 0.0;
        ret.q = Set::Quaternion(0.0,0.0,0.0,0.0);
        //ret.CC = Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor>::Zero();
        return ret;
    }
    static Cubic Random()
    {
        return Random(Util::Random(), Util::Random(), Util::Random());
    }
    static Cubic Random(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44)
    {
        Cubic ret;
        ret.C11 = C11;
        ret.C12 = C12;
        ret.C44 = C44;
        ret.q = Set::Quaternion::UnitRandom();
        return ret;
    }
    static void Parse(Cubic & value, IO::ParmParse & pp)
    {
        // Elastic constant $C_{11}$
        pp.query_default("C11",value.C11, "1.68", Unit::Pressure()); 
        // Elastic constant $C_{12}$
        pp.query_default("C12",value.C12, "1.21", Unit::Pressure()); 
        // Elastic constant $C_{44}$
        pp.query_default("C44",value.C44, "0.75", Unit::Pressure()); 

        if (pp.contains("random"))
        {
            value = Cubic::Random(value.C11,value.C12,value.C44);
            return;
        }

        if (pp.contains_any({"phi1","Phi","phi2"}))
        {
            Set::Scalar phi1 = 0.0, Phi = 0.0, phi2 = 0.0;
            // Bunge Euler angle :math:`\phi_1`
            pp.query_default("phi1",phi1,"0.0_rad", Unit::Angle());  
            // Bunge Euler angle :math:`\Phi`
            pp.query_default("Phi",Phi,"0.0_rad", Unit::Angle());    
            // Bunge Euler angle :math:`\phi_2`
            pp.query_default("phi2",phi2,"0.0_rad", Unit::Angle());  
            Eigen::Matrix3d R;
            R = Eigen::AngleAxisd(phi2, Eigen::Vector3d::UnitX()) *
                Eigen::AngleAxisd(Phi,  Eigen::Vector3d::UnitZ()) *
                Eigen::AngleAxisd(phi1, Eigen::Vector3d::UnitX());
            value.q = Set::Quaternion(R);
            pp.forbid({"X","Y","Z"},"Cannot specify both euler angles and crystallographic axes");
        }
        else if (pp.contains_any({"X","Y","Z"}))
        {
            if (pp.contains_all({"X","Y","Z"}))
            {
                Util::Warning(INFO,"We are not checking to make sure your axes are orthogonal");
            }
            else if (pp.contains_all({"X","Y"}))
            {
                Eigen::Vector3d X, Y;
                pp.queryarr("X",X,Unit::Less());
                pp.queryarr("Y",Y,Unit::Less());
                X.normalize();
                Y.normalize();
                Util::Assert(INFO, TEST(X.dot(Y) < 1E-8), "X and Y must be orthogonal");
                Eigen::Vector3d Z = X.cross(Y);
                
                Eigen::Matrix3d R;
                R.col(0) = X; R.col(1) = Y; R.col(2) = Z;
                
                value.q = Set::Quaternion(R);
            }
            else
            {
                Util::Message(INFO,"Haven't implemented this yet; for now just use X and Y");
            }
        }
    }

#define OP_CLASS Cubic
#define OP_VARS X(C11) X(C12) X(C44) X(q)
#include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"

}
}
}
}

#endif
