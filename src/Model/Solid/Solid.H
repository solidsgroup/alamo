//
// Solid models are used with the :ref:`Integrator::Mechanics` integrator, which
// implements the :ref:`Solver::Nonlocal::Newton` elasticity solver.
// All solid models inherit from the :code:`Model::Solid` base class, which requires
// all of the necessary components to be used in a mechanics solve.
// Model classes have basically two functions:
// 
// #. Provide energy (W), stress (DW), and modulus (DDW) based on a kinematic variable
// #. Evolve internal variables in time.
//

#ifndef MODEL_SOLID_H_
#define MODEL_SOLID_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>
#include <type_traits>

#include "Set/Set.H"
#include "Set/Matrix4.H"
#include "Set/Matrix4_Diagonal.H"
#include "Set/Matrix4_Full.H"
#include "Set/Matrix4_Isotropic.H"
#include "Set/Matrix4_Major.H"
#include "Set/Matrix4_MajorMinor.H"
#include "Util/Util.H"

namespace Model
{
namespace Solid
{

enum KinematicVariable{gradu,epsilon,F};

template<Set::Sym SYM>
class Solid
{
public:
    AMREX_GPU_HOST_DEVICE Solid() = default;
    AMREX_GPU_HOST_DEVICE Solid(const Solid & base) = default;
    AMREX_GPU_HOST_DEVICE virtual ~Solid() = default;

    static constexpr Set::Sym sym = SYM;

    AMREX_GPU_HOST_DEVICE virtual Set::Scalar W(const Set::Matrix &) const      
    {Util::Abort(INFO,"W not implemented"); return 0.0;};
    AMREX_GPU_HOST_DEVICE virtual Set::Matrix DW(const Set::Matrix &) const
    {
        Util::Abort(INFO, "DW not implemented");
        return Set::Matrix::Zero();
    };
    AMREX_GPU_HOST_DEVICE virtual Set::Matrix4<AMREX_SPACEDIM,SYM> DDW(const Set::Matrix &) const
    {
        Util::Abort(INFO, "DDW not implemented");
        return Set::Matrix4<AMREX_SPACEDIM,SYM>::Zero();
    };
    
    AMREX_GPU_HOST_DEVICE
    virtual void Advance(Set::Scalar /*dt*/, Set::Matrix /*strain*/, Set::Matrix /*stress*/, Set::Scalar /*time*/) {}

    // Todo: implement for all models and make this pure virtual
    AMREX_GPU_HOST_DEVICE
    virtual bool ContainsNan() {return true;};


public:
    static const KinematicVariable kinvar = KinematicVariable::F;


    friend std::ostream& operator<<(std::ostream &out, const Solid &a)
    {
        a.Print(out);
        return out;
    }
    AMREX_GPU_HOST_DEVICE
    virtual void Print(std::ostream &out) const
    {
        AMREX_IF_ON_HOST( out << "No print function written for this model.";)
    }

public: 
    template <class T>
    static int PODTest(int verbose = 0)
    {
        int ret = 0;
        if (!std::is_trivially_copyable<T>::value)
        {
            if (verbose) Util::Test::SubWarning("Model is not trivially copyable");
        }
        if (std::has_virtual_destructor<T>::value)
        {
            if (verbose) Util::Test::SubWarning("Model has a virtual destructor");
        }
        if (std::is_polymorphic<T>::value)
        {
            if (verbose) Util::Test::SubWarning("Model is polymorphic");
        }
        if (!std::is_standard_layout<T>::value)
        {
            if (verbose) Util::Test::SubWarning("Model does not have standard layout.");
        }
        return ret;
    }


    template <class T>
    static int ArithmeticTest(int verbose = 0)
    {
        const T z = T::Zero();
        if (!(z==z) && verbose) 
        {
            Util::Message(INFO,"Zero() failed: z = \n",z); return 1;
        }
        
        const T a = T::Random();
        if (!(a==a) && verbose) 
        {
            Util::Message(INFO,"Random() failed: a = \n",a); return 1;
        }

        T c = a;
        if (!(c==a) && verbose) {Util::Message(INFO); return 1;}
        c = 1.0 * a;
        if (!(c==a) && verbose) {Util::Message(INFO); return 1;}
        c = 0.0 * a;
        if (!(c==z) && verbose) 
        {
            Util::Message(INFO,"0.0*a = \n",c); 
            Util::Message(INFO,"Zero() returns \n",z); 
            return 1;
        }

        // Todo: add some more arithmetic checks in here

        return 0;
    }

    template <class T>
    static int DerivativeTest1(int verbose = 0)
    {
        for (int iter = 0; iter < 10; iter++)
        {
            T model = T::Random();

            Set::Scalar dx = 1E-8, tol = 1E-6;

            Set::Matrix F = Set::Matrix::Random();
            while (F.determinant() < 0.1) F = Set::Matrix::Random(); // Ensure that F in GL(3)

            Set::Matrix dw_exact = model.DW(F);
            Set::Matrix dw_numeric = Set::Matrix::Zero();
            for (int i = 0; i < AMREX_SPACEDIM; i++)
                for (int j = 0; j < AMREX_SPACEDIM; j++)
                {
                    Set::Matrix dF = Set::Matrix::Zero();
                    dF(i,j) = dx;
                    dw_numeric(i,j) = (model.W(F+dF) - model.W(F-dF)) / (2.0 * dx);
                }
            Set::Scalar relnorm = (dw_numeric-dw_exact).norm()/(dw_numeric.norm());
            if (relnorm > tol || std::isnan(relnorm) || std::isinf(relnorm))
            {
                if (verbose)
                {
                    Util::Message(INFO,"F \n",F);
                    Util::Message(INFO,"det(F) = ",F.determinant());
                    Util::Message(INFO,"dw exact \n",dw_exact);
                    Util::Message(INFO,"dw numeric \n",dw_numeric);
                    Util::Message(INFO,"error norm \n",relnorm);
                }
                return 1;
            }   
        }
        return 0;
    }

    template <class T>
    static int DerivativeTest2(int verbose = 0)
    {
        for (int iter = 0; iter < 10; iter++)
        {
            T model = T::Random();

            Set::Scalar dx = 1E-8, tol = 1E-4;
            Set::Matrix F = Set::Matrix::Random();
            while (F.determinant() < 0.1) F = Set::Matrix::Random(); // Ensure that F in GL(3)

            Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Major> ddw_exact = model.DDW(F);
            Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Major> ddw_numeric = Set::Matrix4<AMREX_SPACEDIM,T::sym>::Zero();
            for (int i = 0; i < AMREX_SPACEDIM; i++)
                for (int j = 0; j < AMREX_SPACEDIM; j++)
                    for (int k = 0; k < AMREX_SPACEDIM; k++)
                        for (int l = 0; l < AMREX_SPACEDIM; l++)
                        {
                            Set::Matrix dF = Set::Matrix::Zero();
                            dF(k,l) = dx;
                            ddw_numeric(i,j,k,l) = (model.DW(F+dF) - model.DW(F-dF))(i,j) / (2.0 * dx);
                        }
            Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Major> error = ddw_numeric-ddw_exact;
            Set::Scalar relnorm = error.Norm()/ddw_numeric.Norm();
            if (relnorm > tol || std::isnan(relnorm) || std::isinf(relnorm))
            {
                if (verbose)
                {
                    Util::Message(INFO,"F \n",F);
                    Util::Message(INFO,"det(F) = ",F.determinant());
                    Util::Message(INFO,"ddw exact \n",ddw_exact);
                    Util::Message(INFO,"ddw numeric \n",ddw_numeric);
                    Util::Message(INFO,"ddw difference \n",ddw_exact - ddw_numeric);
                    Util::Message(INFO,"error norm \n",relnorm);
                }
                return 1;
            }   
        }
        return 0;
    }

    template <class T>
    static int MaterialFrameIndifference(int verbose = 0)
    {
        if (T::kinvar != Model::Solid::KinematicVariable::F)
            Util::Abort(INFO,"Attempting to test material frame indifference in a non-finite-kinematics model will fail");

        for (int iter = 0; iter < 10; iter++)
        {
            T model = T::Random();

            Set::Scalar tol = 1E-4;
            Set::Matrix F = Set::Matrix::Random();
            while (
                F.determinant() < 0.5 || // skip if determinant negative or too small
                F.determinant() > 2 // also skip if determinant too large
                )
                F = Set::Matrix::Random(); // Ensure that F is reasonably well-behaved.

            Set::Scalar W_unrotated = model.W(F);

            for (int jter = 0; jter < 10; jter++)
            {
                #if AMREX_SPACEDIM == 2
                Set::Scalar theta = 2.0 * Set::Constant::Pi * Util::Random();
                Set::Matrix R;
                R(0,0) =  cos(theta);
                R(0,1) =  sin(theta);
                R(1,0) = -sin(theta);
                R(1,1) =  cos(theta);
                #elif AMREX_SPACEDIM == 3
                Set::Scalar phi1 = 2.0 * Set::Constant::Pi *(Util::Random() - 0.5);
                Set::Scalar Phi  = Set::Constant::Pi * Util::Random();
                Set::Scalar phi2 = 2.0 * Set::Constant::Pi *(Util::Random() - 0.5);
                Eigen::Matrix3d R;
                R = Eigen::AngleAxisd(phi2, Eigen::Vector3d::UnitX()) *
                    Eigen::AngleAxisd(Phi,  Eigen::Vector3d::UnitZ()) *
                    Eigen::AngleAxisd(phi1, Eigen::Vector3d::UnitX());
                #endif

                Util::Assert(INFO,TEST( fabs(R.determinant() - 1.0) < tol ));
                Util::Assert(INFO,TEST( fabs((R.transpose() * R - Set::Matrix::Identity()).norm()) < tol ));
                
                Set::Scalar W_rotated = model.W(R*F);
                Set::Scalar error = fabs(W_rotated - W_unrotated) / fabs(W_rotated + W_unrotated + tol);

                if (error > tol)
                {
                    if (verbose)
                    {
                        Util::Message(INFO,"F = \n",F);
                        Util::Message(INFO,"R = \n",R);
                        Util::Message(INFO,"W unrotated = ",W_unrotated);
                        Util::Message(INFO,"W rotated = ",W_rotated);
                    }
                    return 1;
                }
            }
        }
        return 0;
    }

};

}
}


#endif
