#ifndef MODEL_SOLID_H_
#define MODEL_SOLID_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>

#include "Set/Set.H"

namespace Model
{
namespace Solid
{

enum KinematicVariable{gradu,epsilon,F};

template<int Sym>
class Solid
{
public:
	Solid() {} ;
	virtual ~Solid() {};

    AMREX_FORCE_INLINE
	virtual Set::Matrix operator () (Set::Matrix &F, bool /*a_homogeneous=true*/=true) const
    {
        return ddw*F;
    };

    AMREX_FORCE_INLINE
	virtual Set::Vector operator () (Set::Matrix3 &gradF, bool /*a_homogeneous=true*/=true)
    {
        return ddw*gradF;
    };

    virtual Set::Scalar W(Set::Matrix)           {Util::Abort(INFO,"W not implemented"); return 0.0;};
    virtual Set::Matrix DW(Set::Matrix)          {Util::Abort(INFO,"DW not implemented"); return Set::Matrix::Zero();};
    virtual Set::Matrix4<AMREX_SPACEDIM,Sym> DDW(Set::Matrix) {Util::Abort(INFO,"DDW not implemented"); return ddw;};
	
public:
    mutable Set::Matrix4<AMREX_SPACEDIM,Sym> ddw;

    static const KinematicVariable kinvar = KinematicVariable::F;

   	Solid operator + (const Solid &rhs) const
	{
        Solid ret;
        ret.ddw   = ddw; ret.ddw += rhs.ddw;
        return ret;
	}
   	virtual Solid operator * (const Set::Scalar alpha) const
	{
        Solid ret;
        ret.ddw   = ddw;   ret.ddw   *= alpha;
        return ret;
	}
	virtual Solid operator / (const Set::Scalar alpha) const
	{
        Solid ret;
        ret.ddw   = ddw;   ret.ddw   /= alpha;
        return ret;
	}
	void operator += (const Solid &rhs)
	{
        ddw += rhs.ddw;
	}

	Solid operator - (const Solid &rhs) const
	{
        Solid ret;
        ret.ddw = ddw; ret.ddw -= rhs.ddw;
        return ret;
	}

   	friend std::ostream& operator<<(std::ostream &out, const Solid &a)
    {
        a.Print(out);
        return out;
    }
    virtual void Print(std::ostream &out) const
    {
        out << "No print function written for this model.";
    }


public: 
    template <class T>
    static int DerivativeTest1(int verbose = 0)
    {
        for (int iter = 0; iter < 10; iter++)
        {
            T model = T::Random();

            Set::Scalar dx = 1E-8, tol = 1E-6;

            Set::Matrix F = Set::Matrix::Random();
            while (fabs(F.determinant()) < 0.10) F = Set::Matrix::Random(); // Ensure that F in GL(3)

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

            Set::Scalar dx = 1E-8, tol = 1E-6;
            Set::Matrix F = Set::Matrix::Random();
            while (fabs(F.determinant()) < 1.0) F = Set::Matrix::Random(); // Ensure that F in GL(3)

            Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Major> ddw_exact = model.DDW(F);
            Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Major> ddw_numeric = Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Major>::Zero();
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
            Set::Scalar relnorm = error.norm()/ddw_numeric.norm();
            if (relnorm > tol || std::isnan(relnorm) || std::isinf(relnorm))
            {
                if (verbose)
                {
                    Util::Message(INFO,"F \n",F);
                    Util::Message(INFO,"det(F) = ",F.determinant());
                    Util::Message(INFO,"ddw exact \n",ddw_exact);
                    Util::Message(INFO,"ddw numeric \n",ddw_numeric);
                    Util::Message(INFO,"error norm \n",relnorm);
                }
                return 1;
            }   
        }
        return 0;
    }

};

}
}


#endif


