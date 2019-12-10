#ifndef MODEL_SOLID_ELASTIC_H_
#define MODEL_SOLID_ELASTIC_H_

#include "Model/Solid/Solid.H"
#include "Set/Set.H"

namespace Model
{
namespace Solid
{
namespace Elastic
{
class Elastic : public Solid
{
public:
	Elastic() {};
	virtual ~Elastic() {};
    
    AMREX_FORCE_INLINE
	virtual Set::Matrix operator () (Set::Matrix &F, bool a_homogeneous=true) const
    {
        //return mu*(F + F.transpose()) + kappa*Set::Matrix::Identity()*F.trace();
        //Util::Message(INFO,"begin");
        Set::Matrix P;
        for (int i = 0; i < 3; i++)
            for (int J = 0; J < 3; J++)
            {
                P(i,J) = 0.0;
                for (int k = 0; k < 3; k++)
                    for (int L = 0; L < 3; L++)
                        P(i,J) += ddw(i,J,k,L) * F(k,L);
                if (std::isnan(P(i,J)))
                {
                    Util::Message(INFO,ddw);
                    Util::Message(INFO,F);
                    Util::Abort(INFO);
                }

            }

        if (!a_homogeneous) P += dw;
        //Util::Message(INFO,"end");
        return P;

    };

    AMREX_FORCE_INLINE
	virtual Set::Vector operator () (Set::Matrix3 &gradF, bool a_homogeneous=true)
    {
//        //Util::Message(INFO,"begin");
//		Set::Vector _f = Set::Vector::Zero();
//		for (int i = 0; i < AMREX_SPACEDIM; i++)
//			for (int j=0; j < AMREX_SPACEDIM; j++)
//				{
//					_f(i) += mu*(gradF(i,j,j) + gradF(j,i,j))  + kappa*gradF(j,j,i);
//				}
		//return _f;

        Set::Vector f;
        for (int i = 0; i < 3; i++)
        {
            f(i) = 0.0;
            for (int J = 0; J < 3; J++)
                for (int k = 0; k < 3; k++)
                    for (int L = 0; L < 3; L++)
                    {
                        //if (i==J && k==L) f(i) += kappa * gradF(k,L,J);
                        //if (i==k && J==L) f(i) += mu * gradF(k,L,J);
                        //if (i==L && J==k) f(i) += mu * gradF(k,L,J);
                        f(i) += ddw(i,J,k,L) * gradF(k,L,J);
                    }
            if (std::isnan(f(i)))
            {
                Util::Message(INFO,ddw);
                Util::Abort(INFO,f.transpose());
            }
        }    

        if (!a_homogeneous) f += divdw;
//      if ((f-_f).norm() > 1E-8) Util::Abort(INFO,"_f = ",_f.transpose(),"\nf = ",f.transpose());
        //Util::Message(INFO,"end");
        return f;
    };

    virtual Set::Scalar W(Set::Matrix F) = 0;
    virtual Set::Matrix DW(Set::Matrix F) = 0;
    virtual Set::Matrix4<3,Set::Sym::Major> DDW(Set::Matrix F) = 0;
	
public:
    Set::Scalar                               w;
    Set::Matrix                              dw;
    Set::Vector                           divdw;
    mutable Set::Matrix4<3,Set::Sym::Major> ddw;

public: 
    template <class T>
    static int DerivativeTest1(int verbose = 0)
    {
        for (int iter = 0; iter < 10; iter++)
        {
            T model;
            model.mu = Util::Random();
            model.kappa = Util::Random();

            Set::Scalar dx = 1E-8, tol = 1E-6;

            Set::Matrix F = Set::Matrix::Random();
            while (fabs(F.determinant()) < 0.10) F = Set::Matrix::Random(); // Ensure that F in GL(3)

            Set::Matrix dw_exact = model.DW(F);
            Set::Matrix dw_numeric = Set::Matrix::Zero();
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
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
                    Util::Message(INFO,"mu    = ",model.mu);
                    Util::Message(INFO,"kappa = ",model.kappa);
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
            T model;
            model.mu = Util::Random();
            model.kappa = Util::Random();

            Set::Scalar dx = 1E-8, tol = 1E-6;
            Set::Matrix F = Set::Matrix::Random();
            while (fabs(F.determinant()) < 1.0) F = Set::Matrix::Random(); // Ensure that F in GL(3)

            Set::Matrix4<3,Set::Sym::Major> ddw_exact = model.DDW(F);
            Set::Matrix4<3,Set::Sym::Major> ddw_numeric = Set::Matrix4<3,Set::Sym::Major>::Zero();
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                        for (int l = 0; l < 3; l++)
                        {
                            Set::Matrix dF = Set::Matrix::Zero();
                            dF(k,l) = dx;
                            ddw_numeric(i,j,k,l) = (model.DW(F+dF) - model.DW(F-dF))(i,j) / (2.0 * dx);
                        }
            Set::Matrix4<3,Set::Sym::Major> error = ddw_numeric-ddw_exact;
            Set::Scalar relnorm = error.norm()/ddw_numeric.norm();
            if (relnorm > tol || std::isnan(relnorm) || std::isinf(relnorm))
            {
                if (verbose)
                {
                    Util::Message(INFO,"mu    = ",model.mu);
                    Util::Message(INFO,"kappa = ",model.kappa);
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



//template<Set::Sym sym>
//std::ostream& operator<< (std::ostream& os, const LinearElastic<sym>& b);

}
}
}

#endif

