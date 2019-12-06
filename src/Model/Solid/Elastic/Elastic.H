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
                    Util::Abort(INFO);
                }

            }
        return P;

    };

    AMREX_FORCE_INLINE
	virtual Set::Vector operator () (Set::Matrix3 &gradF, bool a_homogeneous=true)
    {
		//Set::Vector _f = Set::Vector::Zero();
		//for (int i = 0; i < AMREX_SPACEDIM; i++)
		//	for (int j=0; j < AMREX_SPACEDIM; j++)
		//		{
		//			_f(i) += mu*(gradF(i,j,j) + gradF(j,i,j))  + kappa*gradF(j,j,i);
		//		}
		//return _f;

        Set::Vector f;
        for (int i = 0; i < 3; i++)
        {
            f(i) = 0.0;
            for (int J = 0; J < 3; J++)
                for (int k = 0; k < 3; k++)
                    for (int L = 0; L < 3; L++)
                        f(i) += ddw(i,J,k,L) * gradF(i,J,L);
            if (std::isnan(f(i)))
            {
                Util::Message(INFO,ddw);
                Util::Abort(INFO,f.transpose());
            }
        }    
        return f;

    };

    Set::Scalar W(Set::Matrix F)
    {
        return 0.0;

        //Set::Scalar J = F.determinant();
        //Set::Scalar J23 = std::pow(fabs(J),2./3.);
        //Set::Scalar w = 0.0;
        //w += 0.5 * mu * ((F*F.transpose()).trace() / J23 - 3);
        //w += 0.5 * kappa * (J - 1.0) * (J - 1.0);
        //return w;
    }
    Set::Matrix DW(Set::Matrix F)
    {

        

//        Set::Scalar J = F.determinant();
//        Set::Scalar J23 = std::pow(fabs(J),2./3.);
//        Set::Matrix FinvT = F.inverse().transpose();
//
//        dw = Set::Matrix::Zero();
//
//        dw += mu * (F/J23 - (F*F.transpose()).trace()*FinvT / (3.*J23));
//        dw += kappa*(J-1)*J*FinvT;
//
//        return dw;
    }
    Set::Matrix4<3,Set::Sym::Major> DDW(Set::Matrix F)
    {
        Set::Matrix4<3,Set::Sym::Major> ddw;
        Set::Scalar J = F.determinant();
        Set::Scalar J23 = std::pow(fabs(J),2./3.);
        Set::Matrix FinvT = F.inverse().transpose();
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                    for (int l = 0; l < 3; l++)
                    {
//                        Set::Scalar t1 = 0.0, t2 = 0.0;
//
//                        if (i==k && j==l) t1 += 1.0;
//                        t1 -= (2./3.) * F(i,j)*FinvT(k,l);
//                        t1 -= (2./3.) * FinvT(i,j)*F(k,l);
//                        t1 += (2./9.) * (F*F.transpose()).trace() * FinvT(i,j) * FinvT(k,l);
//                        t1 += (1./3.) * (F*F.transpose()).trace() * FinvT(i,l) * FinvT(k,j);
//
//                        t2 += (2.*J - 1.) * FinvT(i,j)*FinvT(k,l);
//                        t2 += (1. - J) * FinvT(i,l) * FinvT(k,j);
//
//                        ddw(i,j,k,l) = (mu/J23)*t1 + kappa*J*t2;
                    }
        return ddw;
    }
	
public:
    Set::Scalar mu = 6.0, kappa = 2.6;

    Set::Scalar                               w;
    Set::Matrix                              dw;
    Set::Vector                           divdw;
    mutable Set::Matrix4<3,Set::Sym::Major> ddw;

public: 
    static int DerivativeTest1(int verbose = 0)
    {
        for (int iter = 0; iter < 10; iter++)
        {
            Elastic model;

            Set::Scalar dx = 1E-8, tol = 1E-6;

            Set::Matrix F = Set::Matrix::Random();
            while (fabs(F.determinant()) < 1.0) F = Set::Matrix::Random(); // Ensure that F in GL(3)

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
            if (relnorm > tol)
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
    static int DerivativeTest2(int verbose = 0)
    {
        for (int iter = 0; iter < 10; iter++)
        {
            Elastic model;
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
            if (relnorm > tol)
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


	Elastic operator + (const Elastic &rhs) const
	{
        Elastic ret;
        ret.mu    = mu + rhs.mu;
        ret.kappa = kappa + rhs.kappa;
        ret.w     = w + rhs.w;
        ret.dw    = dw + rhs.dw;
        ret.divdw = divdw + rhs.divdw;
        ret.ddw   = ddw; ret.ddw += rhs.ddw;
        return ret;
	}

	void operator += (const Elastic &rhs)
	{
        mu += rhs.mu;
        kappa += rhs.kappa;
        w += rhs.w;
        dw += rhs.dw;
        divdw += rhs.divdw;
        ddw += rhs.ddw;
        
		//mu+=rhs.mu; lambda+=rhs.lambda; Fgb += rhs.Fgb;
		//gradFgb += rhs.gradFgb;
	}

	Elastic operator - (const Elastic &rhs) const
	{
        Elastic ret;
        ret.mu = mu - rhs.mu;
        ret.kappa = kappa - rhs.kappa;
        ret.w = w - rhs.w;
        ret.dw = dw - rhs.dw;
        ret.divdw = divdw - rhs.divdw;
        ret.ddw = ddw; ret.ddw -= rhs.ddw;
        return ret;
	}

	Elastic & operator = (const Elastic &rhs)
	{
		mu    = rhs.mu; 
        kappa = rhs.kappa;
        w     = rhs.w;
        dw    = rhs.dw;
        divdw = rhs.divdw;
        ddw   = rhs.ddw;
		return *this;
	}

	virtual Elastic operator / (const Set::Scalar alpha) const
	{
        Elastic ret;
        ret.mu = mu/alpha;
        ret.kappa = kappa/alpha;
        ret.w     = w;     ret.w     /= alpha;
        ret.dw    = dw;    ret.dw    /= alpha;
        ret.divdw = divdw; ret.divdw /= alpha;
        ret.ddw   = ddw;   ret.ddw   /= alpha;
        return ret;
	}

	virtual Elastic operator * (const Set::Scalar alpha) const
	{
        Elastic ret;
        ret.mu = mu*alpha;
        ret.kappa = kappa*alpha;
        ret.w     = w;     ret.w     *= alpha;
        ret.dw    = dw;    ret.dw    *= alpha;
        ret.divdw = divdw; ret.divdw *= alpha;
        ret.ddw   = ddw;   ret.ddw   *= alpha;
        return ret;
	}


};



//template<Set::Sym sym>
//std::ostream& operator<< (std::ostream& os, const LinearElastic<sym>& b);

}
}
}

#endif

