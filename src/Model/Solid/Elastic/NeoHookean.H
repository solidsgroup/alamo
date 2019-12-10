#ifndef MODEL_SOLID_ELASTIC_NEOHOOKEAN_H_
#define MODEL_SOLID_ELASTIC_NEOHOOKEAN_H_

#include "Model/Solid/Elastic/Elastic.H"

namespace Model
{
namespace Solid
{
namespace Elastic
{
class NeoHookean : public Elastic
{
public:
	NeoHookean() {};
	virtual ~NeoHookean() {};

    Set::Scalar W(Set::Matrix F)
    {
        Set::Scalar J = F.determinant();
        Set::Scalar J23 = std::pow(fabs(J),2./3.);
        Set::Scalar w = 0.0;
        w += 0.5 * mu * ((F*F.transpose()).trace() / J23 - 3);
        w += 0.5 * kappa * (J - 1.0) * (J - 1.0);
        return w;
    }
    Set::Matrix DW(Set::Matrix F)
    {
        Set::Scalar J = F.determinant();
        Set::Scalar J23 = std::pow(fabs(J),2./3.);
        Set::Matrix FinvT = F.inverse().transpose();

        dw = Set::Matrix::Zero();

        dw += mu * (F/J23 - (F*F.transpose()).trace()*FinvT / (3.*J23));
        dw += kappa*(J-1)*J*FinvT;

        return dw;
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
                        ddw(i,j,k,l) = 0.0;
//                      if (i==j && k==l) ddw(i,j,k,l) += kappa;
//                      if (i==k && j==l) ddw(i,j,k,l) += mu;
//                      if (i==l && j==k) ddw(i,j,k,l) += mu;

                        Set::Scalar t1 = 0.0, t2 = 0.0;

                        if (i==k && j==l) t1 += 1.0;
                        t1 -= (2./3.) * F(i,j)*FinvT(k,l);
                        t1 -= (2./3.) * FinvT(i,j)*F(k,l);
                        t1 += (2./9.) * (F*F.transpose()).trace() * FinvT(i,j) * FinvT(k,l);
                        t1 += (1./3.) * (F*F.transpose()).trace() * FinvT(i,l) * FinvT(k,j);

                        t2 += (2.*J - 1.) * FinvT(i,j)*FinvT(k,l);
                        t2 += (1. - J) * FinvT(i,l) * FinvT(k,j);

                        ddw(i,j,k,l) = (mu/J23)*t1 + kappa*J*t2;

                        if (std::isnan(ddw(i,j,k,l)) && amrex::ParallelDescriptor::IOProcessor())
                        {
                            Util::Message(INFO,"t1=",t1);
                            Util::Message(INFO,"t2=",t2);
                            Util::Message(INFO,"FinvT=",FinvT);
                            Util::Message(INFO,"J=",J);
                            Util::Message(INFO,"J23=",J23);
                            Util::Message(INFO,"kappa=",kappa);
                            Util::Message(INFO,"mu=",mu);
                            Util::Abort(INFO);
                        }
                    }
        return ddw;
    }
	
public:
    //Set::Scalar mu = 6.0, kappa = 2.6; // works
    
    //Set::Scalar mu = 0.793, kappa = 1.60; // from steel, doesn't work
    Set::Scalar mu = NAN, kappa = NAN;

	NeoHookean operator + (const NeoHookean &rhs) const
	{
        NeoHookean ret;
        ret.mu    = mu + rhs.mu;
        ret.kappa = kappa + rhs.kappa;
        ret.w     = w + rhs.w;
        ret.dw    = dw + rhs.dw;
        ret.divdw = divdw + rhs.divdw;
        ret.ddw   = ddw; ret.ddw += rhs.ddw;
        return ret;
	}

	void operator += (const NeoHookean &rhs)
	{
        mu += rhs.mu;
        kappa += rhs.kappa;
        w += rhs.w;
        dw += rhs.dw;
        divdw += rhs.divdw;
        ddw += rhs.ddw;
	}

	NeoHookean operator - (const NeoHookean &rhs) const
	{
        NeoHookean ret;
        ret.mu = mu - rhs.mu;
        ret.kappa = kappa - rhs.kappa;
        ret.w = w - rhs.w;
        ret.dw = dw - rhs.dw;
        ret.divdw = divdw - rhs.divdw;
        ret.ddw = ddw; ret.ddw -= rhs.ddw;
        return ret;
	}

	NeoHookean & operator = (const NeoHookean &rhs)
	{
		mu    = rhs.mu; 
        kappa = rhs.kappa;
        w     = rhs.w;
        dw    = rhs.dw;
        divdw = rhs.divdw;
        ddw   = rhs.ddw;
		return *this;
	}

	virtual NeoHookean operator / (const Set::Scalar alpha) const
	{
        NeoHookean ret;
        ret.mu = mu/alpha;
        ret.kappa = kappa/alpha;
        ret.w     = w;     ret.w     /= alpha;
        ret.dw    = dw;    ret.dw    /= alpha;
        ret.divdw = divdw; ret.divdw /= alpha;
        ret.ddw   = ddw;   ret.ddw   /= alpha;
        return ret;
	}

	virtual NeoHookean operator * (const Set::Scalar alpha) const
	{
        NeoHookean ret;
        ret.mu = mu*alpha;
        ret.kappa = kappa*alpha;
        ret.w     = w;     ret.w     *= alpha;
        ret.dw    = dw;    ret.dw    *= alpha;
        ret.divdw = divdw; ret.divdw *= alpha;
        ret.ddw   = ddw;   ret.ddw   *= alpha;
        return ret;
	}

public:
    static NeoHookean Random()
    {
        NeoHookean ret;
        ret.mu = Util::Random();
        ret.kappa = Util::Random();
        return ret;
    }


};



//template<Set::Sym sym>
//std::ostream& operator<< (std::ostream& os, const LinearElastic<sym>& b);

}
}
}

#endif

