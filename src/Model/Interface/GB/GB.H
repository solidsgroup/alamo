#ifndef MODEL_INTERFACE_GB_H
#define MODEL_INTERFACE_GB_H

#include <AMReX.H>
#include <AMReX_AmrCore.H>

#include <eigen3/Eigen/Eigenvalues>

#include <iostream>
#include <fstream>

namespace Model
{
namespace Interface
{
namespace GB
{
class GB
{
    public:
    GB() {};
    virtual Set::Scalar W(const Set::Scalar theta) const = 0;
    virtual Set::Scalar DW(const Set::Scalar theta) const = 0;
    virtual Set::Scalar DDW(const Set::Scalar theta) const = 0;
    virtual Set::Scalar W(const Set::Vector &a_n) const {return W(atan2(a_n(1),a_n(0)));};
    virtual Set::Scalar DW(const Set::Vector &, const Set::Vector &) const {return NAN;};
    virtual Set::Scalar DDW(const Set::Vector &, const Set::Vector &) const {return NAN;};

    void ExportToFile(std::string filename, amrex::Real dTheta)
    {
        std::ofstream outFile;
        outFile.open(filename);
    
        for(amrex::Real theta=0; theta<2*pi; theta=theta+dTheta)
            {
                outFile << theta << " " << W(theta) << std::endl;   
            }
        outFile.close();
    }

    enum Regularization{Wilhelm,K23};

    // Return anisotropic driving force, regularization
    std::tuple<Set::Scalar, Set::Scalar>
    DrivingForce(const Set::Vector &Deta, const Set::Matrix &DDeta, const Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Full> &DDDDeta)
    {
        #if AMREX_SPACEDIM == 2
            Set::Scalar Theta = atan2(Deta(1),Deta(0));
            Set::Scalar sigma = /*pf.l_gb*0.75**/W(Theta);
            Set::Scalar Dsigma = /*pf.l_gb*0.75**/DW(Theta);
            Set::Scalar DDsigma = /*pf.l_gb*0.75**/DDW(Theta);
            Set::Scalar sinTheta = sin(Theta);
            Set::Scalar cosTheta = cos(Theta);
            Set::Scalar sin2Theta = sinTheta*sinTheta;
            Set::Scalar cos2Theta = cosTheta*cosTheta;
            Set::Scalar cosThetasinTheta = cosTheta*sinTheta;

            Set::Scalar boundary_term = 
                - (sigma*DDeta.trace() +
                   Dsigma*(cos(2.0*Theta)*DDeta(0,1) + 0.5*sin(2.0*Theta)*(DDeta(1,1) - DDeta(0,0))) + 
                   0.5*DDsigma*(sin2Theta*DDeta(0,0) - 2.*cosThetasinTheta*DDeta(0,1) + cos2Theta*DDeta(1,1)));

            Set::Scalar curvature_term =
                DDDDeta(0,0,0,0)*(    sin2Theta*sin2Theta) +
                DDDDeta(0,0,0,1)*(4.0*sin2Theta*cosThetasinTheta) +
                DDDDeta(0,0,1,1)*(6.0*sin2Theta*cos2Theta) +
                DDDDeta(0,1,1,1)*(4.0*cosThetasinTheta*cos2Theta) +
                DDDDeta(1,1,1,1)*(    cos2Theta*cos2Theta);
            return std::make_tuple(boundary_term, curvature_term);

        #elif AMREX_SPACEDIM == 3

            Set::Vector normal = Deta / Deta.lpNorm<2>();

            // GRAHM-SCHMIDT PROCESS to get orthonormal basis
            const Set::Vector e1(1,0,0), e2(0,1,0), e3(0,0,1);
            Set::Vector _t2, _t3;
            if      (fabs(normal(0)) > fabs(normal(1)) && fabs(normal(0)) > fabs(normal(2)))
            {
                _t2 = e2 - normal.dot(e2)*normal; _t2 /= _t2.lpNorm<2>();
                _t3 = e3 - normal.dot(e3)*normal - _t2.dot(e3)*_t2; _t3 /= _t3.lpNorm<2>();
            }
            else if (fabs(normal(1)) > fabs(normal(0)) && fabs(normal(1)) > fabs(normal(2)))
            {
                _t2 = e1 - normal.dot(e1)*normal; _t2 /= _t2.lpNorm<2>();
                _t3 = e3 - normal.dot(e3)*normal - _t2.dot(e3)*_t2; _t3 /= _t3.lpNorm<2>();
            }
            else
            {
                _t2 = e1 - normal.dot(e1)*normal; _t2 /= _t2.lpNorm<2>();
                _t3 = e2 - normal.dot(e2)*normal - _t2.dot(e2)*_t2; _t3 /= _t3.lpNorm<2>();
            }

            // Compute Hessian projected into tangent space (spanned by _t1,_t2)
            Eigen::Matrix2d DDeta2D;
            DDeta2D <<
                _t2.dot(DDeta*_t2) , _t2.dot(DDeta*_t3),
                _t3.dot(DDeta*_t2) , _t3.dot(DDeta*_t3);
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver(2);
            eigensolver.computeDirect(DDeta2D);
            Eigen::Matrix2d eigenvecs = eigensolver.eigenvectors();

            // Compute tangent vectors embedded in R^3
            Set::Vector t2 = _t2*eigenvecs(0,0) + _t3*eigenvecs(0,1),
                t3 = _t2*eigenvecs(1,0) + _t3*eigenvecs(1,1);

            // Compute components of second Hessian in t2,t3 directions
            Set::Scalar DH2 = 0.0, DH3 = 0.0;
            Set::Scalar DH23 = 0.0;
            for (int p = 0; p < 3; p++)
                for (int q = 0; q < 3; q++)
                    for (int r = 0; r < 3; r++)
                        for (int s = 0; s < 3; s++)
                        {
                            DH2 += DDDDeta(p,q,r,s)*t2(p)*t2(q)*t2(r)*t2(s);
                            DH3 += DDDDeta(p,q,r,s)*t3(p)*t3(q)*t3(r)*t3(s);
                            DH23 += DDDDeta(p,q,r,s)*t2(p)*t2(q)*t3(r)*t3(s);
                        }

            Set::Scalar sigma = W(normal);
            Set::Scalar DDK2 = DDW(normal,_t2);
            Set::Scalar DDK3 = DDW(normal,_t3);

            // GB energy anisotropy term
            Set::Scalar gbenergy_df = - sigma*DDeta.trace() - DDK2*DDeta2D(0,0) - DDK3*DDeta2D(1,1);
                                  
            // Second order curvature term
            Set::Scalar reg_df = NAN;
            if (regularization == Regularization::Wilhelm) reg_df = DH2 + DH3 + 2.0*DH23;
            else if (regularization == Regularization::K23) reg_df = DH2 + DH3;

            return std::make_tuple(gbenergy_df, reg_df);
        
        #endif
    }


protected:
    static constexpr amrex::Real pi = 3.14159265359;
    Regularization regularization = Regularization::Wilhelm;
};
}
}
}

#endif
