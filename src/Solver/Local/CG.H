#include <eigen3/Eigen/Core>
#include "Set/Set.H"

/// A bunch of solvers
namespace Solver
{
/// Local solvers
namespace Local
{
/// This is a basic implementation of the conjugate gradient method to solve
/// systems of the form
/// \f[ \mathbb{A}\mathbf{x} = \mathbf{b}\f]
/// where \f$\mathbf{x},\mathbf{b}\in \text{Sym}(GL(3))\f$, i.e. the set of symmetric
/// 3x3 matrices, and \f$\mathbb{A}\f$ is a 4th order tensor with major and minor
/// symmetries.
Set::Matrix CG (Set::Matrix4<3,Set::Sym::MajorMinor> A, 
                Set::Matrix b, 
                Set::Matrix x = Set::Matrix::Zero(),
                Set::iMatrix mask = Set::iMatrix::Zero(),
                int verbose = false)
{
    Set::Matrix r = b - A*x;
    for (int i = 0; i < AMREX_SPACEDIM; i++) 
        for (int j = 0; j < AMREX_SPACEDIM; j++) 
            if (mask(i,j)) r(i,j) = 0.0;
    if (std::sqrt(r.lpNorm<2>()) < 1E-8) 
    {
        if (verbose) Util::Message(INFO,"No iterations needed");
        return x;
    }

    Set::Matrix p = r;
    for (int k = 0; k < 12; k++)
    {
        Set::Scalar alpha = (r.transpose()*r).trace() / (p.transpose() * (A*p)).trace();
        x = x + (p * alpha);

        Set::Matrix rnew = r - (A*p)*alpha;
        for (int i = 0; i < AMREX_SPACEDIM; i++) 
            for (int j = 0; j < AMREX_SPACEDIM; j++) 
                if (mask(i,j)) rnew(i,j) = 0.0;

        if (verbose) Util::Message(INFO,"Iteration ", k, ": Resid=" , std::sqrt(rnew.lpNorm<2>()));
        if (std::sqrt(rnew.lpNorm<2>()) < 1E-8) return x;

        Set::Scalar beta = (rnew.transpose()*rnew).trace() / (r.transpose()*r).trace();

        p = rnew + p*beta;
        r = rnew;
    }
    Util::Abort(INFO,"Solver failed to converge");
}
}
}
