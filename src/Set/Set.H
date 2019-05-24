#ifndef SET_SET_
#define SET_SET_

#include <iomanip>

#include "AMReX.H"
#include "AMReX_REAL.H"
#include "AMReX_SPACE.H"
#include "AMReX_Vector.H"
#include "AMReX_BLassert.H"
//#define EIGEN_MATRIX_PLUGIN "Set/Eigen_Amrex.H"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/SVD>
#include <AMReX_REAL.H>

/// \brief A collection of data types and simple data structures
namespace Set
{
using Scalar   = amrex::Real;
using Vector   = Eigen::Matrix<amrex::Real,AMREX_SPACEDIM,1>;
using Covector = Eigen::Matrix<amrex::Real,1,AMREX_SPACEDIM>;
using Matrix   = Eigen::Matrix<amrex::Real,AMREX_SPACEDIM,AMREX_SPACEDIM>;

enum Sym {None,Major,Minor,MajorMinor,Full};
template<int dim,int sym>
class Matrix4{};

/// \brief Data structure for a symmetrix 4th order 3D tensor
///
/// Let the tensor
/// \f$\mathbb{C}\in\mathbb{R}^3\times\mathbb{R}^3\times\mathbb{R}^3\times\mathbb{R}^3\f$
/// be fully symmetrix such that \f$\mathbb{C}_{ijkl}=\mathbb{C}_{\sigma(i,j,k,l)\f$ where
/// \f$\sigma\f$ is any permutation.
/// Then there are only 15 unique elements (rather than 81).
///
/// This object acts like a 4D array such that `C(i,j,k,l)` returns the corresponding
/// element, but symmetry is always obeyed. This allows the user code to be much prettier.
///
/// The way this is implemented is quite crude, because I don't know of a nice way to 
/// compute unique indices. Improvements are needed...
template<>
class Matrix4<3,Sym::Full>
{
    Scalar data[15];
public:
    Matrix4() {};
    AMREX_FORCE_INLINE
    Scalar & operator () (const int i, const int j, const int k, const int l)
    {
        int uid = i + 3*j + 9*k + 27*l;
        // [0, 0, 0, 0]
        if (uid==0 ) return data[0];
        // [0, 0, 0, 1]
        if (uid==27 || uid==9 || uid==3 || uid==1 ) return data[1];
        // [0, 0, 0, 2]
        if (uid==54 || uid==18 || uid==6 || uid==2 ) return data[2];
        // [0, 0, 1, 1]
        if (uid==36 || uid==30 || uid==12 || uid==28 || uid==10 || uid==4 ) return data[3];
        // [0, 0, 1, 2]
        if (uid==63 || uid==45 || uid==57 || uid==21 || uid==33 || uid==15 || uid==55 || uid==19 || uid==7 || uid==29 || uid==11 || uid==5 ) return data[4];
        // [0, 0, 2, 2]
        if (uid==72 || uid==60 || uid==24 || uid==56 || uid==20 || uid==8 ) return data[5];
        // [0, 1, 1, 1]
        if (uid==39 || uid==37 || uid==31 || uid==13 ) return data[6];
        // [0, 1, 1, 2]
        if (uid==66 || uid==48 || uid==42 || uid==64 || uid==46 || uid==58 || uid==22 || uid==34 || uid==16 || uid==38 || uid==32 || uid==14 ) return data[7];
        // [0, 1, 2, 2]
        if (uid==75 || uid==69 || uid==51 || uid==73 || uid==61 || uid==25 || uid==65 || uid==47 || uid==59 || uid==23 || uid==35 || uid==17 ) return data[8];
        // [0, 2, 2, 2]
        if (uid==78 || uid==74 || uid==62 || uid==26 ) return data[9];
        // [1, 1, 1, 1]
        if (uid==40 ) return data[10];
        // [1, 1, 1, 2]
        if (uid==67 || uid==49 || uid==43 || uid==41 ) return data[11];
        // [1, 1, 2, 2]
        if (uid==76 || uid==70 || uid==52 || uid==68 || uid==50 || uid==44 ) return data[12];
        // [1, 2, 2, 2]
        if (uid==79 || uid==77 || uid==71 || uid==53 ) return data[13];
        // [2, 2, 2, 2]
        if (uid==80 ) return data[14];
    }
    void Print (std::ostream& os)
    {
        os.precision(4);
        for (int k = 0; k < 3; k++)
        {
            for (int i = -1; i < 4; i++)
            {
                for (int l = 0; l < 3; l++)
                {
                    if (i==-1)       os << "┌                                     ┐ ";
                    else if (i == 3) os << "└                                     ┘ ";
                    else 
                    {
                        os << "│ ";
                        for (int j = 0; j < 3; j++)
                        {
                            os << std::scientific << std::setw(11) << (fabs((*this)(i,j,k,l))>1E-10 ? (*this)(i,j,k,l) : 0);
                            os << " "; 
                        }
                        os << "│ ";
                    }
                }
                os<<std::endl;
            }
        }
    }
};

namespace Constant
{
static const Set::Scalar Pi = 3.14159265359;
}
}


#endif
