#ifndef SET_SET_
#define SET_SET_

#include <iomanip>

#include "AMReX.H"
#include "AMReX_REAL.H"
#include "AMReX_SPACE.H"
#include "AMReX_Vector.H"
#include "AMReX_BLassert.H"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/SVD>
#include <AMReX_REAL.H>

#include "Util/Util.H"

/// \brief A collection of data types and symmetry-reduced data structures
namespace Set
{
using Scalar   = amrex::Real;
using Vector   = Eigen::Matrix<amrex::Real,AMREX_SPACEDIM,1>;
using Covector = Eigen::Matrix<amrex::Real,1,AMREX_SPACEDIM>;
using Matrix   = Eigen::Matrix<amrex::Real,AMREX_SPACEDIM,AMREX_SPACEDIM>;
using iMatrix  = Eigen::Matrix<int,AMREX_SPACEDIM,AMREX_SPACEDIM>;

enum Sym {None,Major,Minor,MajorMinor,Full,Isotropic};
template<int dim,int sym> class Matrix4{};
template<int dim, int sym>
std::ostream&
operator<< (std::ostream& os, const Matrix4<dim,sym>& b)
{
	Matrix4<dim,sym> bcopy = b;
   	bcopy.Print(os);
	return os;
}

}
namespace Util
{
Set::Scalar Random();
}

namespace Set
{
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
template<>
class Matrix4<3,Sym::Full>
{
    Scalar data[15] = { NAN,NAN,NAN,NAN,NAN,
                        NAN,NAN,NAN,NAN,NAN,
                        NAN,NAN,NAN,NAN,NAN};

public:
    Matrix4() {};
    AMREX_FORCE_INLINE
    Scalar & operator () (const int i, const int j, const int k, const int l)
    {
        int uid = i + 3*j + 9*k + 27*l;
        // [0, 0, 0, 0]
        if (uid==0 ) return data[0];
        // [0, 0, 0, 1]
        else if (uid==27 || uid==9 || uid==3 || uid==1 ) return data[1];
        // [0, 0, 0, 2]
        else if (uid==54 || uid==18 || uid==6 || uid==2 ) return data[2];
        // [0, 0, 1, 1]
        else if (uid==36 || uid==30 || uid==12 || uid==28 || uid==10 || uid==4 ) return data[3];
        // [0, 0, 1, 2]
        else if (uid==63 || uid==45 || uid==57 || uid==21 || uid==33 || uid==15 || uid==55 || uid==19 || uid==7 || uid==29 || uid==11 || uid==5 ) return data[4];
        // [0, 0, 2, 2]
        else if (uid==72 || uid==60 || uid==24 || uid==56 || uid==20 || uid==8 ) return data[5];
        // [0, 1, 1, 1]
        else if (uid==39 || uid==37 || uid==31 || uid==13 ) return data[6];
        // [0, 1, 1, 2]
        else if (uid==66 || uid==48 || uid==42 || uid==64 || uid==46 || uid==58 || uid==22 || uid==34 || uid==16 || uid==38 || uid==32 || uid==14 ) return data[7];
        // [0, 1, 2, 2]
        else if (uid==75 || uid==69 || uid==51 || uid==73 || uid==61 || uid==25 || uid==65 || uid==47 || uid==59 || uid==23 || uid==35 || uid==17 ) return data[8];
        // [0, 2, 2, 2]
        else if (uid==78 || uid==74 || uid==62 || uid==26 ) return data[9];
        // [1, 1, 1, 1]
        else if (uid==40 ) return data[10];
        // [1, 1, 1, 2]
        else if (uid==67 || uid==49 || uid==43 || uid==41 ) return data[11];
        // [1, 1, 2, 2]
        else if (uid==76 || uid==70 || uid==52 || uid==68 || uid==50 || uid==44 ) return data[12];
        // [1, 2, 2, 2]
        else if (uid==79 || uid==77 || uid==71 || uid==53 ) return data[13];
        // [2, 2, 2, 2]
        else if (uid==80 ) return data[14];
        else Util::Abort(INFO,"uid not in range (uid=",uid,")");        
        return data[-1];
    }
    void Print (std::ostream& os)
    {
        for (int i = 0; i < 14; i++)
            os << "i = " << i << " " << data[i] << std::endl;
        
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
                            const Set::Scalar &val = (*this)(i,j,k,l);
                            os << std::scientific << std::setw(11) << val ; //(fabs(val)>1E-10 ? val : 0);
                            os << " "; 
                        }
                        os << "│ ";
                    }
                }
                os<<std::endl;
            }
        }
    }
    static Matrix4<3,Sym::Full> Increment()
    {
        Matrix4<3,Sym::Full> ret;
        for (int i = 0 ; i < 15; i++) ret.data[i] = (Set::Scalar)i;
        return ret;
    }
    static Matrix4<3,Sym::Full> Randomize()
    {
        Matrix4<3,Sym::Full> ret;
        for (int i = 0 ; i < 15; i++) ret.data[i] = Util::Random();
        return ret;
    }
    static Matrix4<3,Sym::Full> Zero()
    {
        Matrix4<3,Sym::Full> ret;
        for (int i = 0 ; i < 15; i++) ret.data[i] = 0.0;
        return ret;
     }
};
std::ostream&
operator<< (std::ostream& os, const Matrix4<3,Sym::Full>& b);

template<>
class Matrix4<2,Sym::Full>
{
    Scalar data[5] = {NAN,NAN,NAN,NAN,NAN};
public:
    Matrix4() {};
    AMREX_FORCE_INLINE
    Scalar & operator () (const int i, const int j, const int k, const int l)
    {
        int uid = i + 2*j + 4*k + 8*l;
        // [0, 0, 0, 0]
        if (uid==0 ) return data[0];
        // [0, 0, 0, 1]
        else if (uid==8 || uid==4 || uid==2 || uid==1 ) return data[1];
        // [0, 0, 1, 1]
        else if (uid==12 || uid==10 || uid==6 || uid==9 || uid==5 || uid==3 ) return data[2];
        // [0, 1, 1, 1]
        else if (uid==14 || uid==13 || uid==11 || uid==7 ) return data[3];
        // [1, 1, 1, 1]
        else if (uid==15 ) return data[4];
        else Util::Abort(INFO,"Index out of range");
        return data[-1];
    }
    static Matrix4<2,Sym::Full> Randomize()
    {
        Matrix4<2,Sym::Full> ret;
        for (int i = 0 ; i < 5; i++) ret.data[i] = Util::Random();
        return ret;
    }
    void Print (std::ostream& /* os */)
    {
        Util::Abort(INFO,"Not yet implemented");
    }
    static Matrix4<2,Sym::Full> Zero()
    {
        Matrix4<2,Sym::Full> zero;
        zero.data[0] = 0.0;
        zero.data[1] = 0.0;
        zero.data[2] = 0.0;
        zero.data[3] = 0.0;
        zero.data[4] = 0.0;
        return zero;
    }
};

/// \brief Data structure for a 4th order 3D tensor with major and minor symmetry
///
/// Let the tensor
/// \f$\mathbb{C}\in\mathbb{R}^3\times\mathbb{R}^3\times\mathbb{R}^3\times\mathbb{R}^3\f$
/// be have major symmetry (\f$\mathbb{C}_{ijkl}=\mathbb{C}_{klij}\f$)
/// and minor symmetry (\f$\mathbb{C}_{ijkl}=\mathbb{C}_{jikl}=\mathbb{C}_{ijlk}\f$)
/// Then there are only 21 unique elements (rather than 81).
///
/// This object acts like a 4D array such that `C(i,j,k,l)` returns the corresponding
/// element, but symmetry is always obeyed. This allows the user code to be much prettier, 
/// while maintaining a relatively small object size.
///
template<>
class Matrix4<3,Sym::MajorMinor>
{
    Scalar data[21] = { NAN,NAN,NAN,NAN,NAN,NAN,NAN,
                        NAN,NAN,NAN,NAN,NAN,NAN,NAN,
                        NAN,NAN,NAN,NAN,NAN,NAN,NAN };

public:
    Matrix4() {};
    AMREX_FORCE_INLINE
    Scalar & operator () (const int i, const int j, const int k, const int l)
    {
        int uid = i + 3*j + 9*k + 27*l;
        // [0, 0, 0, 0]
        if (0 || uid==0 ) return data[0];
        // [0, 0, 0, 1]
        else if (0 || uid==27 || uid==9 || uid==3 || uid==1 ) return data[1];
        // [0, 0, 0, 2]
        else if (0 || uid==54 || uid==18 || uid==6 || uid==2 ) return data[2];
        // [0, 0, 1, 1]
        else if (0 || uid==36 || uid==4 ) return data[3];
        // [0, 0, 1, 2]
        else if (0 || uid==63 || uid==45 || uid==7 || uid==5 ) return data[4];
        // [0, 0, 2, 2]
        else if (0 || uid==72 || uid==8 ) return data[5];
        // [0, 1, 0, 1]
        else if (0 || uid==30 || uid==12 || uid==28 || uid==10 ) return data[6];
        // [0, 1, 0, 2]
        else if (0 || uid==57 || uid==21 || uid==33 || uid==15 || uid==55 || uid==19 || uid==29 || uid==11 ) return data[7];
        // [0, 1, 1, 1]
        else if (0 || uid==39 || uid==37 || uid==31 || uid==13 ) return data[8];
        // [0, 1, 1, 2]
        else if (0 || uid==66 || uid==48 || uid==64 || uid==46 || uid==34 || uid==16 || uid==32 || uid==14 ) return data[9];
        // [0, 1, 2, 2]
        else if (0 || uid==75 || uid==73 || uid==35 || uid==17 ) return data[10];
        // [0, 2, 0, 2]
        else if (0 || uid==60 || uid==24 || uid==56 || uid==20 ) return data[11];
        // [0, 2, 1, 1]
        else if (0 || uid==42 || uid==58 || uid==22 || uid==38 ) return data[12];
        // [0, 2, 1, 2]
        else if (0 || uid==69 || uid==51 || uid==61 || uid==25 || uid==65 || uid==47 || uid==59 || uid==23 ) return data[13];
        // [0, 2, 2, 2]
        else if (0 || uid==78 || uid==74 || uid==62 || uid==26 ) return data[14];
        // [1, 1, 1, 1]
        else if (0 || uid==40 ) return data[15];
        // [1, 1, 1, 2]
        else if (0 || uid==67 || uid==49 || uid==43 || uid==41 ) return data[16];
        // [1, 1, 2, 2]
        else if (0 || uid==76 || uid==44 ) return data[17];
        // [1, 2, 1, 2]
        else if (0 || uid==70 || uid==52 || uid==68 || uid==50 ) return data[18];
        // [1, 2, 2, 2]
        else if (0 || uid==79 || uid==77 || uid==71 || uid==53 ) return data[19];
        // [2, 2, 2, 2]
        else if (0 || uid==80 ) return data[20];
        else Util::Abort(INFO,"Index out of range");
        return data[-1]; 
    }
    void Print (std::ostream& os)
    {
        for (int i = 0; i < 14; i++)
            os << "i = " << i << " " << data[i] << std::endl;
        
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
                            const Set::Scalar &val = (*this)(i,j,k,l);
                            os << std::scientific << std::setw(11) << val ; //(fabs(val)>1E-10 ? val : 0);
                            os << " "; 
                        }
                        os << "│ ";
                    }
                }
                os<<std::endl;
            }
        }
    }
    void operator += (Matrix4<3,Sym::MajorMinor> a) {for (int i = 0; i < 21; i++) data[i] += a.data[i];}
    void operator -= (Matrix4<3,Sym::MajorMinor> a) {for (int i = 0; i < 21; i++) data[i] -= a.data[i];}
    void operator *= (Matrix4<3,Sym::MajorMinor> a) {for (int i = 0; i < 21; i++) data[i] *= a.data[i];}
    void operator /= (Matrix4<3,Sym::MajorMinor> a) {for (int i = 0; i < 21; i++) data[i] /= a.data[i];}
    void operator *= (Set::Scalar alpha) {for (int i = 0; i < 21; i++) data[i] *= alpha;}
    void operator /= (Set::Scalar alpha) {for (int i = 0; i < 21; i++) data[i] /= alpha;}

    static Matrix4<3,Sym::MajorMinor> Increment()
    {
        Matrix4<3,Sym::MajorMinor> ret;
        for (int i = 0 ; i < 21; i++) ret.data[i] = (Set::Scalar)i;
        return ret;
     }
    static Matrix4<3,Sym::MajorMinor> Randomize()
    {
        Matrix4<3,Sym::MajorMinor> ret;
        for (int i = 0 ; i < 21; i++) ret.data[i] = (Util::Random());
        return ret;
     }
    static Matrix4<3,Sym::MajorMinor> Zero()
    {
        Matrix4<3,Sym::MajorMinor> ret;
        for (int i = 0 ; i < 21; i++) ret.data[i] = 0.0;
        return ret;
     }
};
Eigen::Matrix<amrex::Real,3,3> operator * (Matrix4<3,Sym::MajorMinor> a, Eigen::Matrix<amrex::Real,3,3> b);
Eigen::Matrix<amrex::Real,2,2> operator * (Matrix4<3,Sym::MajorMinor> a, Eigen::Matrix<amrex::Real,2,2> b);
Eigen::Matrix<amrex::Real,3,1> operator * (Matrix4<3,Sym::MajorMinor> a, std::array<Eigen::Matrix<amrex::Real,3,3>,3> b);
Eigen::Matrix<amrex::Real,2,1> operator * (Matrix4<3,Sym::MajorMinor> a, std::array<Eigen::Matrix<amrex::Real,2,2>,2> b);
 
namespace Constant
{
static const Set::Scalar Pi = 3.14159265359;
}
}

#endif