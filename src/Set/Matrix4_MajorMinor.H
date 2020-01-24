#ifndef SET_MATRIX4_MAJORMINOR_H
#define SET_MATRIX4_MAJORMINOR_H

#include "Base.H"

namespace Set
{
/// \brief Data structure for a 4th order 3D tensor with major and minor symmetry
///
/// 2D version. See full explanation below.
///
template<>
class Matrix4<2,Sym::MajorMinor>
{
    Scalar data[6] = {NAN,NAN,NAN,NAN,NAN,NAN};

public:
    AMREX_GPU_HOST_DEVICE Matrix4() {};
    AMREX_FORCE_INLINE
    Scalar & operator () (const int i, const int j, const int k, const int l)
    {
        int uid = i + 2*j + 4*k + 8*l;
        if      (uid==0 )                                  return data[0]; // [0000]
        else if (uid==8  || uid==4 || uid==2 || uid==1 )   return data[1]; // [0001] [0010] [0100] [1000]
        else if (uid==12 || uid==3 )                       return data[2]; // [0011] [1100]
        else if (uid==10 || uid==6 || uid==9 || uid==5 )   return data[3]; // [0101] [1001] [0110] [1010]
        else if (uid==14 || uid==13 || uid==11 || uid==7 ) return data[4]; // [0111] [1011] [1101] [1110]
        else if (uid==15 )                                 return data[5]; // [1111]
        else Util::Abort(INFO,"Index out of range");
        return data[-1]; 
    }
    void Print (std::ostream& os)
    {
        os.precision(4);
        for (int k = 0; k < 2; k++)
        {
            for (int i = -1; i < 3; i++)
            {
                for (int l = 0; l < 2; l++)
                {
                    if (i==-1)       os << "┌                         ┐ ";
                    else if (i == 2) os << "└                         ┘ ";
                    else 
                    {
                        os << "│ ";
                        for (int j = 0; j < 2; j++)
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
    AMREX_GPU_HOST_DEVICE void operator  = (Matrix4<2,Sym::MajorMinor> a) {for (int i = 0; i < 6; i++) data[i] =  a.data[i];}
    AMREX_GPU_HOST_DEVICE void operator += (Matrix4<2,Sym::MajorMinor> a) {for (int i = 0; i < 6; i++) data[i] += a.data[i];}
    AMREX_GPU_HOST_DEVICE void operator -= (Matrix4<2,Sym::MajorMinor> a) {for (int i = 0; i < 6; i++) data[i] -= a.data[i];}
    AMREX_GPU_HOST_DEVICE void operator *= (Matrix4<2,Sym::MajorMinor> a) {for (int i = 0; i < 6; i++) data[i] *= a.data[i];}
    AMREX_GPU_HOST_DEVICE void operator /= (Matrix4<2,Sym::MajorMinor> a) {for (int i = 0; i < 6; i++) data[i] /= a.data[i];}
    AMREX_GPU_HOST_DEVICE void operator *= (Set::Scalar alpha)            {for (int i = 0; i < 6; i++) data[i] *= alpha;}
    AMREX_GPU_HOST_DEVICE void operator /= (Set::Scalar alpha)            {for (int i = 0; i < 6; i++) data[i] /= alpha;}

    static Matrix4<2,Sym::MajorMinor> Increment()
    {
        Matrix4<2,Sym::MajorMinor> ret;
        for (int i = 0 ; i < 6; i++) ret.data[i] = (Set::Scalar)i;
        return ret;
     }
    static Matrix4<2,Sym::MajorMinor> Randomize()
    {
        Matrix4<2,Sym::MajorMinor> ret;
        for (int i = 0 ; i < 6; i++) ret.data[i] = (Util::Random());
        return ret;
     }
    static Matrix4<2,Sym::MajorMinor> Zero()
    {
        Matrix4<2,Sym::MajorMinor> ret;
        for (int i = 0 ; i < 6; i++) ret.data[i] = 0.0;
        return ret;
    }
    friend Eigen::Matrix<Set::Scalar,2,2> operator * (Matrix4<2,Sym::MajorMinor> a, Eigen::Matrix<Set::Scalar,2,2> b);
    friend Set::Vector operator * (Matrix4<2,Sym::MajorMinor> a, Set::Matrix3 b);
};
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Eigen::Matrix<Set::Scalar,2,2> operator * (Matrix4<2,Sym::MajorMinor> a, Eigen::Matrix<Set::Scalar,2,2> b)
{
   	Eigen::Matrix<Set::Scalar,2,2> ret = Eigen::Matrix<Set::Scalar,2,2>::Zero();
    ret(0,0) = a.data[0]*b(0,0) + a.data[1]*(b(0,1) + b(1,0)) + a.data[2]*b(1,1);
    ret(0,1) = a.data[1]*b(0,0) + a.data[3]*(b(0,1) + b(1,0)) + a.data[4]*b(1,1);
    ret(1,1) = a.data[2]*b(0,0) + a.data[4]*(b(0,1) + b(1,0)) + a.data[5]*b(1,1);
    ret(1,0) = ret(0,1);
	return ret;
}
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Eigen::Matrix<amrex::Real,2,1> operator * (Matrix4<2,Sym::MajorMinor> a, std::array<Eigen::Matrix<amrex::Real,2,2>,2> b);

AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Set::Vector operator * (Matrix4<2,Sym::MajorMinor> a, Set::Matrix3 b)
{
    Set::Vector ret = Set::Vector::Zero();
    ret(0) = a.data[0]*b(0,0,0) + a.data[1]*(b(0,1,0) + b(1,0,0) + b(0,0,1)) + a.data[2]*b(1,1,0)  + a.data[3]*(b(0,1,1) + b(1,0,1)) + a.data[4]*b(1,1,1);
    ret(1) = a.data[1]*b(0,0,0) + a.data[2]*b(0,0,1) + a.data[3]*(b(0,1,0) + b(1,0,0)) + a.data[4]*(b(1,1,0)  + b(0,1,1) + b(1,0,1)) + a.data[5]*b(1,1,1);
    return ret;
}
    
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
    AMREX_GPU_HOST_DEVICE Matrix4() {};
    AMREX_FORCE_INLINE
    Scalar & operator () (const int i, const int j, const int k, const int l)
    {
        int uid = i + 3*j + 9*k + 27*l;
        // [0000]
        if (uid==0 ) return data[0];
        // [0001] [0010] [0100] [1000]
        else if (uid==27 || uid==9 || uid==3 || uid==1 ) return data[1];
        // [0002] [0020] [0200] [2000]
        else if (uid==54 || uid==18 || uid==6 || uid==2 ) return data[2];
        // [0011] [1100]
        else if (uid==36 || uid==4 ) return data[3];
        // [0012] [0021] [1200] [2100]
        else if (uid==63 || uid==45 || uid==7 || uid==5 ) return data[4];
        // [0022] [2200]
        else if (uid==72 || uid==8 ) return data[5];
        // [0101] [0110] [1001] [1010]
        else if (uid==30 || uid==12 || uid==28 || uid==10 ) return data[6];
        // [0102] [0120] [1002] [1020] [0201] [2001] [0210] [2010]
        else if (uid==57 || uid==21 || uid==33 || uid==15 || uid==55 || uid==19 || uid==29 || uid==11 ) return data[7];
        // [0111] [1011] [1101] [1110]
        else if (uid==39 || uid==37 || uid==31 || uid==13 ) return data[8];
        // [0112] [1012] [0121] [1021] [1201] [1210] [2101] [2110]
        else if (uid==66 || uid==48 || uid==64 || uid==46 || uid==34 || uid==16 || uid==32 || uid==14 ) return data[9];
        // [0122] [1022] [2201] [2210]
        else if (uid==75 || uid==73 || uid==35 || uid==17 ) return data[10];
        // [0202] [2002] [0220] [2020]
        else if (uid==60 || uid==24 || uid==56 || uid==20 ) return data[11];
        // [0211] [2011] [1102] [1120]
        else if (uid==42 || uid==58 || uid==22 || uid==38 ) return data[12];
        // [0212] [2012] [0221] [2021] [1202] [1220] [2102] [2120]
        else if (uid==69 || uid==51 || uid==61 || uid==25 || uid==65 || uid==47 || uid==59 || uid==23 ) return data[13];
        // [0222] [2022] [2202] [2220]
        else if (uid==78 || uid==74 || uid==62 || uid==26 ) return data[14];
        // [1111]
        else if (uid==40 ) return data[15];
        // [1112] [1121] [1211] [2111]
        else if (uid==67 || uid==49 || uid==43 || uid==41 ) return data[16];
        // [1122] [2211]
        else if (uid==76 || uid==44 ) return data[17];
        // [1212] [2112] [1221] [2121]
        else if (uid==70 || uid==52 || uid==68 || uid==50 ) return data[18];
        // [1222] [2122] [2212] [2221]
        else if (uid==79 || uid==77 || uid==71 || uid==53 ) return data[19];
        // [2222]
        else if (uid==80 ) return data[20];
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
    AMREX_GPU_HOST_DEVICE void operator  = (Matrix4<3,Sym::MajorMinor> a) {for (int i = 0; i < 21; i++) data[i] = a.data[i];}
    AMREX_GPU_HOST_DEVICE void operator += (Matrix4<3,Sym::MajorMinor> a) {for (int i = 0; i < 21; i++) data[i] += a.data[i];}
    AMREX_GPU_HOST_DEVICE void operator -= (Matrix4<3,Sym::MajorMinor> a) {for (int i = 0; i < 21; i++) data[i] -= a.data[i];}
    AMREX_GPU_HOST_DEVICE void operator *= (Matrix4<3,Sym::MajorMinor> a) {for (int i = 0; i < 21; i++) data[i] *= a.data[i];}
    AMREX_GPU_HOST_DEVICE void operator /= (Matrix4<3,Sym::MajorMinor> a) {for (int i = 0; i < 21; i++) data[i] /= a.data[i];}
    AMREX_GPU_HOST_DEVICE void operator *= (Set::Scalar alpha) {for (int i = 0; i < 21; i++) data[i] *= alpha;}
    AMREX_GPU_HOST_DEVICE void operator /= (Set::Scalar alpha) {for (int i = 0; i < 21; i++) data[i] /= alpha;}

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
    friend Set::Vector operator * (Matrix4<3,Sym::MajorMinor> a, Set::Matrix3 b);
    friend Eigen::Matrix<amrex::Real,3,3> operator * (Matrix4<3,Sym::MajorMinor> a, Eigen::Matrix<amrex::Real,3,3> b);
};
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Eigen::Matrix<amrex::Real,3,3> operator * (Matrix4<3,Sym::MajorMinor> a, Eigen::Matrix<amrex::Real,3,3> b)
{
    Eigen::Matrix<amrex::Real,3,3> ret = Eigen::Matrix<amrex::Real,3,3>::Zero();
	ret(0,0) += a.data[   0]*b(0,0) + a.data[   1]*(b(0,1) + b(1,0)) + a.data[   2]*(b(0,2) + b(2,0)) + a.data[   3]*b(1,1) + a.data[   4]*(b(1,2) + b(2,1)) + a.data[   5]*b(2,2);
   	ret(1,1) += a.data[   3]*b(0,0) + a.data[   8]*(b(0,1) + b(1,0)) + a.data[  12]*(b(0,2) + b(2,0)) + a.data[  15]*b(1,1) + a.data[  16]*(b(1,2) + b(2,1)) + a.data[  17]*b(2,2);
  	ret(2,2) += a.data[   5]*b(0,0) + a.data[  10]*(b(0,1) + b(1,0)) + a.data[  14]*(b(0,2) + b(2,0)) + a.data[  17]*b(1,1) + a.data[  19]*(b(1,2) + b(2,1)) + a.data[  20]*b(2,2);
  	ret(1,2) += a.data[   4]*b(0,0) + a.data[   9]*(b(0,1) + b(1,0)) + a.data[  13]*(b(0,2) + b(2,0)) + a.data[  16]*b(1,1) + a.data[  18]*(b(1,2) + b(2,1)) + a.data[  19]*b(2,2);
    ret(2,0) += a.data[   2]*b(0,0) + a.data[   7]*(b(0,1) + b(1,0)) + a.data[  11]*(b(0,2) + b(2,0)) + a.data[  12]*b(1,1) + a.data[  13]*(b(1,2) + b(2,1)) + a.data[  14]*b(2,2);
  	ret(0,1) += a.data[   1]*b(0,0) + a.data[   6]*(b(0,1) + b(1,0)) + a.data[   7]*(b(0,2) + b(2,0)) + a.data[   8]*b(1,1) + a.data[   9]*(b(1,2) + b(2,1)) + a.data[  10]*b(2,2);

  	ret(2,1) += ret(1,2);
  	ret(0,2) += ret(2,0);
	ret(1,0) += ret(0,1);

	return ret;
}
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Eigen::Matrix<amrex::Real,2,2> operator * (Matrix4<3,Sym::MajorMinor> a, Eigen::Matrix<amrex::Real,2,2> b)
{
   	Eigen::Matrix<amrex::Real,2,2> ret = Eigen::Matrix<amrex::Real,2,2>::Zero();
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			for (int k = 0; k < 2; k++)
				for (int l = 0; l < 2; l++)
					ret(i,j) += a(i,j,k,l)*b(k,l);
	return ret;
}
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Eigen::Matrix<amrex::Real,3,1> operator * (Matrix4<3,Sym::MajorMinor> a, std::array<Eigen::Matrix<amrex::Real,3,3>,3> b);

AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Eigen::Matrix<amrex::Real,2,1> operator * (Matrix4<3,Sym::MajorMinor> a, std::array<Eigen::Matrix<amrex::Real,2,2>,2> b);

AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Set::Vector operator * (Matrix4<3,Sym::MajorMinor> a, Set::Matrix3 b)
{
    // TODO: improve efficiency of this method
    Set::Vector ret = Set::Vector::Zero();
    
    ret(0) += 
    a.data[   0]*b(0,0,0) + 
    a.data[   1]*(b(0,1,0) + b(1,0,0) + b(0,0,1)) + 
    a.data[   2]*(b(0,2,0) + b(2,0,0) + b(0,0,2)) + 
    a.data[   3]*b(1,1,0) + 
    a.data[   4]*(b(1,2,0) + b(2,1,0)) + 
    a.data[   5]*b(2,2,0) + 
    a.data[   6]*(b(0,1,1) + b(1,0,1)) + 
    a.data[   7]*(b(0,2,1) + b(2,0,1) + b(0,1,2) + b(1,0,2)) + 
    a.data[   8]*b(1,1,1) + 
    a.data[   9]*(b(1,2,1) + b(2,1,1)) + 
    a.data[  10]*b(2,2,1) + 
    a.data[  11]*(b(0,2,2) + b(2,0,2)) + 
    a.data[  12]*b(1,1,2) + 
    a.data[  13]*(b(1,2,2) + b(2,1,2)) + 
    a.data[  14]*b(2,2,2);
    
    ret(1) += 
    a.data[   1]*b(0,0,0) + 
    a.data[   4]*b(0,0,2) + 
    a.data[   3]*b(0,0,1) + 
    a.data[   6]*(b(0,1,0) + b(1,0,0)) + 
    a.data[   7]*(b(0,2,0) + b(2,0,0)) + 
    a.data[   8]*(b(1,1,0) + b(0,1,1) + b(1,0,1)) + 
    a.data[   9]*(b(1,2,0) + b(2,1,0) + b(0,1,2) + b(1,0,2)) + 
    a.data[  10]*b(2,2,0) + 
    a.data[  12]*(b(0,2,1) + b(2,0,1)) + 
    a.data[  13]*(b(0,2,2) + b(2,0,2)) + 
    a.data[  15]*b(1,1,1) + 
    a.data[  16]*(b(1,2,1) + b(2,1,1) + b(1,1,2)) + 
    a.data[  17]*b(2,2,1) + 
    a.data[  18]*(b(1,2,2) + b(2,1,2)) + 
    a.data[  19]*b(2,2,2);

    ret(2) += 
    a.data[   2]*b(0,0,0) + 
    a.data[   4]*b(0,0,1) + 
    a.data[   5]*b(0,0,2) + 
    a.data[   7]*(b(0,1,0) + b(1,0,0)) + 
    a.data[   9]*(b(0,1,1) + b(1,0,1)) + 
    a.data[  10]*(b(0,1,2) + b(1,0,2)) + 
    a.data[  11]*(b(0,2,0) + b(2,0,0)) + 
    a.data[  12]*b(1,1,0) + 
    a.data[  13]*(b(1,2,0) + b(2,1,0) + b(0,2,1) + b(2,0,1)) + 
    a.data[  14]*(b(2,2,0) + b(0,2,2) + b(2,0,2)) + 
    a.data[  16]*b(1,1,1) + 
    a.data[  17]*b(1,1,2) + 
    a.data[  18]*(b(1,2,1) + b(2,1,1)) + 
    a.data[  19]*(b(2,2,1) + b(1,2,2) + b(2,1,2)) + 
    a.data[  20]*b(2,2,2);

    return ret;
}

}
#endif