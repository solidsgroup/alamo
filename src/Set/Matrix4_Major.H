#ifndef SET_MATRIX4_MAJOR_H
#define SET_MATRIX4_MAJOR_H

#include "Base.H"

namespace Set
{
    
template <>
class Matrix4<2, Sym::Major>
{
    Scalar data[10] = {NAN, NAN, NAN, NAN, NAN,
                        NAN, NAN, NAN, NAN, NAN};

public:
    AMREX_GPU_HOST_DEVICE Matrix4(){};

    #if AMREX_SPACEDIM == 2
    Matrix4(const Matrix4<2,Sym::Isotropic> &in)
    {
        // TODO: very inefficient, need to optimize
        for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
        for (int k = 0; k < 2; k++)
        for (int l = 0; l < 2; l++)
            (*this)(i,j,k,l) = in(i,j,k,l);
    }
    Matrix4(const Matrix4<2,Sym::Diagonal> &in)
    {
        // TODO: very inefficient, need to optimize
        for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
        for (int k = 0; k < 2; k++)
        for (int l = 0; l < 2; l++)
            (*this)(i,j,k,l) = in(i,j,k,l);
    }
    Matrix4(const Matrix4<2,Sym::MajorMinor> &in)
    {
        // TODO: very inefficient, need to optimize
        for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
        for (int k = 0; k < 2; k++)
        for (int l = 0; l < 2; l++)
            (*this)(i,j,k,l) = in(i,j,k,l);
    }
    #endif
    
    AMREX_FORCE_INLINE
    const Scalar &operator()(const int i, const int j, const int k, const int l) const
    {
        int uid = i + 2 * j + 4 * k + 8 * l;
        if      (uid==0 )             return data[0]; // [0000]
        else if (uid==8  || uid==2 )  return data[1]; // [0001] [0100]
        else if (uid==4  || uid==1 )  return data[2]; // [0010] [1000]
        else if (uid==12 || uid==3 )  return data[3]; // [0011] [1100]
        else if (uid==10 )            return data[4]; // [0101] 
        else if (uid==6  || uid==9 )  return data[5]; // [0110] [1001]
        else if (uid==14 || uid==11 ) return data[6]; // [0111] [1101]
        else if (uid==5 )             return data[7]; // [1010] 
        else if (uid==13 || uid==7 )  return data[8]; // [1011] [1110]
        else if (uid==15 )            return data[9]; // [1111] 
        else
            Util::Abort(INFO, "Index out of range");
        return data[-1];
    }

    AMREX_FORCE_INLINE
    Scalar &operator()(const int i, const int j, const int k, const int l)
    {
        int uid = i + 2 * j + 4 * k + 8 * l;
        if      (uid==0 )             return data[0]; // [0000]
        else if (uid==8  || uid==2 )  return data[1]; // [0001] [0100]
        else if (uid==4  || uid==1 )  return data[2]; // [0010] [1000]
        else if (uid==12 || uid==3 )  return data[3]; // [0011] [1100]
        else if (uid==10 )            return data[4]; // [0101] 
        else if (uid==6  || uid==9 )  return data[5]; // [0110] [1001]
        else if (uid==14 || uid==11 ) return data[6]; // [0111] [1101]
        else if (uid==5 )             return data[7]; // [1010] 
        else if (uid==13 || uid==7 )  return data[8]; // [1011] [1110]
        else if (uid==15 )            return data[9]; // [1111] 
        else
            Util::Abort(INFO, "Index out of range");
        return data[-1];
    }
    void Print(std::ostream &os)
    {
        os.precision(4);
        for (int k = 0; k < 2; k++)
        {
            for (int i = -1; i < 3; i++)
            {
                for (int l = 0; l < 2; l++)
                {
                    if (i == -1)
                        os << "┌                         ┐ ";
                    else if (i == 2)
                        os << "└                         ┘ ";
                    else
                    {
                        os << "│ ";
                        for (int j = 0; j < 2; j++)
                        {
                            const Set::Scalar &val = (*this)(i, j, k, l);
                            os << std::scientific << std::setw(11) << val; //(fabs(val)>1E-10 ? val : 0);
                            os << " ";
                        }
                        os << "│ ";
                    }
                }
                os << std::endl;
            }
        }
    }
    //AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
    //void operator=(const Matrix4<2, Sym::Major> &a)  { for (int i = 0; i < 10; i++) data[i]  = a.data[i]; }
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
    void operator+=(const Matrix4<2, Sym::Major> &a) { for (int i = 0; i < 10; i++) data[i] += a.data[i]; }
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
    void operator-=(const Matrix4<2, Sym::Major> &a) { for (int i = 0; i < 10; i++) data[i] -= a.data[i]; }
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
    void operator*=(const Matrix4<2, Sym::Major> &a) { for (int i = 0; i < 10; i++) data[i] *= a.data[i]; }
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
    void operator/=(const Matrix4<2, Sym::Major> &a) { for (int i = 0; i < 10; i++) data[i] /= a.data[i]; }
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
    void operator*=(const Set::Scalar &alpha)        { for (int i = 0; i < 10; i++) data[i] *= alpha; }
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
    void operator/=(const Set::Scalar &alpha)        { for (int i = 0; i < 10; i++) data[i] /= alpha; }

    static Matrix4<2, Sym::Major> Increment()
    {
        Matrix4<2, Sym::Major> ret;
        for (int i = 0; i < 10; i++) ret.data[i] = (Set::Scalar)i;
        return ret;
    }
    static Matrix4<2, Sym::Major> Randomize()
    {
        Matrix4<2, Sym::Major> ret;
        for (int i = 0; i < 10; i++) ret.data[i] = (Util::Random());
        return ret;
    }
    static Matrix4<2, Sym::Major> Zero()
    {
        Matrix4<2, Sym::Major> ret;
        for (int i = 0; i < 10; i++) ret.data[i] = 0.0;
        return ret;
    }
    Set::Scalar Norm()
    {
        Set::Scalar normsq = 0.0;
        for (int i = 0; i < 10; i++) normsq += data[i]*data[i];
        return std::sqrt(normsq);
    }

    friend Matrix4<2, Sym::Major> operator+(const Matrix4<2, Sym::Major> &a, const Matrix4<2, Sym::Major> &b);
    friend Matrix4<2, Sym::Major> operator-(const Matrix4<2, Sym::Major> &a, const Matrix4<2, Sym::Major> &b);
    friend Matrix4<2, Sym::Major> operator*(const Matrix4<2, Sym::Major> &a, const Set::Scalar &b);
    friend Matrix4<2, Sym::Major> operator/(const Matrix4<2, Sym::Major> &a, const Set::Scalar &b);
    friend Set::Matrix operator*(const Matrix4<2, Sym::Major> &a, const Set::Matrix &b);
};
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Matrix4<2, Sym::Major> operator+(const Matrix4<2, Sym::Major> &a, const Matrix4<2, Sym::Major> &b)
{
    Matrix4<2,Sym::Major> ret;
    for (int i = 0; i < 10; i++) ret.data[i] = a.data[i] + b.data[i];
    return ret;
}
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Matrix4<2, Sym::Major> operator-(const Matrix4<2, Sym::Major> &a, const Matrix4<2, Sym::Major> &b)
{
    Matrix4<2,Sym::Major> ret;
    for (int i = 0; i < 10; i++) ret.data[i] = a.data[i] - b.data[i];
    return ret;
}
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Matrix4<2, Sym::Major> operator*(const Matrix4<2, Sym::Major> &a, const Set::Scalar &b)
{
    Matrix4<2,Sym::Major> ret;
    for (int i = 0; i < 10; i++) ret.data[i] = a.data[i] * b;
    return ret;
}
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Matrix4<2, Sym::Major> operator/(const Matrix4<2, Sym::Major> &a, const Set::Scalar &b)
{
    Matrix4<2,Sym::Major> ret;
    for (int i = 0; i < 10; i++) ret.data[i] = a.data[i] / b;
    return ret;
}
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Set::Matrix operator*(const Matrix4<2, Sym::Major> &a, const Set::Matrix &b)
{
    Set::Matrix ret;
    ret(0,0) = a.data[0]*b(0,0) + a.data[1]*b(0,1) + a.data[2]*b(1,0) + a.data[3]*b(1,1);
    ret(0,1) = a.data[1]*b(0,0) + a.data[4]*b(0,1) + a.data[5]*b(1,0) + a.data[6]*b(1,1);
    ret(1,0) = a.data[2]*b(0,0) + a.data[5]*b(0,1) + a.data[7]*b(1,0) + a.data[8]*b(1,1);
    ret(1,1) = a.data[3]*b(0,0) + a.data[6]*b(0,1) + a.data[8]*b(1,0) + a.data[9]*b(1,1);
    return ret;
}    
    
    
    
template <>
class Matrix4<3, Sym::Major>
{
    Scalar data[45] = {NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN,
                        NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN,
                        NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN,
                        NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN,
                        NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN};

public:
    AMREX_GPU_HOST_DEVICE Matrix4(){};
    
    #if AMREX_SPACEDIM == 3
    Matrix4(const Matrix4<3,Sym::Isotropic> &in)
    {
        // TODO: very inefficient, need to optimize
        for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++)
            (*this)(i,j,k,l) = in(i,j,k,l);
    }
    Matrix4(const Matrix4<3,Sym::Diagonal> &in)
    {
        // TODO: very inefficient, need to optimize
        for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++)
            (*this)(i,j,k,l) = in(i,j,k,l);
    }
    Matrix4(const Matrix4<3,Sym::MajorMinor> &in)
    {
        // TODO: very inefficient, need to optimize
        for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++)
            (*this)(i,j,k,l) = in(i,j,k,l);
    }
    #endif

    AMREX_FORCE_INLINE
    Scalar &operator()(const int i, const int j, const int k, const int l)
    {
        int uid = i + 3 * j + 9 * k + 27 * l;
        
        if (uid == 0)                    return data[0];  // [0000]
        else if (uid == 27 || uid == 3)  return data[1];  // [0001] [0100]
        else if (uid == 54 || uid == 6)  return data[2];  // [0002] [0200]
        else if (uid == 9 || uid == 1)   return data[3];  // [0010] [1000]
        else if (uid == 36 || uid == 4)  return data[4];  // [0011] [1100]
        else if (uid == 63 || uid == 7)  return data[5];  // [0012] [1200]
        else if (uid == 18 || uid == 2)  return data[6];  // [0020] [2000]
        else if (uid == 45 || uid == 5)  return data[7];  // [0021] [2100]
        else if (uid == 72 || uid == 8)  return data[8];  // [0022] [2200]
        else if (uid == 30)              return data[9];  // [0101]
        else if (uid == 57 || uid == 33) return data[10]; // [0102] [0201]
        else if (uid == 12 || uid == 28) return data[11]; // [0110] [1001]
        else if (uid == 39 || uid == 31) return data[12]; // [0111] [1101]
        else if (uid == 66 || uid == 34) return data[13]; // [0112] [1201]
        else if (uid == 21 || uid == 29) return data[14]; // [0120] [2001]
        else if (uid == 48 || uid == 32) return data[15]; // [0121] [2101]
        else if (uid == 75 || uid == 35) return data[16]; // [0122] [2201]
        else if (uid == 60)              return data[17]; // [0202]
        else if (uid == 15 || uid == 55) return data[18]; // [0210] [1002]
        else if (uid == 42 || uid == 58) return data[19]; // [0211] [1102]
        else if (uid == 69 || uid == 61) return data[20]; // [0212] [1202]
        else if (uid == 24 || uid == 56) return data[21]; // [0220] [2002]
        else if (uid == 51 || uid == 59) return data[22]; // [0221] [2102]
        else if (uid == 78 || uid == 62) return data[23]; // [0222] [2202]
        else if (uid == 10)              return data[24]; // [1010]
        else if (uid == 37 || uid == 13) return data[25]; // [1011] [1110]
        else if (uid == 64 || uid == 16) return data[26]; // [1012] [1210]
        else if (uid == 19 || uid == 11) return data[27]; // [1020] [2010]
        else if (uid == 46 || uid == 14) return data[28]; // [1021] [2110]
        else if (uid == 73 || uid == 17) return data[29]; // [1022] [2210]
        else if (uid == 40)              return data[30]; // [1111]
        else if (uid == 67 || uid == 43) return data[31]; // [1112] [1211]
        else if (uid == 22 || uid == 38) return data[32]; // [1120] [2011]
        else if (uid == 49 || uid == 41) return data[33]; // [1121] [2111]
        else if (uid == 76 || uid == 44) return data[34]; // [1122] [2211]
        else if (uid == 70)              return data[35]; // [1212]
        else if (uid == 25 || uid == 65) return data[36]; // [1220] [2012]
        else if (uid == 52 || uid == 68) return data[37]; // [1221] [2112]
        else if (uid == 79 || uid == 71) return data[38]; // [1222] [2212]
        else if (uid == 20)              return data[39]; // [2020]
        else if (uid == 47 || uid == 23) return data[40]; // [2021] [2120]
        else if (uid == 74 || uid == 26) return data[41]; // [2022] [2220]
        else if (uid == 50)              return data[42]; // [2121]
        else if (uid == 77 || uid == 53) return data[43]; // [2122] [2221]
        else if (uid == 80)              return data[44]; // [2222]
        else
            Util::Abort(INFO, "Index out of range");
        return data[-1];
    }
    

    AMREX_FORCE_INLINE
    const Scalar &operator()(const int i, const int j, const int k, const int l) const
    {
        int uid = i + 3 * j + 9 * k + 27 * l;
        
        if (uid == 0)                    return data[0];  // [0000]
        else if (uid == 27 || uid == 3)  return data[1];  // [0001] [0100]
        else if (uid == 54 || uid == 6)  return data[2];  // [0002] [0200]
        else if (uid == 9 || uid == 1)   return data[3];  // [0010] [1000]
        else if (uid == 36 || uid == 4)  return data[4];  // [0011] [1100]
        else if (uid == 63 || uid == 7)  return data[5];  // [0012] [1200]
        else if (uid == 18 || uid == 2)  return data[6];  // [0020] [2000]
        else if (uid == 45 || uid == 5)  return data[7];  // [0021] [2100]
        else if (uid == 72 || uid == 8)  return data[8];  // [0022] [2200]
        else if (uid == 30)              return data[9];  // [0101]
        else if (uid == 57 || uid == 33) return data[10]; // [0102] [0201]
        else if (uid == 12 || uid == 28) return data[11]; // [0110] [1001]
        else if (uid == 39 || uid == 31) return data[12]; // [0111] [1101]
        else if (uid == 66 || uid == 34) return data[13]; // [0112] [1201]
        else if (uid == 21 || uid == 29) return data[14]; // [0120] [2001]
        else if (uid == 48 || uid == 32) return data[15]; // [0121] [2101]
        else if (uid == 75 || uid == 35) return data[16]; // [0122] [2201]
        else if (uid == 60)              return data[17]; // [0202]
        else if (uid == 15 || uid == 55) return data[18]; // [0210] [1002]
        else if (uid == 42 || uid == 58) return data[19]; // [0211] [1102]
        else if (uid == 69 || uid == 61) return data[20]; // [0212] [1202]
        else if (uid == 24 || uid == 56) return data[21]; // [0220] [2002]
        else if (uid == 51 || uid == 59) return data[22]; // [0221] [2102]
        else if (uid == 78 || uid == 62) return data[23]; // [0222] [2202]
        else if (uid == 10)              return data[24]; // [1010]
        else if (uid == 37 || uid == 13) return data[25]; // [1011] [1110]
        else if (uid == 64 || uid == 16) return data[26]; // [1012] [1210]
        else if (uid == 19 || uid == 11) return data[27]; // [1020] [2010]
        else if (uid == 46 || uid == 14) return data[28]; // [1021] [2110]
        else if (uid == 73 || uid == 17) return data[29]; // [1022] [2210]
        else if (uid == 40)              return data[30]; // [1111]
        else if (uid == 67 || uid == 43) return data[31]; // [1112] [1211]
        else if (uid == 22 || uid == 38) return data[32]; // [1120] [2011]
        else if (uid == 49 || uid == 41) return data[33]; // [1121] [2111]
        else if (uid == 76 || uid == 44) return data[34]; // [1122] [2211]
        else if (uid == 70)              return data[35]; // [1212]
        else if (uid == 25 || uid == 65) return data[36]; // [1220] [2012]
        else if (uid == 52 || uid == 68) return data[37]; // [1221] [2112]
        else if (uid == 79 || uid == 71) return data[38]; // [1222] [2212]
        else if (uid == 20)              return data[39]; // [2020]
        else if (uid == 47 || uid == 23) return data[40]; // [2021] [2120]
        else if (uid == 74 || uid == 26) return data[41]; // [2022] [2220]
        else if (uid == 50)              return data[42]; // [2121]
        else if (uid == 77 || uid == 53) return data[43]; // [2122] [2221]
        else if (uid == 80)              return data[44]; // [2222]
        else
            Util::Abort(INFO, "Index out of range");
        return data[-1];
    }

    Set::Scalar Norm()
    {
        Set::Scalar retsq = 0;
        for (int i = 0; i < 45; i++) retsq += data[i]*data[i];
        return std::sqrt(retsq);
    }
    
    
    void Print(std::ostream &os)
    {
        for (int i = 0; i < 45; i++)
            os << "i = " << i << " " << data[i] << std::endl;

        os.precision(4);
        for (int k = 0; k < 3; k++)
        {
            for (int i = -1; i < 4; i++)
            {
                for (int l = 0; l < 3; l++)
                {
                    if (i == -1)
                        os << "┌                                     ┐ ";
                    else if (i == 3)
                        os << "└                                     ┘ ";
                    else
                    {
                        os << "│ ";
                        for (int j = 0; j < 3; j++)
                        {
                            const Set::Scalar &val = (*this)(i, j, k, l);
                            os << std::scientific << std::setw(11) << val; //(fabs(val)>1E-10 ? val : 0);
                            os << " ";
                        }
                        os << "│ ";
                    }
                }
                os << std::endl;
            }
        }
    }
    //AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
    //void operator=(const Matrix4<3, Sym::Major> &a)  { for (int i = 0; i < 45; i++) data[i]  = a.data[i]; }
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
    void operator+=(const Matrix4<3, Sym::Major> &a) { for (int i = 0; i < 45; i++) data[i] += a.data[i]; }
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
    void operator-=(const Matrix4<3, Sym::Major> &a) { for (int i = 0; i < 45; i++) data[i] -= a.data[i]; }
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
    void operator*=(const Matrix4<3, Sym::Major> &a) { for (int i = 0; i < 45; i++) data[i] *= a.data[i]; }
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
    void operator/=(const Matrix4<3, Sym::Major> &a) { for (int i = 0; i < 45; i++) data[i] /= a.data[i]; }
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
    void operator*=(const Set::Scalar &alpha)        { for (int i = 0; i < 45; i++) data[i] *= alpha; }
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
    void operator/=(const Set::Scalar &alpha)        { for (int i = 0; i < 45; i++) data[i] /= alpha; }

    static Matrix4<3, Sym::Major> Increment()
    {
        Matrix4<3, Sym::Major> ret;
        for (int i = 0; i < 45; i++) ret.data[i] = (Set::Scalar)i;
        return ret;
    }
    static Matrix4<3, Sym::Major> Randomize()
    {
        Matrix4<3, Sym::Major> ret;
        for (int i = 0; i < 45; i++) ret.data[i] = (Util::Random());
        return ret;
    }
    static Matrix4<3, Sym::Major> Zero()
    {
        Matrix4<3, Sym::Major> ret;
        for (int i = 0; i < 45; i++) ret.data[i] = 0.0;
        return ret;
    }
    Set::Scalar norm()
    {
        Set::Scalar normsq = 0.0;
        for (int i = 0; i < 45; i++) normsq += data[i]*data[i];
        return std::sqrt(normsq);
    }

    friend Matrix4<3, Sym::Major> operator-(const Matrix4<3, Sym::Major> &a, const Matrix4<3, Sym::Major> &b);
    friend Matrix4<3, Sym::Major> operator+(const Matrix4<3, Sym::Major> &a, const Matrix4<3, Sym::Major> &b);
    friend Set::Matrix operator*(const Matrix4<3, Sym::Major> &a, const Set::Matrix &b);
    friend Matrix4<3, Sym::Major> operator*(const Matrix4<3, Sym::Major> &a, const Set::Scalar &b);
    friend Matrix4<3, Sym::Major> operator/(const Matrix4<3, Sym::Major> &a, const Set::Scalar &b);
};

AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Matrix4<3, Sym::Major> operator+(const Matrix4<3, Sym::Major> &a, const Matrix4<3, Sym::Major> &b)
{
    Matrix4<3,Sym::Major> ret;
    for (int i = 0; i < 45; i++) ret.data[i] = a.data[i] + b.data[i];
    return ret;
}
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Matrix4<3, Sym::Major> operator-(const Matrix4<3, Sym::Major> &a, const Matrix4<3, Sym::Major> &b)
{
    Matrix4<3,Sym::Major> ret;
    for (int i = 0; i < 45; i++) ret.data[i] = a.data[i] - b.data[i];
    return ret;
}
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Matrix4<3, Sym::Major> operator*(const Matrix4<3, Sym::Major> &a, const Set::Scalar &b)
{
    Matrix4<3,Sym::Major> ret;
    for (int i = 0; i < 45; i++) ret.data[i] = a.data[i] * b;
    return ret;
}
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Matrix4<3, Sym::Major> operator/(const Matrix4<3, Sym::Major> &a, const Set::Scalar &b)
{
    Matrix4<3,Sym::Major> ret;
    for (int i = 0; i < 45; i++) ret.data[i] = a.data[i] / b;
    return ret;
}

AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Set::Matrix operator*(const Matrix4<3, Sym::Major> &a, const Set::Matrix &b)
{
    Set::Matrix ret;
    ret(0,0) = a.data[ 0]*b(0,0) + a.data[ 1]*b(0,1) + a.data[ 2]*b(0,2) + a.data[ 3]*b(1,0) + a.data[ 4]*b(1,1) + a.data[ 5]*b(1,2) + a.data[ 6]*b(2,0) + a.data[ 7]*b(2,1) + a.data[ 8]*b(2,2);
    ret(0,1) = a.data[ 1]*b(0,0) + a.data[ 9]*b(0,1) + a.data[10]*b(0,2) + a.data[11]*b(1,0) + a.data[12]*b(1,1) + a.data[13]*b(1,2) + a.data[14]*b(2,0) + a.data[15]*b(2,1) + a.data[16]*b(2,2);
    ret(0,2) = a.data[ 2]*b(0,0) + a.data[10]*b(0,1) + a.data[17]*b(0,2) + a.data[18]*b(1,0) + a.data[19]*b(1,1) + a.data[20]*b(1,2) + a.data[21]*b(2,0) + a.data[22]*b(2,1) + a.data[23]*b(2,2);
    ret(1,0) = a.data[ 3]*b(0,0) + a.data[11]*b(0,1) + a.data[18]*b(0,2) + a.data[24]*b(1,0) + a.data[25]*b(1,1) + a.data[26]*b(1,2) + a.data[27]*b(2,0) + a.data[28]*b(2,1) + a.data[29]*b(2,2);
    ret(1,1) = a.data[ 4]*b(0,0) + a.data[12]*b(0,1) + a.data[19]*b(0,2) + a.data[25]*b(1,0) + a.data[30]*b(1,1) + a.data[31]*b(1,2) + a.data[32]*b(2,0) + a.data[33]*b(2,1) + a.data[34]*b(2,2);
    ret(1,2) = a.data[ 5]*b(0,0) + a.data[13]*b(0,1) + a.data[20]*b(0,2) + a.data[26]*b(1,0) + a.data[31]*b(1,1) + a.data[35]*b(1,2) + a.data[36]*b(2,0) + a.data[37]*b(2,1) + a.data[38]*b(2,2);
    ret(2,0) = a.data[ 6]*b(0,0) + a.data[14]*b(0,1) + a.data[21]*b(0,2) + a.data[27]*b(1,0) + a.data[32]*b(1,1) + a.data[36]*b(1,2) + a.data[39]*b(2,0) + a.data[40]*b(2,1) + a.data[41]*b(2,2);
    ret(2,1) = a.data[ 7]*b(0,0) + a.data[15]*b(0,1) + a.data[22]*b(0,2) + a.data[28]*b(1,0) + a.data[33]*b(1,1) + a.data[37]*b(1,2) + a.data[40]*b(2,0) + a.data[42]*b(2,1) + a.data[43]*b(2,2);
    ret(2,2) = a.data[ 8]*b(0,0) + a.data[16]*b(0,1) + a.data[23]*b(0,2) + a.data[29]*b(1,0) + a.data[34]*b(1,1) + a.data[38]*b(1,2) + a.data[41]*b(2,0) + a.data[43]*b(2,1) + a.data[44]*b(2,2);
    return ret;
}

AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Set::Vector operator * (const Matrix4<AMREX_SPACEDIM,Sym::Major> &a, const Set::Matrix3 &b)
{
    // TODO: improve efficiency of this method
    Set::Vector ret = Set::Vector::Zero();
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        for (int J = 0; J < AMREX_SPACEDIM; J++)
            for (int k = 0; k < AMREX_SPACEDIM; k++)
                for (int L = 0; L < AMREX_SPACEDIM; L++)
                    ret(i) += a(i,J,k,L) * b(k,L,J);
    }    
    return ret;
}

} 
#endif
