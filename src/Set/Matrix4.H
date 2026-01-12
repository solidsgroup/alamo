#ifndef SET_MATRIX4_H
#define SET_MATRIX4_H

#include <iostream>

#include "Util/Util.H"
#include "Set/Base.H"

namespace Set
{

enum Sym {None,Major,Minor,MajorMinor,Diagonal,Full,Isotropic};

template<int DIM, int SYM>
class Matrix4
{
};

template<int dim, int sym>
std::ostream&
operator<< (std::ostream& os, const Matrix4<dim,sym>& b)
{
    Matrix4<dim,sym> bcopy = b;
        bcopy.Print(os);
    return os;
}

}



#endif
