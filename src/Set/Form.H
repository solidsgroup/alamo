#ifndef SET_FORM_H
#define SET_FORM_H

#include "Field.H"

namespace Set
{
    
template<int D, class T>
class Form {};

template<class T>
class Form<0,T> : public Field<T>
{
    public:        
    class ConstPatch : public amrex::Array4<const T>
    {
        public:
        ConstPatch(amrex::Array4<T> &rhs) : 
        amrex::Array4<T>(rhs) {}
    };
};

    
}
#endif