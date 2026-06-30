#ifndef MODEL_GAS_THERMO_H_
#define MODEL_GAS_THERMO_H_

#include <vector>
#include <memory>
#include "Set/Set.H"
#include <cmath>
#include "IO/ParmParse.H"
#include "Set/Set.H"

namespace Model {
namespace Gas {
namespace Thermo {

template<typename... OBJ>
class Thermo {
private:
    friend IO::ParmParse;
    std::tuple<OBJ...> obj;
    static constexpr const char* names[] = { OBJ::name... };
    int selected = -1;


public:
    template<int I = 0>
    AMREX_FORCE_INLINE
    double cp_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).cp_mol(T,X,i,j,k);
            else
                return cp_mol<I+1>(T,X,i,j,k);
        }
        else {
            Util::Abort(INFO, selected);
            return NAN;
        }
    }

    template<int I = 0>
    AMREX_FORCE_INLINE
    double enthalpy_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).enthalpy_mol(T,X,i,j,k);
            else
                return enthalpy_mol<I+1>(T,X,i,j,k);
        }
        else {
            Util::Abort(INFO, selected);
            return NAN;
        }
    }

    template<int I = 0>
    AMREX_FORCE_INLINE
    double entropy_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).entropy_mol(T,X,i,j,k);
            else
                return entropy_mol<I+1>(T,X,i,j,k);
        }
        else {
            Util::Abort(INFO, selected);
            return NAN;
        }
    }

    template<int I = 0>
    AMREX_FORCE_INLINE
    double cp_mol_species(double T, int n) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).cp_mol_species(T,n);
            else
                return cp_mol_species<I+1>(T,n);
        }
        else {
            Util::Abort(INFO, selected);
            return NAN;
        }
    }

    template<int I = 0>
    AMREX_FORCE_INLINE
    double enthalpy_mol_species(double T, int n) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).enthalpy_mol_species(T,n);
            else
                return enthalpy_mol_species<I+1>(T,n);
        }
        else {
            Util::Abort(INFO, selected);
            return NAN;
        }
    }

    template<int I = 0>
    AMREX_FORCE_INLINE
    double entropy_mol_species(double T, int n) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).entropy_mol_species(T,n);
            else
                return entropy_mol_species<I+1>(T,n);
        }
        else {
            Util::Abort(INFO, selected);
            return NAN;
        }
    }

    ~Thermo() = default;
}; // class Thermo

} // namespace Thermo
} // namespace Gas
} // namespace Model

#endif
