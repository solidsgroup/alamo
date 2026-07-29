//
// This header defines the :code:Model::Gas::Thermo::Thermo:code: static dispatch
// wrapper for thermodynamic property models.
//
// The :code:Thermo:code: class selects one thermodynamic model at input parsing time
// and dispatches mixture and species calls for molar heat capacity, enthalpy, and
// entropy. This keeps thermodynamic property evaluation lightweight in per-cell gas,
// chemistry, and equation-of-state calculations.
//

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
    template<typename T>
    const T& Get() const
    {
        return std::get<T>(obj);
    }

    const char* model_name() const
    {
        return selected >= 0 ? names[selected] : "";
    }

    template<class Composition, int I = 0>
    AMREX_FORCE_INLINE
    Set::Scalar cp_mol(Set::Scalar T, const Composition& X, int i, int j, int k) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).cp_mol(T,X,i,j,k);
            else
                return cp_mol<Composition,I+1>(T,X,i,j,k);
        }
        else {
            Util::Abort(INFO, selected);
            return NAN;
        }
    }

    template<class Composition, int I = 0>
    AMREX_FORCE_INLINE
    Set::Scalar enthalpy_mol(Set::Scalar T, const Composition& X, int i, int j, int k) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).enthalpy_mol(T,X,i,j,k);
            else
                return enthalpy_mol<Composition,I+1>(T,X,i,j,k);
        }
        else {
            Util::Abort(INFO, selected);
            return NAN;
        }
    }

    template<class Composition, int I = 0>
    AMREX_FORCE_INLINE
    Set::Scalar entropy_mol(Set::Scalar T, const Composition& X, int i, int j, int k) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).entropy_mol(T,X,i,j,k);
            else
                return entropy_mol<Composition,I+1>(T,X,i,j,k);
        }
        else {
            Util::Abort(INFO, selected);
            return NAN;
        }
    }

    template<int I = 0>
    AMREX_FORCE_INLINE
    Set::Scalar cp_mol_species(Set::Scalar T, int n) const
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
    Set::Scalar dcp_mol_species_dT(Set::Scalar T, int n) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).dcp_mol_species_dT(T,n);
            else
                return dcp_mol_species_dT<I+1>(T,n);
        }
        else {
            Util::Abort(INFO, selected);
            return NAN;
        }
    }

    template<int I = 0>
    AMREX_FORCE_INLINE
    Set::Scalar enthalpy_mol_species(Set::Scalar T, int n) const
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
    Set::Scalar entropy_mol_species(Set::Scalar T, int n) const
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
