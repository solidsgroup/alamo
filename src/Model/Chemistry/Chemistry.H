// 
// This header defines the :code:Model::Chemistry::Chemistry class, which introduces a static 
// dispatch mechanism to replicate traditional dynamic polymorphism without incurring the performance 
// cost of virtual table lookups. 
// 
// The :code:Chemistry class acts as the static equivalent of a pure abstract interface for 
// chemical reaction models. The function :code:compute statically dispatches to the 
// appropriate underlying implementation, encapsulating key properties for reaction rates 
// and heat release. The :code:compute function is resolved at compile-time for GPU 
// compatibility and runtime performance.
// 
// This system enables polymorphic behavior without relying on runtime virtual functions, making it 
// suitable for high-performance computing environments.
// 

#ifndef MODEL_CHEMISTRY_H
#define MODEL_CHEMISTRY_H

#include "AMReX_Extension.H"
#include "IO/ParmParse.H"
#include "Model/Gas/Gas.H"

namespace Model
{
namespace Chemistry
{
template<typename... OBJ>
class Chemistry
{
private:
    friend IO::ParmParse;
    std::tuple<OBJ...> obj;
    static constexpr const char* names[] = { OBJ::name... };
    int selected = -1;

public:

    template<int I = 0>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    std::pair<std::array<double,NSPECIES>,Set::Scalar>
    compute (const Set::Scalar P, const Set::Scalar T, const std::array<double,NSPECIES>& rhoY,
                const Set::Scalar dt, const Model::Gas::Gas* gas = nullptr)
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I) return std::get<I>(obj).compute(P,T,rhoY,dt,gas);
            else return compute<I+1>(P,T,rhoY,dt,gas);
        }
        else 
        {
            Util::Abort(INFO,selected);
            return {};
        }
    }
};
}
}


#endif
