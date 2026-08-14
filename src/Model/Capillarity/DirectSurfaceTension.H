#ifndef MODEL_CAPILLARITY_DIRECTSURFACETENSION_H
#define MODEL_CAPILLARITY_DIRECTSURFACETENSION_H

#include "IO/ParmParse.H"

namespace Model
{
namespace Capillarity
{
// Conservative Korteweg momentum coupling with purely advected phase
// fractions.  This is the least expensive capillarity model and introduces
// no separate phase-field kinetics.
class DirectSurfaceTension
{
public:
    static constexpr const char* name = "direct_surface_tension";
    static void Parse(DirectSurfaceTension&, IO::ParmParse&) {}
};
}
}

#endif
