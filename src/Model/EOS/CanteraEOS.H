#ifndef MODEL_EOS_CANTERA_H_
#define MODEL_EOS_CANTERA_H_

#include "IO/ParmParse.H"
#include "Model/EOS/EOS.H"

namespace Model
{
namespace EOS
{
class ThermallyPerfect : public EOS
{
public:
    ThermallyPerfect() {};
    virtual ~ThermallyPerfect() {};

    Set::Scalar calcPressure(const Set::Scalar density, const Set::Scalar energy, const Cantera::AnyMap Y) const override
    {

#ifdef ALAMO_CANTERA
    static void Parse(ThermallyPerfect& value, IO::ParmParse& pp)
    {
        if (pp.contains("cantera.yaml")) {
            pp_query("yaml", value.yaml);
            
        }

#endif



#endif
