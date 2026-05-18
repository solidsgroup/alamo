#ifndef MODEL_GAS_THERMO_ROCFIRE_H_
#define MODEL_GAS_THERMO_ROCFIRE_H_

#include <vector>
#include <memory>
#include <cmath>
#include "IO/ParmParse.H"
#include "Model/Gas/Thermo/Thermo.H"

namespace Model {
namespace Gas {
namespace Thermo {

class Rocfire : public Thermo {
public:
    static constexpr const char* name = "rocfire";
    const char* model_name() const override { return name; }

private:
    double mw;
    double cp_mass;
    
public:
    Rocfire() {} ;
    Rocfire(IO::ParmParse& pp, std::string name) : Rocfire()
    {
        pp.queryclass(name, *this);
    }

    ~Rocfire() override = default;

    static void Parse(Rocfire & value, IO::ParmParse & pp)
    {
        pp.query_default("mw", value.mw, "26.0_g/mol", Unit::Mass()/Unit::Amount());
        pp.query_default("cp_mass", value.cp_mass, "0.3_cal/g/K", Unit::SpecificHeatCapacity());

        Util::Assert(INFO, TEST(value.mw > 0.0), "gas.thermo.rocfire.mw must be positive.");
        Util::Assert(INFO, TEST(value.cp_mass > 0.0), "gas.thermo.rocfire.cp_mass must be positive.");
    }

public:
    double cp_mol(double /*T*/, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const override {
        double cp = 0.0;
        for (int n = 0; n < NSPECIES; ++n)
        {
            cp += X(i,j,k,n) * cp_mol_species(0.0, n);
        }
        return cp;
    }        
    double enthalpy_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const override {
        double h = 0.0;
        for (int n = 0; n < NSPECIES; ++n)
        {
            h += X(i,j,k,n) * enthalpy_mol_species(T, n);
        }
        return h;
    }
    double entropy_mol(double /*T*/, Set::Patch<const Set::Scalar>& /*X*/, int /*i*/, int /*j*/, int /*k*/) const override {
        return 0.0;
    }
    double cp_mol_species(double /*T*/, int /*n*/) const override {
        return cp_mass * mw;
    }
    double enthalpy_mol_species(double T, int /*n*/) const override {
        return cp_mass * mw * T;
    }
    double entropy_mol_species(double /*T*/, int /*n*/) const override {
        return 0.0;
    }

}; // class Rocfire

} // namespace Thermo
} // namespace Gas
} // namespace Model

#endif
