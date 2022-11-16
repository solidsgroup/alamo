#ifndef MODEL_SOLID_ELASTIC_MYMODEL_H_
#define MODEL_SOLID_ELASTIC_MYMODEL_H_

#include "IO/ParmParse.H"
#include "Model/Solid/Solid.H"

namespace Model
{
namespace Solid
{
namespace Elastic
{
class MyModel : public Solid<Set::Sym::Major>
{
public:
    MyModel() {};
    MyModel(Solid<Set::Sym::Major> base) : Solid<Set::Sym::Major>(base) {};
    virtual ~MyModel() {};

    Set::Scalar W(const Set::Matrix & F) const override
    {
        Set::Scalar w = 0.0;

        //
        // TODO: set the value for w
        // 

        return w;
    }
    Set::Matrix DW(const Set::Matrix & F) const override
    {
        Set::Matrix dw = Set::Matrix::Zero();

        //
        // TODO: compute dw
        //

        return dw;
    }
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Major> DDW(const Set::Matrix & F) const override
    {
        Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Major> ddw;

        //
        // TODO: compute ddw
        //

        return ddw;
    }
    
    void Advance(Set::Scalar dt, Set::Matrix eps, Set::Matrix sig) override
    {
        //
        // TODO: (optional)
        //       Implement kinetic equations to evolve internal
        //       variables. If you are working with a hyperelastic
        //       model, you can leave this blank.
        //
        
    }

public:

    //
    // TODO: replace "prop" with your own material properties
    //
    Set::Scalar prop;

    static constexpr KinematicVariable kinvar = KinematicVariable::F;

public:
    static MyModel Zero()
    {
        MyModel ret;

        //
        // TODO: Set all material properties in "ret" to zero
        //
        ret.prop = 0.0;

        return ret;
    }
    static MyModel Random()
    {
        MyModel ret;

        //
        // TODO: Assign randomized properties in "ret"
        //
        ret.prop = Util::Random();

        return ret;
    }
    static void Parse(MyModel & value, IO::ParmParse & pp)
    {

        //
        // TODO: Read in material properties from input file.
        //       This is optional, you can also hard-code them.
        //
        pp.query("prop",value.prop);

    }

#define OP_CLASS MyModel
//
// TODO: Replace X(prop) with X(prop1) X(prop2)
//       etc, where prop1 and prop2 are your own material properties.
//
#define OP_VARS X(prop)


#include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"

}
}
}

#endif
