//
// Set each point in the field to a random value between 0 and 1
//

#ifndef IC_RANDOM_H_
#define IC_RANDOM_H_

#include "AMReX_Vector.H"

#include "IO/ParmParse.H"
#include "IC/IC.H"
#include "Util/Util.H"
#include "Set/Set.H"
#include "IO/ParmParse.H"

namespace IC
{
/// Set each point to a random value.
class Random : public IC<Set::Scalar>
{
public:
    static constexpr const char* name = "random";

    Random (amrex::Vector<amrex::Geometry> &_geom, Set::Scalar a_mult) :
        IC(_geom), mult(a_mult)
    {}
    Random(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name) : IC(_geom)
    {
        pp.queryclass(name, *this);
    }
  
    void Add(const int &lev, Set::Field<Set::Scalar> &a_field, Set::Scalar)
    {
        for (amrex::MFIter mfi(*a_field[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx;
            amrex::IndexType type = a_field[lev]->ixType();
            if (type == amrex::IndexType::TheCellType())      bx = mfi.growntilebox();
            else if (type == amrex::IndexType::TheNodeType()) bx = mfi.grownnodaltilebox();
            else Util::Abort(INFO, "Unkonwn index type");
            
            Set::Patch<Set::Scalar> field = a_field.Patch(lev,mfi);
            for (int n = 0; n < a_field[lev]->nComp(); n++)
            {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    field(i,j,k,n) += offset + mult * Util::Random();
                });
            }
        }
    }

    using IC::Add;

    static void Parse(Random & value, IO::ParmParse & pp)
    {
        // offset from the [0,1] random number range
        pp.query_default("offset",value.offset,0.0); 
        // multiplier for the [0,1] random number range
        pp.query_default("mult",value.mult,1.0);
    }
    
private:
    Set::Scalar mult = NAN; 
    Set::Scalar offset = NAN;

};
}
#endif

