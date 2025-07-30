#ifndef IC_LS_EXPRESSION_H_
#define IC_LS_EXPRESSION_H_
#include "IC/IC.H"
#include "Util/Util.H"
#include "IO/ParmParse.H"
#include "AMReX_Parser.H"

struct LevelSetData;
struct NarrowBandTubeType;
struct Neighbors;

namespace IC {
namespace LS {

class Expression : public IC<Set::Scalar>
{
public:
    static constexpr const char* name = "expression";

    Expression(amrex::Vector<amrex::Geometry>& _geom);

    Expression(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name);

    void SetTubeField(Set::Field<Set::Scalar>* tube_field_ptr);

    void SetZeroField(Set::Field<Set::Scalar>* zero_field_ptr);

    void SetCPTField(Set::Field<Set::Scalar>* cpt_field_ptr);

    void SetLSoldField(Set::Field<Set::Scalar>* LSold_field_ptr);

    void SetLSData(LevelSetData* data_ptr);

    virtual void Add(const int& lev, Set::Field<Set::Scalar>& a_field, Set::Scalar a_time = 0.0) override;

    static void Parse(Expression& value, IO::ParmParse& pp);

private:
    enum CoordSys { Cartesian, Polar, Spherical };
    std::vector<amrex::Parser> parser;
    std::vector<amrex::ParserExecutor<4>> f;
    CoordSys coord = Cartesian;
    
    Set::MultiVector levelset_id;
    Set::Scalar INNERBAND_Width = 4.0;
    Set::Scalar BAND_WIDTH = 6.0;
    Set::Field<Set::Scalar>* Tube_mf = nullptr;
    Set::Field<Set::Scalar>* Zero_mf = nullptr;
    Set::Field<Set::Scalar>* CPT_mf = nullptr;
    Set::Field<Set::Scalar>* LSold_mf = nullptr;
    LevelSetData* LSData = nullptr;
};

}} // namespace IC::LS

#endif
