#ifndef IC_SHOCK_H_
#define IC_SHOCK_H_

#include "Set/Set.H"
#include "IC/IC.H"
#include "IO/ParmParse.H"
#include "Util/ScimitarX_Util.H"
#include "Numeric/IntegratorVariableAccessLayer.H"
#include "AMReX_Parser.H"
#include <vector> 

namespace IC {
class Shock : public IC<Set::Scalar> {
public:
    // Constructor requiring variable indices
    Shock(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name,
            const Numeric::GenericVariableAccessor::VariableIndices& precomputed_indices);

    // Constructor for pressure (no variable indices required)
    Shock(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name);

    // Add field variables or pressure
    virtual void Add(const int& lev, Set::Field<Set::Scalar>& a_phi, Set::Scalar time);

protected:

    // Implementation for constant values
    void AddConstant(const int& lev, Set::Field<Set::Scalar>& a_phi, Set::Scalar time);
    
    // Implementation for expressions
    void AddExpression(const int& lev, Set::Field<Set::Scalar>& a_phi, Set::Scalar time);

    // Initialize the Shock properties
    void initialize(IO::ParmParse& pp, const std::string& name);

    // MultiFab type name
    std::string mf_name;

    // Variable indices (only required for qvec)
    const Numeric::GenericVariableAccessor::VariableIndices* variable_indices;
    bool requires_variable_indices;

    // Number of shocks (zones)
    int number_of_shocks;   

    // Shock positions in x-direction
    std::vector<Set::Scalar> shock_positions;

    // Shock directions (xdir, ydir and zdir)
    std::string shock_direction;

    int direction_index;  // Direction indices (0=x, 1=y, 2=z)
    std::vector<Set::Scalar> zone_lower_bounds;  // Pre-computed lower bounds for each zone
    std::vector<Set::Scalar> zone_upper_bounds;  // Pre-computed upper bounds for each zone


private:

    // Flag to determine whether to use expressions for initialization
    bool use_expressions = false;

    // pvec initialization variables
    std::vector<Set::Scalar> density_zone;
    std::vector<Set::Scalar> uvel_zone;
    std::vector<Set::Scalar> vvel_zone;
    std::vector<Set::Scalar> wvel_zone;
    // pressure initialization variables
    std::vector<Set::Scalar> pressure_zone;
    //ElastoPlastic
    std::vector<Set::Scalar> sxx_zone;
    std::vector<Set::Scalar> syy_zone;
    std::vector<Set::Scalar> szz_zone;
    std::vector<Set::Scalar> sxy_zone;
    std::vector<Set::Scalar> sxz_zone;
    std::vector<Set::Scalar> syz_zone;
    std::vector<Set::Scalar> ie_elastic_zone;
    std::vector<Set::Scalar> epsbar_zone;

    // Expression strings - to store directly from input, similar to density_zone format
    std::vector<std::string> density_expr;
    std::vector<std::string> pressure_expr;
    std::vector<std::string> uvel_expr;
    std::vector<std::string> vvel_expr;
    std::vector<std::string> wvel_expr;    
    
    // Parsers and executors for expressions (similar to Expression class)
    std::vector<amrex::Parser> density_parser;
    std::vector<amrex::Parser> pressure_parser;
    std::vector<amrex::Parser> uvel_parser;
    std::vector<amrex::Parser> vvel_parser;
    std::vector<amrex::Parser> wvel_parser;
    //ElastoPlastic
    std::vector<amrex::Parser> sxx_parser;
    std::vector<amrex::Parser> syy_parser;
    std::vector<amrex::Parser> szz_parser;
    std::vector<amrex::Parser> sxy_parser;
    std::vector<amrex::Parser> sxz_parser;
    std::vector<amrex::Parser> syz_parser;
    std::vector<amrex::Parser> ie_elastic_parser;
    std::vector<amrex::Parser> epsbar_parser;
    
    std::vector<amrex::ParserExecutor<4>> density_func;
    std::vector<amrex::ParserExecutor<4>> pressure_func;
    std::vector<amrex::ParserExecutor<4>> uvel_func;
    std::vector<amrex::ParserExecutor<4>> vvel_func;
    std::vector<amrex::ParserExecutor<4>> wvel_func;
    //ElastoPlastic
    std::vector<amrex::ParserExecutor<4>> sxx_func;
    std::vector<amrex::ParserExecutor<4>> syy_func;
    std::vector<amrex::ParserExecutor<4>> szz_func;
    std::vector<amrex::ParserExecutor<4>> sxy_func;
    std::vector<amrex::ParserExecutor<4>> sxz_func;
    std::vector<amrex::ParserExecutor<4>> syz_func;
    std::vector<amrex::ParserExecutor<4>> epsbar_func;
    std::vector<amrex::ParserExecutor<4>> ie_elastic_func;

    // Helper method to compile a parser from an expression string
    void SetupParser(std::vector<amrex::Parser>& parser, 
                    std::vector<amrex::ParserExecutor<4>>& func,
                    const std::vector<std::string>& expr_list,
                    const std::set<std::pair<std::string, Set::Scalar>>& constants);
    

};

} // namespace IC

#endif // IC_SHOCK_H_
