#include <AMReX_MLMG.H>

namespace Solver
{
/// \brief Multigrid Linear solver for multicomponent, multi-level operators
/// 
/// This class is a thin wrapper for the `amrex::MLMG` solver.
/// It exists to set a range of default MLMG settings automatically, for instance,
/// `setCFStrategy`, which may not be obvious to the user.
///
/// It also exists as a compatibility layer so that future fixes for compatibility
/// with AMReX can be implemented here.
class Linear : public amrex::MLMG
{
    public:
    /// This constructor sets several defaults for NODAL solves
    Linear (Operator::Operator<Grid::Node>& a_lp) : MLMG(a_lp)
    {
        MLMG::setBottomSolver(MLMG::BottomSolver::bicgstab);
	    MLMG::setCFStrategy(MLMG::CFStrategy::ghostnodes);
        MLMG::setFinalFillBC(false);	
    }
    Linear (Operator::Operator<Grid::Cell>& a_lp): MLMG(a_lp)
    {
        MLMG::setBottomSolver(MLMG::BottomSolver::bicgstab);
	    MLMG::setCFStrategy(MLMG::CFStrategy::none);
        MLMG::setFinalFillBC(false);	
    }
    void setVerbose(int verbosity)
    {
        MLMG::setVerbose(verbosity);
		if (verbosity > 4) MLMG::setCGVerbose(verbosity);
     	else
 	    {
 		    MLMG::setVerbose(0);
 		    MLMG::setCGVerbose(0);
	    }
    }

};
}