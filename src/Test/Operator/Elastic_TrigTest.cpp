#include "Elastic.H"
#include "Set/Set.H"
#include "IC/Trig.H"
#include "IC/Affine.H"
#include "IC/Random.H"
#include "Operator/Elastic.H"
#include "Model/Solid/Linear/Laplacian.H"
#include "BC/Operator/Elastic/Constant.H"
//#include "Solver/Nonlocal/Linear.H"
#include "Solver/Nonlocal/Newton.H"

namespace Test
{
namespace Operator
{
int
Elastic::TrigTest(int verbose, int component, int n, std::string plotfile)
{
    Generate();
    
    Set::Scalar tolerance = 0.02;

    int failed = 0;

    // Define the "model" fab to be a Laplacian, so that this
    // elastic operator acts as a Laplacian on the "component-th" component of the fab.
    //using model_type = Model::Solid::Linear::Laplacian; model_type model;
    using model_type = Model::Solid::Linear::Laplacian; model_type model;
    Set::Field<model_type> modelfab(nlevels,ngrids,dmap,1,2);
    for (int ilev = 0; ilev < nlevels; ++ilev) modelfab[ilev]->setVal(model);

    // Initialize: set the rhs_prescribed to sin(n pi x1 / L) * sin(n pi x2 / L), so that
    // the exact solution is sin(n pi x1 / L) * sin(n pi x2 / L) / pi / 2.
    // Set everything else to zero.
    std::complex<int> i(0,1);
    IC::Trig icrhs(geom,1.0,AMREX_D_DECL(n*i,n*i,n*i),dim);
    icrhs.SetComp(component);
    IC::Trig icexact(geom,-(1./dim/Set::Constant::Pi/Set::Constant::Pi/n/n),AMREX_D_DECL(n*i,n*i,n*i),dim);
    icexact.SetComp(component);
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        icrhs.Initialize(ilev,rhs_prescribed);
        icexact.Initialize(ilev,solution_exact);
    }

    amrex::LPInfo info;
    info.setAgglomeration(m_agglomeration);
    info.setConsolidation(m_consolidation);
    if (m_maxCoarseningLevel > -1) info.setMaxCoarseningLevel(m_maxCoarseningLevel);
    nlevels = geom.size();

    ::Operator::Elastic<model_type::sym> elastic;
    elastic.SetUniform(false);
    elastic.define(geom, cgrids, dmap, info);

    // Set up boundary conditions, and 
    // configure the problem so that it is 1D, 2D, or 3D
    BC::Operator::Elastic::Constant bc;
    if (dim == 1)
    {
        AMREX_D_TERM(,// nothing to do in 1D case
                bc.Set(bc.Face::XLO, bc.Direction::Y, bc.Type::Traction, 0.0, rhs_prescribed, geom);
                bc.Set(bc.Face::YLO, bc.Direction::X, bc.Type::Traction, 0.0, rhs_prescribed, geom);
                bc.Set(bc.Face::XHI, bc.Direction::Y, bc.Type::Traction, 0.0, rhs_prescribed, geom);
                bc.Set(bc.Face::YHI, bc.Direction::X, bc.Type::Traction, 0.0, rhs_prescribed, geom);
                ,
                bc.Set(bc.Face::XLO, bc.Direction::Z, bc.Type::Traction, 0.0, rhs_prescribed, geom);
                bc.Set(bc.Face::YLO, bc.Direction::Z, bc.Type::Traction, 0.0, rhs_prescribed, geom);
                bc.Set(bc.Face::ZLO, bc.Direction::X, bc.Type::Traction, 0.0, rhs_prescribed, geom);
                bc.Set(bc.Face::ZLO, bc.Direction::Y, bc.Type::Traction, 0.0, rhs_prescribed, geom);
                bc.Set(bc.Face::XHI, bc.Direction::Z, bc.Type::Traction, 0.0, rhs_prescribed, geom);
                bc.Set(bc.Face::YHI, bc.Direction::Z, bc.Type::Traction, 0.0, rhs_prescribed, geom);
                bc.Set(bc.Face::ZHI, bc.Direction::X, bc.Type::Traction, 0.0, rhs_prescribed, geom);
                bc.Set(bc.Face::ZHI, bc.Direction::Y, bc.Type::Traction, 0.0, rhs_prescribed, geom););
    }
    if (dim == 2)
    {
        AMREX_D_TERM(, // nothing to do in 1D case
                , // nothing to do in 2D case
                bc.Set(bc.Face::XLO, bc.Direction::Z, bc.Type::Traction, 0.0, rhs_prescribed, geom);
                bc.Set(bc.Face::YLO, bc.Direction::Z, bc.Type::Traction, 0.0, rhs_prescribed, geom);
                bc.Set(bc.Face::XHI, bc.Direction::Z, bc.Type::Traction, 0.0, rhs_prescribed, geom);
                bc.Set(bc.Face::YHI, bc.Direction::Z, bc.Type::Traction, 0.0, rhs_prescribed, geom););
    }
    if (dim == 3)
    {
        // nothing to do - displacement BC is default
    }
    bc.Init(rhs_prescribed,geom);
    elastic.SetBC(&bc);


    // Create MLMG solver and solve
    //amrex::MLMG mlmg(elastic);
    //Solver::Nonlocal::Linear mlmg(elastic);
    Solver::Nonlocal::Newton<model_type> mlmg(elastic);
    if (m_fixedIter > -1)     mlmg.setFixedIter(m_fixedIter);
    if (m_maxIter > -1 )      mlmg.setMaxIter(m_maxIter);
    if (m_maxFmgIter > -1)    mlmg.setMaxFmgIter(m_maxFmgIter);
    mlmg.setVerbose(verbose);
    if (m_bottomMaxIter > -1) mlmg.setBottomMaxIter(m_bottomMaxIter);

    mlmg.solve(solution_numeric, rhs_prescribed, modelfab, m_tol_rel,m_tol_abs);

    // Compute solution error
    for (int i = 0; i < nlevels; i++)
    {
        amrex::MultiFab::Copy(*solution_error[i],*solution_numeric[i],component,component,1,2);
        amrex::MultiFab::Subtract(*solution_error[i],*solution_exact[i],component,component,1,2);
    }

    //Compute numerical right hand side
    mlmg.apply(GetVecOfPtrs(rhs_numeric),GetVecOfPtrs(solution_numeric));

    // Compute exact right hand side
    mlmg.apply(GetVecOfPtrs(rhs_exact),GetVecOfPtrs(solution_exact));

    // Compute numerical residual
    mlmg.compResidual(res_numeric,solution_numeric,rhs_prescribed, modelfab);

    // Compute exact residual
    mlmg.compResidual(res_exact,solution_exact,rhs_prescribed, modelfab);

    // If specified, output plot file
    if (plotfile != "")
    {
        Util::Message(INFO,"Printing plot file to ",plotfile);
        WritePlotFile(plotfile);
    }

    // Find maximum solution error
    std::vector<Set::Scalar> error_norm(nlevels);
    for (int i = 0; i < nlevels; i++) error_norm[i] = solution_error[0]->norm0(component,0,false) / solution_exact[0]->norm0(component,0,false);
    Set::Scalar maxnorm = fabs(*std::max_element(error_norm.begin(),error_norm.end()));

    if (verbose) Util::Message(INFO,"relative error = ", 100*maxnorm, " %");
    if (maxnorm > tolerance) failed += 1;
    return failed;
}
}
}
