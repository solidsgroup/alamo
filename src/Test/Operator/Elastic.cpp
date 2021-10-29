#include "Elastic.H"
#include "Set/Set.H"
#include "IC/Trig.H"
#include "IC/Affine.H"
#include "IC/Random.H"
#include "Operator/Elastic.H"

namespace Test
{
namespace Operator
{
void Elastic::Define(const amrex::IntVect _ncells,
            const int _nlevels,
            const int _dim,
            const Grid _config)
{
    dim = _dim;
    ncells = _ncells;
    nlevels = _nlevels;
    m_config = _config;
}

void Elastic::Define(const int _ncells,
            const int _nlevels,
            const int _dim,
            const Grid _config)
{
    Define(amrex::IntVect(AMREX_D_DECL(_ncells,_ncells,_ncells)),_nlevels,_dim,_config);
}

void Elastic::Generate()
{
    //int max_grid_size = 100000;
    int max_grid_size = ncells[0]/4;
    //std::string orientation = "h";
    geom.resize(nlevels); 
    cgrids.resize(nlevels);
    ngrids.resize(nlevels);
    dmap.resize(nlevels);

    solution_exact.resize(nlevels);   solution_exact.finest_level = nlevels-1;
    solution_numeric.resize(nlevels); solution_numeric.finest_level = nlevels-1;
    solution_error.resize(nlevels);   solution_error.finest_level = nlevels-1;
    rhs_prescribed.resize(nlevels);   rhs_prescribed.finest_level = nlevels-1;
    rhs_numeric.resize(nlevels);      rhs_numeric.finest_level = nlevels-1;
    rhs_exact.resize(nlevels);        rhs_exact.finest_level = nlevels-1;
    res_numeric.resize(nlevels);      res_numeric.finest_level = nlevels-1;
    res_exact.resize(nlevels);        res_exact.finest_level = nlevels-1;
    ghost_force.resize(nlevels);      ghost_force.finest_level = nlevels-1;

    amrex::RealBox rb({AMREX_D_DECL(0.,0.,0.)},
            {AMREX_D_DECL(m_bounds[0],m_bounds[1],m_bounds[2])});
    amrex::Geometry::Setup(&rb, 0);

    amrex::Box NDomain(amrex::IntVect{AMREX_D_DECL(0,0,0)}, ncells,
                amrex::IntVect::TheNodeVector());
    amrex::Box CDomain = amrex::convert(NDomain, amrex::IntVect::TheCellVector());

    amrex::Box domain = CDomain;
    for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            geom[ilev].define(domain);
            domain.refine(ref_ratio);
        }
    amrex::Box cdomain = CDomain;

    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        cgrids[ilev].define(cdomain);
        // if (ilev == 0) cgrids[ilev].maxSize(10000000);
        // if (ilev == 1) cgrids[ilev].maxSize(max_grid_size);
        // if (ilev == 2) cgrids[ilev].maxSize(10000000);
        cgrids[ilev].maxSize(max_grid_size);
        if (m_config == Grid::XYZ)
            cdomain.grow(amrex::IntVect(AMREX_D_DECL(-ncells[0]/4,-ncells[1]/4,-ncells[2]/4))); 
        else if (m_config == Grid::X)
            cdomain.grow(amrex::IntVect(AMREX_D_DECL(-ncells[0]/4,0,0)));
        else if (m_config == Grid::Y)
            cdomain.grow(amrex::IntVect(AMREX_D_DECL(0,-ncells[1]/4,0)));
        else if (m_config == Grid::Z)
            cdomain.grow(amrex::IntVect(AMREX_D_DECL(0,0,-ncells[2]/4)));
        else if (m_config == Grid::YZ)
            cdomain.grow(amrex::IntVect(AMREX_D_DECL(0,-ncells[1]/4,-ncells[2]/4)));
        else if (m_config == Grid::ZX)
            cdomain.grow(amrex::IntVect(AMREX_D_DECL(-ncells[0]/4,0,-ncells[2]/4)));
        else if (m_config == Grid::XY)
            cdomain.grow(amrex::IntVect(AMREX_D_DECL(-ncells[0]/4,-ncells[1]/4,0)));
    
        cdomain.refine(ref_ratio); 
        ngrids[ilev] = cgrids[ilev];
        ngrids[ilev].convert(amrex::IntVect::TheNodeVector());
    }

    int number_of_components = AMREX_SPACEDIM;
    for (int ilev = 0; ilev < nlevels; ++ilev) dmap   [ilev].define(cgrids[ilev]);
    solution_numeric.Define(nlevels, ngrids, dmap, number_of_components, 2); 
    solution_exact  .Define(nlevels, ngrids, dmap, number_of_components, 2);
    solution_error  .Define(nlevels, ngrids, dmap, number_of_components, 2);
    rhs_prescribed  .Define(nlevels, ngrids, dmap, number_of_components, 2);
    rhs_numeric     .Define(nlevels, ngrids, dmap, number_of_components, 2);
    rhs_exact       .Define(nlevels, ngrids, dmap, number_of_components, 2);
    res_numeric     .Define(nlevels, ngrids, dmap, number_of_components, 2); 
    res_exact       .Define(nlevels, ngrids, dmap, number_of_components, 2); 
    ghost_force     .Define(nlevels, ngrids, dmap, number_of_components, 2); 

    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        solution_exact  [ilev]->setVal(0.0);
        solution_numeric[ilev]->setVal(0.0);
        solution_error  [ilev]->setVal(0.0);
        rhs_prescribed  [ilev]->setVal(0.0);
        rhs_exact       [ilev]->setVal(0.0);
        rhs_numeric     [ilev]->setVal(0.0);
        res_exact       [ilev]->setVal(0.0);
        res_numeric     [ilev]->setVal(0.0);
        ghost_force     [ilev]->setVal(0.0);
    }

}
}
}
         
