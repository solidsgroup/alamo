#ifndef TEST_NUMERIC_STENCIL
#define TEST_NUMERIC_STENCIL

#include <AMReX.H>

#include "Set/Set.H"
#include "IC/Trig.H"
#include "IC/Trig2.H"
#include "Numeric/Stencil.H"

namespace Test
{
/// Tests for the Numeric namespace classes
namespace Numeric
{
class Stencil
{

public:
    Stencil() {};
    ~Stencil() {};

    void Define(amrex::IntVect _ncells)
    {
        ncells = _ncells;
        amrex::RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(L,L,L)});
        amrex::Geometry::Setup(&rb, 0);
        
        amrex::Box domain(amrex::IntVect{AMREX_D_DECL(0,0,0)}, ncells,
                amrex::IntVect::TheCellVector());
        geom.resize(1);
        grids.resize(1);
        dmap.resize(1);
        phi.resize(1);
        DphiNumeric.resize(1);
        DphiExact.resize(1);

        geom[0].define(domain);
        grids[0].define(domain);
        grids[0].maxSize(ncells[0]/2);
        dmap[0].define(grids[0]);

        phi        .Define(0,grids[0],dmap[0],1,nghost);
        DphiNumeric.Define(0,grids[0],dmap[0],1,nghost);
        DphiExact  .Define(0,grids[0],dmap[0],1,nghost);

        phi[0]->setVal(0.0);
        // Initial function is
        //   cos(2*pi*x/L) * cos(2*pi*y/L) * cos(2*pi*z / L)
        IC::Trig2 ic(geom,1.0,AMREX_D_DECL(0.0,0.0,0.0),AMREX_D_DECL(1,1,1));
        ic.Add(0,phi);
    }
    void Define(int _ncells)
    {
        Define(amrex::IntVect{AMREX_D_DECL(_ncells,_ncells,_ncells)});
    }

    template<int i,int j,int k>
    bool Derivative(int verbose,std::string plotfile = "")
    {
        const Set::Scalar tolerance = 1E-2;
        
        DphiExact[0]->setVal(0.0);
        DphiNumeric[0]->setVal(0.0);

        //
        // Compute exact derivative
        //
        Set::Scalar AMREX_D_DECL(theta1 = 0.0, theta2 = 0.0, theta3 = 0.0);
        Set::Scalar fac = 1.0;
        AMREX_D_TERM( for (int p = 0; p < i; p++) theta1 += 0.5*Set::Constant::Pi;,
                for (int q = 0; q < j; q++) theta2 += 0.5*Set::Constant::Pi;,
                for (int r = 0; r < k; r++) theta3 += 0.5*Set::Constant::Pi; );
                  
        AMREX_D_TERM( for (int p = 0; p < i; p++) fac *= 1.0*Set::Constant::Pi / L;,
                for (int q = 0; q < j; q++) fac *= 1.0*Set::Constant::Pi / L;,
                for (int r = 0; r < k; r++) fac *= 1.0*Set::Constant::Pi / L;);
        IC::Trig2 ic(geom,fac,AMREX_D_DECL(theta1,theta2,theta3),AMREX_D_DECL(1,1,1));
        ic.Add(0,DphiExact);


        //
        // Compute numeric derivative
        //
        const amrex::Real* DX = geom[0].CellSize();
        for ( amrex::MFIter mfi(*phi[0],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<const amrex::Real> const& Phi = phi[0]->array(mfi);
            amrex::Array4<amrex::Real> const& DPhi    = DphiNumeric[0]->array(mfi);
            amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int p, int q, int r){ DPhi(p,q,r) = ::Numeric::Stencil<Set::Scalar,i,j,k>::D(Phi,p,q,r,0,DX); });
        }

        amrex::MultiFab diff(grids[0],dmap[0],ncomp,nghost);
        amrex::MultiFab::Copy    (diff,*DphiExact[0]  ,0,0,ncomp,nghost);
        amrex::MultiFab::Subtract(diff,*DphiNumeric[0],0,0,ncomp,nghost);

        Set::Scalar error = diff.norm0(0,0,false) / DphiExact[0]->norm0(0,0,false);

        if (verbose) Util::Message(INFO,"Infinity norm of error = ", 100*error,"%");

        if (plotfile != "")
            WritePlotFile(plotfile);

        if (error < tolerance) return 0;
        else return 1;
    }


    void WritePlotFile(std::string plotfile)
    {
        amrex::Vector<amrex::MultiFab> plotmf(1);
        plotmf[0].define(grids[0], dmap[0], 3, 0);
        for (amrex::MFIter mfi(plotmf[0],false); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.validbox();
            plotmf[0][mfi].copy((*phi[0])[mfi]        ,bx, 0, bx, 0, 1); 
            plotmf[0][mfi].copy((*DphiExact[0])[mfi]  ,bx, 0, bx, 1, 1); 
            plotmf[0][mfi].copy((*DphiNumeric[0])[mfi],bx, 0, bx, 2, 1); 
        }
        amrex::WriteSingleLevelPlotfile(plotfile,plotmf[0],varnames,geom[0],0.0,1);
    }        


private:
    amrex::IntVect ncells;
    const Set::Scalar L = 1.0;
    const int nlevels = 1, ref_ratio = 1, ncomp = 1, nghost = 2;
    amrex::Vector<amrex::Geometry> geom;
    amrex::Vector<amrex::BoxArray> grids;
    amrex::Vector<amrex::DistributionMapping> dmap;
    
    Set::Field<Set::Scalar> phi;
    Set::Field<Set::Scalar> DphiExact;
    Set::Field<Set::Scalar> DphiNumeric;
    amrex::Vector<std::string> varnames = {"phi","dphi_exact","dphi_numeric"};
};
}
}



#endif
