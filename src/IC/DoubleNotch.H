#ifndef IC_DOUBLENOTCH_H_
#define IC_DOUBLENOTCH_H_

#include "Set/Set.H"
#include "IC/IC.H"
#include "IO/ParmParse.H"

namespace IC
{
class DoubleNotch : public IC
{
public:
    enum Mollifier {Dirac, Gaussian};

	DoubleNotch (amrex::Vector<amrex::Geometry> &_geom) : IC(_geom) 
    {
    }

	
	void Add(const int lev, amrex::Vector<amrex::MultiFab * > &a_field)
	{
        Set::Vector DX(geom[lev].CellSize());

		for (amrex::MFIter mfi(*a_field[lev],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
		{
			amrex::Box bx = mfi.tilebox();
		    bx.grow(a_field[lev]->nGrow());
			amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {

				Set::Vector x;
                // NODE
                if (a_field[0]->ixType() == amrex::IndexType::TheNodeType())
                {
				    AMREX_D_TERM(x(0) = geom[lev].ProbLo()[0] + ((amrex::Real)(i)) * geom[lev].CellSize()[0];,
				    			 x(1) = geom[lev].ProbLo()[1] + ((amrex::Real)(j)) * geom[lev].CellSize()[1];,
				    			 x(2) = geom[lev].ProbLo()[2] + ((amrex::Real)(k)) * geom[lev].CellSize()[2];);
                }
                else if (a_field[0]->ixType() == amrex::IndexType::TheCellType())
                {
				    AMREX_D_TERM(x(0) = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];,
				    			 x(1) = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];,
				    			 x(2) = geom[lev].ProbLo()[2] + ((amrex::Real)(k) + 0.5) * geom[lev].CellSize()[2];);
                }

                field(i,j,k) = L;
                if ( (x - x0)(0) > 0.5*width || (x - x0)(0) < -0.5*width)
                    if ( (x - x0)(1) < 0.5*thickness && (x-x0)(1) > -0.5*thickness)
                        field(i,j,k) = 0.0;
			});
		}
        a_field[lev]->FillBoundary();
	}
	
private:
    Set::Scalar thickness = 0.05;
    Set::Scalar width = 0.1;
    Set::Vector x0 = Set::Vector(AMREX_D_DECL(1.0,0.5,0.0));
    Set::Scalar L = 1.0;
public:
    static void Parse(DoubleNotch & value, IO::ParmParse & pp)
    {
        pp.query("thickness",value.thickness);
        pp.query("width",value.width);
        pp.queryarr("x0",value.x0);
        pp.query("L",value.L);
        //pp.queryarr("x0",value.x0);
        //pp.query("eps",value.eps);
        //if (pp.contains("A"))
        //{
        //    pp.queryarr("A",value.A);
        //}
        //else if (pp.contains("a"))
        //{
        //    Set::Vector a = Set::Vector::Ones();
        //    pp.queryarr("a",a);
        //    for (int d = 0; d < AMREX_SPACEDIM; d++) value.A(d,d) = 1./a(d)/a(d);
        //}
    }
};
}
#endif
