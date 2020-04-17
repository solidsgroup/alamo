#ifndef IC_BMP_H
#define IC_BMP_H
#include <cmath>

#include "IC/IC.H"
#include "Util/Util.H"
#include "Util/BMP.H"

namespace IC
{
class BMP : public IC
{
public:
	enum Type {XYZ, XY, YZ, XZ};

	BMP (amrex::Vector<amrex::Geometry> &_geom) : IC(_geom) 
    {
        Define();
    }
	void Define()
	{
        bmp.Define("Interface_Blur2.bmp");
	}
  
	void Add(const int lev,
			amrex::Vector<amrex::MultiFab * > &a_field)
	{
        Set::Vector DX(geom[lev].CellSize());

        Set::Scalar width  = geom[lev].ProbHi()[0] - geom[lev].ProbLo()[0];
        Set::Scalar height = geom[lev].ProbHi()[1] - geom[lev].ProbLo()[1];

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

                int I = std::min((int)((Set::Scalar)bmp.width  * ((x(0)-geom[lev].ProbLo()[0]) / width)), bmp.width-1);
                int J = std::min((int)((Set::Scalar)bmp.height * ((x(1)-geom[lev].ProbLo()[1]) / height)), bmp.height-1);
                    
                //Util::Message(INFO,I, " ", J, " " , (Set::Scalar)(bmp(I,J)[1]));
				field(i,j,k) = ((Set::Scalar)(bmp(I,J)[1]) - 136.0) / (255.0 - 136.0);//alpha_out;     
			});
		}
        a_field[lev]->FillBoundary();
	};
  
private:
    Util::BMP bmp;
	int dim = AMREX_SPACEDIM;
	amrex::Vector<Set::Scalar> center;
	Set::Scalar radius,alpha_in,alpha_out;
	Type type;

public:
    static void Parse(BMP & value, IO::ParmParse & pp)
    {
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
