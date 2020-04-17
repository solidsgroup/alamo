#ifndef IC_BMP_H
#define IC_BMP_H
#include <cmath>

#include "IC/IC.H"
#include "Util/Util.H"
#include "Util/BMP.H"
#include "Set/Set.H"

namespace IC
{
class BMP : public IC
{
public:
	//enum Type {XYZ, XY, YZ, XZ};
    enum Fit {Stretch,FitWidth};

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
        amrex::Box domain = geom[lev].Domain();

        Set::Scalar width  = geom[lev].ProbHi()[0] - geom[lev].ProbLo()[0];
        Set::Scalar height = geom[lev].ProbHi()[1] - geom[lev].ProbLo()[1];

        Set::Scalar img_width = (Set::Scalar)(bmp.nx-1);
        Set::Scalar img_height = (Set::Scalar)(bmp.ny-1);
        Set::Scalar img_dx = 1.0;
        Set::Scalar img_dy = 1.0;


		for (amrex::MFIter mfi(*a_field[lev],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
		{
			amrex::Box bx = mfi.tilebox();
		    bx.grow(a_field[lev]->nGrow());
			amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {

				Set::Scalar x, y;

                x = (Set::Scalar)i / (Set::Scalar)(domain.hiVect()[0]);
                y = (Set::Scalar)j / (Set::Scalar)(domain.hiVect()[1]);

                
                Set::Scalar img_x = (img_width)  * x;
                Set::Scalar img_y = (img_height) * y;

                //Util::Message(INFO,"imgx = ", img_x);
                //Util::Message(INFO,"imgy = ", img_y);

                int I = (int)(img_x);
                int J = (int)(img_y);

                Set::Scalar x1 = I*img_dx, x2 = (I+1)*img_dx;
                Set::Scalar y1 = J*img_dy, y2 = (J+1)*img_dy;

                if (I < bmp.nx-1 && J < bmp.ny - 1)
                {
                    Set::Scalar fQ11 = ((Set::Scalar)(bmp(I,  J  )[1]) - 136.0) / (255.0 - 136.0);
                    Set::Scalar fQ12 = ((Set::Scalar)(bmp(I,  J+1)[1]) - 136.0) / (255.0 - 136.0);
                    Set::Scalar fQ21 = ((Set::Scalar)(bmp(I+1,J  )[1]) - 136.0) / (255.0 - 136.0);
                    Set::Scalar fQ22 = ((Set::Scalar)(bmp(I+1,J+1)[1]) - 136.0) / (255.0 - 136.0);

                    field(i,j,k) =
                        (fQ11*(x2-img_x)*(y2-img_y) + fQ21*(img_x-x1)*(y2-img_y) + fQ12*(x2-img_x)*(img_y-y1) + fQ22*(img_x-x1)*(img_y-y1)) / (img_dx * img_dy);
                }
                else if (I == bmp.nx-1 && J <  bmp.ny - 1)
                {
                    Set::Scalar fQ11 = ((Set::Scalar)(bmp(I  ,J  )[1]) - 136.0) / (255.0 - 136.0);
                    Set::Scalar fQ12 = ((Set::Scalar)(bmp(I  ,J+1)[1]) - 136.0) / (255.0 - 136.0);
                    field(i,j,k) = fQ11 + (fQ12-fQ11) * (img_y - y1);
                }
                else if (I <  bmp.nx-1 && J == bmp.ny - 1)
                {
                    Set::Scalar fQ11 = ((Set::Scalar)(bmp(I  ,J  )[1]) - 136.0) / (255.0 - 136.0);
                    Set::Scalar fQ21 = ((Set::Scalar)(bmp(I+1,J  )[1]) - 136.0) / (255.0 - 136.0);
                    field(i,j,k) = fQ11 + (fQ21-fQ11) * (img_x - x1);
                }
                else if (I ==  bmp.nx-1 && J == bmp.ny - 1)
                {
                    Set::Scalar fQ11 = ((Set::Scalar)(bmp(I  ,J  )[1]) - 136.0) / (255.0 - 136.0);
                    field(i,j,k) = fQ11;
                }
                else
                {
                    field(i,j,k) = 0.0;
                }

			});
		}
        a_field[lev]->FillBoundary();
	};
  
private:
    Util::BMP bmp;
    Fit fit = Fit::Stretch;

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
