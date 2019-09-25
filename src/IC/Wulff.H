#ifndef IC_WULFF_H_
#define IC_WULFF_H_

#include "IC/IC.H"
#include "Util/Util.H"
#include "Model/Interface/GB/SH.H"

/// \class Wulff
/// \brief Initialize a spherical inclusion
namespace IC
{
class Wulff : public IC
{
public:
	Wulff (amrex::Vector<amrex::Geometry> &_geom) :
		IC(_geom)
	{
        model.Define(0,0,0.5,0.5);
	}
  
	void Add(const int lev,
			 amrex::Vector<amrex::MultiFab * > &a_field)
	{
        a_field[lev]->setVal(1.0);
		bool cellcentered = (a_field[0]->boxArray().ixType() == amrex::IndexType(amrex::IntVect::TheCellVector()));

		for (amrex::MFIter mfi(*a_field[lev],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
		{
			amrex::Box bx = mfi.tilebox();
			bx.grow(a_field[lev]->nGrow());

            //a_field[lev]->setVal(1.0,1);

			amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
            
            for (Set::Scalar theta = 0.0; theta <= 0.5*Set::Constant::Pi; theta += 0.01)
            {
                for (Set::Scalar phi = 0.0; phi <= 0.5*Set::Constant::Pi; phi += 0.01)
                {
                    Set::Vector n(AMREX_D_DECL(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)));
                    Set::Scalar W = model.W(n);
                    //Util::Message(INFO,"phi=",phi," theta=",theta," W=",W);
            
        			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        Set::Vector x;                        
    				    AMREX_D_TERM(x(0) = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];,
    								 x(1) = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];,
    								 x(2) = geom[lev].ProbLo()[2] + ((amrex::Real)(k) + 0.5) * geom[lev].CellSize()[2];);

                        if (x.dot(n) > W) field(i,j,k) = 0.0;
				    });

                }                        
            }
		}

	};
  
private:
    Model::Interface::GB::SH model;
};
}
#endif
