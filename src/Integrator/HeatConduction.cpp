#include "HeatConduction.H"
#include "BC/Constant.H"
#include "IC/Constant.H"

namespace Integrator
{
/// \fn    HeatConduction::Integrator::Integrator
///
/// Read in the following simulation parameters
///
///     heat.alpha                (default 1.0)
///     heat.refinement_threshold (default 0.01)
///     ic.type
///

/// Initialize initial condition pointer #ic, and register
/// the #Temp, #Temp_old Multifab arrays.

HeatConduction::HeatConduction() :
	Integrator()
{
	amrex::ParmParse pp("heat");
	pp.query("alpha", alpha);
	pp.query("refinement_threshold", refinement_threshold);
	pp.query("ic_type", ic_type);

	// // Determine initial condition
	if (ic_type == "cylinder")
	{
		amrex::Real Tin, Tout;
		pp.query("Tin", Tin);
    	pp.query("Tout", Tout);
		IC = new IC::Cylinder(geom,Tin,Tout);
	}
	else if (ic_type == "constant")
		Util::Abort(INFO,"Need to fix this");
		//ic = new IC::Constant(geom,alpha);
	else
		Util::Abort(INFO,"Need to fix this");
		//IC = new IC::Sphere(radius,center,XY)
    
	{
		amrex::ParmParse pp("bc");
		amrex::Vector<std::string> bc_hi_str(AMREX_SPACEDIM);
		amrex::Vector<std::string> bc_lo_str(AMREX_SPACEDIM);
		pp.queryarr("lo",bc_lo_str,0,BL_SPACEDIM);
		pp.queryarr("hi",bc_hi_str,0,BL_SPACEDIM);
		amrex::Vector<amrex::Real> bc_lo_1, bc_hi_1;
		if (pp.countval("lo_1")) pp.getarr("lo_1",bc_lo_1);
		if (pp.countval("hi_1")) pp.getarr("hi_1",bc_hi_1);
		amrex::Vector<amrex::Real> bc_lo_2, bc_hi_2;
		if (pp.countval("lo_2")) pp.getarr("lo_2",bc_lo_2);
		if (pp.countval("hi_2")) pp.getarr("hi_2",bc_hi_2);
		amrex::Vector<amrex::Real> bc_lo_3, bc_hi_3;
		if (pp.countval("lo_3")) pp.getarr("lo_3",bc_lo_3);
		if (pp.countval("hi_3")) pp.getarr("hi_3",bc_hi_3);

		BC = new BC::Constant(bc_hi_str, bc_lo_str,
					AMREX_D_DECL(bc_lo_1, bc_lo_2, bc_lo_3),
					AMREX_D_DECL(bc_hi_1, bc_hi_2, bc_hi_3));
	}


	RegisterNewFab(TempFab,     BC, number_of_components, number_of_ghost_cells, "Temp",true);
	RegisterNewFab(TempOldFab, BC, number_of_components, number_of_ghost_cells, "Temp old",false);
}

/// \fn HeatConduction::Integrator::Initialize
///
/// Use the #ic object to initialize #Temp
void
HeatConduction::Initialize (int lev)
{
	IC->Initialize(lev,TempFab);
}


/// \fn    Integrator::HeatConduction::Advance
///
/// Integrate the heat diffusion equation
/// \f[\nabla^2T = \alpha \frac{\partial T}{\partial t}\f]
/// using an explicit forward Euler method.
/// \f$\alpha\f$ is stored in #alpha
void
HeatConduction::Advance (int lev, amrex::Real /*time*/, amrex::Real dt)
{
	static amrex::IntVect AMREX_D_DECL(dx(AMREX_D_DECL(1,0,0)),
												  dy(AMREX_D_DECL(0,1,0)),
												  dz(AMREX_D_DECL(0,0,1)));

	std::swap(*TempFab[lev], *TempOldFab[lev]);

	const amrex::Real* DX = geom[lev].CellSize();

	for ( amrex::MFIter mfi(*TempFab[lev],true); mfi.isValid(); ++mfi )
		{
			const amrex::Box& bx = mfi.tilebox();

			amrex::FArrayBox &TempOld = (*TempOldFab[lev])[mfi];
			amrex::FArrayBox &Temp    = (*TempFab[lev])[mfi];

			AMREX_D_TERM(for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++),
							 for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++),
							 for (int k = bx.loVect()[2]; k<=bx.hiVect()[2]; k++))
				{
					amrex::IntVect m(AMREX_D_DECL(i,j,k));
					Temp(m)
						= TempOld(m)
						+ dt * alpha * (AMREX_D_TERM((TempOld(m+dx) + TempOld(m-dx) - 2*TempOld(m)) / DX[0] / DX[0],
															  + (TempOld(m+dy) + TempOld(m-dy) - 2*TempOld(m)) / DX[1] / DX[1],
															  + (TempOld(m+dz) + TempOld(m-dz) - 2*TempOld(m)) / DX[2] / DX[2]));
				}
		}
}


/// \fn    Integrator::HeatConduction::TagCellsForRefinement
///
/// The following criterion is used to determine if a cell should be refined:
/// \f[|\nabla T|\,|\mathbf{r}| > h\f]
/// where
/// \f[\mathbf{r} = \sqrt{\Delta x_1^2 + \Delta x_2^2 + \Delta x_3^2}\f]
/// and \f$h\f$ is stored in #refinement_threshold
void
HeatConduction::TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/)
{
	const amrex::Real* DX      = geom[lev].CellSize();

	static amrex::IntVect AMREX_D_DECL(dx(AMREX_D_DECL(1,0,0)),
												  dy(AMREX_D_DECL(0,1,0)),
												  dz(AMREX_D_DECL(0,0,1)));

	for (amrex::MFIter mfi(*TempFab[lev],true); mfi.isValid(); ++mfi)
		{
			const amrex::Box&  bx  = mfi.tilebox();
			amrex::TagBox&     tag  = tags[mfi];
 	    
			amrex::FArrayBox &Temp = (*TempFab[lev])[mfi];

			AMREX_D_TERM(for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++),
							 for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++),
							 for (int k = bx.loVect()[2]; k<=bx.hiVect()[2]; k++))
						{
							amrex::IntVect m(AMREX_D_DECL(i,j,k));

							AMREX_D_TERM(amrex::Real grad1 = (Temp(m+dx) - Temp(m-dx));,
											 amrex::Real grad2 = (Temp(m+dy) - Temp(m-dy));,
											 amrex::Real grad3 = (Temp(m+dz) - Temp(m-dz));)

							amrex::Real grad = sqrt(AMREX_D_TERM(grad1*grad1, + grad2*grad2, + grad3*grad3));

							amrex::Real dr = sqrt(AMREX_D_TERM(DX[0]*DX[0], + DX[1]*DX[1], + DX[2]*DX[2]));

							if (grad*dr > refinement_threshold)
								tag(amrex::IntVect(AMREX_D_DECL(i,j,k))) = amrex::TagBox::SET;
						}
		}
}

}
