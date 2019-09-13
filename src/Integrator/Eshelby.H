#ifndef INTEGRATOR_ESHELBY_H
#define INTEGRATOR_ESHELBY_H
#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_ParmParse.H"

#include "Integrator/Integrator.H"

#include "IC/IC.H"
#include "BC/BC.H"

#include "IC/HeatConduction/Cylinder.H"
#include "IC/HeatConduction/Constant.H"
#include "IC/Ellipse.H"
#include "Numeric/Stencil.H"

#include "Model/Solid/LinearElastic/MultiWell.H"
#include "Solver/Nonlocal/Linear.H"

namespace Integrator
{
class Eshelby : public Integrator
{
public:

  /// \brief Read in parameters and register field variables
  Eshelby()
  {
      RegisterNodalFab(eta_mf,1,2,"eta");
      RegisterNodalFab(disp_mf,AMREX_SPACEDIM,2,"disp");
      RegisterNodalFab(rhs_mf,AMREX_SPACEDIM,2,"rhs");
      RegisterNodalFab(stress_mf,AMREX_SPACEDIM*AMREX_SPACEDIM,2,"stress");
      Set::Matrix A = Set::Matrix::Identity();
      A(0,0) = 64.0; A(1,1) = 32.0; A(2,2) = 0.0;
      Set::Vector x0 = Set::Vector::Ones()*0.5;
      Set::Scalar eps = 0.01;
      ic = new IC::Ellipse(geom,A,x0,IC::Ellipse::Mollifier::Gaussian,eps);
      //ic = new IC::Ellipse(geom,A,x0,IC::Ellipse::Mollifier::Dirac,eps);
  }

protected:

  /// \brief Use the #ic object to initialize#Temp
  void Initialize (int lev) override
  {
      eta_mf[lev]->setVal(0.0);
      ic->Initialize(lev,eta_mf);
  }

  void TimeStepBegin(Set::Scalar, int) override 
  {
      Util::Message(INFO,"finest_level = ", finest_level);
      Util::Message(INFO,"rhs_mf.size() = ",rhs_mf.size());
      
    if (finest_level != rhs_mf.size()-1)
        Util::Abort(INFO,"amr.max_level is larger than necessary. Set to ", finest_level," or less");
   	for (int lev = 0; lev < rhs_mf.size(); lev++) rhs_mf[lev]->setVal(0.0);
	for (int lev = 0; lev < rhs_mf.size(); lev++) disp_mf[lev]->setVal(0.0);

    Set::Scalar lame = 2.6, shear = 6.0;
	using model_type = Model::Solid::LinearElastic::Multiwell;
	Operator::Elastic<model_type> elastic;
    elastic.SetUniform(false);
    elastic.define(geom,grids,dmap);

    // Set linear elastic model
   	model_type mymodel(lame,shear,Set::Matrix::Zero());
	amrex::Vector<amrex::FabArray<amrex::BaseFab<model_type> > > model_mf;
	model_mf.resize(disp_mf.size());
    Set::Matrix Fmatrix    = Set::Matrix::Zero();
    Set::Matrix Finclusion = Set::Matrix::Zero();
    Finclusion(0,0) = 0.1;
	for (int lev = 0; lev < rhs_mf.size(); ++lev)
	{
		amrex::Box domain(geom[lev].Domain());
		domain.convert(amrex::IntVect::TheNodeVector());
		model_mf[lev].define(disp_mf[lev]->boxArray(), disp_mf[lev]->DistributionMap(), 1, 2);
		model_mf[lev].setVal(mymodel);
		
        eta_mf[lev]->FillBoundary();

        Set::Vector DX(geom[lev].CellSize());

		for (MFIter mfi(model_mf[lev],amrex::TilingIfNotGPU());mfi.isValid();++mfi)
		{
			//amrex::Box bx = mfi.growntilebox(2);
			amrex::Box bx = mfi.tilebox();
			bx.grow(2);

			amrex::Array4<model_type> const & model = model_mf[lev].array(mfi);
			amrex::Array4<const Set::Scalar> const & eta = eta_mf[lev]->array(mfi);

			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) 
			{
                Set::Matrix Fgb = (1.0 - eta(i,j,k))*Fmatrix + eta(i,j,k)*Finclusion;
				model(i,j,k) = model_type(lame,shear,Fgb);
			});
			bx = mfi.tilebox();
			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) 
			{
                for (int p = 0; p < AMREX_SPACEDIM; p++)
                for (int q = 0; q < AMREX_SPACEDIM; q++)
                {
				    AMREX_D_TERM(model(i,j,k).gradFgb[p](q,0) = ((model(i+1,j,k).Fgb - model(i-1,j,k).Fgb)/2./DX(0))(p,q);,
				                 model(i,j,k).gradFgb[p](q,1) = ((model(i,j+1,k).Fgb - model(i,j-1,k).Fgb)/2./DX(1))(p,q);,
				                 model(i,j,k).gradFgb[p](q,2) = ((model(i,j,k+1).Fgb - model(i,j,k-1).Fgb)/2./DX(2))(p,q);)
                }
			});
		}
	}
    elastic.SetModel(model_mf);

    // Set boundary conditions
   	BC::Operator::Elastic<model_type> bc;
	elastic.SetBC(&bc);

   	Set::Scalar tol_rel = 1E-8, tol_abs = 1E-8;
	Solver::Nonlocal::Linear linearsolver(elastic);
    linearsolver.setVerbose(2);
    linearsolver.setFixedIter(10);
    linearsolver.solveaffine(disp_mf,rhs_mf,tol_rel,tol_abs,true);
    
    //elastic.SetHomogeneous(false);
    //linearsolver.apply(GetVecOfPtrs(rhs_mf),GetVecOfPtrs(disp_mf));
    
    for (int lev = 0; lev < disp_mf.size(); lev++)
    {
        elastic.Stress(lev,*stress_mf[lev],*disp_mf[lev]);
    }



  }

  /// \brief Integrate the heat equation
  void Advance (int /*lev*/, Set::Scalar /*time*/, Set::Scalar /*dt*/) override
  {
      // Nothing to do here.
  }

  /// \brief Tag cells for mesh refinement based on temperature gradient
  void TagCellsForRefinement (int lev, amrex::TagBoxArray& a_tags, amrex::Real /*time*/, int /*ngrow*/) override
  {
    Set::Vector DX(geom[lev].CellSize());
	Set::Scalar DXnorm = DX.lpNorm<2>();
    a_tags.setVal(amrex::TagBox::CLEAR);
	for (amrex::MFIter mfi(*eta_mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		amrex::Box bx = mfi.tilebox();
		amrex::Array4<char> const& tags        = a_tags.array(mfi);
		amrex::Array4<Set::Scalar> const & Eta = eta_mf[lev]->array(mfi);
		amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k)
		{
			Set::Vector grad = Numeric::Gradient(Eta,i,j,k,0,DX.data());
			if (grad.lpNorm<2>()*DXnorm > 0.01) tags(i,j,k) = amrex::TagBox::SET;
		});
	}
	

  }

private:
  int number_of_components = 1;					///< Number of components
  int number_of_ghost_cells = 1;				///< Number of ghost cells

  amrex::Vector<std::unique_ptr<amrex::MultiFab> > eta_mf;
  amrex::Vector<std::unique_ptr<amrex::MultiFab> > disp_mf;
  amrex::Vector<std::unique_ptr<amrex::MultiFab> > rhs_mf;
  amrex::Vector<std::unique_ptr<amrex::MultiFab> > stress_mf;

  amrex::Real alpha = 1.0;					///< Thermal diffusivity
  amrex::Real refinement_threshold = 0.01;			///< Criterion for cell refinement

  std::string ic_type;						///< String to determine what kind of initial condition to use
  IC::IC *ic;            					///< Pointer to abstract IC object
  BC::BC *mybc;                   				///< Stock generic boundary condition object
};
}
#endif
