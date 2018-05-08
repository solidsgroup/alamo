#include <AMReX_MultiFabUtil.H>
#include "Operator.H"
#include "eigen3/Eigen/Core"
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>
#include <AMReX_ArrayLim.H>


Operator::Operator::Operator (const Vector<Geometry>& a_geom,
		    const Vector<BoxArray>& a_grids,
		    const Vector<DistributionMapping>& a_dmap,
		    const LPInfo& a_info)
{
  define(a_geom, a_grids, a_dmap, a_info);
}

void
Operator::Operator::define (const Vector<Geometry>& a_geom,
			   const Vector<BoxArray>& a_grids,
			   const Vector<DistributionMapping>& a_dmap,
			   const LPInfo& a_info)
{
  MLCellLinOp::define(a_geom, a_grids, a_dmap, a_info);
}

void
Operator::Operator::RegisterNewFab(amrex::Vector<std::unique_ptr<amrex::MultiFab> > &input,
				   BC::BC &new_bc)
{
  /// \todo assertions here
  m_a_coeffs.resize(m_a_coeffs.size() + 1);
  m_a_coeffs[m_num_a_fabs].resize(m_num_amr_levels);
  for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
	{
		m_a_coeffs[m_num_a_fabs][amrlev].resize(m_num_mg_levels[amrlev]);
		for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
			m_a_coeffs[m_num_a_fabs][amrlev][mglev].define(m_grids[amrlev][mglev],
				m_dmap[amrlev][mglev],
				input[amrlev].get()->nComp(),
				input[amrlev].get()->nGrow());

		MultiFab::Copy(m_a_coeffs[m_num_a_fabs][amrlev][0],
			       *input[amrlev].get(), 0, 0,
			       input[amrlev].get()->nComp(),
			       input[amrlev].get()->nGrow());
	}
  m_num_a_fabs++;

  physbc_array.push_back(&new_bc); 
}

const amrex::FArrayBox &
Operator::Operator::GetFab(const int num, const int amrlev, const int mglev, const amrex::MFIter &mfi) const
{
  return m_a_coeffs[num][amrlev][mglev][mfi];
}

void
Operator::Operator::averageDownCoeffs ()
{
	for (int i = 0; i < m_num_a_fabs; i++)
	{
		for (int amrlev = m_num_amr_levels-1; amrlev > 0; --amrlev)
		{
			auto& fine_a_coeffs = m_a_coeffs[i][amrlev];
			averageDownCoeffsSameAmrLevel(fine_a_coeffs);
		}
		averageDownCoeffsSameAmrLevel(m_a_coeffs[i][0]);

		for (int amrlev = 0; amrlev < m_num_amr_levels; amrlev++)
		{
			physbc_array[i]->SetLevel(amrlev);
			for (int mglev = 0 ; mglev < m_num_mg_levels[amrlev]; mglev++)
			{
				/// \todo The last three arguments of FillBoundary are currently unused,
				///       but we will need to modify them here if we ever use them
				physbc_array[i]->FillBoundary(m_a_coeffs[i][amrlev][mglev],0,0,0.0);
			}
		}
	}
}

void
Operator::Operator::averageDownCoeffsSameAmrLevel (Vector<MultiFab>& a)
{
	int nmglevs = a.size();
	for (int mglev = 1; mglev < nmglevs; ++mglev)
	{
		amrex::average_down(a[mglev-1], a[mglev], 0, a[0].nComp(), mg_coarsen_ratio);
	}
}

void
Operator::Operator::applyMetricTermsCoeffs ()
{
#if (AMREX_SPACEDIM != 3)
	for (int i = 0; i < m_num_a_fabs; i++)
	{
		for (int alev = 0; alev < m_num_amr_levels; ++alev)
		{
			const int mglev = 0;
			applyMetricTerm(alev, mglev, m_a_coeffs[i][alev][mglev]);
		}
	}
#endif
}

Operator::Operator::~Operator ()
{}

/// \fn prepareForSolve
///
/// Relay function and distribute coefficients to all MG levels (?)
void
Operator::Operator::prepareForSolve ()
{
  MLCellLinOp::prepareForSolve();
  applyMetricTermsCoeffs();
  averageDownCoeffs();
}
