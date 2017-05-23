#include "AmrAdv_F.H"

void initdata(const int& level, const amrex::Real& time, 
	      const int* lo, const int* hi,
	      BL_FORT_FAB_ARG_3D(state),
	      const amrex::Real* dx, const amrex::Real* problo)
{
  std::cout << "WARNING - still calling Fortran version of: init data" << std::endl;
}

void get_face_velocity(const int& level, const amrex::Real& time, 
		       D_DECL(BL_FORT_FAB_ARG(xvel),
			      BL_FORT_FAB_ARG(yvel),
			      BL_FORT_FAB_ARG(zvel)),
		       const amrex::Real* dx, const amrex::Real* problo)
{
  std::cout << "WARNING - still calling Fortran version of: get face velocity" << std::endl;
}

void state_error(int* tag, const int* tag_lo, const int* tag_hi,
		 const BL_FORT_FAB_ARG_3D(state),
		 const int* tagval, const int* clearval,
		 const int* lo, const int* hi,
		 const amrex::Real* dx, const amrex::Real* problo,
		 const amrex::Real* time, const amrex::Real* phierr)
{
  std::cout << "WARNING - still calling Fortran version of: state error" << std::endl;

}

void advect(const amrex::Real& time, const int* lo, const int*hi,
	    const BL_FORT_FAB_ARG_3D(statein),
	    BL_FORT_FAB_ARG_3D(stateout),
	    D_DECL(const BL_FORT_FAB_ARG_3D(xvel),
		   const BL_FORT_FAB_ARG_3D(yvel),
		   const BL_FORT_FAB_ARG_3D(zvel)),
	    D_DECL(BL_FORT_FAB_ARG_3D(fx),
		   BL_FORT_FAB_ARG_3D(fy),
		   BL_FORT_FAB_ARG_3D(fz)),
	    const amrex::Real* dx, const amrex::Real& dt)
{
  std::cout << "WARNING - still calling Fortran version of: advect" << std::endl;

}
