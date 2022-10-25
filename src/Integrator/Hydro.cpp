#include "Hydro.H"
#include "IO/ParmParse.H"
#include "BC/Constant.H"
#include "Numeric/Stencil.H"
#include "IC/Laminate.H"
#include "IC/PSRead.H"
#include "IC/Expression.H"

namespace Integrator
{

    void Hydro::Parse(Hydro &Value, IO::ParmParse &pp)
    {
      BL_PROFILE("Integrator::Hydro::Hydro()");
      //General Variables Input Read:
      {
	pp.query("r_refinement_criterion", value.r_refinement_criterion);
	pp.query("e_refinement_criterion", value.e_refinement_criterion);
	pp.query("m_refinement_criterion", value.m_refinement_criterion);
	pp.query("gamma", value.gamma);
	pp.query("cfl", value.cfl);

      }
      // Register FabFields:
      {
	value.RegisterNewFab(value.eta_mf, value.bc_eta, 1, 2, "eta", true);
	value.RegisterNewFab(value.eta_old_mf, value.bc_eta, 1, 2, "eta_old", false);

	value.RegisterNewFab(value.rho_mf, value.bc_rho, 1, 2, "rho", true);
	value.RegisterNewFab(value.rho_old_mf, 1, 2, "rho_old", false);

	value.RegisterNewFab(value.E_mf, value.bc_E, 1, 2, "E", true);
	value.RegisterNewFab(value.E_old_mf, value.bc_E, 1, 2, "E_old", false);

	value.RegisterNewFab(value.Px_mf, value.bc_Px, 1, 2, "Px", true);
	value.RegisterNewFab(value.Py_mf, value.bc_Py, 1, 2, "Py", true);
        value.RegisterNewFab(value.Pz_mf, value.bc_Pz, 1, 2, "Pz", true);


	value.RegisterNewFab(value.Px_old_mf, value.bc_Px, 1, 2, "Px_old", false);
	value.RegisterNewFab(value.Py_old_mf, value.bc_Py, 1, 2, "Py_old", false);
	value.RegisterNewFab(value.Pz_old_mf, value.bc_Pz, 1, 2, "Pz_old", false);

      }


    }


    void Hydro::Initialize(int lev)
    {
      BL_PROFILE("Integrator::Hydro::Initialize");

      ic_eta -> Initialize(lev, eta_mf);
      ic_eta -> Initialize(lev, eta_old_mf);

      Energy_mf[lev] -> setVal(0.0);
      Energy_old_mf[lev] -> setVal(0.0);
      
      Density_mf -> setVal(0.0);
      Density_old_mf -> setVal(0.0);

      Momentum_mf -> setVal(0.0);
      Momentum_old_mf -> setVal(0.0);

      c_max = 0.0;
    }

    void Hydro::TimeStepBegin(Set::Scalar a_time, int a_iter)
    {
      BL_PROFILE("Integrator::Hydro::TimeStepBegin");
    }

  
    void Hydro::TimeStepEnd(Set::Scalar a_time, int a_iter)
    {
      // Syncronize c_max between processors so that they all have the same minimum value
      amrex::ParallelDescriptor::ReduceRealMax(c_max);

      // TODO calculate timestep?

      // This will set the timestep on the coarsest level
      //SetTimestep(new_timestep);
    }


    void Hydro::Advance(int lev, Set::Scalar time, Set::Scalar dt)
    {

      std::swap(eta_old_mf[lev], eta_mf[lev]);
      std::swap(Momentum_old_mf[lev], Momentum_mf[lev]);
      std::swap(Energy_old_mf[lev], Energy_mf[lev]);
      std::swap(Density_old_mf[lev], Density_mf[lev]);

      for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi)
      {
	const amrex::Box &bx = mfi.tilebox();

	amrex::Array4<Set::Scalar> const &eta = (*eta_mf[lev]).array(mfi);
	amrex::Array4<const Set::Scalar> const &etaold = (*eta_old_mf[lev]).array(mfi);
	amrex::Array4<Set::Scalar> const &E = (*Energy_mf[lev]).array(mfi);
	amrex::Array4<const Set::Scalar> const &Eold = (*Energy_old_mf[lev]).array(mfi);
	amrex::Array4<Set::Scalar> const &rho = (*Densitiy_mf[lev]).array(mfi);
	amrex::Array4<const Set::Scalar> const &rhoold = (*Density_old_mf[lev]).array(mfi);
	amrex::Array4<Set::Vector> const &M = (*Momentum_mf[lev]).array(mfi);
	amrex::Array4<const Set::Vector> const &Mold = (*Momentum_old_mf[lev]).array(mfi);
	amrex::Array4<Set::Vector> const &V = (*Velocity_mf[lev]).array(mfi);
	amrex::Array4<Set::Vector> const &p = (*Pressure_mf[lev]).array(mfi);
        //Computes Velocity and Pressure over the domain
 
	amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k){
	    V(i, j, k, 0) = M(i, j, k, 0) / rho(i,j,k);
	    V(i, j, k, 1) = M(i, j, k, 1) / rho(i,j,k);
	    V(i, j, k, 2) = M(i, j, k, 2) / rho(i,j,k);

	    Set::Scalar ke = V(i,j,k, 0);

	    p(i,j,k) = (gamma - 1.0) * rho(i,j,k) * (E(i,j,k) / rho(i,j,k) - 0.5 * ke * ke);

	    Set::Scalar c = gamma * p(i,j,k) / rho(i,j,k);

	    if (c > c_max){ c_max = c;}
 	});

	// Here you can add compute_dt

	
	//this loop will be running the godnov solver over the space
	amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
	{
	  Set::Vector dq;

	  //dq[0] = Numeric::Slopes();
	  //dq[1] = Numeric::Slopes();

	  // first derivative in x direction
	  Numeric::Stencil<Set::Scalar,1,0,0>::D(rho,i,j,k,0,dx);
	  // second derivative in x direction
	  Numeric::Stencil<Set::Scalar,2,0,0>::D(rho,i,j,k,0,dx);
	  // first derivative in y direction
	  Numeric::Stencil<Set::Scalar,0,1,0>::D(rho,i,j,k,0,dx);
 
	  

	});	
      }      
    }//end Advance

  void Hydro::Regrid(int lev, Set::Scalar /* time */)
  {
    

  }//end regrid

  void Hydro::TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, Set::Scalar time, int ngrow)
  {
    BL_PROFILE("Integrator::Hydro::TagCellsForRefinement");
    Base::Mechanics<Model::Solid::Affine::Isotropic>::TagCellsForRefinement(lev, a_tags, time, ngrow);

    const Set::Scalar *DX = geom[lev].CellSize();
    Set::Scalar dr = sqrt(AMREX_D_TERM(DX[0] * DX[0], +DX[1] * DX[1], +DX[2] * DX[2]));

    //Eta Criterion
    for (amrex::MFIter mfi(*eta_mf[lev], true), mfi.isValid(), ++mfi)
    {
      const amrex::Box &bx = mfi.tilebox();
      amrex::Array4<char> const &tags = a_tags.array(mfi);
      amrex::Array4<const Set::Scalar> const &eta = (*eta_mf[lev]).array(mfi);

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k){
	Set::Vector grad_eta = Numeric::Gradient(eta, i, j, k, 0, DX);
	if(grad_eta.lpnorm<2>() * dr * 2 > m_refinement_criterion) tags(i,j,k) = amrex::TagBox::SET;
     });
    }
    // E criterion
    for (amrex::MFIter mfi(*E_mf[lev], true), mfi.isValid(), ++mfi)
    {
      const amrex::Box &bx = mfi.tilebox();
      amrex::Array4<char> const &tags = a_tags.array(mfi);
      amrex::Array4<const Set::Scalar> const &E = (*E_mf[lev]).array(mfi);
      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k){
	Set::Vector grad_E = Numeric::Gradient(E, i, j, k, 0, DX);
	if(grad_E.lpNorm<2>() * dr > e_refinement_criterion) tags(i,j,k) = amrex::TagBox::SET;


	});
    }
    // rho criterion
    for (amrex::MFIter mfi(*rho_mf[lev] , true), mfi.isValid(), ++mfi)
    {
      const amrex::Box &bx = mfi.tilebox();
      amrex::Array4<char> const &tags = a_tags.array(mfi);
      amrex::Array4<const Set::Scalar> const &rho = (*rho_mf[lev]).array(mfi);
      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k){
	Set::Vector grad_rho = Numeric::Gradient(rho, i, j, k, 0, DX);
	if (grad_rho.lpNorm<2>() * dr > r_refinment_criterion) tags(i,j,k) = amrex::TagBox::SET;
      });
    }

    // P criterion


  }//end TagCells

  void Hydro::Integrate(int amrlev, Set::Scalar /*time*/, int /*step*/, const amrex::MFIter &mfi, const amrex::Box &box)
  {

  }//end Integrate

  void Hydro::UpdateModel(int /*a_step*/)
  {


  }//end update

  
}//end code






class hydroUtils(object):
    def __init__(self, param):
        self.gamma = param.gamma
        self.slope_type = param.slope_type
        riemann_solver_str = 'roe'
        if riemann_solver_str == 'diffuse':
            self.riemann_2d = self.riemann_diffuse
        elif riemann_solver_str == 'roe':
            self.riemann_2d = self.riemann_roe
        else:
            self.riemann_2d = self.riemann_hllc

    def cpg(self, rho, eint):
        p = (self.gamma - 1.0) * rho * eint
        c = np.sqrt(self.gamma * p / rho)
        return (p, c)

    def computePrimitives_ij(self, U, i, j):
        rho = U[(i, j, ID)]
        etot = U[(i, j, IP)]
        mom_x = U[(i, j, IU)]
        mom_y = U[(i, j, IV)]
        q = np.zeros(NBVAR)
        q[ID] = rho
        q[IU] = mom_x / q[ID]
        q[IV] = mom_y / q[ID]
        ke = 0.5 * (q[IU] * q[IU] + q[IV] * q[IV])
        e = etot / q[ID] - ke
        if e < 0:
            print('FATAL ERROR at {},{}: hydro eint < 0  : e {} ke {} d {} u {} v {}'.format(i, j, etot, ke, rho, mom_x, mom_y))
            return None
        q[IP], c = self.cpg(q[ID], e)
        return (q, c)

    #compute flux
    def get_dynamic_flux(self, qgdnv):
        ue = qgdnv[IP] / (self.gamma - 1.0)
        ke = 0.5 * qgdnv[ID] * (qgdnv[IU] * qgdnv[IU] + qgdnv[IV] * qgdnv[IV])
        etot = ue + ke
        ue = etot - ke

        fl = np.empty(NBVAR)
        fl[ID] = qgdnv[ID] * qgdnv[IU]
        fl[IU] = fl[ID] * qgdnv[IU] 
        fl[IV] = fl[ID] * qgdnv[IV]
        fl[IP] = qgdnv[IU] * ke
            
        return fl

    def get_static_flux(self, qgdnv):
        ue = qgdnv[IP] / (self.gamma - 1.0)
        ke = 0.5 * qgdnv[ID] * (qgdnv[IU] * qgdnv[IU] + qgdnv[IV] * qgdnv[IV])
        etot = ue + ke
        ue = etot - ke

        fl = np.empty(NBVAR)
        fl[ID] = 0.
        fl[IU] = qgdnv[IP]
        fl[IV] = 0.
        fl[IP] = qgdnv[IU] * (ue + qgdnv[IP])
            
        return fl

    def slope_unsplit_hydro_2d(self, q, qNeighbors):
        #compute primitive variables slope (vector dq) from q and its neighbors.
        qPlusX = qNeighbors[0]
        qMinusX = qNeighbors[1]
        qPlusY = qNeighbors[2]
        qMinusY = qNeighbors[3]
        dqX = np.zeros(NBVAR)
        dqY = np.zeros(NBVAR)
        if self.slope_type == 1 or self.slope_type == 2:
            dlft = self.slope_type * (q - qMinusX)
            drgt = self.slope_type * (qPlusX - q)
            dcen = 0.5 * (qPlusX - qMinusX)
            dsgn = np.sign(dcen)
            slop = np.minimum(abs(dlft), abs(drgt))
            dlim = slop
            dlim *= (np.sign(dlft * drgt) + 1) / 2
            dqX = dsgn * np.minimum(dlim, abs(dcen))
            dlft = self.slope_type * (q - qMinusY)
            drgt = self.slope_type * (qPlusY - q)
            dcen = 0.5 * (qPlusY - qMinusY)
            dsgn = np.sign(dcen)
            slop = np.minimum(abs(dlft), abs(drgt))
            dlim = slop
            dlim *= (np.sign(dlft * drgt) + 1) / 2
            dqY = dsgn * np.minimum(dlim, abs(dcen))
        return (dqX, dqY)

    def trace_unsplit_hydro_2d_by_direction(self, q, dqX, dqY, dtdx, dtdy, faceId):
        #trace computations for unsplit Godunov scheme
        qface = np.zeros(NBVAR)
        r = q[ID]
        p = q[IP]
        u = q[IU]
        v = q[IV]
        drx = dqX[ID]
        drx *= 0.5
        dpx = dqX[IP]
        dpx *= 0.5
        dux = dqX[IU]
        dux *= 0.5
        dvx = dqX[IV]
        dvx *= 0.5
        dry = dqY[ID]
        dry *= 0.5
        dpy = dqY[IP]
        dpy *= 0.5
        duy = dqY[IU]
        duy *= 0.5
        dvy = dqY[IV]
        dvy *= 0.5
        sr0 = (-u * drx - dux * r) * dtdx + (-v * dry - dvy * r) * dtdy
        su0 = (-u * dux - dpx / r) * dtdx + -v * duy * dtdy
        sv0 = -u * dvx * dtdx + (-v * dvy - dpy / r) * dtdy
        sp0 = (-u * dpx - dux * self.gamma * p) * dtdx + (-v * dpy - dvy * self.gamma * p) * dtdy
        r = r + sr0
        u = u + su0
        v = v + sv0
        p = p + sp0
        if faceId == FACE_XMIN:
            qface[ID] = r - drx
            qface[IU] = u - dux
            qface[IV] = v - dvx
            qface[IP] = p - dpx
        if faceId == FACE_XMAX:
            qface[ID] = r + drx
            qface[IU] = u + dux
            qface[IV] = v + dvx
            qface[IP] = p + dpx
        if faceId == FACE_YMIN:
            qface[ID] = r - dry
            qface[IU] = u - duy
            qface[IV] = v - dvy
            qface[IP] = p - dpy
        if faceId == FACE_YMAX:
            qface[ID] = r + dry
            qface[IU] = u + duy
            qface[IV] = v + dvy
            qface[IP] = p + dpy
        return qface

    ###################################
    ### RIEMANN SOLVER - ROE SOLVER ###
    ###################################

    def riemann_roe(self, qleft, qright):
        #Computation Gas Dynamics - Culbert B. Laney
        
        rho_l = qleft[ID]
        vel_l = qleft[IU]
        p_l = qleft[IP]
        ke_l = 0.5 * rho_l * vel_l**2
        ue_l = p_l / (self.gamma - 1.0)
        h_l = (ke_l + ue_l + p_l) / rho_l
        
        rho_r = qright[ID]
        vel_r = qright[IU]
        p_r = qright[IP]
        ke_r = 0.5 * rho_r * vel_r**2
        ue_r = p_r / (self.gamma - 1.0)
        h_r = (ke_r + ue_r + p_r) / rho_r
        
        #roe averages
        srl = np.sqrt(rho_l)
        srr = np.sqrt(rho_r)
        rho_RL = srl * srr
        vel_RL = (srr * vel_r + srl * vel_l) / (srl + srr)
        h_RL = (srr * h_r + srl * h_l) / (srl + srr)
        a_RL2 = (self.gamma - 1) * (h_RL - 0.5 * vel_RL ** 2)
        a_RL = np.sqrt(a_RL2)
        
        #compute eigenvalues
        ev = np.ndarray(shape = 3, dtype=float)
        ev[0] = vel_RL
        ev[1] = vel_RL + a_RL
        ev[2] = vel_RL - a_RL
        
        #compute wave strengths
        drho = rho_r - rho_l
        dp = p_r - p_l
        dvel = vel_r - vel_l

        dv = np.empty_like(ev)
        dv[0] = drho - dp / a_RL2
        dv[1] = dvel + dp / (rho_RL * a_RL)
        dv[2] = dvel - dp / (rho_RL * a_RL)
        
        #construct right characteristic eigenvectors
        rev = np.empty([3, 3])
        rev[0, :] = np.array([1, vel_RL, 0.5 * vel_RL ** 2])
        rev[1, :] = (rho_RL / 2 / a_RL) * np.array([1, vel_RL + a_RL, h_RL + a_RL * vel_RL])
        rev[2, :] = (-rho_RL / 2 / a_RL) * np.array([1, vel_RL - a_RL, h_RL - a_RL * vel_RL])

        #compute flux
        fl = 0.5 * (self.get_dynamic_flux(qleft) + self.get_dynamic_flux(qright)) + 0.5 * (self.get_static_flux(qleft) + self.get_static_flux(qright))

        #upwinding term
        for i in range(0, 3):
            upwind = -0.5 * rev[i, :] * np.abs(ev[i]) * dv[i]
            fl[ID] += upwind[0]
            fl[IU] += upwind[1]
            fl[IP] += upwind[2]
        
        return fl

########################################################
# `hydroRun` class
########################################################
class hydroRun(object):
    def __init__(self):

        self.param = hydroParams()

        self.utils = hydroUtils(self.param)

        self.U  = np.zeros((self.param.isize, self.param.jsize, NBVAR))
        self.U_old = np.zeros((self.param.isize, self.param.jsize, NBVAR))
        self.Q  = np.zeros((self.param.isize, self.param.jsize, NBVAR))

    ########################################################
    def init_condition(self):
        
        gw = self.param.ghostWidth

        density_in = 1.2
        density_out = 1.0
        pressure_in = 5.0
        pressure_out = 0.1

        epsilon = 16
        shift_factor = 1.4

        eta1 = np.empty([self.param.isize, self.param.jsize])
        eta2 = np.empty([self.param.isize, self.param.jsize])
        eta3 = np.empty([self.param.isize, self.param.jsize])
        eta4 = np.empty([self.param.isize, self.param.jsize])

        for j in range(0,self.param.jsize):
            for i in range(0,self.param.isize):
                r1 = np.sqrt((i-self.param.isize/2.+self.param.isize/(6*shift_factor))**2 + (j-self.param.jsize/2.+self.param.jsize/(6*shift_factor))**2)*6
                eta1[i,j] = np.tanh((r1+144)/epsilon)/2 - np.tanh((r1-144)/epsilon)/2

                r2 = np.sqrt((i-self.param.isize/2.+self.param.isize/(6*shift_factor))**2 + (j-self.param.jsize/2.-self.param.jsize/(6*shift_factor))**2)*6
                eta2[i,j] = np.tanh((r2+144)/epsilon)/2 - np.tanh((r2-144)/epsilon)/2

                r3 = np.sqrt((i-self.param.isize/2.-self.param.isize/(6*shift_factor))**2 + (j-self.param.jsize/2.+self.param.jsize/(6*shift_factor))**2)*6
                eta3[i,j] = np.tanh((r3+144)/epsilon)/2 - np.tanh((r3-144)/epsilon)/2
            
                r4 = np.sqrt((i-self.param.isize/2.-self.param.isize/(6*shift_factor))**2 + (j-self.param.jsize/2.-self.param.jsize/(6*shift_factor))**2)*6
                eta4[i,j] = np.tanh((r4+144)/epsilon)/2 - np.tanh((r4-144)/epsilon)/2
        
        eta = eta1 + eta2 + eta3 + eta4
        eta[eta>1] = 1.

        U = self.U

        U[:,:,ID] = density_in*eta + density_out*(1.-eta)
        U[:,:,IP] = (pressure_in*eta + pressure_out*(1.-eta))/(self.param.gamma - 1.0)
        U[:,:,IU] = 0.0
        U[:,:,IV] = 0.0

    ########################################################
    def compute_dt(self,useU):
        if useU == 0:
            U = self.U
        else:
            U = self.U_old

        invDt = 0
        dx = self.param.dx
        dy = self.param.dy

        for j in range(0,self.param.jsize):
            for i in range(0,self.param.isize):
                qLoc, c = self.utils.computePrimitives_ij(U,i,j)
                vx = c + abs(qLoc[IU])
                vy = c + abs(qLoc[IV])

                invDt = max(invDt, vx/dx + vy/dy);

        return self.param.cfl / invDt

    ########################################################
    
    def compute_primitives(self,U):
        for j in range(0,self.param.jsize):
            for i in range(0,self.param.isize):
                qLoc, c = self.utils.computePrimitives_ij(U,i,j)
                self.Q[i,j,:] = qLoc[:]
                
    ########################################################
    
    def make_boundaries(self, useU):
        if useU == 0:
            U = self.U
        else:
            U = self.U_old
        
        gw     = self.param.ghostWidth
        b_xmin = self.param.boundary_type_xmin 
        b_xmax = self.param.boundary_type_xmax
        b_ymin = self.param.boundary_type_ymin 
        b_ymax = self.param.boundary_type_ymax 

        nx = self.param.nx
        ny = self.param.ny

        imin = self.param.imin
        imax = self.param.imax
        jmin = self.param.jmin
        jmax = self.param.jmax

        #boundary xmin
        for iVar in range(NBVAR):
            for i in range(gw):
                sign = 1.0
                if   b_xmin == BC_DIRICHLET:
                    i0 = 2*gw-1-i
                    if iVar==IU:
                        sign = -1.0
                elif b_xmin == BC_NEUMANN:
                    i0 = gw
                else: # periodic
                    i0 = nx+i

                for j in range(jmin+gw, jmax-gw+1):
                    U[i,j,iVar] = U[i0,j,iVar]*sign


        #boundary xmax
        for iVar in range(NBVAR):
            for i in range (nx+gw, nx+2*gw):
                sign = 1.0
                if b_xmax == BC_DIRICHLET:
                    i0 = 2*nx + 2*gw-1-i
                    if iVar==IU:
                        sign = -1.0
                elif b_xmax == BC_NEUMANN:
                    i0 = nx+gw-1
                else:  # periodic
                    i0 = i-nx
                
                for j in range(jmin+gw, jmax-gw+1):
                    U[i,j,iVar] = U[i0,j,iVar]*sign
  
        #boundary ymin
        for iVar in range(NBVAR):
            for j in range(gw):
                sign = 1.0
                if b_ymin == BC_DIRICHLET:
                    j0 = 2*gw-1-j
                    if iVar==IV:
                        sign = -1.0
                elif b_ymin == BC_NEUMANN:
                    j0 = gw
                else:  # periodic
                    j0 = ny+j
      
                for i in range(imin+gw, imax-gw+1):
                    U[i,j,iVar] =  U[i,j0,iVar]*sign
        
        #boundary ymax
        for iVar in range(NBVAR):
            for j in range(ny+gw, ny+2*gw):
                sign = 1.0
                if b_ymax == BC_DIRICHLET:
                    j0 = 2*ny+2*gw-1-j
                    if iVar==IV:
                        sign = -1.0
                elif b_ymax == BC_NEUMANN:
                    j0 = ny+gw-1
                else:  # periodic
                    j0 = j-ny
      
                for i in range(imin+gw, imax-gw+1):
                    U[i,j,iVar] = U[i,j0,iVar]*sign

    ########################################################
    
    def godunov_unsplit(self,nStep,dt):
        if nStep%2 == 0:
            U  = self.U
            U_old = self.U_old
        else:
            U  = self.U_old
            U_old = self.U

        self.godunov_unsplit_cpu(U , U_old, dt, nStep)

    ########################################################
    
    def godunov_unsplit_cpu(self, U, U_old, dt, nStep):
        dtdx  = dt / self.param.dx
        dtdy  = dt / self.param.dy
        isize = self.param.isize
        jsize = self.param.jsize
        gw    = self.param.ghostWidth

        #fill ghost cell in data_in
        self.make_boundaries(nStep%2)

        #copy U to U_old
        U_old[:,:,:] = U[:,:,:]

        #main computation
        #convert to primitive variables
        self.compute_primitives(U)

        for j in range(gw, jsize-gw+1):
            for i in range(gw, isize-gw+1):
	
                #primitive variables in neighborhood
                qLoc = np.zeros(NBVAR)
                qLocN = np.zeros(NBVAR)
                qNeighbors = np.zeros((2*TWO_D,NBVAR))

                #get slopes in current cell
                qLoc[:] = self.Q[i,j,:] 

                qNeighbors[0] = self.Q[i+1,j  ,:]
                qNeighbors[1] = self.Q[i-1,j  ,:]
                qNeighbors[2] = self.Q[i  ,j+1,:]
                qNeighbors[3] = self.Q[i  ,j-1,:]
	  
                #compute slopes in current cell
                dqX, dqY = self.utils.slope_unsplit_hydro_2d(qLoc, qNeighbors)
    
                ##################################
                # left interface along x direction
                ##################################

                #get primitive variables state vector in left neighbor along x
                qLocN[:] = self.Q[i-1,j,:] 

                qNeighbors[0] = self.Q[i  ,j  ,:]
                qNeighbors[1] = self.Q[i-2,j  ,:]
                qNeighbors[2] = self.Q[i-1,j+1,:]
                qNeighbors[3] = self.Q[i-1,j-1,:]

                #compute slopes in left neighbor along x
                dqX_n, dqY_n = self.utils.slope_unsplit_hydro_2d(qLocN, qNeighbors)

                #compute reconstructed states at left interface along x in current cell
                #left interface: right state
                qright = self.utils.trace_unsplit_hydro_2d_by_direction(qLoc, dqX, dqY, dtdx, dtdy, FACE_XMIN)

                #left interface: left state
                qleft = self.utils.trace_unsplit_hydro_2d_by_direction(qLocN, dqX_n, dqY_n, dtdx, dtdy, FACE_XMAX)

                flux_x = self.utils.riemann_2d(qleft, qright)

                ##################################
                # left interface along y direction
                ##################################

                #get primitive variables state vector in left neighbor along y
                qLocN[:] = self.Q[i,j-1,:] 

                qNeighbors[0] = self.Q[i+1,j-1,:]
                qNeighbors[1] = self.Q[i-1,j-1,:]
                qNeighbors[2] = self.Q[i  ,j  ,:]
                qNeighbors[3] = self.Q[i  ,j-2,:]
	  
                #compute slopes in current cell
                dqX_n, dqY_n = self.utils.slope_unsplit_hydro_2d(qLocN, qNeighbors)

                #compute reconstructed states at left interface along y in current cell
                #left interface: right state
                qright = self.utils.trace_unsplit_hydro_2d_by_direction(qLoc, dqX, dqY, dtdx, dtdy, FACE_YMIN)
                
                #left interface: left state
                qleft = self.utils.trace_unsplit_hydro_2d_by_direction(qLocN, dqX_n, dqY_n,dtdx, dtdy, FACE_YMAX)
                
                #IU, IV permutations
                qleft[IU], qleft[IV] = qleft[IV], qleft[IU] 
                qright[IU], qright[IV]  = qright[IV], qright[IU] 
                
                flux_y = self.utils.riemann_2d(qleft, qright)
                
                #swap flux_y components
                flux_y[IU], flux_y[IV] = flux_y[IV], flux_y[IU] 
                
                #update hydro array
                U_old[i-1,j  ,:] += (-flux_x[:]*dtdx)
                U_old[i  ,j  ,:] += ( flux_x[:]*dtdx)
                
                U_old[i  ,j-1,:] += (-flux_y[:]*dtdy)
                U_old[i  ,j  ,:] += ( flux_y[:]*dtdy)

########################################################

def output(U, filename):
    nx = 140
    ny = 140
    gw = 2
    
    xaxis = np.linspace(0, 1, nx + 2*gw)
    yaxis = np.linspace(0, 1, ny + 2*gw)
    
    rho_vec = U[:,:,ID].ravel()
    rho = rho_vec.reshape((nx + 2*gw, ny + 2*gw))
    #print(np.shape(rho))
    
    mx_vec = U[:,:,IU].ravel()
    mx = mx_vec.reshape((nx + 2*gw, ny + 2*gw))
    
    my_vec = U[:,:,IV].ravel()
    my = my_vec.reshape((nx + 2*gw, ny + 2*gw))
    
    E_vec = U[:,:,IP].ravel()
    E = E_vec.reshape((nx + 2*gw, ny + 2*gw))

    p = (E - np.sqrt(mx**2 + my**2)/rho)/(0.4)
    
    np.savetxt(filename, p)
    
# #########################
# MAIN
# #########################
def main():
    t  = 0.0
    dt = 0.0
    nStep = 0

    par = hydroParams()

    par.printSetUp()

    #create run
    hr = hydroRun()

    #initial condition
    hr.init_condition()

    #initialize time step
    dt = hr.compute_dt(0)

    #initialize boundaries
    hr.make_boundaries(0)
    hr.make_boundaries(1)

    #start one integration
    print("Start computation....")
    
    t     = 0.0
    nStep = 0
    while t < hr.param.tEnd and nStep < hr.param.nStepmax:
        #output
        if nStep % hr.param.nOutput == 0:
            print("time t={0:16.13f} step {1:05d} dt={2:13.10f}".format(t, nStep, dt))
            filename = "U_{0:03d}.dat".format(nStep)
            if nStep % 2 == 0:
                output(hr.U,  filename)
            else:
                output(hr.U_old, filename)

        #compute new dt
        dt =  hr.compute_dt(nStep%2)
    
        #perform one step integration
        hr.godunov_unsplit(nStep, dt)

        #increase time
        nStep += 1
        t+=dt
       
    print("time : {0:5.3f} seconds".format(t))
    print("Perf : {0:10.2f} number of cell-updates/s".format(nStep*par.isize*par.jsize/t_tot))
  

# ################################
# EULER2D MAIN
# ################################
if __name__ == '__main__':
    main()
