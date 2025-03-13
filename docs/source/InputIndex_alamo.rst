alamo
--------------------------


.. raw:: html

      <table class='api-inputs-table'>
      <thead><tr>
      <th>Input name</th>
      <th>Description</th>
      </tr></thead>
       <tr>
         <td style='padding-left: 0px'>
           alamo.program
         </td>
         <td>
           The integrator to select 
         </td>
       </tr>
       <tr>
         <td style='padding-left: 0px'>
           alamo.program.microstructure.model
         </td>
         <td>
           which model to use (if using PFM with elasticity) 
         </td>
       </tr>
      </table><br/><br/>
       <h3> Integrator::PhaseFieldMicrostructure&lt;Model::Solid::Affine::Cubic&gt;  </h3>
      <table class='api-inputs-table'>
      <thead><tr>
      <th>Input name</th>
      <th>Description</th>
      </tr></thead>
         <tr>
           <td style='padding-left: 10px'>
             pf.number_of_grains
           </td>
           <td>
             Number of grain fields (may be more if using different IC) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.M
           </td>
           <td>
             Mobility
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.gamma
           </td>
           <td>
             pp_query_required("pf.mu", value.pf.mu);                       // Phase field :math:`\mu` Phase field :math:`\gamma`
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.sigma0
           </td>
           <td>
             Initial GB energy if not using  anisotropy
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.l_gb
           </td>
           <td>
             Mobility
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.elastic_df
           </td>
           <td>
             Determine whether to use elastic driving force
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.elastic_mult
           </td>
           <td>
             Multiplier of elastic energy
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.threshold.value
           </td>
           <td>
             Value used for thresholding kinetic relation
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.threshold.chempot
           </td>
           <td>
             Whether to include chemical potential in threshold
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.threshold.boundary
           </td>
           <td>
             Whether to include boundary energy in threshold
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.threshold.corner
           </td>
           <td>
             Whether to include corner regularization in threshold
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.threshold.lagrange
           </td>
           <td>
             Whether to include lagrange multiplier in threshold
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.threshold.mechanics
           </td>
           <td>
             Whether to include mechanical driving force in threshold
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.threshold.type
           </td>
           <td>
             Type of thresholding to use
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             amr.max_level
           </td>
           <td>
             Maximum AMR level
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             amr.ref_threshold
           </td>
           <td>
             Phase field refinement threshold
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.tstart
           </td>
           <td>
             Elasticity 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.mix_order
           </td>
           <td>
             Mixing order  
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.model_neuman_boundary
           </td>
           <td>
             Force Neumann BCs on the model 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.type
           </td>
           <td>
             Type of mecahnics to use. Static: do full implicit solve. Dynamic: evolve dynamic equations with explicit dynamics Disable: do nothing. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.time_evolving
           </td>
           <td>
             Treat mechanics fields as changing in time. [false] You should use this if you care about other physics driven by the output of this integrator. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.plot_disp
           </td>
           <td>
             Include displacement field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.plot_rhs
           </td>
           <td>
             Include right-hand side in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.plot_psi
           </td>
           <td>
             Include :math:`\psi` field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.plot_stress
           </td>
           <td>
             Include stress in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.plot_strain
           </td>
           <td>
             Include strain in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.viscous.mu_dashpot
           </td>
           <td>
             Dashpot damping (damps velocity)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.viscous.mu_newton
           </td>
           <td>
             Newtonian viscous damping (damps velocity gradient)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.velocity.ic.type
           </td>
           <td>
             Initializer for RHS
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              mechanics.bc.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              mechanics.bc.type
              = tensiontest
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.tensiontest.type
             </td>
             <td>
               Tension test type. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.tensiontest.disp
             </td>
             <td>
               Applied displacement (can be interpolator)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.tensiontest.trac
             </td>
             <td>
               Applied traction (can be interpolator)
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              mechanics.bc.type
              = expression
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.print_model
           </td>
           <td>
             Print out model variables (if enabled by model)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              mechanics.rhs.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.rhs.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              mechanics.rhs.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.rhs.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              mechanics.rhs.type
              = trig
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.rhs.trig.nr
             </td>
             <td>
               Number of real (cosin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.rhs.trig.ni
             </td>
             <td>
               Number of imaginary (sin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.rhs.trig.dim
             </td>
             <td>
               Spatial dimension
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.rhs.trig.alpha
             </td>
             <td>
               Multiplier
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.interval
           </td>
           <td>
             Timestep interval for elastic solves (default - solve every time) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.max_coarsening_level
           </td>
           <td>
             Maximum multigrid coarsening level (default - none, maximum coarsening) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.print_residual
           </td>
           <td>
             Whether to include residual output field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.elastic_ref_threshold
           </td>
           <td>
             Whether to refine based on elastic solution 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.zero_out_displacement
           </td>
           <td>
             Set this to true to zero out the displacement before each solve. (This is a temporary fix - we need to figure out why this is needed.) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.tstart
           </td>
           <td>
             Time to start doing the elastic solve (by default, start immediately) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             lagrange.on
           </td>
           <td>
             Lagrange multiplier method for enforcing volumes 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             lagrange.lambda
           </td>
           <td>
             Lagrange multiplier value
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             lagrange.vol0
           </td>
           <td>
             Prescribed volume
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             lagrange.tstart
           </td>
           <td>
             Lagrange multipler start time
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             sdf.on
           </td>
           <td>
             synthetic driving force (SDF)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             sdf.val
           </td>
           <td>
             value of SDF for each grain
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             sdf.tstart
           </td>
           <td>
             time to begin applying SDF
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropy.on
           </td>
           <td>
             Turn on
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropy.beta
           </td>
           <td>
             Regularization para m 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropy.tstart
           </td>
           <td>
             Time to turn on anisotropy 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropy.timestep
           </td>
           <td>
             Modify timestep when turned on 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropy.plot_int
           </td>
           <td>
             Modify plot_int when turned on 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropy.plot_dt
           </td>
           <td>
             Modify plot_dt when turned on 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropy.thermo_int
           </td>
           <td>
             Modify thermo int when turned on 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropy.thermo_plot_int
           </td>
           <td>
             Modify thermo plot int when turned on 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropy.elastic_int
           </td>
           <td>
             Frequency of elastic calculation 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropy.regularization
           </td>
           <td>
             Type of regularization to use   
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              anisotropy.type
              = abssin
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.abssin.theta0
             </td>
             <td>
               Angle offset (degrees)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.abssin.sigma0
             </td>
             <td>
               Minimum energy
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.abssin.sigma1
             </td>
             <td>
               Energy multiplier
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              anisotropy.type
              = sin
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.sin.theta0
             </td>
             <td>
               Theta offset (degrees)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.sin.sigma0
             </td>
             <td>
               Minimum energy
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.sin.sigma1
             </td>
             <td>
               Energy multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.sin.n
             </td>
             <td>
               Frequency number (integer)
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              anisotropy.type
              = read
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.read.filename
             </td>
             <td>
               Filename containing GB data
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              anisotropy.type
              = sh
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.sh.theta0
             </td>
             <td>
               Theta offset (degrees)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.sh.phi0
             </td>
             <td>
               Phi offset (radians)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.sh.sigma0
             </td>
             <td>
               Minimum energy value
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.sh.sigma1
             </td>
             <td>
               Energy multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.sh.regularization
             </td>
             <td>
               Type of regularization to use: {wilhelm,k23}
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             fluctuation.on
           </td>
           <td>
             Thermal fluctuations 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             fluctuation.amp
           </td>
           <td>
             fluctuation amplitude
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             fluctuation.sd
           </td>
           <td>
             fluctuation stadard deviation
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             fluctuation.tstart
           </td>
           <td>
             time to start applying fluctuation
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.on
           </td>
           <td>
             Disconnection generation 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.tstart
           </td>
           <td>
             time to start applying disconnections 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.nucleation_energy
           </td>
           <td>
             nucleation energy 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.tau_vol
           </td>
           <td>
             characteristic time 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.temp
           </td>
           <td>
             temperature 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.box_size
           </td>
           <td>
             characteristic size 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.interval
           </td>
           <td>
             interval between generation events 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.epsilon
           </td>
           <td>
             regularization epsilon 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.disconnection.fixed.on
           </td>
           <td>
             whether to manually specify disconnection nucleation points 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.fixed.sitex
           </td>
           <td>
             array of x locations 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.fixed.sitey
           </td>
           <td>
             array of y locations 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.fixed.phases
           </td>
           <td>
             array of order parameter number 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.fixed.time
           </td>
           <td>
             time to appear 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.verbose
           </td>
           <td>
             verbosity 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             shearcouple.on
           </td>
           <td>
             Shear coupling matrices 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.eta.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.type.xlo
             </td>
             <td>
               BC type on the lower x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.type.xhi
             </td>
             <td>
               BC type on the upper x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.type.ylo
             </td>
             <td>
               BC type on the lower y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.type.yhi
             </td>
             <td>
               BC type on the upper y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.type.zlo
             </td>
             <td>
               BC type on the lower z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.type.zhi
             </td>
             <td>
               BC type on the upper z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.val.xlo
             </td>
             <td>
               BC value on the lower x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.val.xhi
             </td>
             <td>
               BC value on the upper x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.val.ylo
             </td>
             <td>
               BC value on the lower y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.val.yhi
             </td>
             <td>
               BC value on the upper y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.val.zlo
             </td>
             <td>
               BC value on the lower z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.val.zhi
             </td>
             <td>
               BC value on the upper z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = perturbedinterface
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.perturbedinterface.wave_numbers
             </td>
             <td>
               Wave numbers
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.perturbedinterface.wave_amplitudes
             </td>
             <td>
               Wave amplitudes
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.perturbedinterface.normal
             </td>
             <td>
               Which axis is normal to the interface (x,y,z)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.perturbedinterface.offset
             </td>
             <td>
               Interface offset from origin
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.perturbedinterface.reverse
             </td>
             <td>
               If true, flip the interface (default:false)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.perturbedinterface.mollifier
             </td>
             <td>
               Mollifier (options: dirac, [gaussian])
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.perturbedinterface.eps
             </td>
             <td>
               Magnitude of mollifier
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = voronoi
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.voronoi.number_of_grains
             </td>
             <td>
               Number of grains
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.voronoi.alpha
             </td>
             <td>
               Value to take in the region [1.0]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.voronoi.seed
             </td>
             <td>
               Random seed to use
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = sphere
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.sphere.radius
             </td>
             <td>
               Radius of the sphere
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.sphere.center
             </td>
             <td>
               Vector location of the sphere center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.sphere.inside
             </td>
             <td>
               Value of the field inside the sphere
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.sphere.outside
             </td>
             <td>
               Value of the field outside teh sphere
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.sphere.type
             </td>
             <td>
               Type - can be cylinder oriented along the x, y, z directions or full sphere. 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = ellipse
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.x0
             </td>
             <td>
               Coorinates of ellipse center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.eps
             </td>
             <td>
               Diffuse boundary thickness
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.A
             </td>
             <td>
               DxD square matrix defining an ellipse. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.a
             </td>
             <td>
               If :code:`A` is not defined, then assume a sphere with radius :code:`a`
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.number_of_inclusions
             </td>
             <td>
               Number of ellipses
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.center
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.x0
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.A
             </td>
             <td>
               either a vector containing ellipse radii, or a matrix defining the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.A
             </td>
             <td>
               Same
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.radius
             </td>
             <td>
               Array of radii [depricated]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.eps
             </td>
             <td>
               Regularization for smooth boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.invert
             </td>
             <td>
               Flip the inside and the outside 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = random
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.random.offset
             </td>
             <td>
               offset from the [0,1] random number range 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.random.mult
             </td>
             <td>
               multiplier for the [0,1] random number range 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropic_kinetics.on
           </td>
           <td>
             Anisotropic mobility 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropic_kinetics.tstart
           </td>
           <td>
             simulation time when anisotropic kinetics is activated 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropic_kinetics.mobility
           </td>
           <td>
             file containing anisotropic mobility data 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropic_kinetics.threshold
           </td>
           <td>
             file containing anisotropic mobility data 
           </td>
         </tr>
      </table><br/><br/>
       <h3> Integrator::PhaseFieldMicrostructure&lt;Model::Solid::Affine::Hexagonal&gt;  </h3>
      <table class='api-inputs-table'>
      <thead><tr>
      <th>Input name</th>
      <th>Description</th>
      </tr></thead>
         <tr>
           <td style='padding-left: 10px'>
             pf.number_of_grains
           </td>
           <td>
             Number of grain fields (may be more if using different IC) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.M
           </td>
           <td>
             Mobility
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.gamma
           </td>
           <td>
             pp_query_required("pf.mu", value.pf.mu);                       // Phase field :math:`\mu` Phase field :math:`\gamma`
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.sigma0
           </td>
           <td>
             Initial GB energy if not using  anisotropy
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.l_gb
           </td>
           <td>
             Mobility
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.elastic_df
           </td>
           <td>
             Determine whether to use elastic driving force
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.elastic_mult
           </td>
           <td>
             Multiplier of elastic energy
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.threshold.value
           </td>
           <td>
             Value used for thresholding kinetic relation
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.threshold.chempot
           </td>
           <td>
             Whether to include chemical potential in threshold
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.threshold.boundary
           </td>
           <td>
             Whether to include boundary energy in threshold
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.threshold.corner
           </td>
           <td>
             Whether to include corner regularization in threshold
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.threshold.lagrange
           </td>
           <td>
             Whether to include lagrange multiplier in threshold
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.threshold.mechanics
           </td>
           <td>
             Whether to include mechanical driving force in threshold
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.threshold.type
           </td>
           <td>
             Type of thresholding to use
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             amr.max_level
           </td>
           <td>
             Maximum AMR level
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             amr.ref_threshold
           </td>
           <td>
             Phase field refinement threshold
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.tstart
           </td>
           <td>
             Elasticity 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.mix_order
           </td>
           <td>
             Mixing order  
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.model_neuman_boundary
           </td>
           <td>
             Force Neumann BCs on the model 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.type
           </td>
           <td>
             Type of mecahnics to use. Static: do full implicit solve. Dynamic: evolve dynamic equations with explicit dynamics Disable: do nothing. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.time_evolving
           </td>
           <td>
             Treat mechanics fields as changing in time. [false] You should use this if you care about other physics driven by the output of this integrator. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.plot_disp
           </td>
           <td>
             Include displacement field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.plot_rhs
           </td>
           <td>
             Include right-hand side in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.plot_psi
           </td>
           <td>
             Include :math:`\psi` field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.plot_stress
           </td>
           <td>
             Include stress in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.plot_strain
           </td>
           <td>
             Include strain in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.viscous.mu_dashpot
           </td>
           <td>
             Dashpot damping (damps velocity)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.viscous.mu_newton
           </td>
           <td>
             Newtonian viscous damping (damps velocity gradient)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.velocity.ic.type
           </td>
           <td>
             Initializer for RHS
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              mechanics.bc.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              mechanics.bc.type
              = tensiontest
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.tensiontest.type
             </td>
             <td>
               Tension test type. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.tensiontest.disp
             </td>
             <td>
               Applied displacement (can be interpolator)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.tensiontest.trac
             </td>
             <td>
               Applied traction (can be interpolator)
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              mechanics.bc.type
              = expression
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.print_model
           </td>
           <td>
             Print out model variables (if enabled by model)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              mechanics.rhs.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.rhs.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              mechanics.rhs.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.rhs.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              mechanics.rhs.type
              = trig
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.rhs.trig.nr
             </td>
             <td>
               Number of real (cosin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.rhs.trig.ni
             </td>
             <td>
               Number of imaginary (sin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.rhs.trig.dim
             </td>
             <td>
               Spatial dimension
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.rhs.trig.alpha
             </td>
             <td>
               Multiplier
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.interval
           </td>
           <td>
             Timestep interval for elastic solves (default - solve every time) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.max_coarsening_level
           </td>
           <td>
             Maximum multigrid coarsening level (default - none, maximum coarsening) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.print_residual
           </td>
           <td>
             Whether to include residual output field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.elastic_ref_threshold
           </td>
           <td>
             Whether to refine based on elastic solution 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.zero_out_displacement
           </td>
           <td>
             Set this to true to zero out the displacement before each solve. (This is a temporary fix - we need to figure out why this is needed.) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.tstart
           </td>
           <td>
             Time to start doing the elastic solve (by default, start immediately) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             lagrange.on
           </td>
           <td>
             Lagrange multiplier method for enforcing volumes 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             lagrange.lambda
           </td>
           <td>
             Lagrange multiplier value
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             lagrange.vol0
           </td>
           <td>
             Prescribed volume
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             lagrange.tstart
           </td>
           <td>
             Lagrange multipler start time
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             sdf.on
           </td>
           <td>
             synthetic driving force (SDF)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             sdf.val
           </td>
           <td>
             value of SDF for each grain
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             sdf.tstart
           </td>
           <td>
             time to begin applying SDF
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropy.on
           </td>
           <td>
             Turn on
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropy.beta
           </td>
           <td>
             Regularization para m 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropy.tstart
           </td>
           <td>
             Time to turn on anisotropy 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropy.timestep
           </td>
           <td>
             Modify timestep when turned on 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropy.plot_int
           </td>
           <td>
             Modify plot_int when turned on 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropy.plot_dt
           </td>
           <td>
             Modify plot_dt when turned on 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropy.thermo_int
           </td>
           <td>
             Modify thermo int when turned on 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropy.thermo_plot_int
           </td>
           <td>
             Modify thermo plot int when turned on 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropy.elastic_int
           </td>
           <td>
             Frequency of elastic calculation 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropy.regularization
           </td>
           <td>
             Type of regularization to use   
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              anisotropy.type
              = abssin
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.abssin.theta0
             </td>
             <td>
               Angle offset (degrees)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.abssin.sigma0
             </td>
             <td>
               Minimum energy
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.abssin.sigma1
             </td>
             <td>
               Energy multiplier
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              anisotropy.type
              = sin
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.sin.theta0
             </td>
             <td>
               Theta offset (degrees)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.sin.sigma0
             </td>
             <td>
               Minimum energy
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.sin.sigma1
             </td>
             <td>
               Energy multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.sin.n
             </td>
             <td>
               Frequency number (integer)
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              anisotropy.type
              = read
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.read.filename
             </td>
             <td>
               Filename containing GB data
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              anisotropy.type
              = sh
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.sh.theta0
             </td>
             <td>
               Theta offset (degrees)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.sh.phi0
             </td>
             <td>
               Phi offset (radians)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.sh.sigma0
             </td>
             <td>
               Minimum energy value
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.sh.sigma1
             </td>
             <td>
               Energy multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.sh.regularization
             </td>
             <td>
               Type of regularization to use: {wilhelm,k23}
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             fluctuation.on
           </td>
           <td>
             Thermal fluctuations 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             fluctuation.amp
           </td>
           <td>
             fluctuation amplitude
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             fluctuation.sd
           </td>
           <td>
             fluctuation stadard deviation
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             fluctuation.tstart
           </td>
           <td>
             time to start applying fluctuation
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.on
           </td>
           <td>
             Disconnection generation 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.tstart
           </td>
           <td>
             time to start applying disconnections 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.nucleation_energy
           </td>
           <td>
             nucleation energy 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.tau_vol
           </td>
           <td>
             characteristic time 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.temp
           </td>
           <td>
             temperature 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.box_size
           </td>
           <td>
             characteristic size 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.interval
           </td>
           <td>
             interval between generation events 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.epsilon
           </td>
           <td>
             regularization epsilon 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.disconnection.fixed.on
           </td>
           <td>
             whether to manually specify disconnection nucleation points 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.fixed.sitex
           </td>
           <td>
             array of x locations 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.fixed.sitey
           </td>
           <td>
             array of y locations 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.fixed.phases
           </td>
           <td>
             array of order parameter number 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.fixed.time
           </td>
           <td>
             time to appear 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.verbose
           </td>
           <td>
             verbosity 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             shearcouple.on
           </td>
           <td>
             Shear coupling matrices 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.eta.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.type.xlo
             </td>
             <td>
               BC type on the lower x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.type.xhi
             </td>
             <td>
               BC type on the upper x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.type.ylo
             </td>
             <td>
               BC type on the lower y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.type.yhi
             </td>
             <td>
               BC type on the upper y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.type.zlo
             </td>
             <td>
               BC type on the lower z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.type.zhi
             </td>
             <td>
               BC type on the upper z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.val.xlo
             </td>
             <td>
               BC value on the lower x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.val.xhi
             </td>
             <td>
               BC value on the upper x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.val.ylo
             </td>
             <td>
               BC value on the lower y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.val.yhi
             </td>
             <td>
               BC value on the upper y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.val.zlo
             </td>
             <td>
               BC value on the lower z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.val.zhi
             </td>
             <td>
               BC value on the upper z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = perturbedinterface
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.perturbedinterface.wave_numbers
             </td>
             <td>
               Wave numbers
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.perturbedinterface.wave_amplitudes
             </td>
             <td>
               Wave amplitudes
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.perturbedinterface.normal
             </td>
             <td>
               Which axis is normal to the interface (x,y,z)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.perturbedinterface.offset
             </td>
             <td>
               Interface offset from origin
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.perturbedinterface.reverse
             </td>
             <td>
               If true, flip the interface (default:false)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.perturbedinterface.mollifier
             </td>
             <td>
               Mollifier (options: dirac, [gaussian])
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.perturbedinterface.eps
             </td>
             <td>
               Magnitude of mollifier
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = voronoi
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.voronoi.number_of_grains
             </td>
             <td>
               Number of grains
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.voronoi.alpha
             </td>
             <td>
               Value to take in the region [1.0]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.voronoi.seed
             </td>
             <td>
               Random seed to use
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = sphere
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.sphere.radius
             </td>
             <td>
               Radius of the sphere
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.sphere.center
             </td>
             <td>
               Vector location of the sphere center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.sphere.inside
             </td>
             <td>
               Value of the field inside the sphere
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.sphere.outside
             </td>
             <td>
               Value of the field outside teh sphere
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.sphere.type
             </td>
             <td>
               Type - can be cylinder oriented along the x, y, z directions or full sphere. 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = ellipse
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.x0
             </td>
             <td>
               Coorinates of ellipse center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.eps
             </td>
             <td>
               Diffuse boundary thickness
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.A
             </td>
             <td>
               DxD square matrix defining an ellipse. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.a
             </td>
             <td>
               If :code:`A` is not defined, then assume a sphere with radius :code:`a`
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.number_of_inclusions
             </td>
             <td>
               Number of ellipses
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.center
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.x0
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.A
             </td>
             <td>
               either a vector containing ellipse radii, or a matrix defining the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.A
             </td>
             <td>
               Same
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.radius
             </td>
             <td>
               Array of radii [depricated]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.eps
             </td>
             <td>
               Regularization for smooth boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.invert
             </td>
             <td>
               Flip the inside and the outside 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = random
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.random.offset
             </td>
             <td>
               offset from the [0,1] random number range 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.random.mult
             </td>
             <td>
               multiplier for the [0,1] random number range 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropic_kinetics.on
           </td>
           <td>
             Anisotropic mobility 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropic_kinetics.tstart
           </td>
           <td>
             simulation time when anisotropic kinetics is activated 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropic_kinetics.mobility
           </td>
           <td>
             file containing anisotropic mobility data 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropic_kinetics.threshold
           </td>
           <td>
             file containing anisotropic mobility data 
           </td>
         </tr>
      </table><br/><br/>
       <h3> Integrator::PhaseFieldMicrostructure&lt;Model::Solid::Finite::PseudoAffine::Cubic&gt;  </h3>
      <table class='api-inputs-table'>
      <thead><tr>
      <th>Input name</th>
      <th>Description</th>
      </tr></thead>
         <tr>
           <td style='padding-left: 10px'>
             pf.number_of_grains
           </td>
           <td>
             Number of grain fields (may be more if using different IC) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.M
           </td>
           <td>
             Mobility
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.gamma
           </td>
           <td>
             pp_query_required("pf.mu", value.pf.mu);                       // Phase field :math:`\mu` Phase field :math:`\gamma`
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.sigma0
           </td>
           <td>
             Initial GB energy if not using  anisotropy
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.l_gb
           </td>
           <td>
             Mobility
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.elastic_df
           </td>
           <td>
             Determine whether to use elastic driving force
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.elastic_mult
           </td>
           <td>
             Multiplier of elastic energy
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.threshold.value
           </td>
           <td>
             Value used for thresholding kinetic relation
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.threshold.chempot
           </td>
           <td>
             Whether to include chemical potential in threshold
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.threshold.boundary
           </td>
           <td>
             Whether to include boundary energy in threshold
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.threshold.corner
           </td>
           <td>
             Whether to include corner regularization in threshold
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.threshold.lagrange
           </td>
           <td>
             Whether to include lagrange multiplier in threshold
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.threshold.mechanics
           </td>
           <td>
             Whether to include mechanical driving force in threshold
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.threshold.type
           </td>
           <td>
             Type of thresholding to use
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             amr.max_level
           </td>
           <td>
             Maximum AMR level
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             amr.ref_threshold
           </td>
           <td>
             Phase field refinement threshold
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.tstart
           </td>
           <td>
             Elasticity 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.mix_order
           </td>
           <td>
             Mixing order  
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.model_neuman_boundary
           </td>
           <td>
             Force Neumann BCs on the model 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.type
           </td>
           <td>
             Type of mecahnics to use. Static: do full implicit solve. Dynamic: evolve dynamic equations with explicit dynamics Disable: do nothing. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.time_evolving
           </td>
           <td>
             Treat mechanics fields as changing in time. [false] You should use this if you care about other physics driven by the output of this integrator. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.plot_disp
           </td>
           <td>
             Include displacement field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.plot_rhs
           </td>
           <td>
             Include right-hand side in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.plot_psi
           </td>
           <td>
             Include :math:`\psi` field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.plot_stress
           </td>
           <td>
             Include stress in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.plot_strain
           </td>
           <td>
             Include strain in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.viscous.mu_dashpot
           </td>
           <td>
             Dashpot damping (damps velocity)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.viscous.mu_newton
           </td>
           <td>
             Newtonian viscous damping (damps velocity gradient)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.velocity.ic.type
           </td>
           <td>
             Initializer for RHS
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              mechanics.bc.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.type.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.constant.val.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              mechanics.bc.type
              = tensiontest
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.tensiontest.type
             </td>
             <td>
               Tension test type. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.tensiontest.disp
             </td>
             <td>
               Applied displacement (can be interpolator)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.bc.tensiontest.trac
             </td>
             <td>
               Applied traction (can be interpolator)
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              mechanics.bc.type
              = expression
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.print_model
           </td>
           <td>
             Print out model variables (if enabled by model)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              mechanics.rhs.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.rhs.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              mechanics.rhs.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.rhs.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              mechanics.rhs.type
              = trig
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.rhs.trig.nr
             </td>
             <td>
               Number of real (cosin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.rhs.trig.ni
             </td>
             <td>
               Number of imaginary (sin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.rhs.trig.dim
             </td>
             <td>
               Spatial dimension
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               mechanics.rhs.trig.alpha
             </td>
             <td>
               Multiplier
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.interval
           </td>
           <td>
             Timestep interval for elastic solves (default - solve every time) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.max_coarsening_level
           </td>
           <td>
             Maximum multigrid coarsening level (default - none, maximum coarsening) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.print_residual
           </td>
           <td>
             Whether to include residual output field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.elastic_ref_threshold
           </td>
           <td>
             Whether to refine based on elastic solution 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.zero_out_displacement
           </td>
           <td>
             Set this to true to zero out the displacement before each solve. (This is a temporary fix - we need to figure out why this is needed.) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mechanics.tstart
           </td>
           <td>
             Time to start doing the elastic solve (by default, start immediately) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             lagrange.on
           </td>
           <td>
             Lagrange multiplier method for enforcing volumes 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             lagrange.lambda
           </td>
           <td>
             Lagrange multiplier value
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             lagrange.vol0
           </td>
           <td>
             Prescribed volume
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             lagrange.tstart
           </td>
           <td>
             Lagrange multipler start time
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             sdf.on
           </td>
           <td>
             synthetic driving force (SDF)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             sdf.val
           </td>
           <td>
             value of SDF for each grain
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             sdf.tstart
           </td>
           <td>
             time to begin applying SDF
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropy.on
           </td>
           <td>
             Turn on
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropy.beta
           </td>
           <td>
             Regularization para m 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropy.tstart
           </td>
           <td>
             Time to turn on anisotropy 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropy.timestep
           </td>
           <td>
             Modify timestep when turned on 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropy.plot_int
           </td>
           <td>
             Modify plot_int when turned on 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropy.plot_dt
           </td>
           <td>
             Modify plot_dt when turned on 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropy.thermo_int
           </td>
           <td>
             Modify thermo int when turned on 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropy.thermo_plot_int
           </td>
           <td>
             Modify thermo plot int when turned on 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropy.elastic_int
           </td>
           <td>
             Frequency of elastic calculation 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropy.regularization
           </td>
           <td>
             Type of regularization to use   
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              anisotropy.type
              = abssin
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.abssin.theta0
             </td>
             <td>
               Angle offset (degrees)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.abssin.sigma0
             </td>
             <td>
               Minimum energy
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.abssin.sigma1
             </td>
             <td>
               Energy multiplier
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              anisotropy.type
              = sin
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.sin.theta0
             </td>
             <td>
               Theta offset (degrees)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.sin.sigma0
             </td>
             <td>
               Minimum energy
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.sin.sigma1
             </td>
             <td>
               Energy multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.sin.n
             </td>
             <td>
               Frequency number (integer)
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              anisotropy.type
              = read
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.read.filename
             </td>
             <td>
               Filename containing GB data
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              anisotropy.type
              = sh
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.sh.theta0
             </td>
             <td>
               Theta offset (degrees)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.sh.phi0
             </td>
             <td>
               Phi offset (radians)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.sh.sigma0
             </td>
             <td>
               Minimum energy value
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.sh.sigma1
             </td>
             <td>
               Energy multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               anisotropy.sh.regularization
             </td>
             <td>
               Type of regularization to use: {wilhelm,k23}
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             fluctuation.on
           </td>
           <td>
             Thermal fluctuations 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             fluctuation.amp
           </td>
           <td>
             fluctuation amplitude
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             fluctuation.sd
           </td>
           <td>
             fluctuation stadard deviation
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             fluctuation.tstart
           </td>
           <td>
             time to start applying fluctuation
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.on
           </td>
           <td>
             Disconnection generation 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.tstart
           </td>
           <td>
             time to start applying disconnections 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.nucleation_energy
           </td>
           <td>
             nucleation energy 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.tau_vol
           </td>
           <td>
             characteristic time 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.temp
           </td>
           <td>
             temperature 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.box_size
           </td>
           <td>
             characteristic size 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.interval
           </td>
           <td>
             interval between generation events 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.epsilon
           </td>
           <td>
             regularization epsilon 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.disconnection.fixed.on
           </td>
           <td>
             whether to manually specify disconnection nucleation points 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.fixed.sitex
           </td>
           <td>
             array of x locations 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.fixed.sitey
           </td>
           <td>
             array of y locations 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.fixed.phases
           </td>
           <td>
             array of order parameter number 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.fixed.time
           </td>
           <td>
             time to appear 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             disconnection.verbose
           </td>
           <td>
             verbosity 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             shearcouple.on
           </td>
           <td>
             Shear coupling matrices 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.eta.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.type.xlo
             </td>
             <td>
               BC type on the lower x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.type.xhi
             </td>
             <td>
               BC type on the upper x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.type.ylo
             </td>
             <td>
               BC type on the lower y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.type.yhi
             </td>
             <td>
               BC type on the upper y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.type.zlo
             </td>
             <td>
               BC type on the lower z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.type.zhi
             </td>
             <td>
               BC type on the upper z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.val.xlo
             </td>
             <td>
               BC value on the lower x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.val.xhi
             </td>
             <td>
               BC value on the upper x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.val.ylo
             </td>
             <td>
               BC value on the lower y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.val.yhi
             </td>
             <td>
               BC value on the upper y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.val.zlo
             </td>
             <td>
               BC value on the lower z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.eta.constant.val.zhi
             </td>
             <td>
               BC value on the upper z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = perturbedinterface
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.perturbedinterface.wave_numbers
             </td>
             <td>
               Wave numbers
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.perturbedinterface.wave_amplitudes
             </td>
             <td>
               Wave amplitudes
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.perturbedinterface.normal
             </td>
             <td>
               Which axis is normal to the interface (x,y,z)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.perturbedinterface.offset
             </td>
             <td>
               Interface offset from origin
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.perturbedinterface.reverse
             </td>
             <td>
               If true, flip the interface (default:false)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.perturbedinterface.mollifier
             </td>
             <td>
               Mollifier (options: dirac, [gaussian])
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.perturbedinterface.eps
             </td>
             <td>
               Magnitude of mollifier
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = voronoi
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.voronoi.number_of_grains
             </td>
             <td>
               Number of grains
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.voronoi.alpha
             </td>
             <td>
               Value to take in the region [1.0]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.voronoi.seed
             </td>
             <td>
               Random seed to use
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = sphere
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.sphere.radius
             </td>
             <td>
               Radius of the sphere
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.sphere.center
             </td>
             <td>
               Vector location of the sphere center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.sphere.inside
             </td>
             <td>
               Value of the field inside the sphere
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.sphere.outside
             </td>
             <td>
               Value of the field outside teh sphere
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.sphere.type
             </td>
             <td>
               Type - can be cylinder oriented along the x, y, z directions or full sphere. 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = ellipse
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.x0
             </td>
             <td>
               Coorinates of ellipse center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.eps
             </td>
             <td>
               Diffuse boundary thickness
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.A
             </td>
             <td>
               DxD square matrix defining an ellipse. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.a
             </td>
             <td>
               If :code:`A` is not defined, then assume a sphere with radius :code:`a`
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.number_of_inclusions
             </td>
             <td>
               Number of ellipses
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.center
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.x0
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.A
             </td>
             <td>
               either a vector containing ellipse radii, or a matrix defining the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.A
             </td>
             <td>
               Same
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.radius
             </td>
             <td>
               Array of radii [depricated]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.eps
             </td>
             <td>
               Regularization for smooth boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.ellipse.invert
             </td>
             <td>
               Flip the inside and the outside 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = random
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.random.offset
             </td>
             <td>
               offset from the [0,1] random number range 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.random.mult
             </td>
             <td>
               multiplier for the [0,1] random number range 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropic_kinetics.on
           </td>
           <td>
             Anisotropic mobility 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropic_kinetics.tstart
           </td>
           <td>
             simulation time when anisotropic kinetics is activated 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropic_kinetics.mobility
           </td>
           <td>
             file containing anisotropic mobility data 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             anisotropic_kinetics.threshold
           </td>
           <td>
             file containing anisotropic mobility data 
           </td>
         </tr>
      </table><br/><br/>
       <h3> Integrator::Flame  </h3>
      <table class='api-inputs-table'>
      <thead><tr>
      <th>Input name</th>
      <th>Description</th>
      </tr></thead>
         <tr>
           <td style='padding-left: 10px'>
             plot_field
           </td>
           <td>
             Whether to include extra fields (such as mdot, etc) in the plot output 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             timestep
           </td>
           <td>
             Simulation timestep
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.eps
           </td>
           <td>
             Burn width thickness
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.kappa
           </td>
           <td>
             Interface energy param
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.gamma
           </td>
           <td>
             Scaling factor for mobility
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.lambda
           </td>
           <td>
             Chemical potential multiplier
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.w1
           </td>
           <td>
             Unburned rest energy
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.w12
           </td>
           <td>
             Barrier energy
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pf.w0
           </td>
           <td>
             Burned rest energy
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             amr.ghost_cells
           </td>
           <td>
             number of ghost cells in all fields
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             geometry.x_len
           </td>
           <td>
             Domain x length
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             geometry.y_len
           </td>
           <td>
             Domain y length
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              pf.eta.bc.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.bc.constant.type.xlo
             </td>
             <td>
               BC type on the lower x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.bc.constant.type.xhi
             </td>
             <td>
               BC type on the upper x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.bc.constant.type.ylo
             </td>
             <td>
               BC type on the lower y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.bc.constant.type.yhi
             </td>
             <td>
               BC type on the upper y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.bc.constant.type.zlo
             </td>
             <td>
               BC type on the lower z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.bc.constant.type.zhi
             </td>
             <td>
               BC type on the upper z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.bc.constant.val.xlo
             </td>
             <td>
               BC value on the lower x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.bc.constant.val.xhi
             </td>
             <td>
               BC value on the upper x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.bc.constant.val.ylo
             </td>
             <td>
               BC value on the lower y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.bc.constant.val.yhi
             </td>
             <td>
               BC value on the upper y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.bc.constant.val.zlo
             </td>
             <td>
               BC value on the lower z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.bc.constant.val.zhi
             </td>
             <td>
               BC value on the upper z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              pf.eta.ic.type
              = laminate
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.ic.laminate.number_of_inclusions
             </td>
             <td>
               How many laminates (MUST be greater than or equal to 1). 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.ic.laminate.orientation
             </td>
             <td>
               Vector normal to the interface of the laminate 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.ic.laminate.eps
             </td>
             <td>
               Diffuse thickness 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.ic.laminate.mollifier
             </td>
             <td>
               Type of mollifer to use (options: dirac, [gaussian])
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.ic.laminate.singlefab
             </td>
             <td>
               Switch to mode where only one component is used. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.ic.laminate.invert
             </td>
             <td>
               Take the complement of the laminate 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              pf.eta.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              pf.eta.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              pf.eta.ic.type
              = bmp
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.ic.bmp.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.ic.bmp.fit
             </td>
             <td>
               How to position image in space 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.ic.bmp.coord.lo
             </td>
             <td>
               Location of lower-left corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.ic.bmp.coord.hi
             </td>
             <td>
               Location of upper-right corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.ic.bmp.channel
             </td>
             <td>
               Color channel to use 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.ic.bmp.min
             </td>
             <td>
               Scaling value - minimum
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.ic.bmp.max
             </td>
             <td>
               Scaling value - maximum
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              pf.eta.ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               pf.eta.ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             thermal.on
           </td>
           <td>
             IO::ParmParse pp("thermal"); Whether to use the Thermal Transport Model
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic.on
           </td>
           <td>
             Whether to use Neo-hookean Elastic model
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             thermal.bound
           </td>
           <td>
             System Initial Temperature
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic.traction
           </td>
           <td>
             Body force
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic.phirefinement
           </td>
           <td>
             Phi refinement criteria 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             thermal.rho_ap
           </td>
           <td>
             AP Density
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             thermal.rho_htpb
           </td>
           <td>
             HTPB Density
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             thermal.k_ap
           </td>
           <td>
             AP Thermal Conductivity
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             thermal.k_htpb
           </td>
           <td>
             HTPB Thermal Conductivity
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             thermal.cp_ap
           </td>
           <td>
             AP Specific Heat
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             thermal.cp_htpb
           </td>
           <td>
             HTPB Specific Heat
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             thermal.q0
           </td>
           <td>
             Baseline heat flux
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             thermal.m_ap
           </td>
           <td>
             AP Pre-exponential factor for Arrhenius Law
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             thermal.m_htpb
           </td>
           <td>
             HTPB Pre-exponential factor for Arrhenius Law
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             thermal.E_ap
           </td>
           <td>
             AP Activation Energy for Arrhenius Law
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             thermal.E_htpb
           </td>
           <td>
             HTPB Activation Energy for Arrhenius Law
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             thermal.hc
           </td>
           <td>
             Used to change heat flux units
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             thermal.massfraction
           </td>
           <td>
             Systen AP mass fraction
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             thermal.mlocal_ap
           </td>
           <td>
             AP mass flux reference value 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             thermal.mlocal_htpb
           </td>
           <td>
             HTPB mass flux reference value 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             thermal.mlocal_comb
           </td>
           <td>
             AP/HTPB mass flux reference value 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             thermal.T_fluid
           </td>
           <td>
             Temperature of the Standin Fluid 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             thermal.disperssion1
           </td>
           <td>
             K; dispersion variables are use to set the outter field properties for the void grain case.
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             thermal.disperssion2
           </td>
           <td>
             rho; dispersion variables are use to set the outter field properties for the void grain case.
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             thermal.disperssion3
           </td>
           <td>
             cp; dispersion variables are use to set the outter field properties for the void grain case.
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             thermal.modeling_ap
           </td>
           <td>
             Scaling factor for AP thermal conductivity (default = 1.0)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             thermal.modeling_htpb
           </td>
           <td>
             Scaling factor for HTPB thermal conductivity (default = 1.0)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              thermal.temp.bc.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               thermal.temp.bc.constant.type.xlo
             </td>
             <td>
               BC type on the lower x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               thermal.temp.bc.constant.type.xhi
             </td>
             <td>
               BC type on the upper x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               thermal.temp.bc.constant.type.ylo
             </td>
             <td>
               BC type on the lower y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               thermal.temp.bc.constant.type.yhi
             </td>
             <td>
               BC type on the upper y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               thermal.temp.bc.constant.type.zlo
             </td>
             <td>
               BC type on the lower z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               thermal.temp.bc.constant.type.zhi
             </td>
             <td>
               BC type on the upper z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               thermal.temp.bc.constant.val.xlo
             </td>
             <td>
               BC value on the lower x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               thermal.temp.bc.constant.val.xhi
             </td>
             <td>
               BC value on the upper x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               thermal.temp.bc.constant.val.ylo
             </td>
             <td>
               BC value on the lower y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               thermal.temp.bc.constant.val.yhi
             </td>
             <td>
               BC value on the upper y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               thermal.temp.bc.constant.val.zlo
             </td>
             <td>
               BC value on the lower z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               thermal.temp.bc.constant.val.zhi
             </td>
             <td>
               BC value on the upper z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              laser.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               laser.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              laser.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               laser.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              temp.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               temp.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              temp.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               temp.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              temp.ic.type
              = bmp
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               temp.ic.bmp.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               temp.ic.bmp.fit
             </td>
             <td>
               How to position image in space 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               temp.ic.bmp.coord.lo
             </td>
             <td>
               Location of lower-left corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               temp.ic.bmp.coord.hi
             </td>
             <td>
               Location of upper-right corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               temp.ic.bmp.channel
             </td>
             <td>
               Color channel to use 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               temp.ic.bmp.min
             </td>
             <td>
               Scaling value - minimum
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               temp.ic.bmp.max
             </td>
             <td>
               Scaling value - maximum
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              temp.ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               temp.ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               temp.ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               temp.ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               temp.ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               temp.ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               temp.ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               temp.ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             pressure.P
           </td>
           <td>
             Constant pressure value
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pressure.a1
           </td>
           <td>
             Surgate heat flux model paramater - AP
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pressure.a2
           </td>
           <td>
             Surgate heat flux model paramater - HTPB
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pressure.a3
           </td>
           <td>
             Surgate heat flux model paramater - Total
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pressure.b1
           </td>
           <td>
             Surgate heat flux model paramater - AP
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pressure.b2
           </td>
           <td>
             Surgate heat flux model paramater - HTPB
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pressure.b3
           </td>
           <td>
             Surgate heat flux model paramater - Total
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pressure.c1
           </td>
           <td>
             Surgate heat flux model paramater - Total
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pressure.mob_ap
           </td>
           <td>
             Whether to include pressure to the arrhenius law
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pressure.dependency
           </td>
           <td>
             Whether to use pressure to determined the reference Zeta 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pressure.h1
           </td>
           <td>
             Surgate heat flux model paramater - Homogenized
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pressure.h2
           </td>
           <td>
             Surgate heat flux model paramater - Homogenized
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pressure.r_ap
           </td>
           <td>
             AP power pressure law parameter (r*P^n)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pressure.r_htpb
           </td>
           <td>
             HTPB power pressure law parameter (r*P^n)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pressure.r_comb
           </td>
           <td>
             AP/HTPB power pressure law parameter (r*P^n)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pressure.n_ap
           </td>
           <td>
             AP power pressure law parameter (r*P^n)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pressure.n_htpb
           </td>
           <td>
             HTPB power pressure law parameter (r*P^n)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pressure.n_comb
           </td>
           <td>
             AP/HTPB power pressure law parameter (r*P^n)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             variable_pressure
           </td>
           <td>
             Whether to compute the pressure evolution
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             homogeneousSystem
           </td>
           <td>
             Whether to initialize Phi with homogenized properties
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             amr.refinement_criterion
           </td>
           <td>
             Refinement criterion for eta field   
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             amr.refinement_criterion_temp
           </td>
           <td>
             Refinement criterion for temperature field    
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             amr.refinament_restriction
           </td>
           <td>
             Eta value to restrict the refinament for the temperature field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             amr.phi_refinement_criterion
           </td>
           <td>
             Refinement criterion for phi field [infinity]
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             small
           </td>
           <td>
             Lowest value of Eta.
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             phi.ic.type
           </td>
           <td>
             IC type (psread, laminate, constant)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             phi.ic.psread.eps
           </td>
           <td>
             value.ic_phicell = new IC::PSRead(value.geom, pp, "phi.ic.psread"); AP/HTPB interface length
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             phi.zeta_0
           </td>
           <td>
             Reference interface length for heat integration
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             phi.ic.laminate.eps
           </td>
           <td>
             value.ic_phicell = new IC::Laminate(value.geom, pp, "phi.ic.laminate"); AP/HTPB interface length
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             phi.zeta_0
           </td>
           <td>
             Reference interface length for heat integration
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             phi.zeta_0
           </td>
           <td>
             value.ic_phicell = new IC::Expression(value.geom, pp, "phi.ic.expression"); Reference interface length for heat integration
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             phi.zeta
           </td>
           <td>
             AP/HTPB interface length
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             phi.zeta_0
           </td>
           <td>
             value.ic_phicell = new IC::Constant(value.geom, pp, "phi.ic.constant"); Reference interface length for heat integration
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             phi.zeta
           </td>
           <td>
             AP/HTPB interface length
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             phi.zeta_0
           </td>
           <td>
             value.ic_phicell = new IC::BMP(value.geom, pp, "phi.ic.bmp"); Reference interface length for heat integration
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             phi.zeta
           </td>
           <td>
             AP/HTPB interface length
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             phi.zeta_0
           </td>
           <td>
             Reference interface length for heat integration
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             phi.zeta
           </td>
           <td>
             AP/HTPB interface length
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic.type
           </td>
           <td>
             Type of mecahnics to use. Static: do full implicit solve. Dynamic: evolve dynamic equations with explicit dynamics Disable: do nothing. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic.time_evolving
           </td>
           <td>
             Treat mechanics fields as changing in time. [false] You should use this if you care about other physics driven by the output of this integrator. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic.plot_disp
           </td>
           <td>
             Include displacement field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic.plot_rhs
           </td>
           <td>
             Include right-hand side in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic.plot_psi
           </td>
           <td>
             Include :math:`\psi` field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic.plot_stress
           </td>
           <td>
             Include stress in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic.plot_strain
           </td>
           <td>
             Include strain in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic.viscous.mu_dashpot
           </td>
           <td>
             Dashpot damping (damps velocity)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic.viscous.mu_newton
           </td>
           <td>
             Newtonian viscous damping (damps velocity gradient)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic.velocity.ic.type
           </td>
           <td>
             Initializer for RHS
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              elastic.bc.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.type.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.type.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.type.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.type.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.type.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.type.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.type.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.type.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.type.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.type.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.type.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.type.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.type.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.type.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.type.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.type.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.type.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.type.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.type.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.type.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.type.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.type.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.type.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.type.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.type.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.type.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.val.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.val.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.val.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.val.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.val.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.val.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.val.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.val.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.val.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.val.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.val.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.val.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.val.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.val.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.val.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.val.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.val.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.val.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.val.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.val.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.val.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.val.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.val.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.val.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.val.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.constant.val.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              elastic.bc.type
              = tensiontest
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.tensiontest.type
             </td>
             <td>
               Tension test type. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.tensiontest.disp
             </td>
             <td>
               Applied displacement (can be interpolator)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.bc.tensiontest.trac
             </td>
             <td>
               Applied traction (can be interpolator)
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              elastic.bc.type
              = expression
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic.print_model
           </td>
           <td>
             Print out model variables (if enabled by model)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              elastic.rhs.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.rhs.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              elastic.rhs.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.rhs.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              elastic.rhs.type
              = trig
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.rhs.trig.nr
             </td>
             <td>
               Number of real (cosin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.rhs.trig.ni
             </td>
             <td>
               Number of imaginary (sin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.rhs.trig.dim
             </td>
             <td>
               Spatial dimension
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               elastic.rhs.trig.alpha
             </td>
             <td>
               Multiplier
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic.interval
           </td>
           <td>
             Timestep interval for elastic solves (default - solve every time) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic.max_coarsening_level
           </td>
           <td>
             Maximum multigrid coarsening level (default - none, maximum coarsening) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic.print_residual
           </td>
           <td>
             Whether to include residual output field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic.elastic_ref_threshold
           </td>
           <td>
             Whether to refine based on elastic solution 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic.zero_out_displacement
           </td>
           <td>
             Set this to true to zero out the displacement before each solve. (This is a temporary fix - we need to figure out why this is needed.) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic.tstart
           </td>
           <td>
             Time to start doing the elastic solve (by default, start immediately) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             Tref
           </td>
           <td>
             Initial temperature for thermal expansion computation
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_ap.shear
           </td>
           <td>
             Shear modulus
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_ap.kappa
           </td>
           <td>
             Bulk modulus
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_ap.mu
           </td>
           <td>
             Alternative input for shear modulus
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_ap.kappa
           </td>
           <td>
             Bulk modulus
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_ap.shear
           </td>
           <td>
             Shear modulus
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_ap.lame
           </td>
           <td>
             Lame parameter
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_ap.E
           </td>
           <td>
             Young's modulus
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_ap.nu
           </td>
           <td>
             Poisson's ratio
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_ap.F0
           </td>
           <td>
             Large-deformation eigendeformation (Identity = no deformation)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_ap.eps0
           </td>
           <td>
             Small-deformation eigendeformation (Zero = no deformation)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_htpb.shear
           </td>
           <td>
             Shear modulus
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_htpb.kappa
           </td>
           <td>
             Bulk modulus
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_htpb.mu
           </td>
           <td>
             Alternative input for shear modulus
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_htpb.kappa
           </td>
           <td>
             Bulk modulus
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_htpb.shear
           </td>
           <td>
             Shear modulus
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_htpb.lame
           </td>
           <td>
             Lame parameter
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_htpb.E
           </td>
           <td>
             Young's modulus
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_htpb.nu
           </td>
           <td>
             Poisson's ratio
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_htpb.F0
           </td>
           <td>
             Large-deformation eigendeformation (Identity = no deformation)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_htpb.eps0
           </td>
           <td>
             Small-deformation eigendeformation (Zero = no deformation)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             allow_unused
           </td>
           <td>
             Set this to true to allow unused inputs without error. (Not recommended.) 
           </td>
         </tr>
      </table><br/><br/>
       <h3> Integrator::HeatConduction  </h3>
      <table class='api-inputs-table'>
      <thead><tr>
      <th>Input name</th>
      <th>Description</th>
      </tr></thead>
         <tr>
           <td style='padding-left: 10px'>
             heat.alpha
           </td>
           <td>
             Diffusion coefficient :math:`\alpha`. *This is an example of a required input variable - - program will terminate unless it is provided.* 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             heat.refinement_threshold
           </td>
           <td>
             Criterion for mesh refinement. *This is an example of a default input variable. The default value is provided here, not in the  definition of the variable.* 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = sphere
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.sphere.radius
             </td>
             <td>
               Radius of the sphere
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.sphere.center
             </td>
             <td>
               Vector location of the sphere center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.sphere.inside
             </td>
             <td>
               Value of the field inside the sphere
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.sphere.outside
             </td>
             <td>
               Value of the field outside teh sphere
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.sphere.type
             </td>
             <td>
               Type - can be cylinder oriented along the x, y, z directions or full sphere. 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.temp.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.temp.constant.type.xlo
             </td>
             <td>
               BC type on the lower x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.temp.constant.type.xhi
             </td>
             <td>
               BC type on the upper x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.temp.constant.type.ylo
             </td>
             <td>
               BC type on the lower y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.temp.constant.type.yhi
             </td>
             <td>
               BC type on the upper y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.temp.constant.type.zlo
             </td>
             <td>
               BC type on the lower z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.temp.constant.type.zhi
             </td>
             <td>
               BC type on the upper z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.temp.constant.val.xlo
             </td>
             <td>
               BC value on the lower x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.temp.constant.val.xhi
             </td>
             <td>
               BC value on the upper x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.temp.constant.val.ylo
             </td>
             <td>
               BC value on the lower y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.temp.constant.val.yhi
             </td>
             <td>
               BC value on the upper y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.temp.constant.val.zlo
             </td>
             <td>
               BC value on the lower z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.temp.constant.val.zhi
             </td>
             <td>
               BC value on the upper z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.temp.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.temp.expression.type.xlo
             </td>
             <td>
               BC type on the lower x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.temp.expression.type.xhi
             </td>
             <td>
               BC type on the upper x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.temp.expression.type.ylo
             </td>
             <td>
               BC type on the lower y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.temp.expression.type.yhi
             </td>
             <td>
               BC type on the upper y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.temp.expression.type.zlo
             </td>
             <td>
               BC type on the lower z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.temp.expression.type.zhi
             </td>
             <td>
               BC type on the upper z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
      </table><br/><br/>
       <h3> Integrator::ThermoElastic  </h3>
      <table class='api-inputs-table'>
      <thead><tr>
      <th>Input name</th>
      <th>Description</th>
      </tr></thead>
         <tr>
           <td style='padding-left: 10px'>
             hc.heat.alpha
           </td>
           <td>
             Diffusion coefficient :math:`\alpha`. *This is an example of a required input variable - - program will terminate unless it is provided.* 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             hc.heat.refinement_threshold
           </td>
           <td>
             Criterion for mesh refinement. *This is an example of a default input variable. The default value is provided here, not in the  definition of the variable.* 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hc.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hc.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hc.ic.type
              = sphere
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hc.ic.sphere.radius
             </td>
             <td>
               Radius of the sphere
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hc.ic.sphere.center
             </td>
             <td>
               Vector location of the sphere center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hc.ic.sphere.inside
             </td>
             <td>
               Value of the field inside the sphere
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hc.ic.sphere.outside
             </td>
             <td>
               Value of the field outside teh sphere
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hc.ic.sphere.type
             </td>
             <td>
               Type - can be cylinder oriented along the x, y, z directions or full sphere. 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hc.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hc.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hc.bc.temp.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hc.bc.temp.constant.type.xlo
             </td>
             <td>
               BC type on the lower x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hc.bc.temp.constant.type.xhi
             </td>
             <td>
               BC type on the upper x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hc.bc.temp.constant.type.ylo
             </td>
             <td>
               BC type on the lower y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hc.bc.temp.constant.type.yhi
             </td>
             <td>
               BC type on the upper y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hc.bc.temp.constant.type.zlo
             </td>
             <td>
               BC type on the lower z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hc.bc.temp.constant.type.zhi
             </td>
             <td>
               BC type on the upper z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hc.bc.temp.constant.val.xlo
             </td>
             <td>
               BC value on the lower x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hc.bc.temp.constant.val.xhi
             </td>
             <td>
               BC value on the upper x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hc.bc.temp.constant.val.ylo
             </td>
             <td>
               BC value on the lower y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hc.bc.temp.constant.val.yhi
             </td>
             <td>
               BC value on the upper y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hc.bc.temp.constant.val.zlo
             </td>
             <td>
               BC value on the lower z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hc.bc.temp.constant.val.zhi
             </td>
             <td>
               BC value on the upper z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hc.bc.temp.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hc.bc.temp.expression.type.xlo
             </td>
             <td>
               BC type on the lower x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hc.bc.temp.expression.type.xhi
             </td>
             <td>
               BC type on the upper x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hc.bc.temp.expression.type.ylo
             </td>
             <td>
               BC type on the lower y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hc.bc.temp.expression.type.yhi
             </td>
             <td>
               BC type on the upper y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hc.bc.temp.expression.type.zlo
             </td>
             <td>
               BC type on the lower z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hc.bc.temp.expression.type.zhi
             </td>
             <td>
               BC type on the upper z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             el.type
           </td>
           <td>
             Type of mecahnics to use. Static: do full implicit solve. Dynamic: evolve dynamic equations with explicit dynamics Disable: do nothing. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             el.time_evolving
           </td>
           <td>
             Treat mechanics fields as changing in time. [false] You should use this if you care about other physics driven by the output of this integrator. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             el.plot_disp
           </td>
           <td>
             Include displacement field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             el.plot_rhs
           </td>
           <td>
             Include right-hand side in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             el.plot_psi
           </td>
           <td>
             Include :math:`\psi` field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             el.plot_stress
           </td>
           <td>
             Include stress in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             el.plot_strain
           </td>
           <td>
             Include strain in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             el.viscous.mu_dashpot
           </td>
           <td>
             Dashpot damping (damps velocity)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             el.viscous.mu_newton
           </td>
           <td>
             Newtonian viscous damping (damps velocity gradient)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             el.velocity.ic.type
           </td>
           <td>
             Initializer for RHS
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              el.bc.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.type.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.type.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.type.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.type.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.type.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.type.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.type.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.type.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.type.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.type.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.type.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.type.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.type.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.type.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.type.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.type.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.type.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.type.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.type.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.type.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.type.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.type.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.type.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.type.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.type.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.type.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.val.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.val.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.val.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.val.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.val.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.val.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.val.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.val.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.val.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.val.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.val.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.val.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.val.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.val.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.val.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.val.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.val.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.val.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.val.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.val.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.val.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.val.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.val.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.val.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.val.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.constant.val.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              el.bc.type
              = tensiontest
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.tensiontest.type
             </td>
             <td>
               Tension test type. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.tensiontest.disp
             </td>
             <td>
               Applied displacement (can be interpolator)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.bc.tensiontest.trac
             </td>
             <td>
               Applied traction (can be interpolator)
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              el.bc.type
              = expression
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             el.print_model
           </td>
           <td>
             Print out model variables (if enabled by model)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              el.rhs.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.rhs.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              el.rhs.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.rhs.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              el.rhs.type
              = trig
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.rhs.trig.nr
             </td>
             <td>
               Number of real (cosin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.rhs.trig.ni
             </td>
             <td>
               Number of imaginary (sin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.rhs.trig.dim
             </td>
             <td>
               Spatial dimension
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.rhs.trig.alpha
             </td>
             <td>
               Multiplier
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             el.interval
           </td>
           <td>
             Timestep interval for elastic solves (default - solve every time) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             el.max_coarsening_level
           </td>
           <td>
             Maximum multigrid coarsening level (default - none, maximum coarsening) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             el.print_residual
           </td>
           <td>
             Whether to include residual output field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             el.elastic_ref_threshold
           </td>
           <td>
             Whether to refine based on elastic solution 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             el.zero_out_displacement
           </td>
           <td>
             Set this to true to zero out the displacement before each solve. (This is a temporary fix - we need to figure out why this is needed.) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             el.tstart
           </td>
           <td>
             Time to start doing the elastic solve (by default, start immediately) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             el.nmodels
           </td>
           <td>
             Number of elastic model varieties
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             el.eta_ref_threshold
           </td>
           <td>
             Refinement threshold for eta field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             el.ref_threshold
           </td>
           <td>
             Refinement threshold for strain gradient 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             el.model_neumann_boundary
           </td>
           <td>
             Explicity impose neumann condition on model at domain boundaries (2d only) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              el.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              el.ic.type
              = ellipse
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.ellipse.x0
             </td>
             <td>
               Coorinates of ellipse center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.ellipse.eps
             </td>
             <td>
               Diffuse boundary thickness
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.ellipse.A
             </td>
             <td>
               DxD square matrix defining an ellipse. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.ellipse.a
             </td>
             <td>
               If :code:`A` is not defined, then assume a sphere with radius :code:`a`
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.ellipse.number_of_inclusions
             </td>
             <td>
               Number of ellipses
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.ellipse.center
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.ellipse.x0
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.ellipse.A
             </td>
             <td>
               either a vector containing ellipse radii, or a matrix defining the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.ellipse.A
             </td>
             <td>
               Same
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.ellipse.radius
             </td>
             <td>
               Array of radii [depricated]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.ellipse.eps
             </td>
             <td>
               Regularization for smooth boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.ellipse.invert
             </td>
             <td>
               Flip the inside and the outside 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              el.ic.type
              = voronoi
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.voronoi.number_of_grains
             </td>
             <td>
               Number of grains
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.voronoi.alpha
             </td>
             <td>
               Value to take in the region [1.0]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.voronoi.seed
             </td>
             <td>
               Random seed to use
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              el.ic.type
              = bmp
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.bmp.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.bmp.fit
             </td>
             <td>
               How to position image in space 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.bmp.coord.lo
             </td>
             <td>
               Location of lower-left corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.bmp.coord.hi
             </td>
             <td>
               Location of upper-right corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.bmp.channel
             </td>
             <td>
               Color channel to use 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.bmp.min
             </td>
             <td>
               Scaling value - minimum
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.bmp.max
             </td>
             <td>
               Scaling value - maximum
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              el.ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              el.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              el.ic.type
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             el.eta.reset_on_regrid
           </td>
           <td>
             Whether to re-initialize eta when re-gridding occurs. Default is false unless eta ic is set, then default is. true. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              el.psi.ic.type
              = ellipse
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.psi.ic.ellipse.x0
             </td>
             <td>
               Coorinates of ellipse center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.psi.ic.ellipse.eps
             </td>
             <td>
               Diffuse boundary thickness
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.psi.ic.ellipse.A
             </td>
             <td>
               DxD square matrix defining an ellipse. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.psi.ic.ellipse.a
             </td>
             <td>
               If :code:`A` is not defined, then assume a sphere with radius :code:`a`
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.psi.ic.ellipse.number_of_inclusions
             </td>
             <td>
               Number of ellipses
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.psi.ic.ellipse.center
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.psi.ic.ellipse.x0
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.psi.ic.ellipse.A
             </td>
             <td>
               either a vector containing ellipse radii, or a matrix defining the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.psi.ic.ellipse.A
             </td>
             <td>
               Same
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.psi.ic.ellipse.radius
             </td>
             <td>
               Array of radii [depricated]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.psi.ic.ellipse.eps
             </td>
             <td>
               Regularization for smooth boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.psi.ic.ellipse.invert
             </td>
             <td>
               Flip the inside and the outside 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              el.psi.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.psi.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              el.psi.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.psi.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              el.psi.ic.type
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.psi.ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.psi.ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.psi.ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.psi.ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.psi.ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.psi.ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              el.psi.ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.psi.ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.psi.ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.psi.ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.psi.ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.psi.ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.psi.ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.psi.ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             el.psi.reset_on_regrid
           </td>
           <td>
             Whether to re-initialize psi when re-gridding occurs. Default is false unless a psi ic is set, then default is true. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              el.trac_normal.ic.type
              = ellipse
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.trac_normal.ic.ellipse.x0
             </td>
             <td>
               Coorinates of ellipse center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.trac_normal.ic.ellipse.eps
             </td>
             <td>
               Diffuse boundary thickness
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.trac_normal.ic.ellipse.A
             </td>
             <td>
               DxD square matrix defining an ellipse. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.trac_normal.ic.ellipse.a
             </td>
             <td>
               If :code:`A` is not defined, then assume a sphere with radius :code:`a`
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.trac_normal.ic.ellipse.number_of_inclusions
             </td>
             <td>
               Number of ellipses
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.trac_normal.ic.ellipse.center
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.trac_normal.ic.ellipse.x0
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.trac_normal.ic.ellipse.A
             </td>
             <td>
               either a vector containing ellipse radii, or a matrix defining the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.trac_normal.ic.ellipse.A
             </td>
             <td>
               Same
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.trac_normal.ic.ellipse.radius
             </td>
             <td>
               Array of radii [depricated]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.trac_normal.ic.ellipse.eps
             </td>
             <td>
               Regularization for smooth boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.trac_normal.ic.ellipse.invert
             </td>
             <td>
               Flip the inside and the outside 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              el.trac_normal.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.trac_normal.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              el.trac_normal.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.trac_normal.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              el.trac_normal.ic.type
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.trac_normal.ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.trac_normal.ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.trac_normal.ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.trac_normal.ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.trac_normal.ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               el.trac_normal.ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             alpha
           </td>
           <td>
             Diffusion coefficient
           </td>
         </tr>
      </table><br/><br/>
       <h3> Integrator::Fracture  </h3>
      <table class='api-inputs-table'>
      <thead><tr>
      <th>Input name</th>
      <th>Description</th>
      </tr></thead>
         <tr>
           <td style='padding-left: 10px'>
             type
           </td>
           <td>
             Type of crack {notch,ellipsoid}
           </td>
         </tr>
      </table><br/><br/>
       <h3> Integrator::Dendrite  </h3>
      <table class='api-inputs-table'>
      <thead><tr>
      <th>Input name</th>
      <th>Description</th>
      </tr></thead>
         <tr>
           <td style='padding-left: 10px'>
             alpha
           </td>
           <td>
             Pre-multiplier of "m" barrier height
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             delta
           </td>
           <td>
             Anisotropy factor
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             gamma
           </td>
           <td>
             Anisotropic temperature coupling factor
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             diffusion
           </td>
           <td>
             Thermal constant
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             eps
           </td>
           <td>
             Diffuse boundary width
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             tau
           </td>
           <td>
             Diffusive timescale
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             theta
           </td>
           <td>
             Orientation about z axis (Deg)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             heat.refinement_threshold
           </td>
           <td>
             Refinement criteria for temperature 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             phi.refinement_threshold
           </td>
           <td>
             Refinement criteria for phi 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.temp.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.temp.constant.type.xlo
             </td>
             <td>
               BC type on the lower x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.temp.constant.type.xhi
             </td>
             <td>
               BC type on the upper x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.temp.constant.type.ylo
             </td>
             <td>
               BC type on the lower y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.temp.constant.type.yhi
             </td>
             <td>
               BC type on the upper y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.temp.constant.type.zlo
             </td>
             <td>
               BC type on the lower z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.temp.constant.type.zhi
             </td>
             <td>
               BC type on the upper z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.temp.constant.val.xlo
             </td>
             <td>
               BC value on the lower x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.temp.constant.val.xhi
             </td>
             <td>
               BC value on the upper x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.temp.constant.val.ylo
             </td>
             <td>
               BC value on the lower y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.temp.constant.val.yhi
             </td>
             <td>
               BC value on the upper y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.temp.constant.val.zlo
             </td>
             <td>
               BC value on the lower z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.temp.constant.val.zhi
             </td>
             <td>
               BC value on the upper z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.phi.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.phi.constant.type.xlo
             </td>
             <td>
               BC type on the lower x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.phi.constant.type.xhi
             </td>
             <td>
               BC type on the upper x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.phi.constant.type.ylo
             </td>
             <td>
               BC type on the lower y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.phi.constant.type.yhi
             </td>
             <td>
               BC type on the upper y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.phi.constant.type.zlo
             </td>
             <td>
               BC type on the lower z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.phi.constant.type.zhi
             </td>
             <td>
               BC type on the upper z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.phi.constant.val.xlo
             </td>
             <td>
               BC value on the lower x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.phi.constant.val.xhi
             </td>
             <td>
               BC value on the upper x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.phi.constant.val.ylo
             </td>
             <td>
               BC value on the lower y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.phi.constant.val.yhi
             </td>
             <td>
               BC value on the upper y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.phi.constant.val.zlo
             </td>
             <td>
               BC value on the lower z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.phi.constant.val.zhi
             </td>
             <td>
               BC value on the upper z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
      </table><br/><br/>
       <h3> Integrator::AllenCahn  </h3>
      <table class='api-inputs-table'>
      <thead><tr>
      <th>Input name</th>
      <th>Description</th>
      </tr></thead>
         <tr>
           <td style='padding-left: 10px'>
             refinement_threshold
           </td>
           <td>
             Criterion for mesh refinement [0.01] 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             ch.L
           </td>
           <td>
             Value for :math:`L` (mobility) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             ch.eps
           </td>
           <td>
             Value for :math:`\epsilon` (diffuse boundary width) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             ch.grad
           </td>
           <td>
             Value for :math:`\kappa` (Interface energy parameter) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             ch.chempot
           </td>
           <td>
             Value for :math:`\lambda` (Chemical potential coefficient) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             ch.direction
           </td>
           <td>
             Force directional growth: 0=no growth, 1=only positive, -1=only negative 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             ch.direction_tstart
           </td>
           <td>
             Time to start forcing directional growth 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              alpha.ic.type
              = sphere
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.ic.sphere.radius
             </td>
             <td>
               Radius of the sphere
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.ic.sphere.center
             </td>
             <td>
               Vector location of the sphere center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.ic.sphere.inside
             </td>
             <td>
               Value of the field inside the sphere
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.ic.sphere.outside
             </td>
             <td>
               Value of the field outside teh sphere
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.ic.sphere.type
             </td>
             <td>
               Type - can be cylinder oriented along the x, y, z directions or full sphere. 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              alpha.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              alpha.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              alpha.ic.type
              = bmp
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.ic.bmp.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.ic.bmp.fit
             </td>
             <td>
               How to position image in space 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.ic.bmp.coord.lo
             </td>
             <td>
               Location of lower-left corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.ic.bmp.coord.hi
             </td>
             <td>
               Location of upper-right corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.ic.bmp.channel
             </td>
             <td>
               Color channel to use 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.ic.bmp.min
             </td>
             <td>
               Scaling value - minimum
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.ic.bmp.max
             </td>
             <td>
               Scaling value - maximum
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              alpha.ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              alpha.ic.type
              = random
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.ic.random.offset
             </td>
             <td>
               offset from the [0,1] random number range 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.ic.random.mult
             </td>
             <td>
               multiplier for the [0,1] random number range 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              alpha.ic.type
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              alpha.bc.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.bc.constant.type.xlo
             </td>
             <td>
               BC type on the lower x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.bc.constant.type.xhi
             </td>
             <td>
               BC type on the upper x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.bc.constant.type.ylo
             </td>
             <td>
               BC type on the lower y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.bc.constant.type.yhi
             </td>
             <td>
               BC type on the upper y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.bc.constant.type.zlo
             </td>
             <td>
               BC type on the lower z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.bc.constant.type.zhi
             </td>
             <td>
               BC type on the upper z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.bc.constant.val.xlo
             </td>
             <td>
               BC value on the lower x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.bc.constant.val.xhi
             </td>
             <td>
               BC value on the upper x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.bc.constant.val.ylo
             </td>
             <td>
               BC value on the lower y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.bc.constant.val.yhi
             </td>
             <td>
               BC value on the upper y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.bc.constant.val.zlo
             </td>
             <td>
               BC value on the lower z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               alpha.bc.constant.val.zhi
             </td>
             <td>
               BC value on the upper z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
      </table><br/><br/>
       <h3> Integrator::CahnHilliard  </h3>
      <table class='api-inputs-table'>
      <thead><tr>
      <th>Input name</th>
      <th>Description</th>
      </tr></thead>
         <tr>
           <td style='padding-left: 10px'>
             gamma
           </td>
           <td>
             Interface energy 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             L
           </td>
           <td>
             Mobility 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             refinement_threshold
           </td>
           <td>
             Regridding criterion 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              eta.ic.type
              = random
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.ic.random.offset
             </td>
             <td>
               offset from the [0,1] random number range 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.ic.random.mult
             </td>
             <td>
               multiplier for the [0,1] random number range 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              eta.bc.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.bc.constant.type.xlo
             </td>
             <td>
               BC type on the lower x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.bc.constant.type.xhi
             </td>
             <td>
               BC type on the upper x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.bc.constant.type.ylo
             </td>
             <td>
               BC type on the lower y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.bc.constant.type.yhi
             </td>
             <td>
               BC type on the upper y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.bc.constant.type.zlo
             </td>
             <td>
               BC type on the lower z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.bc.constant.type.zhi
             </td>
             <td>
               BC type on the upper z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.bc.constant.val.xlo
             </td>
             <td>
               BC value on the lower x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.bc.constant.val.xhi
             </td>
             <td>
               BC value on the upper x edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.bc.constant.val.ylo
             </td>
             <td>
               BC value on the lower y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.bc.constant.val.yhi
             </td>
             <td>
               BC value on the upper y edge (2d) face (3d)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.bc.constant.val.zlo
             </td>
             <td>
               BC value on the lower z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.bc.constant.val.zhi
             </td>
             <td>
               BC value on the upper z face (processed but ignored in 2d to prevent unused input errors)
             </td>
           </tr>
      </table><br/><br/>
