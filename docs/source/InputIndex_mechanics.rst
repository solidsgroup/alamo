mechanics
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
           which integrator to use (can only be mechanics) 
         </td>
       </tr>
       <tr>
         <td style='padding-left: 0px'>
           alamo.program.mechanics.model
         </td>
         <td>
           which mechanics model to use 
         </td>
       </tr>
      </table><br/><br/>
       <h3> Integrator::Mechanics&lt;Model::Solid::Linear::Isotropic&gt;  </h3>
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
             Type of mecahnics to use. Static: do full implicit solve. Dynamic: evolve dynamic equations with explicit dynamics Disable: do nothing. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             time_evolving
           </td>
           <td>
             Treat mechanics fields as changing in time. [false] You should use this if you care about other physics driven by the output of this integrator. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_disp
           </td>
           <td>
             Include displacement field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_rhs
           </td>
           <td>
             Include right-hand side in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_psi
           </td>
           <td>
             Include :math:`\psi` field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_stress
           </td>
           <td>
             Include stress in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_strain
           </td>
           <td>
             Include strain in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             viscous.mu_dashpot
           </td>
           <td>
             Dashpot damping (damps velocity)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             viscous.mu_newton
           </td>
           <td>
             Newtonian viscous damping (damps velocity gradient)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             velocity.ic.type
           </td>
           <td>
             Initializer for RHS
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = tensiontest
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.type
             </td>
             <td>
               Tension test type. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.disp
             </td>
             <td>
               Applied displacement (can be interpolator)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.trac
             </td>
             <td>
               Applied traction (can be interpolator)
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = expression
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             print_model
           </td>
           <td>
             Print out model variables (if enabled by model)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = trig
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.nr
             </td>
             <td>
               Number of real (cosin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.ni
             </td>
             <td>
               Number of imaginary (sin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.dim
             </td>
             <td>
               Spatial dimension
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.alpha
             </td>
             <td>
               Multiplier
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             interval
           </td>
           <td>
             Timestep interval for elastic solves (default - solve every time) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             max_coarsening_level
           </td>
           <td>
             Maximum multigrid coarsening level (default - none, maximum coarsening) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             print_residual
           </td>
           <td>
             Whether to include residual output field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic_ref_threshold
           </td>
           <td>
             Whether to refine based on elastic solution 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             zero_out_displacement
           </td>
           <td>
             Set this to true to zero out the displacement before each solve. (This is a temporary fix - we need to figure out why this is needed.) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             tstart
           </td>
           <td>
             Time to start doing the elastic solve (by default, start immediately) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             nmodels
           </td>
           <td>
             Number of elastic model varieties
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             eta_ref_threshold
           </td>
           <td>
             Refinement threshold for eta field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             ref_threshold
           </td>
           <td>
             Refinement threshold for strain gradient 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_neumann_boundary
           </td>
           <td>
             Explicity impose neumann condition on model at domain boundaries (2d only) 
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
              = bmp
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.fit
             </td>
             <td>
               How to position image in space 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.coord.lo
             </td>
             <td>
               Location of lower-left corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.coord.hi
             </td>
             <td>
               Location of upper-right corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.channel
             </td>
             <td>
               Color channel to use 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.min
             </td>
             <td>
               Scaling value - minimum
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.max
             </td>
             <td>
               Scaling value - maximum
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
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
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             eta.reset_on_regrid
           </td>
           <td>
             Whether to re-initialize eta when re-gridding occurs. Default is false unless eta ic is set, then default is. true. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = ellipse
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.x0
             </td>
             <td>
               Coorinates of ellipse center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.eps
             </td>
             <td>
               Diffuse boundary thickness
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               DxD square matrix defining an ellipse. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.a
             </td>
             <td>
               If :code:`A` is not defined, then assume a sphere with radius :code:`a`
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.number_of_inclusions
             </td>
             <td>
               Number of ellipses
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.center
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.x0
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               either a vector containing ellipse radii, or a matrix defining the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               Same
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.radius
             </td>
             <td>
               Array of radii [depricated]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.eps
             </td>
             <td>
               Regularization for smooth boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.invert
             </td>
             <td>
               Flip the inside and the outside 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             psi.reset_on_regrid
           </td>
           <td>
             Whether to re-initialize psi when re-gridding occurs. Default is false unless a psi ic is set, then default is true. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = ellipse
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.x0
             </td>
             <td>
               Coorinates of ellipse center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.eps
             </td>
             <td>
               Diffuse boundary thickness
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               DxD square matrix defining an ellipse. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.a
             </td>
             <td>
               If :code:`A` is not defined, then assume a sphere with radius :code:`a`
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.number_of_inclusions
             </td>
             <td>
               Number of ellipses
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.center
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.x0
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               either a vector containing ellipse radii, or a matrix defining the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               Same
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.radius
             </td>
             <td>
               Array of radii [depricated]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.eps
             </td>
             <td>
               Regularization for smooth boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.invert
             </td>
             <td>
               Flip the inside and the outside 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
      </table><br/><br/>
       <h3> Integrator::Mechanics&lt;Model::Solid::Linear::Cubic&gt;  </h3>
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
             Type of mecahnics to use. Static: do full implicit solve. Dynamic: evolve dynamic equations with explicit dynamics Disable: do nothing. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             time_evolving
           </td>
           <td>
             Treat mechanics fields as changing in time. [false] You should use this if you care about other physics driven by the output of this integrator. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_disp
           </td>
           <td>
             Include displacement field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_rhs
           </td>
           <td>
             Include right-hand side in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_psi
           </td>
           <td>
             Include :math:`\psi` field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_stress
           </td>
           <td>
             Include stress in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_strain
           </td>
           <td>
             Include strain in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             viscous.mu_dashpot
           </td>
           <td>
             Dashpot damping (damps velocity)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             viscous.mu_newton
           </td>
           <td>
             Newtonian viscous damping (damps velocity gradient)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             velocity.ic.type
           </td>
           <td>
             Initializer for RHS
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = tensiontest
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.type
             </td>
             <td>
               Tension test type. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.disp
             </td>
             <td>
               Applied displacement (can be interpolator)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.trac
             </td>
             <td>
               Applied traction (can be interpolator)
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = expression
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             print_model
           </td>
           <td>
             Print out model variables (if enabled by model)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = trig
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.nr
             </td>
             <td>
               Number of real (cosin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.ni
             </td>
             <td>
               Number of imaginary (sin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.dim
             </td>
             <td>
               Spatial dimension
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.alpha
             </td>
             <td>
               Multiplier
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             interval
           </td>
           <td>
             Timestep interval for elastic solves (default - solve every time) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             max_coarsening_level
           </td>
           <td>
             Maximum multigrid coarsening level (default - none, maximum coarsening) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             print_residual
           </td>
           <td>
             Whether to include residual output field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic_ref_threshold
           </td>
           <td>
             Whether to refine based on elastic solution 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             zero_out_displacement
           </td>
           <td>
             Set this to true to zero out the displacement before each solve. (This is a temporary fix - we need to figure out why this is needed.) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             tstart
           </td>
           <td>
             Time to start doing the elastic solve (by default, start immediately) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             nmodels
           </td>
           <td>
             Number of elastic model varieties
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             eta_ref_threshold
           </td>
           <td>
             Refinement threshold for eta field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             ref_threshold
           </td>
           <td>
             Refinement threshold for strain gradient 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_neumann_boundary
           </td>
           <td>
             Explicity impose neumann condition on model at domain boundaries (2d only) 
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
              = bmp
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.fit
             </td>
             <td>
               How to position image in space 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.coord.lo
             </td>
             <td>
               Location of lower-left corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.coord.hi
             </td>
             <td>
               Location of upper-right corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.channel
             </td>
             <td>
               Color channel to use 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.min
             </td>
             <td>
               Scaling value - minimum
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.max
             </td>
             <td>
               Scaling value - maximum
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
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
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             eta.reset_on_regrid
           </td>
           <td>
             Whether to re-initialize eta when re-gridding occurs. Default is false unless eta ic is set, then default is. true. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = ellipse
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.x0
             </td>
             <td>
               Coorinates of ellipse center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.eps
             </td>
             <td>
               Diffuse boundary thickness
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               DxD square matrix defining an ellipse. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.a
             </td>
             <td>
               If :code:`A` is not defined, then assume a sphere with radius :code:`a`
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.number_of_inclusions
             </td>
             <td>
               Number of ellipses
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.center
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.x0
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               either a vector containing ellipse radii, or a matrix defining the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               Same
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.radius
             </td>
             <td>
               Array of radii [depricated]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.eps
             </td>
             <td>
               Regularization for smooth boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.invert
             </td>
             <td>
               Flip the inside and the outside 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             psi.reset_on_regrid
           </td>
           <td>
             Whether to re-initialize psi when re-gridding occurs. Default is false unless a psi ic is set, then default is true. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = ellipse
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.x0
             </td>
             <td>
               Coorinates of ellipse center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.eps
             </td>
             <td>
               Diffuse boundary thickness
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               DxD square matrix defining an ellipse. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.a
             </td>
             <td>
               If :code:`A` is not defined, then assume a sphere with radius :code:`a`
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.number_of_inclusions
             </td>
             <td>
               Number of ellipses
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.center
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.x0
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               either a vector containing ellipse radii, or a matrix defining the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               Same
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.radius
             </td>
             <td>
               Array of radii [depricated]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.eps
             </td>
             <td>
               Regularization for smooth boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.invert
             </td>
             <td>
               Flip the inside and the outside 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
      </table><br/><br/>
       <h3> Integrator::Mechanics&lt;Model::Solid::Affine::Cubic&gt;  </h3>
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
             Type of mecahnics to use. Static: do full implicit solve. Dynamic: evolve dynamic equations with explicit dynamics Disable: do nothing. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             time_evolving
           </td>
           <td>
             Treat mechanics fields as changing in time. [false] You should use this if you care about other physics driven by the output of this integrator. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_disp
           </td>
           <td>
             Include displacement field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_rhs
           </td>
           <td>
             Include right-hand side in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_psi
           </td>
           <td>
             Include :math:`\psi` field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_stress
           </td>
           <td>
             Include stress in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_strain
           </td>
           <td>
             Include strain in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             viscous.mu_dashpot
           </td>
           <td>
             Dashpot damping (damps velocity)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             viscous.mu_newton
           </td>
           <td>
             Newtonian viscous damping (damps velocity gradient)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             velocity.ic.type
           </td>
           <td>
             Initializer for RHS
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = tensiontest
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.type
             </td>
             <td>
               Tension test type. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.disp
             </td>
             <td>
               Applied displacement (can be interpolator)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.trac
             </td>
             <td>
               Applied traction (can be interpolator)
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = expression
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             print_model
           </td>
           <td>
             Print out model variables (if enabled by model)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = trig
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.nr
             </td>
             <td>
               Number of real (cosin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.ni
             </td>
             <td>
               Number of imaginary (sin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.dim
             </td>
             <td>
               Spatial dimension
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.alpha
             </td>
             <td>
               Multiplier
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             interval
           </td>
           <td>
             Timestep interval for elastic solves (default - solve every time) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             max_coarsening_level
           </td>
           <td>
             Maximum multigrid coarsening level (default - none, maximum coarsening) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             print_residual
           </td>
           <td>
             Whether to include residual output field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic_ref_threshold
           </td>
           <td>
             Whether to refine based on elastic solution 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             zero_out_displacement
           </td>
           <td>
             Set this to true to zero out the displacement before each solve. (This is a temporary fix - we need to figure out why this is needed.) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             tstart
           </td>
           <td>
             Time to start doing the elastic solve (by default, start immediately) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             nmodels
           </td>
           <td>
             Number of elastic model varieties
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             eta_ref_threshold
           </td>
           <td>
             Refinement threshold for eta field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             ref_threshold
           </td>
           <td>
             Refinement threshold for strain gradient 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_neumann_boundary
           </td>
           <td>
             Explicity impose neumann condition on model at domain boundaries (2d only) 
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
              = bmp
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.fit
             </td>
             <td>
               How to position image in space 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.coord.lo
             </td>
             <td>
               Location of lower-left corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.coord.hi
             </td>
             <td>
               Location of upper-right corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.channel
             </td>
             <td>
               Color channel to use 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.min
             </td>
             <td>
               Scaling value - minimum
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.max
             </td>
             <td>
               Scaling value - maximum
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
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
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             eta.reset_on_regrid
           </td>
           <td>
             Whether to re-initialize eta when re-gridding occurs. Default is false unless eta ic is set, then default is. true. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = ellipse
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.x0
             </td>
             <td>
               Coorinates of ellipse center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.eps
             </td>
             <td>
               Diffuse boundary thickness
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               DxD square matrix defining an ellipse. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.a
             </td>
             <td>
               If :code:`A` is not defined, then assume a sphere with radius :code:`a`
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.number_of_inclusions
             </td>
             <td>
               Number of ellipses
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.center
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.x0
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               either a vector containing ellipse radii, or a matrix defining the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               Same
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.radius
             </td>
             <td>
               Array of radii [depricated]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.eps
             </td>
             <td>
               Regularization for smooth boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.invert
             </td>
             <td>
               Flip the inside and the outside 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             psi.reset_on_regrid
           </td>
           <td>
             Whether to re-initialize psi when re-gridding occurs. Default is false unless a psi ic is set, then default is true. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = ellipse
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.x0
             </td>
             <td>
               Coorinates of ellipse center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.eps
             </td>
             <td>
               Diffuse boundary thickness
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               DxD square matrix defining an ellipse. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.a
             </td>
             <td>
               If :code:`A` is not defined, then assume a sphere with radius :code:`a`
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.number_of_inclusions
             </td>
             <td>
               Number of ellipses
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.center
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.x0
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               either a vector containing ellipse radii, or a matrix defining the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               Same
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.radius
             </td>
             <td>
               Array of radii [depricated]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.eps
             </td>
             <td>
               Regularization for smooth boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.invert
             </td>
             <td>
               Flip the inside and the outside 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
      </table><br/><br/>
       <h3> Integrator::Mechanics&lt;Model::Solid::Affine::Hexagonal&gt;  </h3>
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
             Type of mecahnics to use. Static: do full implicit solve. Dynamic: evolve dynamic equations with explicit dynamics Disable: do nothing. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             time_evolving
           </td>
           <td>
             Treat mechanics fields as changing in time. [false] You should use this if you care about other physics driven by the output of this integrator. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_disp
           </td>
           <td>
             Include displacement field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_rhs
           </td>
           <td>
             Include right-hand side in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_psi
           </td>
           <td>
             Include :math:`\psi` field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_stress
           </td>
           <td>
             Include stress in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_strain
           </td>
           <td>
             Include strain in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             viscous.mu_dashpot
           </td>
           <td>
             Dashpot damping (damps velocity)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             viscous.mu_newton
           </td>
           <td>
             Newtonian viscous damping (damps velocity gradient)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             velocity.ic.type
           </td>
           <td>
             Initializer for RHS
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = tensiontest
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.type
             </td>
             <td>
               Tension test type. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.disp
             </td>
             <td>
               Applied displacement (can be interpolator)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.trac
             </td>
             <td>
               Applied traction (can be interpolator)
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = expression
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             print_model
           </td>
           <td>
             Print out model variables (if enabled by model)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = trig
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.nr
             </td>
             <td>
               Number of real (cosin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.ni
             </td>
             <td>
               Number of imaginary (sin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.dim
             </td>
             <td>
               Spatial dimension
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.alpha
             </td>
             <td>
               Multiplier
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             interval
           </td>
           <td>
             Timestep interval for elastic solves (default - solve every time) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             max_coarsening_level
           </td>
           <td>
             Maximum multigrid coarsening level (default - none, maximum coarsening) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             print_residual
           </td>
           <td>
             Whether to include residual output field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic_ref_threshold
           </td>
           <td>
             Whether to refine based on elastic solution 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             zero_out_displacement
           </td>
           <td>
             Set this to true to zero out the displacement before each solve. (This is a temporary fix - we need to figure out why this is needed.) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             tstart
           </td>
           <td>
             Time to start doing the elastic solve (by default, start immediately) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             nmodels
           </td>
           <td>
             Number of elastic model varieties
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             eta_ref_threshold
           </td>
           <td>
             Refinement threshold for eta field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             ref_threshold
           </td>
           <td>
             Refinement threshold for strain gradient 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_neumann_boundary
           </td>
           <td>
             Explicity impose neumann condition on model at domain boundaries (2d only) 
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
              = bmp
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.fit
             </td>
             <td>
               How to position image in space 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.coord.lo
             </td>
             <td>
               Location of lower-left corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.coord.hi
             </td>
             <td>
               Location of upper-right corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.channel
             </td>
             <td>
               Color channel to use 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.min
             </td>
             <td>
               Scaling value - minimum
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.max
             </td>
             <td>
               Scaling value - maximum
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
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
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             eta.reset_on_regrid
           </td>
           <td>
             Whether to re-initialize eta when re-gridding occurs. Default is false unless eta ic is set, then default is. true. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = ellipse
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.x0
             </td>
             <td>
               Coorinates of ellipse center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.eps
             </td>
             <td>
               Diffuse boundary thickness
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               DxD square matrix defining an ellipse. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.a
             </td>
             <td>
               If :code:`A` is not defined, then assume a sphere with radius :code:`a`
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.number_of_inclusions
             </td>
             <td>
               Number of ellipses
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.center
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.x0
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               either a vector containing ellipse radii, or a matrix defining the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               Same
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.radius
             </td>
             <td>
               Array of radii [depricated]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.eps
             </td>
             <td>
               Regularization for smooth boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.invert
             </td>
             <td>
               Flip the inside and the outside 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             psi.reset_on_regrid
           </td>
           <td>
             Whether to re-initialize psi when re-gridding occurs. Default is false unless a psi ic is set, then default is true. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = ellipse
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.x0
             </td>
             <td>
               Coorinates of ellipse center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.eps
             </td>
             <td>
               Diffuse boundary thickness
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               DxD square matrix defining an ellipse. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.a
             </td>
             <td>
               If :code:`A` is not defined, then assume a sphere with radius :code:`a`
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.number_of_inclusions
             </td>
             <td>
               Number of ellipses
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.center
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.x0
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               either a vector containing ellipse radii, or a matrix defining the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               Same
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.radius
             </td>
             <td>
               Array of radii [depricated]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.eps
             </td>
             <td>
               Regularization for smooth boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.invert
             </td>
             <td>
               Flip the inside and the outside 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
      </table><br/><br/>
       <h3> Integrator::Mechanics&lt;Model::Solid::Affine::Isotropic&gt;  </h3>
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
             Type of mecahnics to use. Static: do full implicit solve. Dynamic: evolve dynamic equations with explicit dynamics Disable: do nothing. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             time_evolving
           </td>
           <td>
             Treat mechanics fields as changing in time. [false] You should use this if you care about other physics driven by the output of this integrator. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_disp
           </td>
           <td>
             Include displacement field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_rhs
           </td>
           <td>
             Include right-hand side in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_psi
           </td>
           <td>
             Include :math:`\psi` field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_stress
           </td>
           <td>
             Include stress in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_strain
           </td>
           <td>
             Include strain in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             viscous.mu_dashpot
           </td>
           <td>
             Dashpot damping (damps velocity)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             viscous.mu_newton
           </td>
           <td>
             Newtonian viscous damping (damps velocity gradient)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             velocity.ic.type
           </td>
           <td>
             Initializer for RHS
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = tensiontest
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.type
             </td>
             <td>
               Tension test type. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.disp
             </td>
             <td>
               Applied displacement (can be interpolator)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.trac
             </td>
             <td>
               Applied traction (can be interpolator)
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = expression
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             print_model
           </td>
           <td>
             Print out model variables (if enabled by model)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = trig
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.nr
             </td>
             <td>
               Number of real (cosin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.ni
             </td>
             <td>
               Number of imaginary (sin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.dim
             </td>
             <td>
               Spatial dimension
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.alpha
             </td>
             <td>
               Multiplier
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             interval
           </td>
           <td>
             Timestep interval for elastic solves (default - solve every time) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             max_coarsening_level
           </td>
           <td>
             Maximum multigrid coarsening level (default - none, maximum coarsening) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             print_residual
           </td>
           <td>
             Whether to include residual output field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic_ref_threshold
           </td>
           <td>
             Whether to refine based on elastic solution 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             zero_out_displacement
           </td>
           <td>
             Set this to true to zero out the displacement before each solve. (This is a temporary fix - we need to figure out why this is needed.) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             tstart
           </td>
           <td>
             Time to start doing the elastic solve (by default, start immediately) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             nmodels
           </td>
           <td>
             Number of elastic model varieties
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             eta_ref_threshold
           </td>
           <td>
             Refinement threshold for eta field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             ref_threshold
           </td>
           <td>
             Refinement threshold for strain gradient 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_neumann_boundary
           </td>
           <td>
             Explicity impose neumann condition on model at domain boundaries (2d only) 
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
              = bmp
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.fit
             </td>
             <td>
               How to position image in space 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.coord.lo
             </td>
             <td>
               Location of lower-left corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.coord.hi
             </td>
             <td>
               Location of upper-right corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.channel
             </td>
             <td>
               Color channel to use 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.min
             </td>
             <td>
               Scaling value - minimum
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.max
             </td>
             <td>
               Scaling value - maximum
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
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
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             eta.reset_on_regrid
           </td>
           <td>
             Whether to re-initialize eta when re-gridding occurs. Default is false unless eta ic is set, then default is. true. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = ellipse
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.x0
             </td>
             <td>
               Coorinates of ellipse center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.eps
             </td>
             <td>
               Diffuse boundary thickness
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               DxD square matrix defining an ellipse. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.a
             </td>
             <td>
               If :code:`A` is not defined, then assume a sphere with radius :code:`a`
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.number_of_inclusions
             </td>
             <td>
               Number of ellipses
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.center
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.x0
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               either a vector containing ellipse radii, or a matrix defining the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               Same
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.radius
             </td>
             <td>
               Array of radii [depricated]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.eps
             </td>
             <td>
               Regularization for smooth boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.invert
             </td>
             <td>
               Flip the inside and the outside 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             psi.reset_on_regrid
           </td>
           <td>
             Whether to re-initialize psi when re-gridding occurs. Default is false unless a psi ic is set, then default is true. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = ellipse
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.x0
             </td>
             <td>
               Coorinates of ellipse center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.eps
             </td>
             <td>
               Diffuse boundary thickness
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               DxD square matrix defining an ellipse. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.a
             </td>
             <td>
               If :code:`A` is not defined, then assume a sphere with radius :code:`a`
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.number_of_inclusions
             </td>
             <td>
               Number of ellipses
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.center
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.x0
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               either a vector containing ellipse radii, or a matrix defining the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               Same
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.radius
             </td>
             <td>
               Array of radii [depricated]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.eps
             </td>
             <td>
               Regularization for smooth boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.invert
             </td>
             <td>
               Flip the inside and the outside 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
      </table><br/><br/>
       <h3> Integrator::Mechanics&lt;Model::Solid::Linear::Laplacian&gt;  </h3>
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
             Type of mecahnics to use. Static: do full implicit solve. Dynamic: evolve dynamic equations with explicit dynamics Disable: do nothing. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             time_evolving
           </td>
           <td>
             Treat mechanics fields as changing in time. [false] You should use this if you care about other physics driven by the output of this integrator. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_disp
           </td>
           <td>
             Include displacement field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_rhs
           </td>
           <td>
             Include right-hand side in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_psi
           </td>
           <td>
             Include :math:`\psi` field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_stress
           </td>
           <td>
             Include stress in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_strain
           </td>
           <td>
             Include strain in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             viscous.mu_dashpot
           </td>
           <td>
             Dashpot damping (damps velocity)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             viscous.mu_newton
           </td>
           <td>
             Newtonian viscous damping (damps velocity gradient)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             velocity.ic.type
           </td>
           <td>
             Initializer for RHS
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = tensiontest
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.type
             </td>
             <td>
               Tension test type. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.disp
             </td>
             <td>
               Applied displacement (can be interpolator)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.trac
             </td>
             <td>
               Applied traction (can be interpolator)
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = expression
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             print_model
           </td>
           <td>
             Print out model variables (if enabled by model)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = trig
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.nr
             </td>
             <td>
               Number of real (cosin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.ni
             </td>
             <td>
               Number of imaginary (sin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.dim
             </td>
             <td>
               Spatial dimension
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.alpha
             </td>
             <td>
               Multiplier
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             interval
           </td>
           <td>
             Timestep interval for elastic solves (default - solve every time) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             max_coarsening_level
           </td>
           <td>
             Maximum multigrid coarsening level (default - none, maximum coarsening) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             print_residual
           </td>
           <td>
             Whether to include residual output field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic_ref_threshold
           </td>
           <td>
             Whether to refine based on elastic solution 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             zero_out_displacement
           </td>
           <td>
             Set this to true to zero out the displacement before each solve. (This is a temporary fix - we need to figure out why this is needed.) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             tstart
           </td>
           <td>
             Time to start doing the elastic solve (by default, start immediately) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             nmodels
           </td>
           <td>
             Number of elastic model varieties
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             eta_ref_threshold
           </td>
           <td>
             Refinement threshold for eta field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             ref_threshold
           </td>
           <td>
             Refinement threshold for strain gradient 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_neumann_boundary
           </td>
           <td>
             Explicity impose neumann condition on model at domain boundaries (2d only) 
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
              = bmp
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.fit
             </td>
             <td>
               How to position image in space 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.coord.lo
             </td>
             <td>
               Location of lower-left corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.coord.hi
             </td>
             <td>
               Location of upper-right corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.channel
             </td>
             <td>
               Color channel to use 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.min
             </td>
             <td>
               Scaling value - minimum
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.max
             </td>
             <td>
               Scaling value - maximum
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
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
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             eta.reset_on_regrid
           </td>
           <td>
             Whether to re-initialize eta when re-gridding occurs. Default is false unless eta ic is set, then default is. true. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = ellipse
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.x0
             </td>
             <td>
               Coorinates of ellipse center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.eps
             </td>
             <td>
               Diffuse boundary thickness
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               DxD square matrix defining an ellipse. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.a
             </td>
             <td>
               If :code:`A` is not defined, then assume a sphere with radius :code:`a`
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.number_of_inclusions
             </td>
             <td>
               Number of ellipses
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.center
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.x0
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               either a vector containing ellipse radii, or a matrix defining the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               Same
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.radius
             </td>
             <td>
               Array of radii [depricated]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.eps
             </td>
             <td>
               Regularization for smooth boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.invert
             </td>
             <td>
               Flip the inside and the outside 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             psi.reset_on_regrid
           </td>
           <td>
             Whether to re-initialize psi when re-gridding occurs. Default is false unless a psi ic is set, then default is true. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = ellipse
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.x0
             </td>
             <td>
               Coorinates of ellipse center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.eps
             </td>
             <td>
               Diffuse boundary thickness
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               DxD square matrix defining an ellipse. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.a
             </td>
             <td>
               If :code:`A` is not defined, then assume a sphere with radius :code:`a`
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.number_of_inclusions
             </td>
             <td>
               Number of ellipses
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.center
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.x0
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               either a vector containing ellipse radii, or a matrix defining the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               Same
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.radius
             </td>
             <td>
               Array of radii [depricated]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.eps
             </td>
             <td>
               Regularization for smooth boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.invert
             </td>
             <td>
               Flip the inside and the outside 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
      </table><br/><br/>
       <h3> Integrator::Mechanics&lt;Model::Solid::Finite::NeoHookean&gt;  </h3>
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
             Type of mecahnics to use. Static: do full implicit solve. Dynamic: evolve dynamic equations with explicit dynamics Disable: do nothing. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             time_evolving
           </td>
           <td>
             Treat mechanics fields as changing in time. [false] You should use this if you care about other physics driven by the output of this integrator. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_disp
           </td>
           <td>
             Include displacement field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_rhs
           </td>
           <td>
             Include right-hand side in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_psi
           </td>
           <td>
             Include :math:`\psi` field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_stress
           </td>
           <td>
             Include stress in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_strain
           </td>
           <td>
             Include strain in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             viscous.mu_dashpot
           </td>
           <td>
             Dashpot damping (damps velocity)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             viscous.mu_newton
           </td>
           <td>
             Newtonian viscous damping (damps velocity gradient)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             velocity.ic.type
           </td>
           <td>
             Initializer for RHS
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = tensiontest
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.type
             </td>
             <td>
               Tension test type. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.disp
             </td>
             <td>
               Applied displacement (can be interpolator)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.trac
             </td>
             <td>
               Applied traction (can be interpolator)
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = expression
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             print_model
           </td>
           <td>
             Print out model variables (if enabled by model)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = trig
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.nr
             </td>
             <td>
               Number of real (cosin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.ni
             </td>
             <td>
               Number of imaginary (sin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.dim
             </td>
             <td>
               Spatial dimension
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.alpha
             </td>
             <td>
               Multiplier
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             interval
           </td>
           <td>
             Timestep interval for elastic solves (default - solve every time) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             max_coarsening_level
           </td>
           <td>
             Maximum multigrid coarsening level (default - none, maximum coarsening) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             print_residual
           </td>
           <td>
             Whether to include residual output field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic_ref_threshold
           </td>
           <td>
             Whether to refine based on elastic solution 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             zero_out_displacement
           </td>
           <td>
             Set this to true to zero out the displacement before each solve. (This is a temporary fix - we need to figure out why this is needed.) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             tstart
           </td>
           <td>
             Time to start doing the elastic solve (by default, start immediately) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             nmodels
           </td>
           <td>
             Number of elastic model varieties
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             eta_ref_threshold
           </td>
           <td>
             Refinement threshold for eta field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             ref_threshold
           </td>
           <td>
             Refinement threshold for strain gradient 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_neumann_boundary
           </td>
           <td>
             Explicity impose neumann condition on model at domain boundaries (2d only) 
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
              = bmp
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.fit
             </td>
             <td>
               How to position image in space 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.coord.lo
             </td>
             <td>
               Location of lower-left corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.coord.hi
             </td>
             <td>
               Location of upper-right corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.channel
             </td>
             <td>
               Color channel to use 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.min
             </td>
             <td>
               Scaling value - minimum
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.max
             </td>
             <td>
               Scaling value - maximum
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
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
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             eta.reset_on_regrid
           </td>
           <td>
             Whether to re-initialize eta when re-gridding occurs. Default is false unless eta ic is set, then default is. true. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = ellipse
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.x0
             </td>
             <td>
               Coorinates of ellipse center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.eps
             </td>
             <td>
               Diffuse boundary thickness
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               DxD square matrix defining an ellipse. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.a
             </td>
             <td>
               If :code:`A` is not defined, then assume a sphere with radius :code:`a`
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.number_of_inclusions
             </td>
             <td>
               Number of ellipses
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.center
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.x0
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               either a vector containing ellipse radii, or a matrix defining the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               Same
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.radius
             </td>
             <td>
               Array of radii [depricated]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.eps
             </td>
             <td>
               Regularization for smooth boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.invert
             </td>
             <td>
               Flip the inside and the outside 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             psi.reset_on_regrid
           </td>
           <td>
             Whether to re-initialize psi when re-gridding occurs. Default is false unless a psi ic is set, then default is true. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = ellipse
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.x0
             </td>
             <td>
               Coorinates of ellipse center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.eps
             </td>
             <td>
               Diffuse boundary thickness
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               DxD square matrix defining an ellipse. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.a
             </td>
             <td>
               If :code:`A` is not defined, then assume a sphere with radius :code:`a`
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.number_of_inclusions
             </td>
             <td>
               Number of ellipses
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.center
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.x0
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               either a vector containing ellipse radii, or a matrix defining the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               Same
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.radius
             </td>
             <td>
               Array of radii [depricated]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.eps
             </td>
             <td>
               Regularization for smooth boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.invert
             </td>
             <td>
               Flip the inside and the outside 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
      </table><br/><br/>
       <h3> Integrator::Mechanics&lt;Model::Solid::Finite::NeoHookeanPredeformed&gt;  </h3>
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
             Type of mecahnics to use. Static: do full implicit solve. Dynamic: evolve dynamic equations with explicit dynamics Disable: do nothing. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             time_evolving
           </td>
           <td>
             Treat mechanics fields as changing in time. [false] You should use this if you care about other physics driven by the output of this integrator. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_disp
           </td>
           <td>
             Include displacement field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_rhs
           </td>
           <td>
             Include right-hand side in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_psi
           </td>
           <td>
             Include :math:`\psi` field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_stress
           </td>
           <td>
             Include stress in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_strain
           </td>
           <td>
             Include strain in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             viscous.mu_dashpot
           </td>
           <td>
             Dashpot damping (damps velocity)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             viscous.mu_newton
           </td>
           <td>
             Newtonian viscous damping (damps velocity gradient)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             velocity.ic.type
           </td>
           <td>
             Initializer for RHS
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = tensiontest
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.type
             </td>
             <td>
               Tension test type. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.disp
             </td>
             <td>
               Applied displacement (can be interpolator)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.trac
             </td>
             <td>
               Applied traction (can be interpolator)
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = expression
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             print_model
           </td>
           <td>
             Print out model variables (if enabled by model)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = trig
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.nr
             </td>
             <td>
               Number of real (cosin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.ni
             </td>
             <td>
               Number of imaginary (sin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.dim
             </td>
             <td>
               Spatial dimension
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.alpha
             </td>
             <td>
               Multiplier
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             interval
           </td>
           <td>
             Timestep interval for elastic solves (default - solve every time) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             max_coarsening_level
           </td>
           <td>
             Maximum multigrid coarsening level (default - none, maximum coarsening) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             print_residual
           </td>
           <td>
             Whether to include residual output field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic_ref_threshold
           </td>
           <td>
             Whether to refine based on elastic solution 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             zero_out_displacement
           </td>
           <td>
             Set this to true to zero out the displacement before each solve. (This is a temporary fix - we need to figure out why this is needed.) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             tstart
           </td>
           <td>
             Time to start doing the elastic solve (by default, start immediately) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             nmodels
           </td>
           <td>
             Number of elastic model varieties
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             eta_ref_threshold
           </td>
           <td>
             Refinement threshold for eta field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             ref_threshold
           </td>
           <td>
             Refinement threshold for strain gradient 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_neumann_boundary
           </td>
           <td>
             Explicity impose neumann condition on model at domain boundaries (2d only) 
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
              = bmp
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.fit
             </td>
             <td>
               How to position image in space 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.coord.lo
             </td>
             <td>
               Location of lower-left corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.coord.hi
             </td>
             <td>
               Location of upper-right corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.channel
             </td>
             <td>
               Color channel to use 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.min
             </td>
             <td>
               Scaling value - minimum
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.max
             </td>
             <td>
               Scaling value - maximum
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
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
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             eta.reset_on_regrid
           </td>
           <td>
             Whether to re-initialize eta when re-gridding occurs. Default is false unless eta ic is set, then default is. true. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = ellipse
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.x0
             </td>
             <td>
               Coorinates of ellipse center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.eps
             </td>
             <td>
               Diffuse boundary thickness
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               DxD square matrix defining an ellipse. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.a
             </td>
             <td>
               If :code:`A` is not defined, then assume a sphere with radius :code:`a`
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.number_of_inclusions
             </td>
             <td>
               Number of ellipses
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.center
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.x0
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               either a vector containing ellipse radii, or a matrix defining the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               Same
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.radius
             </td>
             <td>
               Array of radii [depricated]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.eps
             </td>
             <td>
               Regularization for smooth boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.invert
             </td>
             <td>
               Flip the inside and the outside 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             psi.reset_on_regrid
           </td>
           <td>
             Whether to re-initialize psi when re-gridding occurs. Default is false unless a psi ic is set, then default is true. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = ellipse
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.x0
             </td>
             <td>
               Coorinates of ellipse center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.eps
             </td>
             <td>
               Diffuse boundary thickness
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               DxD square matrix defining an ellipse. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.a
             </td>
             <td>
               If :code:`A` is not defined, then assume a sphere with radius :code:`a`
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.number_of_inclusions
             </td>
             <td>
               Number of ellipses
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.center
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.x0
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               either a vector containing ellipse radii, or a matrix defining the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               Same
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.radius
             </td>
             <td>
               Array of radii [depricated]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.eps
             </td>
             <td>
               Regularization for smooth boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.invert
             </td>
             <td>
               Flip the inside and the outside 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
      </table><br/><br/>
       <h3> Integrator::Mechanics&lt;Model::Solid::Finite::PseudoLinear::Cubic&gt;  </h3>
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
             Type of mecahnics to use. Static: do full implicit solve. Dynamic: evolve dynamic equations with explicit dynamics Disable: do nothing. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             time_evolving
           </td>
           <td>
             Treat mechanics fields as changing in time. [false] You should use this if you care about other physics driven by the output of this integrator. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_disp
           </td>
           <td>
             Include displacement field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_rhs
           </td>
           <td>
             Include right-hand side in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_psi
           </td>
           <td>
             Include :math:`\psi` field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_stress
           </td>
           <td>
             Include stress in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_strain
           </td>
           <td>
             Include strain in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             viscous.mu_dashpot
           </td>
           <td>
             Dashpot damping (damps velocity)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             viscous.mu_newton
           </td>
           <td>
             Newtonian viscous damping (damps velocity gradient)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             velocity.ic.type
           </td>
           <td>
             Initializer for RHS
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = tensiontest
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.type
             </td>
             <td>
               Tension test type. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.disp
             </td>
             <td>
               Applied displacement (can be interpolator)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.trac
             </td>
             <td>
               Applied traction (can be interpolator)
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = expression
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             print_model
           </td>
           <td>
             Print out model variables (if enabled by model)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = trig
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.nr
             </td>
             <td>
               Number of real (cosin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.ni
             </td>
             <td>
               Number of imaginary (sin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.dim
             </td>
             <td>
               Spatial dimension
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.alpha
             </td>
             <td>
               Multiplier
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             interval
           </td>
           <td>
             Timestep interval for elastic solves (default - solve every time) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             max_coarsening_level
           </td>
           <td>
             Maximum multigrid coarsening level (default - none, maximum coarsening) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             print_residual
           </td>
           <td>
             Whether to include residual output field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic_ref_threshold
           </td>
           <td>
             Whether to refine based on elastic solution 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             zero_out_displacement
           </td>
           <td>
             Set this to true to zero out the displacement before each solve. (This is a temporary fix - we need to figure out why this is needed.) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             tstart
           </td>
           <td>
             Time to start doing the elastic solve (by default, start immediately) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             nmodels
           </td>
           <td>
             Number of elastic model varieties
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             eta_ref_threshold
           </td>
           <td>
             Refinement threshold for eta field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             ref_threshold
           </td>
           <td>
             Refinement threshold for strain gradient 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_neumann_boundary
           </td>
           <td>
             Explicity impose neumann condition on model at domain boundaries (2d only) 
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
              = bmp
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.fit
             </td>
             <td>
               How to position image in space 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.coord.lo
             </td>
             <td>
               Location of lower-left corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.coord.hi
             </td>
             <td>
               Location of upper-right corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.channel
             </td>
             <td>
               Color channel to use 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.min
             </td>
             <td>
               Scaling value - minimum
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.max
             </td>
             <td>
               Scaling value - maximum
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
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
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             eta.reset_on_regrid
           </td>
           <td>
             Whether to re-initialize eta when re-gridding occurs. Default is false unless eta ic is set, then default is. true. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = ellipse
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.x0
             </td>
             <td>
               Coorinates of ellipse center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.eps
             </td>
             <td>
               Diffuse boundary thickness
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               DxD square matrix defining an ellipse. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.a
             </td>
             <td>
               If :code:`A` is not defined, then assume a sphere with radius :code:`a`
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.number_of_inclusions
             </td>
             <td>
               Number of ellipses
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.center
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.x0
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               either a vector containing ellipse radii, or a matrix defining the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               Same
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.radius
             </td>
             <td>
               Array of radii [depricated]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.eps
             </td>
             <td>
               Regularization for smooth boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.invert
             </td>
             <td>
               Flip the inside and the outside 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             psi.reset_on_regrid
           </td>
           <td>
             Whether to re-initialize psi when re-gridding occurs. Default is false unless a psi ic is set, then default is true. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = ellipse
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.x0
             </td>
             <td>
               Coorinates of ellipse center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.eps
             </td>
             <td>
               Diffuse boundary thickness
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               DxD square matrix defining an ellipse. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.a
             </td>
             <td>
               If :code:`A` is not defined, then assume a sphere with radius :code:`a`
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.number_of_inclusions
             </td>
             <td>
               Number of ellipses
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.center
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.x0
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               either a vector containing ellipse radii, or a matrix defining the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               Same
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.radius
             </td>
             <td>
               Array of radii [depricated]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.eps
             </td>
             <td>
               Regularization for smooth boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.invert
             </td>
             <td>
               Flip the inside and the outside 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
      </table><br/><br/>
       <h3> Integrator::Mechanics&lt;Model::Solid::Finite::PseudoAffine::Cubic&gt;  </h3>
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
             Type of mecahnics to use. Static: do full implicit solve. Dynamic: evolve dynamic equations with explicit dynamics Disable: do nothing. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             time_evolving
           </td>
           <td>
             Treat mechanics fields as changing in time. [false] You should use this if you care about other physics driven by the output of this integrator. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_disp
           </td>
           <td>
             Include displacement field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_rhs
           </td>
           <td>
             Include right-hand side in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_psi
           </td>
           <td>
             Include :math:`\psi` field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_stress
           </td>
           <td>
             Include stress in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_strain
           </td>
           <td>
             Include strain in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             viscous.mu_dashpot
           </td>
           <td>
             Dashpot damping (damps velocity)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             viscous.mu_newton
           </td>
           <td>
             Newtonian viscous damping (damps velocity gradient)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             velocity.ic.type
           </td>
           <td>
             Initializer for RHS
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = tensiontest
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.type
             </td>
             <td>
               Tension test type. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.disp
             </td>
             <td>
               Applied displacement (can be interpolator)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.trac
             </td>
             <td>
               Applied traction (can be interpolator)
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = expression
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             print_model
           </td>
           <td>
             Print out model variables (if enabled by model)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = trig
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.nr
             </td>
             <td>
               Number of real (cosin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.ni
             </td>
             <td>
               Number of imaginary (sin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.dim
             </td>
             <td>
               Spatial dimension
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.alpha
             </td>
             <td>
               Multiplier
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             interval
           </td>
           <td>
             Timestep interval for elastic solves (default - solve every time) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             max_coarsening_level
           </td>
           <td>
             Maximum multigrid coarsening level (default - none, maximum coarsening) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             print_residual
           </td>
           <td>
             Whether to include residual output field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic_ref_threshold
           </td>
           <td>
             Whether to refine based on elastic solution 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             zero_out_displacement
           </td>
           <td>
             Set this to true to zero out the displacement before each solve. (This is a temporary fix - we need to figure out why this is needed.) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             tstart
           </td>
           <td>
             Time to start doing the elastic solve (by default, start immediately) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             nmodels
           </td>
           <td>
             Number of elastic model varieties
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             eta_ref_threshold
           </td>
           <td>
             Refinement threshold for eta field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             ref_threshold
           </td>
           <td>
             Refinement threshold for strain gradient 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_neumann_boundary
           </td>
           <td>
             Explicity impose neumann condition on model at domain boundaries (2d only) 
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
              = bmp
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.fit
             </td>
             <td>
               How to position image in space 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.coord.lo
             </td>
             <td>
               Location of lower-left corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.coord.hi
             </td>
             <td>
               Location of upper-right corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.channel
             </td>
             <td>
               Color channel to use 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.min
             </td>
             <td>
               Scaling value - minimum
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.max
             </td>
             <td>
               Scaling value - maximum
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
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
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             eta.reset_on_regrid
           </td>
           <td>
             Whether to re-initialize eta when re-gridding occurs. Default is false unless eta ic is set, then default is. true. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = ellipse
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.x0
             </td>
             <td>
               Coorinates of ellipse center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.eps
             </td>
             <td>
               Diffuse boundary thickness
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               DxD square matrix defining an ellipse. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.a
             </td>
             <td>
               If :code:`A` is not defined, then assume a sphere with radius :code:`a`
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.number_of_inclusions
             </td>
             <td>
               Number of ellipses
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.center
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.x0
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               either a vector containing ellipse radii, or a matrix defining the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               Same
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.radius
             </td>
             <td>
               Array of radii [depricated]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.eps
             </td>
             <td>
               Regularization for smooth boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.invert
             </td>
             <td>
               Flip the inside and the outside 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             psi.reset_on_regrid
           </td>
           <td>
             Whether to re-initialize psi when re-gridding occurs. Default is false unless a psi ic is set, then default is true. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = ellipse
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.x0
             </td>
             <td>
               Coorinates of ellipse center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.eps
             </td>
             <td>
               Diffuse boundary thickness
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               DxD square matrix defining an ellipse. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.a
             </td>
             <td>
               If :code:`A` is not defined, then assume a sphere with radius :code:`a`
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.number_of_inclusions
             </td>
             <td>
               Number of ellipses
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.center
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.x0
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               either a vector containing ellipse radii, or a matrix defining the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               Same
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.radius
             </td>
             <td>
               Array of radii [depricated]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.eps
             </td>
             <td>
               Regularization for smooth boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.invert
             </td>
             <td>
               Flip the inside and the outside 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
      </table><br/><br/>
       <h3> Integrator::Mechanics&lt;Model::Solid::Affine::J2&gt;  </h3>
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
             Type of mecahnics to use. Static: do full implicit solve. Dynamic: evolve dynamic equations with explicit dynamics Disable: do nothing. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             time_evolving
           </td>
           <td>
             Treat mechanics fields as changing in time. [false] You should use this if you care about other physics driven by the output of this integrator. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_disp
           </td>
           <td>
             Include displacement field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_rhs
           </td>
           <td>
             Include right-hand side in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_psi
           </td>
           <td>
             Include :math:`\psi` field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_stress
           </td>
           <td>
             Include stress in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_strain
           </td>
           <td>
             Include strain in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             viscous.mu_dashpot
           </td>
           <td>
             Dashpot damping (damps velocity)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             viscous.mu_newton
           </td>
           <td>
             Newtonian viscous damping (damps velocity gradient)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             velocity.ic.type
           </td>
           <td>
             Initializer for RHS
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = tensiontest
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.type
             </td>
             <td>
               Tension test type. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.disp
             </td>
             <td>
               Applied displacement (can be interpolator)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.trac
             </td>
             <td>
               Applied traction (can be interpolator)
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = expression
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             print_model
           </td>
           <td>
             Print out model variables (if enabled by model)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = trig
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.nr
             </td>
             <td>
               Number of real (cosin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.ni
             </td>
             <td>
               Number of imaginary (sin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.dim
             </td>
             <td>
               Spatial dimension
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.alpha
             </td>
             <td>
               Multiplier
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             interval
           </td>
           <td>
             Timestep interval for elastic solves (default - solve every time) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             max_coarsening_level
           </td>
           <td>
             Maximum multigrid coarsening level (default - none, maximum coarsening) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             print_residual
           </td>
           <td>
             Whether to include residual output field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic_ref_threshold
           </td>
           <td>
             Whether to refine based on elastic solution 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             zero_out_displacement
           </td>
           <td>
             Set this to true to zero out the displacement before each solve. (This is a temporary fix - we need to figure out why this is needed.) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             tstart
           </td>
           <td>
             Time to start doing the elastic solve (by default, start immediately) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             nmodels
           </td>
           <td>
             Number of elastic model varieties
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             eta_ref_threshold
           </td>
           <td>
             Refinement threshold for eta field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             ref_threshold
           </td>
           <td>
             Refinement threshold for strain gradient 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_neumann_boundary
           </td>
           <td>
             Explicity impose neumann condition on model at domain boundaries (2d only) 
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
              = bmp
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.fit
             </td>
             <td>
               How to position image in space 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.coord.lo
             </td>
             <td>
               Location of lower-left corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.coord.hi
             </td>
             <td>
               Location of upper-right corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.channel
             </td>
             <td>
               Color channel to use 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.min
             </td>
             <td>
               Scaling value - minimum
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.max
             </td>
             <td>
               Scaling value - maximum
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
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
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             eta.reset_on_regrid
           </td>
           <td>
             Whether to re-initialize eta when re-gridding occurs. Default is false unless eta ic is set, then default is. true. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = ellipse
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.x0
             </td>
             <td>
               Coorinates of ellipse center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.eps
             </td>
             <td>
               Diffuse boundary thickness
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               DxD square matrix defining an ellipse. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.a
             </td>
             <td>
               If :code:`A` is not defined, then assume a sphere with radius :code:`a`
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.number_of_inclusions
             </td>
             <td>
               Number of ellipses
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.center
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.x0
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               either a vector containing ellipse radii, or a matrix defining the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               Same
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.radius
             </td>
             <td>
               Array of radii [depricated]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.eps
             </td>
             <td>
               Regularization for smooth boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.invert
             </td>
             <td>
               Flip the inside and the outside 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             psi.reset_on_regrid
           </td>
           <td>
             Whether to re-initialize psi when re-gridding occurs. Default is false unless a psi ic is set, then default is true. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = ellipse
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.x0
             </td>
             <td>
               Coorinates of ellipse center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.eps
             </td>
             <td>
               Diffuse boundary thickness
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               DxD square matrix defining an ellipse. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.a
             </td>
             <td>
               If :code:`A` is not defined, then assume a sphere with radius :code:`a`
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.number_of_inclusions
             </td>
             <td>
               Number of ellipses
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.center
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.x0
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               either a vector containing ellipse radii, or a matrix defining the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               Same
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.radius
             </td>
             <td>
               Array of radii [depricated]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.eps
             </td>
             <td>
               Regularization for smooth boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.invert
             </td>
             <td>
               Flip the inside and the outside 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
      </table><br/><br/>
       <h3> Integrator::Mechanics&lt;Model::Solid::Finite::CrystalPlastic&gt;  </h3>
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
             Type of mecahnics to use. Static: do full implicit solve. Dynamic: evolve dynamic equations with explicit dynamics Disable: do nothing. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             time_evolving
           </td>
           <td>
             Treat mechanics fields as changing in time. [false] You should use this if you care about other physics driven by the output of this integrator. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_disp
           </td>
           <td>
             Include displacement field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_rhs
           </td>
           <td>
             Include right-hand side in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_psi
           </td>
           <td>
             Include :math:`\psi` field in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_stress
           </td>
           <td>
             Include stress in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             plot_strain
           </td>
           <td>
             Include strain in output
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             viscous.mu_dashpot
           </td>
           <td>
             Dashpot damping (damps velocity)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             viscous.mu_newton
           </td>
           <td>
             Newtonian viscous damping (damps velocity gradient)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             velocity.ic.type
           </td>
           <td>
             Initializer for RHS
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.type.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylozlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylozhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhizlo
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhizhi
             </td>
             <td>
               3D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylozlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylozhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhizlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhizhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zloxlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zloxhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhixlo
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhixhi
             </td>
             <td>
               3D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xloyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiylo
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhiyhi
             </td>
             <td>
               3D Edge / 2D Corner
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xlo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.xhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.ylo
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.yhi
             </td>
             <td>
               3D Face / 2D Edge
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zlo
             </td>
             <td>
               3D Face
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.constant.val.zhi
             </td>
             <td>
               3D Face
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = tensiontest
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.type
             </td>
             <td>
               Tension test type. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.disp
             </td>
             <td>
               Applied displacement (can be interpolator)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               bc.tensiontest.trac
             </td>
             <td>
               Applied traction (can be interpolator)
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              bc.type
              = expression
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             print_model
           </td>
           <td>
             Print out model variables (if enabled by model)
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              rhs.type
              = trig
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.nr
             </td>
             <td>
               Number of real (cosin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.ni
             </td>
             <td>
               Number of imaginary (sin) waves
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.dim
             </td>
             <td>
               Spatial dimension
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               rhs.trig.alpha
             </td>
             <td>
               Multiplier
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             interval
           </td>
           <td>
             Timestep interval for elastic solves (default - solve every time) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             max_coarsening_level
           </td>
           <td>
             Maximum multigrid coarsening level (default - none, maximum coarsening) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             print_residual
           </td>
           <td>
             Whether to include residual output field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             elastic_ref_threshold
           </td>
           <td>
             Whether to refine based on elastic solution 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             zero_out_displacement
           </td>
           <td>
             Set this to true to zero out the displacement before each solve. (This is a temporary fix - we need to figure out why this is needed.) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             tstart
           </td>
           <td>
             Time to start doing the elastic solve (by default, start immediately) 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             nmodels
           </td>
           <td>
             Number of elastic model varieties
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             eta_ref_threshold
           </td>
           <td>
             Refinement threshold for eta field 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             ref_threshold
           </td>
           <td>
             Refinement threshold for strain gradient 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             model_neumann_boundary
           </td>
           <td>
             Explicity impose neumann condition on model at domain boundaries (2d only) 
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
              = bmp
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.fit
             </td>
             <td>
               How to position image in space 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.coord.lo
             </td>
             <td>
               Location of lower-left corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.coord.hi
             </td>
             <td>
               Location of upper-right corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.channel
             </td>
             <td>
               Color channel to use 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.min
             </td>
             <td>
               Scaling value - minimum
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.bmp.max
             </td>
             <td>
               Scaling value - maximum
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
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
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             eta.reset_on_regrid
           </td>
           <td>
             Whether to re-initialize eta when re-gridding occurs. Default is false unless eta ic is set, then default is. true. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = ellipse
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.x0
             </td>
             <td>
               Coorinates of ellipse center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.eps
             </td>
             <td>
               Diffuse boundary thickness
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               DxD square matrix defining an ellipse. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.a
             </td>
             <td>
               If :code:`A` is not defined, then assume a sphere with radius :code:`a`
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.number_of_inclusions
             </td>
             <td>
               Number of ellipses
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.center
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.x0
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               either a vector containing ellipse radii, or a matrix defining the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.A
             </td>
             <td>
               Same
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.radius
             </td>
             <td>
               Array of radii [depricated]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.eps
             </td>
             <td>
               Regularization for smooth boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.ellipse.invert
             </td>
             <td>
               Flip the inside and the outside 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              psi.ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               psi.ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             psi.reset_on_regrid
           </td>
           <td>
             Whether to re-initialize psi when re-gridding occurs. Default is false unless a psi ic is set, then default is true. 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = ellipse
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.x0
             </td>
             <td>
               Coorinates of ellipse center
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.eps
             </td>
             <td>
               Diffuse boundary thickness
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               DxD square matrix defining an ellipse. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.a
             </td>
             <td>
               If :code:`A` is not defined, then assume a sphere with radius :code:`a`
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.number_of_inclusions
             </td>
             <td>
               Number of ellipses
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.center
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.x0
             </td>
             <td>
               center of the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               either a vector containing ellipse radii, or a matrix defining the ellipse
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.A
             </td>
             <td>
               Same
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.radius
             </td>
             <td>
               Array of radii [depricated]
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.eps
             </td>
             <td>
               Regularization for smooth boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.ellipse.invert
             </td>
             <td>
               Flip the inside and the outside 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              trac_normal.ic.type
              = psread
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.eps
             </td>
             <td>
               Diffuseness of the sphere boundary
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.filename
             </td>
             <td>
               Location of .xyzr file
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.verbose
             </td>
             <td>
               Verbosity (used in parser only)
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.mult
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.invert
             </td>
             <td>
               Coordinate multiplier
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               trac_normal.ic.psread.x0
             </td>
             <td>
               Coordinate offset
             </td>
           </tr>
      </table><br/><br/>
