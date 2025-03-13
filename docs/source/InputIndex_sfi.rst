sfi
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
           which integrator to use with SFI 
         </td>
       </tr>
      </table><br/><br/>
       <h3> Integrator::SFI&lt;Integrator::AllenCahn&gt;  </h3>
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
         <tr>
           <td style='padding-left: 10px'>
             hydro.eta_refinement_criterion
           </td>
           <td>
             eta-based refinement 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             hydro.omega_refinement_criterion
           </td>
           <td>
             vorticity-based refinement 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             hydro.gradu_refinement_criterion
           </td>
           <td>
             velocity gradient-based refinement 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             hydro.p_refinement_criterion
           </td>
           <td>
             pressure-based refinement 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             hydro.rho_refinement_criterion
           </td>
           <td>
             density-based refinement 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             hydro.gamma
           </td>
           <td>
             gamma for gamma law
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             hydro.cfl
           </td>
           <td>
             cfl condition
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             hydro.cfl_v
           </td>
           <td>
             cfl condition
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             hydro.mu
           </td>
           <td>
             linear viscosity coefficient
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             hydro.pref
           </td>
           <td>
             pp_query_default("Pfactor", value.Pfactor,1.0); // (to be removed) test factor for viscous source reference pressure for Roe solver
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             hydro.small
           </td>
           <td>
             small regularization value
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             hydro.cutoff
           </td>
           <td>
             cutoff value
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             hydro.lagrange
           </td>
           <td>
             lagrange no-penetration factor
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.eta.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.eta.ic.type
              = laminate
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.laminate.number_of_inclusions
             </td>
             <td>
               How many laminates (MUST be greater than or equal to 1). 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.laminate.orientation
             </td>
             <td>
               Vector normal to the interface of the laminate 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.laminate.eps
             </td>
             <td>
               Diffuse thickness 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.laminate.mollifier
             </td>
             <td>
               Type of mollifer to use (options: dirac, [gaussian])
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.laminate.singlefab
             </td>
             <td>
               Switch to mode where only one component is used. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.laminate.invert
             </td>
             <td>
               Take the complement of the laminate 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.eta.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.eta.ic.type
              = bmp
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.bmp.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.bmp.fit
             </td>
             <td>
               How to position image in space 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.bmp.coord.lo
             </td>
             <td>
               Location of lower-left corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.bmp.coord.hi
             </td>
             <td>
               Location of upper-right corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.bmp.channel
             </td>
             <td>
               Color channel to use 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.bmp.min
             </td>
             <td>
               Scaling value - minimum
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.bmp.max
             </td>
             <td>
               Scaling value - maximum
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.eta.ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.velocity.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.velocity.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.velocity.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.velocity.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.pressure.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.pressure.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.pressure.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.pressure.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.density.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.density.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.density.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.density.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.solid.momentum.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.solid.momentum.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.solid.momentum.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.solid.momentum.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.solid.density.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.solid.density.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.solid.density.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.solid.density.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.solid.energy.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.solid.energy.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.solid.energy.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.solid.energy.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.m0.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.m0.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.m0.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.m0.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.u0.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.u0.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.u0.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.u0.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.q.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.q.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.q.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.q.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.solver.type
              = roe
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.solver.roe.verbose
             </td>
             <td>
               enable to dump diagnostic data if the roe solver fails 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.solver.roe.entropy_fix
             </td>
             <td>
               apply entropy fix if tru 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             tstart
           </td>
           <td>
             time to activate hydro integrator 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             invert
           </td>
           <td>
             If true, set hydro_eta to 1-pf_eta 
           </td>
         </tr>
      </table><br/><br/>
       <h3> Integrator::SFI&lt;Integrator::Dendrite&gt;  </h3>
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
         <tr>
           <td style='padding-left: 10px'>
             hydro.eta_refinement_criterion
           </td>
           <td>
             eta-based refinement 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             hydro.omega_refinement_criterion
           </td>
           <td>
             vorticity-based refinement 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             hydro.gradu_refinement_criterion
           </td>
           <td>
             velocity gradient-based refinement 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             hydro.p_refinement_criterion
           </td>
           <td>
             pressure-based refinement 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             hydro.rho_refinement_criterion
           </td>
           <td>
             density-based refinement 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             hydro.gamma
           </td>
           <td>
             gamma for gamma law
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             hydro.cfl
           </td>
           <td>
             cfl condition
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             hydro.cfl_v
           </td>
           <td>
             cfl condition
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             hydro.mu
           </td>
           <td>
             linear viscosity coefficient
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             hydro.pref
           </td>
           <td>
             pp_query_default("Pfactor", value.Pfactor,1.0); // (to be removed) test factor for viscous source reference pressure for Roe solver
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             hydro.small
           </td>
           <td>
             small regularization value
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             hydro.cutoff
           </td>
           <td>
             cutoff value
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             hydro.lagrange
           </td>
           <td>
             lagrange no-penetration factor
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.eta.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.eta.ic.type
              = laminate
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.laminate.number_of_inclusions
             </td>
             <td>
               How many laminates (MUST be greater than or equal to 1). 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.laminate.orientation
             </td>
             <td>
               Vector normal to the interface of the laminate 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.laminate.eps
             </td>
             <td>
               Diffuse thickness 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.laminate.mollifier
             </td>
             <td>
               Type of mollifer to use (options: dirac, [gaussian])
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.laminate.singlefab
             </td>
             <td>
               Switch to mode where only one component is used. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.laminate.invert
             </td>
             <td>
               Take the complement of the laminate 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.eta.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.eta.ic.type
              = bmp
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.bmp.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.bmp.fit
             </td>
             <td>
               How to position image in space 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.bmp.coord.lo
             </td>
             <td>
               Location of lower-left corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.bmp.coord.hi
             </td>
             <td>
               Location of upper-right corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.bmp.channel
             </td>
             <td>
               Color channel to use 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.bmp.min
             </td>
             <td>
               Scaling value - minimum
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.bmp.max
             </td>
             <td>
               Scaling value - maximum
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.eta.ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.eta.ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.velocity.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.velocity.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.velocity.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.velocity.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.pressure.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.pressure.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.pressure.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.pressure.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.density.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.density.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.density.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.density.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.solid.momentum.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.solid.momentum.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.solid.momentum.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.solid.momentum.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.solid.density.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.solid.density.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.solid.density.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.solid.density.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.solid.energy.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.solid.energy.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.solid.energy.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.solid.energy.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.m0.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.m0.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.m0.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.m0.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.u0.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.u0.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.u0.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.u0.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.q.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.q.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.q.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.q.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              hydro.solver.type
              = roe
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.solver.roe.verbose
             </td>
             <td>
               enable to dump diagnostic data if the roe solver fails 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               hydro.solver.roe.entropy_fix
             </td>
             <td>
               apply entropy fix if tru 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px'>
             tstart
           </td>
           <td>
             time to activate hydro integrator 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             invert
           </td>
           <td>
             If true, set hydro_eta to 1-pf_eta 
           </td>
         </tr>
      </table><br/><br/>
