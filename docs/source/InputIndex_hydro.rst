hydro
--------------------------


.. raw:: html

      <table class='api-inputs-table'>
      <thead><tr>
      <th>Input name</th>
      <th>Description</th>
      </tr></thead>
      </table><br/><br/>
       <h3> Integrator::Hydro  </h3>
      <table class='api-inputs-table'>
      <thead><tr>
      <th>Input name</th>
      <th>Description</th>
      </tr></thead>
         <tr>
           <td style='padding-left: 10px'>
             eta_refinement_criterion
           </td>
           <td>
             eta-based refinement 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             omega_refinement_criterion
           </td>
           <td>
             vorticity-based refinement 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             gradu_refinement_criterion
           </td>
           <td>
             velocity gradient-based refinement 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             p_refinement_criterion
           </td>
           <td>
             pressure-based refinement 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             rho_refinement_criterion
           </td>
           <td>
             density-based refinement 
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             gamma
           </td>
           <td>
             gamma for gamma law
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             cfl
           </td>
           <td>
             cfl condition
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             cfl_v
           </td>
           <td>
             cfl condition
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             mu
           </td>
           <td>
             linear viscosity coefficient
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             pref
           </td>
           <td>
             pp_query_default("Pfactor", value.Pfactor,1.0); // (to be removed) test factor for viscous source reference pressure for Roe solver
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             small
           </td>
           <td>
             small regularization value
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             cutoff
           </td>
           <td>
             cutoff value
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             lagrange
           </td>
           <td>
             lagrange no-penetration factor
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              eta.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              eta.ic.type
              = laminate
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.ic.laminate.number_of_inclusions
             </td>
             <td>
               How many laminates (MUST be greater than or equal to 1). 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.ic.laminate.orientation
             </td>
             <td>
               Vector normal to the interface of the laminate 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.ic.laminate.eps
             </td>
             <td>
               Diffuse thickness 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.ic.laminate.mollifier
             </td>
             <td>
               Type of mollifer to use (options: dirac, [gaussian])
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.ic.laminate.singlefab
             </td>
             <td>
               Switch to mode where only one component is used. 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.ic.laminate.invert
             </td>
             <td>
               Take the complement of the laminate 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              eta.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              eta.ic.type
              = bmp
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.ic.bmp.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.ic.bmp.fit
             </td>
             <td>
               How to position image in space 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.ic.bmp.coord.lo
             </td>
             <td>
               Location of lower-left corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.ic.bmp.coord.hi
             </td>
             <td>
               Location of upper-right corner in the domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.ic.bmp.channel
             </td>
             <td>
               Color channel to use 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.ic.bmp.min
             </td>
             <td>
               Scaling value - minimum
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.ic.bmp.max
             </td>
             <td>
               Scaling value - maximum
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              eta.ic.type
              = png
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.ic.png.channel
             </td>
             <td>
               Color channel to use (options: r, R, g, G, b, B, a, A)         
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.ic.png.filename
             </td>
             <td>
               BMP filename.
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.ic.png.fit
             </td>
             <td>
               how to position the image 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.ic.png.coord.lo
             </td>
             <td>
               Lower-left coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.ic.png.coord.hi
             </td>
             <td>
               Upper-right coordinates of image in domain
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.ic.png.min
             </td>
             <td>
               Desired minimum value to scale pixels by 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               eta.ic.png.max
             </td>
             <td>
               Desired maximum value to scale pixels by 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              velocity.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               velocity.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              velocity.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               velocity.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              pressure.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               pressure.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              pressure.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               pressure.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              density.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               density.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              density.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               density.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              solid.momentum.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               solid.momentum.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              solid.momentum.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               solid.momentum.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              solid.density.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               solid.density.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              solid.density.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               solid.density.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              solid.energy.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               solid.energy.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              solid.energy.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               solid.energy.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              m0.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               m0.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              m0.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               m0.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              u0.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               u0.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              u0.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               u0.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              q.ic.type
              = constant
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               q.ic.constant.value
             </td>
             <td>
               Array of constant values. The number of values should equal either 1 or N where N is the number of fab components 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              q.ic.type
              = expression
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               q.ic.expression.coord
             </td>
             <td>
               coordinate system to use: "cartesian" (for x,y,z,t) and  "polar" (for r, theta, z, t) 
             </td>
           </tr>
         <tr>
           <td style='padding-left: 10px' colspan=2>
              <b> if </b>
              solver.type
              = roe
           </td>
         </tr>
           <tr>
             <td style='padding-left: 20px'>
               solver.roe.verbose
             </td>
             <td>
               enable to dump diagnostic data if the roe solver fails 
             </td>
           </tr>
           <tr>
             <td style='padding-left: 20px'>
               solver.roe.entropy_fix
             </td>
             <td>
               apply entropy fix if tru 
             </td>
           </tr>
      </table><br/><br/>
