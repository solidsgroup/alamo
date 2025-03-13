thermoelastic
--------------------------


.. raw:: html

      <table class='api-inputs-table'>
      <thead><tr>
      <th>Input name</th>
      <th>Description</th>
      </tr></thead>
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
