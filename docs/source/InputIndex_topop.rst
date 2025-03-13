topop
--------------------------


.. raw:: html

      <table class='api-inputs-table'>
      <thead><tr>
      <th>Input name</th>
      <th>Description</th>
      </tr></thead>
      </table><br/><br/>
       <h3> Integrator::TopOp&lt;Model::Solid::Linear::Isotropic&gt;  </h3>
      <table class='api-inputs-table'>
      <thead><tr>
      <th>Input name</th>
      <th>Description</th>
      </tr></thead>
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
           <td style='padding-left: 10px'>
             eta_ref_threshold
           </td>
           <td>
             Refinement threshold based on eta
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             alpha
           </td>
           <td>
             :math:`\alpha` parameter
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             beta
           </td>
           <td>
             :math:`\beta` parameter
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             gamma
           </td>
           <td>
             :math:`\gamma` parameter
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             volume0frac
           </td>
           <td>
             Prescribed volume fraction
           </td>
         </tr>
         <tr>
           <td style='padding-left: 10px'>
             volume0
           </td>
           <td>
             Prescribed total vlume
           </td>
         </tr>
      </table><br/><br/>
