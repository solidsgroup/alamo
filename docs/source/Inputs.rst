
======
Inputs
======





BC
==


BC::Constant
************

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].type.xlo`
      - queryarr
      - 
    * - :code:`[prefix].type.xhi`
      - queryarr
      - 
    * - :code:`[prefix].type.ylo`
      - queryarr
      - 
    * - :code:`[prefix].type.yhi`
      - queryarr
      - 
    * - :code:`[prefix].type.zlo`
      - queryarr
      - 
    * - :code:`[prefix].type.zhi`
      - queryarr
      - 
    * - :code:`[prefix].val.xlo`
      - queryarr
      - 
    * - :code:`[prefix].val.xhi`
      - queryarr
      - 
    * - :code:`[prefix].val.ylo`
      - queryarr
      - 
    * - :code:`[prefix].val.yhi`
      - queryarr
      - 
    * - :code:`[prefix].val.zlo`
      - queryarr
      - 
    * - :code:`[prefix].val.zhi`
      - queryarr
      - 


BC::Step
********

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].h1`
      - query
      - 
    * - :code:`[prefix].h2`
      - query
      - 


BC::Operator
************


BC::Operator::Elastic
---------------------


BC::Operator::Elastic::Constant
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - 
         This BC defines boundary conditions that are constant with respect to space.
         However they may change in time using the :ref:`Numeric::Interpolator::Linear` method.
         
         In 2D, BCs are prescribed for each edge (4) and corner (4) on the boundary, for a total of 8 possible regions.
         In 3D, BCs are prescribed for each face (6), each edge (12), and each corner (8), for a total of 26 possible regions.
         Care must be taken to define BCs for edges/corners consistently with the faces/edges, or you will get poor 
         convergence and inaccurate behavior.
         (See :ref:`BC::Operator::Elastic::TensionTest` for a reduced BC with streamlined load options.)
        
         To define a boundary condition, you must define both a "type" and a "value" in each direction.
         The type can be "displacement", "disp", "neumann", "traction", "trac".
         The value is the corresponding value.
         All BCs are, by default, dirichlet (displacement) with a value of zero.
        
          
    * - :code:`[prefix].type.xloylozlo`
      - queryarr
      -  3D Corner
        
    * - :code:`[prefix].type.xloylozhi`
      - queryarr
      -  3D Corner
        
    * - :code:`[prefix].type.xloyhizlo`
      - queryarr
      -  3D Corner
        
    * - :code:`[prefix].type.xloyhizhi`
      - queryarr
      -  3D Corner
        
    * - :code:`[prefix].type.xhiylozlo`
      - queryarr
      -  3D Corner
        
    * - :code:`[prefix].type.xhiylozhi`
      - queryarr
      -  3D Corner
        
    * - :code:`[prefix].type.xhiyhizlo`
      - queryarr
      -  3D Corner
        
    * - :code:`[prefix].type.xhiyhizhi`
      - queryarr
      -  3D Corner
        
    * - :code:`[prefix].type.ylozlo`
      - queryarr
      -  3D Edge
        
    * - :code:`[prefix].type.ylozhi`
      - queryarr
      -  3D Edge
        
    * - :code:`[prefix].type.yhizlo`
      - queryarr
      -  3D Edge
        
    * - :code:`[prefix].type.yhizhi`
      - queryarr
      -  3D Edge
        
    * - :code:`[prefix].type.zloxlo`
      - queryarr
      -  3D Edge
        
    * - :code:`[prefix].type.zloxhi`
      - queryarr
      -  3D Edge
        
    * - :code:`[prefix].type.zhixlo`
      - queryarr
      -  3D Edge
        
    * - :code:`[prefix].type.zhixhi`
      - queryarr
      -  3D Edge
        
    * - :code:`[prefix].type.xloylo`
      - queryarr
      -  3D Edge / 2D Corner
        
    * - :code:`[prefix].type.xloyhi`
      - queryarr
      -  3D Edge / 2D Corner
        
    * - :code:`[prefix].type.xhiylo`
      - queryarr
      -  3D Edge / 2D Corner
        
    * - :code:`[prefix].type.xhiyhi`
      - queryarr
      -  3D Edge / 2D Corner
        
    * - :code:`[prefix].type.xlo`
      - queryarr
      -  3D Face  2D Edge
        
    * - :code:`[prefix].type.xhi`
      - queryarr
      -  3D Face  2D Edge
        
    * - :code:`[prefix].type.ylo`
      - queryarr
      -  3D Face  2D Edge
        
    * - :code:`[prefix].type.yhi`
      - queryarr
      -  3D Face  2D Edge
        
    * - :code:`[prefix].type.zlo`
      - queryarr
      -  3D Face
        
    * - :code:`[prefix].type.zhi`
      - queryarr
      -  3D Face
        
    * - :code:`[prefix].val.xloylozlo`
      - queryarr
      -  3D Corner
        
    * - :code:`[prefix].val.xloylozhi`
      - queryarr
      -  3D Corner
        
    * - :code:`[prefix].val.xloyhizlo`
      - queryarr
      -  3D Corner
        
    * - :code:`[prefix].val.xloyhizhi`
      - queryarr
      -  3D Corner
        
    * - :code:`[prefix].val.xhiylozlo`
      - queryarr
      -  3D Corner
        
    * - :code:`[prefix].val.xhiylozhi`
      - queryarr
      -  3D Corner
        
    * - :code:`[prefix].val.xhiyhizlo`
      - queryarr
      -  3D Corner
        
    * - :code:`[prefix].val.xhiyhizhi`
      - queryarr
      -  3D Corner
        
    * - :code:`[prefix].val.ylozlo`
      - queryarr
      -  3D Edge
        
    * - :code:`[prefix].val.ylozhi`
      - queryarr
      -  3D Edge
        
    * - :code:`[prefix].val.yhizlo`
      - queryarr
      -  3D Edge
        
    * - :code:`[prefix].val.yhizhi`
      - queryarr
      -  3D Edge
        
    * - :code:`[prefix].val.zloxlo`
      - queryarr
      -  3D Edge
        
    * - :code:`[prefix].val.zloxhi`
      - queryarr
      -  3D Edge
        
    * - :code:`[prefix].val.zhixlo`
      - queryarr
      -  3D Edge
        
    * - :code:`[prefix].val.zhixhi`
      - queryarr
      -  3D Edge
        
    * - :code:`[prefix].val.xloylo`
      - queryarr
      -  3D Edge / 2D Corner
        
    * - :code:`[prefix].val.xloyhi`
      - queryarr
      -  3D Edge / 2D Corner
        
    * - :code:`[prefix].val.xhiylo`
      - queryarr
      -  3D Edge / 2D Corner
        
    * - :code:`[prefix].val.xhiyhi`
      - queryarr
      -  3D Edge / 2D Corner
        
    * - :code:`[prefix].val.xlo`
      - queryarr
      -  3D Face / 2D Edge
        
    * - :code:`[prefix].val.xhi`
      - queryarr
      -  3D Face / 2D Edge
        
    * - :code:`[prefix].val.ylo`
      - queryarr
      -  3D Face / 2D Edge
        
    * - :code:`[prefix].val.yhi`
      - queryarr
      -  3D Face / 2D Edge
        
    * - :code:`[prefix].val.zlo`
      - queryarr
      -  3D Face
        
    * - :code:`[prefix].val.zhi`
      - queryarr
      -  3D Face
        


BC::Operator::Elastic::Expression
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - 
         This is basically the same as :ref:`BC::Operator::Elastic::Constant` except that
         you can use dynamically compiled expressions in space and time to define values.
        
         Usage is exactly the same except that the "val" inputs can depend on x, y, z, and t.
        
          
    * - :code:`[prefix].type.xloylozlo`
      - queryarr
      -  3D Corners
        
    * - :code:`[prefix].type.xloylozhi`
      - queryarr
      -  3D Corners
        
    * - :code:`[prefix].type.xloyhizlo`
      - queryarr
      -  3D Corners
        
    * - :code:`[prefix].type.xloyhizhi`
      - queryarr
      -  3D Corners
        
    * - :code:`[prefix].type.xhiylozlo`
      - queryarr
      -  3D Corners
        
    * - :code:`[prefix].type.xhiylozhi`
      - queryarr
      -  3D Corners
        
    * - :code:`[prefix].type.xhiyhizlo`
      - queryarr
      -  3D Corners
        
    * - :code:`[prefix].type.xhiyhizhi`
      - queryarr
      -  3D Corners
        
    * - :code:`[prefix].type.ylozlo`
      - queryarr
      -  3D Edges
        
    * - :code:`[prefix].type.ylozhi`
      - queryarr
      -  3D Edges
        
    * - :code:`[prefix].type.yhizlo`
      - queryarr
      -  3D Edges
        
    * - :code:`[prefix].type.yhizhi`
      - queryarr
      -  3D Edges
        
    * - :code:`[prefix].type.zloxlo`
      - queryarr
      -  3D Edges
        
    * - :code:`[prefix].type.zloxhi`
      - queryarr
      -  3D Edges
        
    * - :code:`[prefix].type.zhixlo`
      - queryarr
      -  3D Edges
        
    * - :code:`[prefix].type.zhixhi`
      - queryarr
      -  3D Edges
        
    * - :code:`[prefix].type.xloylo`
      - queryarr
      -  3D Edges / 2D Corners
        
    * - :code:`[prefix].type.xloyhi`
      - queryarr
      -  3D Edges / 2D Corners
        
    * - :code:`[prefix].type.xhiylo`
      - queryarr
      -  3D Edges / 2D Corners
        
    * - :code:`[prefix].type.xhiyhi`
      - queryarr
      -  3D Edges / 2D Corners
        
    * - :code:`[prefix].type.xlo`
      - queryarr
      -  3D Faces / 2D Edges
        
    * - :code:`[prefix].type.xhi`
      - queryarr
      -  3D Faces / 2D Edges
        
    * - :code:`[prefix].type.ylo`
      - queryarr
      -  3D Faces / 2D Edges
        
    * - :code:`[prefix].type.yhi`
      - queryarr
      -  3D Faces / 2D Edges
        
    * - :code:`[prefix].type.zlo`
      - queryarr
      -  3D Faces
        
    * - :code:`[prefix].type.zhi`
      - queryarr
      -  3D Faces
        


BC::Operator::Elastic::TensionTest
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].type`
      - query
      -  Tension test type. (Only option right now is :code:`uniaxial_stress_clamp`)
        
    * - :code:`[prefix].disp`
      - query
      -  Displacement of the test
        
    * - :code:`[prefix].trac`
      - query
      - 


IC
==


IC::Affine
**********

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`ic.n`
      - queryarr
      - 
    * - :code:`ic.alpha`
      - query
      - 


IC::BMP
*******

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].filename`
      - query
      -  BMP filename.
         Note that in GIMP, you must select "do not write color space information"
         and "24 bit R8 G8 B8" when exporting the BMP file.
        
    * - :code:`[prefix].fit`
      - query
      - 
    * - :code:`[prefix].channel`
      - query
      - 
    * - :code:`[prefix].min`
      - query
      - 
    * - :code:`[prefix].max`
      - query
      - 


IC::Constant
************

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].value`
      - queryarr
      -  Array of constant values. The number of values should equal either 1 or N where N is the number of fab components
        


IC::Cuboid
**********

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`ic.center`
      - queryarr
      -  Coordinates (X Y Z) of the center of the square/cube. 
        
    * - :code:`ic.length`
      - queryarr
      -  Lenth of the square/cube edges
        


IC::DoubleNotch
***************

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].thickness`
      - query
      - 
    * - :code:`[prefix].width`
      - query
      - 
    * - :code:`[prefix].x0`
      - queryarr
      - 
    * - :code:`[prefix].L`
      - query
      - 


IC::Ellipse
***********

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * -  If :code:`number_of_inclusions` is specified, then multiple ellipses are specified.
         In this case, each parameter must have number_of_inclusion*M values, where M is the
         number of values specified for the single ellipse case.
          
    * - :code:`[prefix].x0`
      - queryarr
      -  Coorinates of ellipse center
        
    * - :code:`[prefix].eps`
      - query
      -  Diffuse boundary thickness
        
    * - :code:`[prefix].A`
      - queryarr
      -  DxD square matrix defining an ellipse. 
        
    * - :code:`[prefix].a`
      - queryarr
      -  If :code:`A` is not defined, then assume a sphere with radius :code:`a`
        
    * - :code:`[prefix].number_of_inclusions`
      - query
      -  Number of ellipses
        
    * - :code:`[prefix].center`
      - queryarr
      - 
    * - :code:`[prefix].A`
      - queryarr
      - 
    * - :code:`[prefix].A`
      - queryarr
      - 
    * - :code:`[prefix].radius`
      - queryarr
      - 
    * - :code:`[prefix].eps`
      - queryarr
      - 
    * - :code:`[prefix].invert`
      - query
      -  Flip the inside and the outside
        


IC::Ellipsoid
*************

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * -  This IC initializes a single-component fab using the formula
         
         :math:`\phi = \begin{cases} \alpha_{in}  & (\mathbf{x}-\mathbf{x}_0)^T\mathbf{A}(\mathbf{x}-\mathbf{x}_0) \le r /                          \\ \alpha_{out} & \text{else} \end{cases}`
        
         The boundary can be mollified using an error function with parameter :math:`\varepsilon`
          
    * - :code:`[prefix].center`
      - queryarr
      -  Center of the ellipse :math:`\mathbf{x}_0`
        
    * - :code:`[prefix].A`
      - queryarr
      - 
    * - :code:`[prefix].radius`
      - queryarr
      - 
    * - :code:`[prefix].eps`
      - queryarr
      - 
    * - :code:`[prefix].in_value`
      - query
      - 
    * - :code:`[prefix].out_value`
      - query
      - 
    * - :code:`[prefix].mollifier`
      - query
      - 


IC::Expression
**************

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * -  This is an extremely general IC that uses JIT fparser to compile an
         IC based on a mathematical expression at runtime. Eventually we will
         probably replace all ICs with this one. 
         
         This can be used with an arbitrary number of components, named `region1`,
         `region2`, etc. You do not need to specify a number ahead of time.
         Each region should be a string that represents a function in terms of 
         x, y, z, and t. 
         It can be a boolean expression (returning 1 or 0) or it can return a 
         value. 
          
    * - :code:`[prefix].coord`
      - query
      - 


IC::Laminate
************

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].number_of_inclusions`
      - query
      -  How many laminates (MUST be greater than or equal to 1). Default = 1
        
    * - :code:`[prefix].center`
      - queryarr
      -  (x,y,[z]) values for the center point of the laminate
        
    * - :code:`[prefix].thickness`
      - queryarr
      -  thickness of the laminate
        
    * - :code:`[prefix].orientation`
      - queryarr
      -  Vector normal to the interface of the laminate
        
    * - :code:`[prefix].eps`
      - queryarr
      -  Diffuse thickness
        
    * - :code:`[prefix].mollifier`
      - query
      - 
    * - :code:`[prefix].singlefab`
      - query
      -  Switch to mode where only one component is used.
        
    * - :code:`[prefix].invert`
      - query
      -  Take the complement of the laminate
        


IC::Notch
*********

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].center`
      - queryarr
      - 
    * - :code:`[prefix].orientation`
      - queryarr
      - 
    * - :code:`[prefix].thickness`
      - queryarr
      - 
    * - :code:`[prefix].length`
      - queryarr
      - 
    * - :code:`[prefix].radius`
      - queryarr
      - 
    * - :code:`[prefix].eps`
      - query
      - 
    * - :code:`[prefix].mollifier`
      - query
      - 


IC::PS
******

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].nspheres`
      - query
      - 
    * - :code:`[prefix].matrix`
      - query
      - 
    * - :code:`[prefix].inclusion`
      - query
      - 


IC::PSRead
**********

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].eps`
      - query
      - 
    * - :code:`[prefix].filename`
      - query
      - 


IC::PerturbedInterface
**********************

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].wave_numbers`
      - queryarr
      - 
    * - :code:`[prefix].wave_amplitudes`
      - queryarr
      - 
    * - :code:`[prefix].normal`
      - query
      - 
    * - :code:`[prefix].offset`
      - query
      - 
    * - :code:`[prefix].mollifier`
      - query
      - 
    * - :code:`[prefix].eps`
      - query
      - 


IC::Sphere
**********

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * -  This is a somewhat antiquated IC that will eventually be replaced
         with the Expression IC.
          
    * - :code:`[prefix].radius`
      - query
      -  Radius of the sphere
        
    * - :code:`[prefix].center`
      - queryarr
      -  Vector location of the sphere center
        
    * - :code:`[prefix].inside`
      - query
      -  Value of the field inside the sphere
        
    * - :code:`[prefix].outside`
      - query
      -  Value of the field outside teh sphere
        
    * - :code:`[prefix].type`
      - query
      -  Type - can be cylinder oriented along the x, y, z directions or full sphere.
        


IC::TabulatedInterface
**********************

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].xs`
      - queryarr
      - 
    * - :code:`[prefix].ys`
      - queryarr
      - 


IC::Trig
********

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].nr`
      - queryarr
      - 
    * - :code:`[prefix].ni`
      - queryarr
      - 
    * - :code:`[prefix].dim`
      - query
      - 
    * - :code:`[prefix].alpha`
      - query
      - 


IC::Voronoi
***********

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].number_of_grains`
      - query
      - 
    * - :code:`[prefix].alpha`
      - queryarr
      - 
    * - :code:`[prefix].seed`
      - query
      - 


Integrator
==========


Integrator::CahnHilliard
************************


 \file CahnHilliard.H

Integrator::Flame
*****************

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description

    * - :code:`pf.eps`
      - query
      -  Burn width thickness
        
    * - :code:`pf.kappa`
      - query
      -  Interface energy param
        
    * - :code:`pf.gamma`
      - query
      -  Scaling factor for mobility
        
    * - :code:`pf.lambda`
      - query
      -  Chemical potential multiplier
        
    * - :code:`pf.w1`
      - query
      -  Unburned rest energy
        
    * - :code:`pf.w12`
      - query
      -  Barrier energy
        
    * - :code:`pf.w0`
      - query
      -  Burned rest energy
        
    * - :code:`pf.P`
      - query
      -  Pressure [UNITS?]
        
    * - :code:`pf.r_ap`
      - query
      -  AP Power law multiplier
        
    * - :code:`pf.n_ap`
      - query
      -  AP Power law exponent
        
    * - :code:`pf.r_htpb`
      - query
      -  HTPB Power law multiplier
        
    * - :code:`pf.n_htpb`
      - query
      -  HTPB Power law exponent
        
    * - :code:`pf.r_comb`
      - query
      -  Combination power law multiplier
        
    * - :code:`pf.n_comb`
      - query
      -  Combination power law exponent
        
    * - :code:`pf.eta.bc`
      - queryclass
      -  See :ref:`BC::Constant`
        
    * - :code:`eta.ic.type`
      - query
      -  IC type - [packedspheres,laminate] - see classes for more information
        
    * - :code:`thermal.on`
      - query
      -  These parameters are for the **Thermal transport model**
         Whether to use the thermal model
        
    * - :code:`thermal.rho1`
      - query
      -  Density (before)
        
    * - :code:`thermal.rho0`
      - query
      -  Density (after)
        
    * - :code:`thermal.ka`
      - query
      -  Thermal conductivity (before and after)
        
    * - :code:`thermal.kh`
      - query
      -  Thermal conductivity (before and after)
        
    * - :code:`thermal.k0`
      - query
      -  Thermal conductivity (before and after)
        
    * - :code:`thermal.cp1`
      - query
      -  Specific heat (before and after)
        
    * - :code:`thermal.cp0`
      - query
      -  Specific heat (before and after)
        
    * - :code:`thermal.delA`
      - query
      -  Thermal flux of each material
        
    * - :code:`thermal.delH`
      - query
      -  Thermal flux of each material
        
    * - :code:`amr.refinement_criterion`
      - query
      -  Refinement criterion for eta field
        
    * - :code:`amr.refinement_criterion_temp`
      - query
      -  Refinement criterion for temperature field
        
    * - :code:`amr.refinament_restriction`
      - query
      -  Eta value to restrict the refinament for the temperature field
        
    * - :code:`phi.ic.type`
      - query
      -  IC type (psread, laminate, constant)
        
    * - :code:`model_ap`
      - queryclass
      - 
    * - :code:`model_htpb`
      - queryclass
      - 


Integrator::Fracture
********************

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description

    * -  IC for crack field
          
    * - :code:`crack.ic.type`
      - queryarr
      - 
    * - :code:`crack.ic.notch`
      - queryclass
      - 
    * - :code:`crack.ic.ellipsoid`
      - queryclass
      - 

    * - :code:`elastic.solver`
      - queryclass
      - Solver::Nonlocal::Linear<brittle_fracture_model_type_test>  solver(op_b);
        


Integrator::HeatConduction
**************************

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].heat.alpha`
      - query
      - 
    * - :code:`[prefix].heat.refinement_threshold`
      - query
      - 
    * - :code:`[prefix].ic.type`
      - query
      - 
    * - :code:`[prefix].bc.temp`
      - queryclass
      - 


Integrator::Integrator
**********************

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * -  These are basic parameters that are, in 
         general, common to all Alamo simulations.
          
    * - :code:`max_step`
      - query
      -  Number of iterations before ending
        
    * - :code:`stop_time`
      - query
      -  Simulation time before ending
        
    * - :code:`timestep`
      - query
      -  Nominal timestep on amrlev = 0
        
    * - :code:`restart`
      - query
      -  Name of restart file to READ from
        
    * - :code:`restart_cell`
      - query
      -  Name of cell-fab restart file to read from
        
    * - :code:`restart_node`
      - query
      -  Name of node-fab restart file to read from
        
    * - :code:`ignore`
      - queryarr
      -  Space-separated list of entries to ignore
        

    * -  These are parameters that are specific to
         the AMR/regridding part of the code.
          
    * - :code:`amr.regrid_int`
      - query
      -  Regridding interval in step numbers
        
    * - :code:`amr.base_regrid_int`
      - query
      -  Regridding interval based on coarse level only
        
    * - :code:`amr.plot_int`
      - query
      -  Interval (in timesteps) between plotfiles
        
    * - :code:`amr.plot_dt`
      - query
      -  Interval (in simulation time) between plotfiles
        
    * - :code:`amr.plot_file`
      - query
      -  Output file
        
    * - :code:`amr.cell.all`
      - query
      -  Turn on to write all output in cell fabs (default: off)
        
    * - :code:`amr.cell.any`
      - query
      -  Turn off to prevent any cell based output (default: on)
        
    * - :code:`amr.node.all`
      - query
      -  Turn on to write all output in node fabs (default: off)
        
    * - :code:`amr.node.any`
      - query
      -  Turn off to prevent any node based output (default: on)
        
    * - :code:`amr.max_plot_level`
      - query
      -  Specify a maximum level of refinement for output files
        
    * - :code:`amr.nsubsteps`
      - queryarr
      -  Number of substeps to take on each level (default: 2)
        
    * - :code:`amr.nsubsteps`
      - query
      -  Number of substeps to take on each level (set all levels to this value)
        

    * -  Information on how to generate thermodynamic
         data (to show up in thermo.dat)
          
    * - :code:`amr.thermo.int`
      - query
      -  Default: integrate every time.
         Integration interval (1)
        
    * - :code:`amr.thermo.plot_int`
      - query
      -  Interval (in timesteps) between writing
        
    * - :code:`amr.thermo.plot_dt`
      - query
      -  Interval (in simulation time) between writing
        

    * -  Instead of using AMR, prescribe an explicit, user-defined
         set of grids to work on. This is pretty much always used
         for testing purposes only.
          
    * - :code:`explicitmesh.on`
      - query
      -  Use explicit mesh instead of AMR
        


Integrator::Mechanics
*********************

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * -  The mechanics integrator manages the solution of an elastic 
         solve using the MLMG solver. 
          
    * - :code:`[prefix].nmodels`
      - query
      -  Number of elastic model varieties
        
    * - :code:`[prefix].eta_ref_threshold`
      - query
      -  Refinement threshold for eta field
        
    * - :code:`[prefix].ref_threshold`
      - query
      -  Refinement threshold for strain gradient
        
    * - :code:`[prefix].ic.type`
      - query
      -  Read IC type for the eta field
        
    * - :code:`[prefix].psi.ic.type`
      - query
      -  Read IC type for the eta field
        


Integrator::PhaseFieldMicrostructure
************************************


 \file PhaseFieldMicrostructure.H

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].pf.number_of_grains`
      - query
      -  Number of grain fields (may be more if using different IC)
        
    * - :code:`[prefix].pf.M`
      - query
      -  Mobility
        
    * - :code:`[prefix].pf.mu`
      - query
      -  Phase field :math:`\mu`
        
    * - :code:`[prefix].pf.gamma`
      - query
      -  Phase field :math:`\gamma`
        
    * - :code:`[prefix].pf.sigma0`
      - query
      -  Initial GB energy if not using GB anisotropy
        
    * - :code:`[prefix].pf.l_gb`
      - query
      -  Mobility
        
    * - :code:`[prefix].pf.elastic_mult`
      - query
      -  Multiplier of elastic energy
        
    * - :code:`[prefix].pf.elastic_threshold`
      - query
      -  Elastic threshold (:math:`\phi_0`)
        
    * - :code:`[prefix].amr.max_level`
      - query
      -  Maximum AMR level
        
    * - :code:`[prefix].amr.ref_threshold`
      - query
      -  Phase field refinement threshold
        
    * - :code:`[prefix].mechanics.tstart`
      - query
      -  Elasticity
        
    * - :code:`[prefix].mechanics.model`
      - queryclass
      -  Type of model to use
        
    * - :code:`[prefix].lagrange.on`
      - query
      -  Lagrange multiplier method for enforcing volumes
         Turn on
        
    * - :code:`[prefix].lagrange.lambda`
      - query
      -  Lagrange multiplier value
        
    * - :code:`[prefix].lagrange.vol0`
      - query
      -  Prescribed volume
        
    * - :code:`[prefix].lagrange.tstart`
      - query
      -  Time to start enforcing Lagrange multipler
        
    * - :code:`[prefix].anisotropy.on`
      - query
      -  Anisotropic grain boundary energy parameters
         Turn on
        
    * - :code:`[prefix].anisotropy.beta`
      - query
      -  Regularization param 
        
    * - :code:`[prefix].anisotropy.tstart`
      - query
      -  Time to turn on anisotropy
        
    * - :code:`[prefix].anisotropy.timestep`
      - query
      -  Modify timestep when turned on
        
    * - :code:`[prefix].anisotropy.plot_int`
      - query
      -  Modify plot_int when turned on
        
    * - :code:`[prefix].anisotropy.plot_dt`
      - query
      -  Modify plot_dt when turned on
        
    * - :code:`[prefix].anisotropy.thermo_int`
      - query
      -  Modify thermo int when turned on
        
    * - :code:`[prefix].anisotropy.thermo_plot_int`
      - query
      -  Modify thermo plot int when turned on
        
    * - :code:`[prefix].anisotropy.elastic_int`
      - query
      -  Frequency of elastic calculation
        
    * - :code:`[prefix].anisotropy.regularization`
      - query
      -  Type of regularization to use  
        
    * - :code:`[prefix].anisotropy.gb_type`
      - query
      -  Set the anisotropic GB model
         Type of GB energy to use
        
    * - :code:`[prefix].bc.eta.type`
      - query
      -  Type (constnat)
        
    * - :code:`[prefix].ic.type`
      - query
      -  IC Type
        


Integrator::PolymerDegradation
******************************


 \file PolymerDegradation.H

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * -  Water diffusion
          
    * - :code:`water.on`
      - query
      - 
    * - :code:`water.diffusivity`
      - query
      -  Diffusivity 
        
    * - :code:`water.refinement_threshold`
      - query
      -  AMR refinement criterion
        
    * - :code:`water.ic_type`
      - query
      -  
        

    * - :code:`water.ic.value`
      - queryarr
      - 
    * - :code:`water.ic.bc`
      - queryclass
      - 



Integrator::SutureCrack
***********************

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description

    * -  IC for crack field
          
    * - :code:`crack.ic.type`
      - query
      - 
    * - :code:`crack.ic.notch`
      - queryclass
      - 


Integrator::ThermoElastic
*************************

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].hc`
      - queryclass
      - 
    * - :code:`[prefix].el`
      - queryclass
      - 
    * - :code:`[prefix].alpha`
      - queryarr
      - 


Integrator::TopOp
*****************

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * -  The mechanics integrator manages the solution of an elastic 
         solve using the MLMG solver. 
          
    * - :code:`[prefix].model`
      - queryclass
      - 
    * - :code:`[prefix].psi.ic.type`
      - query
      -  Read IC type for the eta field
        
    * - :code:`[prefix].eta_ref_threshold`
      - query
      - 
    * - :code:`[prefix].alpha`
      - query
      - 
    * - :code:`[prefix].beta`
      - query
      - 
    * - :code:`[prefix].gamma`
      - query
      - 
    * - :code:`[prefix].L`
      - queryclass
      - 
    * - :code:`[prefix].volume0`
      - query
      - 
    * - :code:`[prefix].lambda`
      - query
      - 


Integrator::Base
****************


Integrator::Base::Mechanics
---------------------------

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * -  The mechanics integrator manages the solution of an elastic 
         solve using the MLMG solver. 
          
    * - :code:`[prefix].type`
      - query
      - 
    * - :code:`[prefix].time_evolving`
      - query
      - 
    * - :code:`[prefix].solver`
      - queryclass
      -  Read parameters for :ref:`Solver::Nonlocal::Newton` solver
        
    * - :code:`[prefix].viscous.mu`
      - query
      - 
    * - :code:`[prefix].bc.type`
      - query
      -  Determine the boundary condition type (contant, tension_test, expression)
        
    * - :code:`[prefix].print_model`
      - query
      - 
    * - :code:`[prefix].rhs.type`
      - query
      -  Initializer for RHS
        
    * - :code:`[prefix].interval`
      - query
      -  Timestep interval for elastic solves (default - solve every time)
        
    * - :code:`[prefix].max_coarsening_level`
      - query
      -  Maximum multigrid coarsening level (default - none, maximum coarsening)
        
    * - :code:`[prefix].print_residual`
      - query
      - 
    * - :code:`[prefix].elastic_ref_threshold`
      - query
      -  Whether to refine based on elastic solution
        
    * - :code:`[prefix].zero_out_displacement`
      - query
      -  Set this to true to zero out the displacement before each solve.
         (This is a temporary fix - we need to figure out why this is needed.)
        



Model
=====


Model::Model
************


In Alamo, any type of constitutive behavior is encapsulated in a "Model".
There is a broad range of model types, and so there is no abstract Model class.
Rather, the base Model classes are stored in subsets of Model.

Model::Solid
************


Model::Solid::Solid
-------------------


Solid models are used with the :ref:`Integrator::Mechanics` integrator, which
implements the :ref:`Solver::Nonlocal::Newton` elasticity solver.
All solid models inherit from the :code:`Model::Solid` base class, which requires
all of the necessary components to be used in a mechanics solve.
Model classes have basically two functions:

#. Provide energy (W), stress (DW), and modulus (DDW) based on a kinematic variable
#. Evolve internal variables in time.

Model::Solid::Linear
--------------------


Model::Solid::Linear::Cubic
~~~~~~~~~~~~~~~~~~~~~~~~~~~


This class implements basic cubic elasticity.
For a discussion on cubic elasticity, `please see this link <http://solidmechanics.org/text/Chapter3_2/Chapter3_2.htm#Sect3_2_16>`_.

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].C11`
      - query
      -  Elastic constant (default: 1.68)
        
    * - :code:`[prefix].C12`
      - query
      -  Elastic constant (default: 1.21)
        
    * - :code:`[prefix].C44`
      - query
      -  Elastic constant (default: 0.75)
        
    * - :code:`[prefix].phi1`
      - query
      -  Bunge Euler angle :math:`\phi_1`
        
    * - :code:`[prefix].Phi`
      - query
      -  Bunge Euler angle :math:`\Phi`
        
    * - :code:`[prefix].phi2`
      - query
      -  Bunge Euler angle :math:`\phi_2`
        


Model::Solid::Linear::Isotropic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].planestress`
      - query
      - 
    * - :code:`[prefix].lame`
      - query
      - 
    * - :code:`[prefix].shear`
      - query
      - 
    * - :code:`[prefix].lambda`
      - query
      - 
    * - :code:`[prefix].mu`
      - query
      - 
    * - :code:`[prefix].E`
      - query
      - 
    * - :code:`[prefix].nu`
      - query
      - 


Model::Solid::Linear::IsotropicDegradable
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].lambda`
      - query
      - 
    * - :code:`[prefix].mu`
      - query
      - 
    * - :code:`[prefix].E`
      - query
      - 
    * - :code:`[prefix].nu`
      - query
      - 


Model::Solid::Linear::IsotropicDegradableTanh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].E1`
      - query
      - 
    * - :code:`[prefix].E2`
      - query
      - 
    * - :code:`[prefix].Tg`
      - query
      - 
    * - :code:`[prefix].Ts`
      - query
      - 
    * - :code:`[prefix].nu`
      - query
      - 
    * - :code:`[prefix].temp`
      - query
      - 


Model::Solid::Linear::Laplacian
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description


Model::Solid::Affine
--------------------


Model::Solid::Affine::Affine
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


"Affine" generally means "linear with an offset". Here we use "affine" to
refer to models that are elastic with an eigenstrain, i.e.

.. math::

   \sigma = \mathbb{C}(\varepsilon - \varepsilon_0)

The quantity :math:`\varepsilon_0` is any kind of eigenstrain - examples 
include thermal strain, plastic strain, or viscous strain.
This class can be used directly, where the eigenstrain is read in or
set by the integrator.
There are also classes (particularly visco/plastic models) that inherit
from this type of model.

Model::Solid::Affine::Cubic
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * -  This class extends :ref:`Model::Solid::Linear::Cubic` by adding
         an eigenstrain. (See the Linear::Cubic class for more inputs for this model)
          
    * - :code:`[prefix].F0`
      - queryarr
      -  Eigenstrain
        


Model::Solid::Affine::CubicDegradable
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * -  This class inherits from :ref:`Model::Solid::Affine::Cubic`.
         It provides the ability to "degrade" while retaining information
         about its original, pristine state
          
    * - :code:`[prefix].C11`
      - query
      -  Original, undegraded :math:`\mathbb{C}_{11}`
        
    * - :code:`[prefix].C12`
      - query
      -  Original, undegraded :math:`\mathbb{C}_{12}`
        
    * - :code:`[prefix].C44`
      - query
      -  Original, undegraded :math:`\mathbb{C}_{44}`
        
    * - :code:`[prefix].phi1`
      - query
      -  Bunge Euler angles :math:`\phi_1`
        
    * - :code:`[prefix].Phi`
      - query
      -  Bunge Euler angles :math:`\Phi`
        
    * - :code:`[prefix].phi2`
      - query
      -  Bunge Euler angles :math:`\phi_2`
        


Model::Solid::Affine::Isotropic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].lame`
      - query
      - 
    * - :code:`[prefix].shear`
      - query
      - 
    * - :code:`[prefix].E`
      - query
      - 
    * - :code:`[prefix].nu`
      - query
      - 
    * - :code:`[prefix].F0`
      - queryarr
      - 


Model::Solid::Affine::IsotropicDegradable
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].lame`
      - query
      - 
    * - :code:`[prefix].shear`
      - query
      - 
    * - :code:`[prefix].E`
      - query
      - 
    * - :code:`[prefix].nu`
      - query
      - 


Model::Solid::Affine::J2
~~~~~~~~~~~~~~~~~~~~~~~~


This models an isotropic elastic-perfectly-plastic, non-time-dependent solid model.

The energy and derivatives are:

.. math::
   :nowrap:

   \begin{gather}
   W = \frac{1}{2}(\varepsilon - \varepsilon_p):\mathbb{C}(\varepsilon-\varepsilon_p) \\
   DW = \mathbb{C}(\varepsilon-\varepsilon_p) \\
   DDW = \mathbb{C}
   \end{gather}

where :math:`\mathbb{C}` is an isotropic :code:`Set::Matrix4` and :math:`\varepsilon_p` is 
is stored in the :code:`F0` eigenstrain.

The plastic strain is evolved according to the following:

#. Calculate the deviatoric stress :math:`\sigma_v=\sigma - \frac{1}{3}tr(\sigma)\mathbf{I}`
#. Calculate :math:`J_2=\sqrt{\frac{3}{2}\sigma_v:\sigma_v}`
#. If :math:`J_2<\sigma_0` then quit - no plasticity occurs
#. Calculate :math:`\Delta\sigma = (1-\frac{\sigma_0}{J_2})`, which projects the stress 
   back on the yield surface.
#. Convert to change in plastic strain, :math:`\Delta\varepsilon=\mathbb{C}^{-1}\Delta\sigma`
#. Update plastic strain: :math:`\varepsilon_p += \Delta\varepsilon`

Notes:

* This does not implement any kind of hardening model. Rate hardening, isotropic hardening,
  and kinematic hardening have yet to be implemneted.

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].sigma0`
      - query
      -  J2 Yield criterion
        


Model::Solid::Affine::J2Plastic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].lambda`
      - query
      - 
    * - :code:`[prefix].mu`
      - query
      - 
    * - :code:`[prefix].E`
      - query
      - 
    * - :code:`[prefix].nu`
      - query
      - 
    * - :code:`[prefix].yield`
      - query
      - 
    * - :code:`[prefix].hardening`
      - query
      - 
    * - :code:`[prefix].theta`
      - query
      - 


Model::Solid::Affine::J2PlasticDegradable
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].lambda`
      - query
      - 
    * - :code:`[prefix].mu`
      - query
      - 
    * - :code:`[prefix].E`
      - query
      - 
    * - :code:`[prefix].nu`
      - query
      - 
    * - :code:`[prefix].yield`
      - query
      - 
    * - :code:`[prefix].hardening`
      - query
      - 
    * - :code:`[prefix].theta`
      - query
      - 


Model::Solid::Elastic
---------------------


Model::Solid::Elastic::NeoHookean
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].mu`
      - query
      - 
    * - :code:`[prefix].kappa`
      - query
      - 


Model::Interface
****************


Model::Interface::GB
--------------------


Model::Interface::GB::AbsSin
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].theta0`
      - query
      - 
    * - :code:`[prefix].sigma0`
      - query
      -  convert degrees into radians
        
    * - :code:`[prefix].sigma1`
      - query
      - 


Model::Interface::GB::Read
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].filename`
      - query
      - 


Model::Interface::GB::SH
~~~~~~~~~~~~~~~~~~~~~~~~

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].theta0`
      - query
      - 
    * - :code:`[prefix].phi0`
      - query
      -  convert degrees into radians
        
    * - :code:`[prefix].sigma0`
      - query
      -  convert degrees into radians
        
    * - :code:`[prefix].sigma1`
      - query
      - 
    * - :code:`[prefix].regularization`
      - query
      - 


Model::Interface::GB::Sin
~~~~~~~~~~~~~~~~~~~~~~~~~

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].theta0`
      - query
      - 
    * - :code:`[prefix].sigma0`
      - query
      -  convert degrees into radians
        
    * - :code:`[prefix].sigma1`
      - query
      - 


Model::Interface::Crack
-----------------------


Model::Interface::Crack::Constant
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].G_c`
      - query
      - 
    * - :code:`[prefix].zeta`
      - query
      - 
    * - :code:`[prefix].mobility`
      - query
      - 
    * - :code:`[prefix].threshold`
      - query
      - 
    * - :code:`[prefix].gtype`
      - query
      - 
    * - :code:`[prefix].wtype`
      - query
      - 
    * - :code:`[prefix].exponent`
      - query
      - 


Model::Interface::Crack::Sin
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].Gc0`
      - query
      - 
    * - :code:`[prefix].Gc1`
      - query
      - 
    * - :code:`[prefix].theta0`
      - query
      - 
    * - :code:`[prefix].zeta`
      - query
      - 
    * - :code:`[prefix].mobility`
      - query
      - 
    * - :code:`[prefix].threshold`
      - query
      - 
    * - :code:`[prefix].gtype`
      - query
      - 
    * - :code:`[prefix].wtype`
      - query
      - 
    * - :code:`[prefix].exponent`
      - query
      - 


Numeric
=======


Numeric::Interpolator
*********************


Numeric::Interpolator::Linear
-----------------------------

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].str`
      - query
      - 


Util
====


Util::Util
**********

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`plot_file`
      - query
      - 


Solver
======


Solver::Nonlocal
****************


Solver::Nonlocal::Linear
------------------------

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * -  These are the parameters that are read in for a standard 
         multigrid linear solve.
          
    * - :code:`[prefix].max_iter`
      - query
      -  Max number of iterations to perform before erroring out
        
    * - :code:`[prefix].bottom_max_iter`
      - query
      -  Max number of iterations on the bottom solver
        
    * - :code:`[prefix].max_fmg_iter`
      - query
      -  Max number of F-cycle iterations to perform
        
    * - :code:`[prefix].fixed_iter`
      - query
      -  Number of fixed iterations to perform before exiting gracefully
        
    * - :code:`[prefix].verbose`
      - query
      -  Verbosity of the solver (1-5)
        
    * - :code:`[prefix].pre_smooth`
      - query
      -  Number of smoothing operations before bottom solve (2)
        
    * - :code:`[prefix].post_smooth`
      - query
      -  Number of smoothing operations after bottom solve (2)
        
    * - :code:`[prefix].final_smooth`
      - query
      -  Number of final smoothing operations when smoother is used as bottom solver (8)
        
    * - :code:`[prefix].bottom_smooth`
      - query
      -  Additional smoothing after bottom CG solver (0)
        
    * - :code:`[prefix].bottom_solver`
      - query
      -  The method that is used for the multigrid bottom solve (cg, bicgstab, smoother)
        
    * - :code:`[prefix].bottom_tol_rel`
      - query
      -  Relative tolerance on bottom solver
        
    * - :code:`[prefix].bottom_tol_abs`
      - query
      -  Absolute tolerance on bottom solver
        
    * - :code:`[prefix].tol_rel`
      - query
      -  Relative tolerance
        
    * - :code:`[prefix].tol_abs`
      - query
      -  Absolute tolerance
        
    * - :code:`[prefix].omega`
      - query
      -  Omega (used in gauss-seidel solver)
        
    * - :code:`[prefix].average_down_coeffs`
      - query
      -  Whether to average down coefficients or use the ones given.
         (Setting this to true is important for fracture.)
        
    * - :code:`[prefix].normalize_ddw`
      - query
      -  Whether to normalize DDW when calculating the diagonal.
         This is primarily used when DDW is near-singular - like when there
         is a "void" region or when doing phase field fracture.
        


Solver::Nonlocal::Newton
------------------------

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * -  These paramters control a standard Newton-Raphson solve.
         
         **Note**: 
         This class inherits all of the linear solve paramters
         from its parent class, :ref:`Solver::Nonlocal::Linear`
          
    * - :code:`[prefix].nriters`
      - query
      -  Number of newton-raphson iterations.
        
    * - :code:`[prefix].nrtolerance`
      - query
      - 


IO
==


IO::ParmParse
*************

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description


IO::WriteMetaData
*****************

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description


Operator
========


Operator::Elastic
*****************

.. flat-table:: 
    :widths: 20 10 70
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`[prefix].small`
      - query
      - 


