---
title: 'The Alamo multiphysics solver for phase field simulations with block structured adaptive mesh refinement'
tags:
  - c++
  - computational physics
  - materials science
  - phase field
  - adaptive mesh refinement
authors:
  - name: Brandon Runnels
    orcid: 0000-0003-3043-5227
affiliations:
 - name: Department of Aerospace Engineering, Iowa State University, Ames, IA, USA
date: 2025
bibliography: paper.bib
---

# Summary




# Statement of need

The phase field (PF) method is a powerful theoretical framework that enables the systematic description of complex physical systems.
PF methods have been successfully used to describe phenomena such as solidification, microstructure evolution, fracture, damage, dislocations, and many more.
Beyond materials science, PF methods have also enjoyed great success in other applications ranging from deflagration of solid rocket propellant to topology optimization.

The success of the PF method is derived from its implicit, diffuse representation of boundaries and surfaces, which avoids the need for cumbersome interface tracking.
However, the PF method also can incur great computational expense, due to the need for high grid resolution across the diffuse boundary.
In order for the PF method to be feasible, strategic algorithms are necessary in order to provide sufficient boundary resolution without wasting grid points on uninteresting regions.
Such algirthms fall typically into two main categories.
(1) Spectral methods solve the phase field equations in the frequency domain, e.g. [@kochmann2015phase], and 
(2) Real-space methods employing adaptive mesh refinement (AMR).
Spectral methods offer a number of performance advantages, especially when coupling to global mechanical solvers.
However, they can be limited in their ability to resolve fine-scale features, and can be very cumbersome to use when implementing novel types of models.
On the other hand, real-space methods with AMR are often able to attain very good performance, can be easily suited to the domain of interest, and provide an attractive platform for prototyping new physical models.

A number of open-source real-space PF codes exist and have enjoyed significant popularity.
Some of the most widely known are Moose, [@giudicelli2024moose], Fenics [@baratta2023dolfinx], and Prisms-PF [@dewitt2020prisms], which employ octree style AMR.
Though effective, octree-AMR can result on complex and expensive mesh management.
It is can also be challenging to achieve optimal load balancing, due to the high degree of unpredicitable connectivity within the octree mesh.

Block-structured AMR (BSAMR) is an alternative AMR strategy.
BSAMR divides the domains into distinct levels, with each level usually consisting of a collection of Cartesian grid regions (patches), that effectively evolve independently.
Communication between patches and levels is then handled through ghost cells, interpolation, and restriction.
This datastructure is extremely efficient and scalable, while also being highly amenable to efficient code prototyping.
Importantly, BSAMR also acts as a seamless extension to geometric multigrid, making naturally efficient at performing global implicit solves.
The AMReX framework [@zhang2019amrex] provides a powerful platform for development of BSAMR codes.
However, the use of AMReX (and BSAMR in general) has been limited in PF and solid mechanics, due to the inherent challenges of solving implicit mechanical equilibrium equations on a patch-based mesh.

The Alamo multiphysics solver leverages the power of BSAMR for phase-field problems.
Alamo provides a unique, strong-form finite-deformation, matrix-free mechanics solver, enabling the efficient solution of implicit solid mechanics calculations.
It also provides a set of numerical integration routines, myriad material models, and numerous examples covering a broad cross-section of PF modeling interests.


# Methods


## Mechanical solver

The Alamo mechanical solver extends the multi-level multi-grid (MLMG) solver to address the problem of quasi-static mechanical equilibrium, that is,
\begin{equation} \label{eq:mech_equil}
\operatorname{Div}\big( \mathrm{DW}(\mathbf{F})\big) + \mathbf{B} = \mathbf{0},
\end{equation}
where $\mathbf{F}$ is the deformation gradient, $\mathbf{B}$ is a body force, and $\mathrm{W}$ is an arbitrary Helmholtz free energy with derivatives $\mathbf{P}=\mathrm{DW}(\mathbf{F}) = d\mathrm{W}/d\mathbf{F}$ the Piola-Kirchhoff stress tensor, and $\mathbb{C}=\mathrm{DDW}(\mathbf{F})=d^2\mathrm{W}/d\mathbf{F}/d\mathbf{F}$ the tangent modulus.
Usually, the solution of \autoref{eq:mech_equil} is achieved using the finite element method (FEM).
However, FEM is not readily amenable to the BSAMR structure, due to the difficulties in achieving consistent shape functions between levels, and problems with achieiving consistent derivatives at the coarse/fine boundary without a global matrix.

The Alamo solver is developed to be native to the BSAMR framework, and takes full advantage of the integration with geometric multigrid.
It is matrix-free, which is necessary in order to avoid additional communication overhead.
It is strong-form, using finite differences instead of shape functions to calculate derivatives, ensuring consistency between levels and compatibility with restriction/prolongation operations.
It handles coarse-fine boundaries using a novel "reflux-free" method, which avoids the special treatment of boundaries by including an extra layer of smoothed ghost nodes.
Details on these aspects of the solver are available in [@runnels2021massively].

The Alamo mechanical solver is versatile, allowing any type of mechanical model to be used though the abstract solid model interface (per the norm for most FEM codes).
Alamo uses templated AMReX BaseFab structures to encapsulate model parameters, which enables uses to implement sophisticated models without any knowledge of the external Alamo/AMReX infrastructure.
An optional interface is provided that allows users to communicate model information to the Alamo I/O routines, in order to include internal solid mechanics variables (such as accumulated plastic slip) with the field output.
Because mechanical models are instantiated as field variables, they require arithmetic operations in order to be interpolated/averaged as the mesh is adapated.
Alamo implements a "model vector space" framework that provides an intuitive and efficient framework for automatic creation of arithmetic operations.
This allows model developers to implement all necessary arithmetic operators with only three additional lines of code, ensuring systematic compliance of the model and increased readability of the code.
The flexibility of this framework is exemplified through Alamo's library of solid models, which include implementations of solids ranging from linear elastic isotropic materials to finite deformation crystal plasticity.

Because the mechanical solver is coupled to problems defined with diffuse boundaries, often involving "void" regions in which there is no mechanical strength, additional steps are necessary to avoid convergence issues.
Alamo contains methods for accounting for diffuse boundary conditions, and uses a joint cell/node based interpolation scheme to ensure good convergence even when the solution is not uniquely defined.
Details on the near-singular solver capability, as well as the model vector space implementation, have been documented in [@agrawal2023robust].

## Multiple inheritance polymorphic integrators 

In Alamo, each type of physical behavior is encapsulated by an "Integrator" class.
Each integrator is responsible for solving a certain set of equations: the HeatConduction integrator solves the thermal diffusion equation; the AllenCahn integrator solves the Allen-Cahn equations; the Mechanics integrator solves the equations of mechanical equilibrium; etc.
All integrators inherit from the base integrator, which interfaces with AMReX and manages the creation/deletion/evolution of field variables.
This partitions the code so that each integrator is responsible only for its own physics, without excessive bookkeeping infrastructure.

Most phase field problems of interest feature the non-trivial interaction of multiple disparate physical behaviors.
For example, microstructure evolution in metals is strongly coupled to mechanical loading.
Or, some phase field models may require a fully resolved flow-field.
Such couplings can rapidly increase the complexity of the code base.

The Alamo solution is **multiple-inheritance polymorphic integrators** (MIPI). 
The MIPI schema encapsulates each physical system into its own integrator.
Single-application integrators inherit directly from the base integrator, whereas multi-physics application link other integrators together using multiple inheritance.

![The MIPI schema for solving thermoelasticity.\label{fig:example}](mipi.pdf){width=60%}

The MIPI method is exemplified through the case of coupled thermal evolution and mechanical equilibrium (\autoref{fig:example}).
In this example, the Mechanics integrator is only responsible for solving mechanical equilibrium and can be run independently.
Similarly, the Heat Conduction integrator is concerned only with thermal evolution; all physical aspects related to heat conduction are encapsulated here.
Both Mechanics and Heat Conduction inherit from the integrator base class, which tracks and updates all field variables on a shared grid.
To link these together, the ThermalElastic integrator inherits jointly from Mechancis and HeatConduction.
Its only function is to link the Mechanics strain field to the Heat Conduction temperature field; everything else is managed by the respective integrator.
This allows for integrators to be combined in arbitrary ways, with minimal (or possible no) instrumentation required within the linked integrator itself.

Numerous alamo integrators have been developed and used for scientific applications.
A brief summary is included here:

- Microstructure evolution is simulated using the multi-phase field method [@eren2022comparison] combined with the strong-form mechanics solver to simulate grain boundary anisotropy [@ribot2019new], phase field disconnections [@runnels2020phase,gokuli2021multiphase], and twin growth in magnesium [@hu2024atomistic].

- The phase field fracture mechanics model couples crack evolution and mechanics to capture crack propagation [@agrawal2021block,@agrawal2023robust].

- Deflagration of solid rocket propellant is captured using a phase field method [@kanagarajan2022diffuse], coupled to heat transfer [@meier2024diffuse] and the hyperelastic solver [@meier2024finite] to accurately predict ignition, deflagration, and mechanical response.

- Recent work demonstrates the ability of alamo to simulate a hydrodynamic compressible flow through a domain governed by phase field equations including the Allen Cahn equation (for porous media) and the dendrite growth equation [@boyd2025diffuse].

# Infrastructure

Alamo is intended to bridge the gap between *prototype/research codes* (often written with *ad-hoc* structure in a non-performant language with minimal documentation) and *production codes* (which are well-organized and documented, but have a steep learning curve and restrictive contribution requirements).
In other words, it is designed to be easy for a graduate student to learn, contribute to, and run without excessive time spent learning the infrastructure.
Automation is therefore the key to ensure that Alamo retains production-level quality without requiring excessive effort on the part of the developers.

## Parameter parsing and automatic documentation

Alamo employs a recursive object-oriented parameter parsing system.
Each class contains a standard `Parse` function to read its member variables, and is located as close as possible to the location in the code as where the variables are used.
Parser commands self-document the nature of each variable: for instance, `query_default` provides a default value, `query_required` triggers an error if the variable is not provided, `query_validate` provides a list of acceptable values, and so on.
This is effective at eliminating uninitialized variables, and localizes the variable's usage with its parsing statement.
All parameters are read in hierarchically, systematically using prefixes to avoid naming conflicts (especially in the case of MIPI integrators).

In addition, Alamo contains a set of python-based code scrapers to scan the source code for all input parameters.
The inputs are catalogued, along with their source code location and comment string, and formatted into the automatically-generated documentation.
This allows users to browse the inputs and link directly to their usage in the code (a common difficulty when encountering a new code), or to search an input index to determine how the inputs are used. 
It is also used to automatically generate input file builders that are guaranteed to be accurate and consistent with the current source code.

Continuous integration is used to require that all inputs be adequately documented.
Then, the documentation is automatically generated and posted online with every addition to the main development branches.
This ensures that all documentaiton is kept current with the source code, without requiring anything of the developer beyond a single comment string for each input.

## Automatic regression testing framework


coverage

## Guaranteed reproducibility

# Acknowledgements

(citations) 0

# References
