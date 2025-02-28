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
    affiliation: 1 #"1, 2" # (Multiple affiliations must be quoted)
affiliations:
 - name: Department of Aerospace Engineering, Iowa State University, Ames, IA, USA
   index: 1
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


# Overview of methods



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




## Multiphysics integrators


# Code infrastructure

Easy for a graduate student to program in.


## Parameter parsing

## Automatic documentation

## Automatic regression testing framework

## Guaranteed reproducibility

# References
