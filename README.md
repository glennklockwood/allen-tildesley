Allen & Tildesley Sample Programs
================================================================================

These are the sample codes that supplemented "Computer Simulation of Liquids" by
M. P. Allen and D. Tildesley, originally published by Oxford University Press in
1987.  They demonstrate many key algorithms that are foundational to molecular
dynamics simulation of periodic systems.

A&T Code | Source Name         |  Descriptions
---------|---------------------|----------------------------------------
F.01     | pbc.f               | Periodic boundary conditions in various geometries.
F.02     | gear.f              | 5-value Gear predictor-corrector algorithm.
F.03     | leapfrog.f          | Low-storage MD programs using leapfrog Verlet algorithm
F.04     | velocityverlet.f    | Velocity version of Verlet algorithm
F.05     | quaternion.f        | Quaternion parameter predictor-corrector algorithm
F.06     | leapfrog-rotation.f | Leapfrog algorithms for rotational motion
F.07     | constraint.f        | Constraint dynamics for a nonlinear triatomic molecule
F.08     | shake.f             | Shake algorithm for constraint dynamics of a chain molecule
F.09     | rattle.f            | Rattle algorithm for constraint dynamics of a chain molecule
F.10     | hardsphere.f        | Hard sphere molecular dynamics program
F.11     | mc-nvt.f            | Constant-NVT Monte Carlo for Lennard-Jones atoms
F.12     | mc-npt.f            | Constant-NPT Monte Carlo algorithm
F.13     | mc-muvt.f           | The heart of a constant muVT Monte Carlo program
F.14     | mc-muvt-indices.f   | Algorithm to handle indices in constant muVT Monte Carlo
F.15     | random-rotate.f     | Routines to randomly rotate molecules
F.16     | mc-dumbbell.f       | Hard dumb-bell Monte Carlo program
F.17     | lj.f                | A simple Lennard-Jones force routine
F.18     | avoid-sqrt.f        | Algorithm for avoiding the square root operation
F.19     | verlet-list.f       | The Verlet neighbour list
F.20     | link-cell.f         | Routines to construct and use cell linked-list method
F.21     | md-multitimestep.f  | Multiple timestep molecular dynamics
F.22     | ewald.f             | Routines to perform the Ewald sum
F.23     | alpha-fcc.f         | Routine to set up alpha fcc lattice of linear molecules
F.24     | init-velocity.f     | Initial velocity distribution
F.25     | order-parameter.f   | Routine to calculate translational order parameter
F.26     | unfold-pbc.f        | Routines to fold unfold trajectories in periodic boundaries
F.27     | time-correlation.f  | Program to compute time correlation functions
F.28     | md-nvt-extended.f   | Constant-NVT molecular dynamics - extended system method
F.29     | md-nvt-constraint.f | Constant-NVT molecular dynamics - constraint method
F.30     | md-nph-extended.f   | Constant-NPH molecular dynamics - extended system method
F.31     | md-npt-constraint.f | Constant-NPT molecular dynamics - constraint method
F.32     | link-cell-sheared.f | Cell linked-lists in sheared boundaries
F.33     | brownian.f          | Brownian dynamics for a Lennard-Jones fluid
F.34     | clustering.f        | An efficient clustering routine
F.35     | voronai.f           | The Voronoi construction in 2d and 3d
F.36     | mc-hardlines.f      | Monte Carlo simulation of hard lines in 2d
F.37     | fft.f               | Routines to calculate Fourier transforms
