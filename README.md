# Advection --> 3D Pressure Poisson Equation --> Channel Flow Solver

This project develops a fully spectral incompressible channel-flow solver by progressing from scalar advection, to a three-dimensional pressure Poisson formulation, and finally to a coupled pressure–momentum channel-flow DNS framework. The solver uses Fourier discretization in the periodic streamwise direction and Chebyshev polynomial expansion in the wall-normal direction to resolve velocity and pressure fields with high spectral accuracy.

The final channel-flow formulation advances the incompressible Navier–Stokes equations using explicit treatment of nonlinear terms and implicit treatment of viscous terms, while enforcing continuity and wall boundary conditions through a coupled pressure–momentum approach.

All 3 different solvers have been documented with code implemented in MATLAB.
