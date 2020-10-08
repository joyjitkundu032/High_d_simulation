# High_d_simulation
Simulates systems of hard spheres in arbitrary dimensions efficiently using the Wigner-Seitz cell of the checkerboard lattice  

This prog simulates a monodisperse system of hard spheres in arbitrary dimension d using standard Monte Carlo using a noncubic periodic boundary condition. Noncubic periodic boundary conditions are extremely efficient as they require a smaller number of particles to probe similar length scales if standard cubic boundary condition is used.

The code spits out particle trajectories at a fixed time interval.
