# Verification plots
First of all watch out for the timescale!

The [combined] folder contains plots of the numerical solution for different initial conditions and coefficients as well as the plots of the half-analytical solutions, which were calculated on the same time and space resolution as their numerical counterparts with the use of a numerical Dirac-delta function. The [t-cuts] folder contains 2D slices of both the analytical and the numerical solutions at a given timestep.

We can see a significant numerical diffusion in the convection example with lower resolutions, but around the 10^4 by 10^4 resolution, the difference between the numerical and analytical solution starts to diminish.

# The verification process
The notebook used for the verifications can be found in the main directory. It contains the analytical solutions and the commented code which calculates the analytical solutions.


