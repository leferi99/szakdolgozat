# Verification plots

This folder contains plots of the numerical solution for different initial conditions and coefficients as well as the plots of the half-analytical solutions, which were calculated on the same time and space resolution as their numerical counterparts with the use of a numerical Dirac-delta function. 
Given the not fully analytical solution, the time and space resolution have a significant effect on both the numerical and the half-analytical solutions. 
In the lower resolution numerical solutions we can see a significant numerical diffusion as well, but around 10^4 by 10^4 resolutions, the difference between the numerical and analytical solution diminishes.

# Naming
If the name contains r and t, the following numbers note the spatial resolution and the time resolution respectively. In every other case the time and space resolution are both the same.

- conv/diff - convection or diffusion is the examined transport
- number after conv/diff - value of the corresponding coefficient
- e4/e6/e8 - resolution e. g. e4 means 1e4 = 100 * 100 (time * space)
- A/N - Analytical or Numerical solution


