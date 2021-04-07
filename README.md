# szakdolgozat - Bachelor's thesis

1D convection-diffusion equation solver (in progress)
- with time and space dependent coefficinets
- with coefficient reading from CSV files
- with Dirac-delta and uniform density initial value conditions
- with 2D time vs. space contour plotting
- with solving on cylindrical mesh
- with partial HDF5 export and import capabilities

Setup:
1. New Anaconda enviromnet -> Anaconda CMD Prompt:
2. `>>> conda create --name <MYFIPYENV> --channel conda-forge python=3.9.2 numpy scipy matplotlib future`
3. `>>> pip install fipy`
4. Run Jupyter Notebook 6.1.4 or higher in MYFIPYENV
5. Run scripts

Goals:
- 1D modelling of runaway electron transport outside of a cylindrical magnetic field in the ergodic zone along the radius. (Approximating tokamak fusion reactors' toroidal vacuum chambers)
- Using runaway electron generation code from OSREP/Runaphys (https://github.com/osrep/Runaphys) as input for the solver

Using FiPy (https://www.ctcms.nist.gov/fipy/index.html)
