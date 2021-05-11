# szakdolgozat - Bachelor's thesis

1D convection-diffusion equation solver (in progress)
- time and space dependent coefficinets
- coefficient reading from CSV files
- Dirac-delta and uniform density initial value conditions
- 2D time vs. space contour plotting
- solving on 1D cylindrical mesh or simple 1D mesh
- partial HDF5 export and import capabilities
- convection verified with analytical solution
- diffusion verified on slab geometry as well as cylindrical geometry
- runaway generation code implementation
- in specific functions runaway electrons in magnetic islands are treated differently

# Setup:
1. New Anaconda enviromnet -> Anaconda CMD Prompt:
2. `>>> conda create --name <MYFIPYENV> --channel conda-forge python=3.9.4 numpy scipy matplotlib future h5py`
3. `>>> pip install fipy`
4. Run Jupyter Notebook 6.1.4 or higher in MYFIPYENV
5. Run scripts

# Goals:
- 1D modelling of runaway electron transport outside of a cylindrical magnetic field in the ergodic zone along the radius. (Approximating tokamak fusion reactors' toroidal vacuum chambers)
- Using runaway electron generation code from OSREP/Runaphys (https://github.com/osrep/Runaphys) as input for the solver

# Some obtained results
Dreicer generation 1000x1000 resolution
![Dreicer](https://github.com/leferi99/szakdolgozat/blob/main/results/Dreicer_1000x1000.jpg)

Avalanche generation 1000x1000 resolution with islands treated differently - This is the closest to the reference material. (The avalanche generation is basically none because of the collision time.)
![avalanche_withI](https://github.com/leferi99/szakdolgozat/blob/main/results/avalanche_withI_1000x1000.jpg)

# Reference:
Model based on the following paper: K.  Särkimäki,  E.  Hirvijoki,  J.  Decker,  J.  Varje  and  T.  Kurki-Suonio,  „An  advection–diffusion model for cross-field runaway electron transport in perturbed magnetic fields”, Plasma Physics and Controlled Fusion, vol. 58., no. 12., p. 125017., 2016., ISSN: 1361-6587. (https://iopscience.iop.org/article/10.1088/0741-3335/58/12/125017)

Using FiPy (https://www.ctcms.nist.gov/fipy/index.html)

