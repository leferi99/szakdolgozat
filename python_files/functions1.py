import fipy as fp  ## finite volume PDE solver
from fipy.tools import numerix  ## requirement for FiPy, in practice same as numpy
import copy  ## we need the deepcopy() function because some FiPy objects are mutable
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
import h5py
import datetime
import ctypes as ct
from os import path
from scipy.signal import savgol_filter

plt.rcParams.update({'font.size': 18})

class MODULE(ct.Structure):
    _fields_ = [("dreicer_formula", ct.c_char_p),
                ("dreicer_toroidicity", ct.c_bool),
                ("avalanche_formula", ct.c_char_p),
                ("avalanche_toroidicity", ct.c_bool),
                ("hdf5_output", ct.c_bool),
                ("warning_percentage_limit", ct.c_double),
                ("rho_edge_calculation_limit", ct.c_double)]

class PLASMA(ct.Structure):
    _fields_ = [("rho", ct.c_double),
                ("electron_density", ct.c_double),
                ("electron_temperature", ct.c_double),
                ("effective_charge", ct.c_double),
                ("electric_field", ct.c_double),
                ("magnetic_field", ct.c_double),
                ("runaway_density", ct.c_double)]

#basepath = path.dirname("functions1.py")
#if platform.system() == 'Windows':
#	library_path = path.abspath(path.join(basepath, "..", "build/src/libRunaphys.dll"))
#else:
#	library_path = path.abspath(path.join(basepath, "..", "build/src/libRunaphys.so"))
	
library_path=("/pfs/work/g2fleng/Runaphys/build/src/libRunaphys.so")
lib_Runaphys = ct.CDLL(library_path)
adv_RE_pop = lib_Runaphys.runaphys_advance_runaway_population
adv_RE_pop.argtypes = ct.POINTER(PLASMA),ct.c_double,ct.c_double,ct.c_double,ct.POINTER(MODULE),ct.POINTER(ct.c_double)
adv_RE_pop.restype = ct.c_double
rate_values_type = (ct.c_double * 4)

def delta_func(x, epsilon, coeff):
    return ((x < epsilon) & (x > -epsilon)) * \
        (coeff * (1 + numerix.cos(numerix.pi * x / epsilon)) / (2 * epsilon))

def plot_solution(solution_array,ticks=None,levels=300,logdiff=5,figsize=(10,5),duration=0.001,nt=1000, saveplot=False):
    if ticks == None:
        ticks = logdiff
    else:
        ticks = ticks
    
    dt = duration / nt
    sol_min = solution_array[:,:,1].min()
    sol_max = solution_array[:,:,1].max()
    logmax = math.ceil(np.log10(sol_max))
    logmin = logmax - logdiff
    numofticks = ticks
    div = logdiff // numofticks
    power = np.arange((logmax - (numofticks * div)), logmax, div)
    array1 = np.zeros(len(power)) + 10.
    ticks1 = np.power(array1, power)
    levels1 = np.logspace(logmin,logmax,levels, base=10.0)
    norm = matplotlib.colors.LogNorm(vmin=(10 ** logmin),vmax=sol_max)
    formatter = matplotlib.ticker.LogFormatter(10, labelOnlyBase=False)
    fig = plt.figure(figsize=figsize)
    plot = plt.contourf(solution_array[0,:,0], np.linspace(dt,duration,nt), solution_array[:,:,1],
                        levels=levels1, norm=norm, cmap=plt.cm.jet)
    axes = plt.gca()
    axes.set_facecolor((0.,0.,0.25)) ## dark blue will display where the values are too small
    axes.set_yscale('log')
    axes.yaxis.set_major_formatter(matplotlib.ticker.LogFormatterExponent(base=10))
    cbar = plt.colorbar(ticks=ticks1, format=formatter)
    cbar.set_label(r'log$_{10}$(Density (m$^{-3}$)/1m$^{-3}$)')
    cbar.formatter = matplotlib.ticker.LogFormatterExponent(base=10) # 10 is the default
    cbar.update_ticks()
    plt.xlabel(r'$\rho$')
    plt.ylabel(r'log$_{10}$[time(s)/1s]')
    plt.title(r'Density of runaway electrons')
    plt.gcf().subplots_adjust(bottom=0.2)
    if saveplot == True:
        now = datetime.datetime.now()
        date_time = now.strftime("%Y%m%d_%H.%M.%S")
        plt.savefig("solution_" + date_time + ".jpg", dpi=150)
    else:
        1
    return

def coeff_plot(conv_i, diff_i):
    fig_conv = plt.figure(figsize=(8.5,3))
    plt.plot(*conv_i.T)
    plt.title('Convection coefficient')
    fig_diff = plt.figure(figsize=(8.5,3))
    plt.plot(*diff_i.T)
    plt.title('Diffusion coefficient')
    return
    
def hdf5_save(fname, solution, diff, conv, duration):
    now = datetime.datetime.now()
    date_time = now.strftime("%Y%m%d_%H.%M.%S")
    f_Out = h5py.File(fname + "_" + date_time + ".hdf5", "w")
    h5solution = f_Out.create_dataset("01_solution", data=solution)
    h5duration = f_Out.create_dataset("02_duration", data=duration)
    h5conv = f_Out.create_dataset("03_convcoeff", data=conv)
    h5diff = f_Out.create_dataset("04_diffcoeff", data=diff) 
    return
    
def plot_hdf5(fname, plotcoeff=False, logdiff=5):
    f_In = h5py.File(fname, "r")
    h5solution = f_In["01_solution"]
    h5duration = f_In["02_duration"]
    nt = h5solution.shape[0]
    plot_solution(h5solution[()], logdiff=logdiff, duration=h5duration[()], nt=nt)
    if plotcoeff == True:
        h5conv = f_In["03_convcoeff"]
        h5diff = f_In["04_diffcoeff"]
        coeff_plot(conv_i=h5conv[()], diff_i=h5diff[()])
    else:
        1
        
    return

def solve_DiracIC(saveplot = False, R_from = 0.7, R_to = 1.0, nr = 1000, duration = 0.001, nt = 1000,
                  diracLoc = 0.85, diracCoeff = 1., diracPercentage = 2,
                  conv_file = 'convC.txt', diff_file = 'diffC.txt',  plotcoeff = False,
                  levels = 300, logdiff = 6, ticks = None, figsize=(10,5), hdf5 = False):
    
    dr = (R_to - R_from) / nr  ## distance between the centers of the mesh cells
    dt = duration / nt  ## length of one timestep
    solution = np.zeros((nt,nr,2))
    for j in range(nr):
        solution[:,j,0] = (j * dr) + (dr / 2) + R_from

    mesh = fp.CylindricalGrid1D(dx=dr, nx=nr)  ## 1D mesh based on the radial coordinates 
    mesh = mesh + (R_from,)  ## translation of the mesh to R_from
    n = fp.CellVariable(mesh=mesh)  ## fipy.CellVariable for the density solution in each timestep
    diracWidth = int((nr / 100) * diracPercentage)
    n.setValue(delta_func(mesh.x - diracLoc, diracWidth * dr, diracCoeff))
    conv_data = np.genfromtxt(conv_file, delimiter=',')
    diff_data = np.genfromtxt(diff_file, delimiter=',')
    conv_i = np.zeros((nr, 2))
    diff_i = np.zeros((nr, 2))
    for i in range(conv_i.shape[0]):
        conv_i[i, 0] = R_from + (i * dr) + (dr / 2)

    for i in range(diff_i.shape[0]):
        diff_i[i, 0] = R_from + (i * dr) + (dr / 2)

    conv_i[:,1] = np.interp(conv_i[:,0],conv_data[:,0],conv_data[:,1])
    diff_i[:,1] = np.interp(diff_i[:,0],diff_data[:,0],diff_data[:,1])
    dC = diff_i[:,1]
    diffCoeff = fp.CellVariable(mesh=mesh, value=dC)
    cC = conv_i[:,1]
    convCoeff = fp.CellVariable(mesh=mesh, value=[cC])
    gradLeft = (0.,)  ## density gradient (at the "left side of the radius") - must be a vector
    valueRight = 0.  ## density value (at the "right end of the radius")
    n.faceGrad.constrain(gradLeft, where=mesh.facesLeft)  ## applying Neumann boundary condition
    n.constrain(valueRight, mesh.facesRight)  ## applying Dirichlet boundary condition
    convCoeff.setValue(0, where=mesh.x<(R_from + dr))  ## convection coefficient 0 at the inner edge
    diffCoeff.setValue(0.001, where=mesh.x<(R_from + dr))  ## diffusion coefficient almost 0 at inner edge
    eq = (fp.TransientTerm() == fp.DiffusionTerm(coeff=diffCoeff)
          - fp.ConvectionTerm(coeff=convCoeff))
    for i in range(nt):
        eq.solve(var=n, dt=dt)
        solution[i,0:nr,1]=copy.deepcopy(n.value)

    plot_solution(solution,ticks=ticks,levels=levels,logdiff=logdiff,figsize=figsize,
                  duration=duration, nt=nt, saveplot=saveplot)
    if plotcoeff == True:
        coeff_plot(conv_i=conv_i, diff_i=diff_i)
    else:
        1
    
    if hdf5 == True:
        hdf5_save(fname="DiracIC",solution=solution, conv=conv_i, diff=diff_i, duration=duration)
    else:
        1
    
    return solution

def solve_uniformIC(saveplot = False, R_from = 0.7, R_to = 1.0, nr = 1000, duration = 0.001, nt = 1000,
                    conv_file = 'convC.txt', diff_file = 'diffC.txt', plotcoeff = False,
                    levels = 300, logdiff = 5, ticks = None, figsize=(10,5), hdf5 = False):
    
    dr = (R_to - R_from) / nr  ## distance between the centers of the mesh cells
    dt = duration / nt  ## length of one timestep
    solution = np.zeros((nt,nr,2))
    for j in range(nr):
        solution[:,j,0] = (j * dr) + (dr / 2) + R_from

    mesh = fp.CylindricalGrid1D(dx=dr, nx=nr)  ## 1D mesh based on the radial coordinates 
    mesh = mesh + (R_from,)  ## translation of the mesh to R_from
    n = fp.CellVariable(mesh=mesh)  ## fipy.CellVariable for the density solution in each timestep
    conv_data = np.genfromtxt(conv_file, delimiter=',')
    diff_data = np.genfromtxt(diff_file, delimiter=',')
    conv_i = np.zeros((nr, 2))
    diff_i = np.zeros((nr, 2))
    for i in range(conv_i.shape[0]):
        conv_i[i, 0] = R_from + (i * dr) + (dr / 2)

    for i in range(diff_i.shape[0]):
        diff_i[i, 0] = R_from + (i * dr) + (dr / 2)

    conv_i[:,1] = np.interp(conv_i[:,0],conv_data[:,0],conv_data[:,1])
    diff_i[:,1] = np.interp(diff_i[:,0],diff_data[:,0],diff_data[:,1])
    dC = diff_i[:,1]
    diffCoeff = fp.CellVariable(mesh=mesh, value=dC)
    cC = conv_i[:,1]
    convCoeff = fp.CellVariable(mesh=mesh, value=[cC])
    n.setValue(1.0)
    gradLeft = (0.,)  ## density gradient (at the "left side of the radius") - must be a vector
    valueRight = 0.  ## density value (at the "right end of the radius")
    n.faceGrad.constrain(gradLeft, where=mesh.facesLeft)  ## applying Neumann boundary condition
    n.constrain(valueRight, mesh.facesRight)  ## applying Dirichlet boundary condition
    convCoeff.setValue(0, where=mesh.x<(R_from + dr))  ## convection coefficient 0 at the inner edge
    diffCoeff.setValue(0.001, where=mesh.x<(R_from + dr))  ## diffusion coefficient almost 0 at inner edge
    eq = (fp.TransientTerm() == fp.DiffusionTerm(coeff=diffCoeff)
          - fp.ConvectionTerm(coeff=convCoeff))
    for i in range(nt):
        eq.solve(var=n, dt=dt)
        solution[i,0:nr,1]=copy.deepcopy(n.value)

    plot_solution(solution,ticks=ticks,levels=levels,logdiff=logdiff,figsize=figsize,
                  duration=duration, nt=nt, saveplot=saveplot)
    if plotcoeff == True:
        coeff_plot(conv_i=conv_i, diff_i=diff_i)
    else:
        1
    
    if hdf5 == True:
        hdf5_save(fname="uniformIC", solution=solution, conv=conv_i, diff=diff_i, duration=duration)
    else:
        1
    
    return solution
    
def solve_Dreicer(saveplot = False, R_from = 0.7, R_to = 1.0, nr = 1000, duration = 0.001, nt = 1000,
                    conv_file = 'convC.txt', diff_file = 'diffC.txt', plotcoeff = False,
                    levels = 300, logdiff = 5, ticks = None, figsize=(10,5), hdf5 = False):
    
    dr = (R_to - R_from) / nr  ## distance between the centers of the mesh cells
    dt = duration / nt  ## length of one timestep
    solution = np.zeros((nt,nr,2))
    for j in range(nr):
        solution[:,j,0] = (j * dr) + (dr / 2) + R_from

    mesh = fp.CylindricalGrid1D(dx=dr, nx=nr)  ## 1D mesh based on the radial coordinates 
    mesh = mesh + (R_from,)  ## translation of the mesh to R_from
    n = fp.CellVariable(mesh=mesh)  ## fipy.CellVariable for the density solution in each timestep
    conv_data = np.genfromtxt(conv_file, delimiter=',')
    diff_data = np.genfromtxt(diff_file, delimiter=',')
    conv_i = np.zeros((nr, 2))
    diff_i = np.zeros((nr, 2))
    for i in range(conv_i.shape[0]):
        conv_i[i, 0] = R_from + (i * dr) + (dr / 2)

    for i in range(diff_i.shape[0]):
        diff_i[i, 0] = R_from + (i * dr) + (dr / 2)

    conv_i[:,1] = np.interp(conv_i[:,0],conv_data[:,0],conv_data[:,1])
    diff_i[:,1] = np.interp(diff_i[:,0],diff_data[:,0],diff_data[:,1])
    dC = diff_i[:,1]
    diffCoeff = fp.CellVariable(mesh=mesh, value=dC)
    cC = conv_i[:,1]
    convCoeff = fp.CellVariable(mesh=mesh, value=[cC])
    n.setValue(0.0)
    gradLeft = (0.,)  ## density gradient (at the "left side of the radius") - must be a vector
    valueRight = 0.  ## density value (at the "right end of the radius")
    n.faceGrad.constrain(gradLeft, where=mesh.facesLeft)  ## applying Neumann boundary condition
    n.constrain(valueRight, mesh.facesRight)  ## applying Dirichlet boundary condition
    convCoeff.setValue(0, where=mesh.x<(R_from + dr))  ## convection coefficient 0 at the inner edge
    diffCoeff.setValue(0.001, where=mesh.x<(R_from + dr))  ## diffusion coefficient almost 0 at inner edge
    
    modules = MODULE(b"hc_formula_63",False,b"",False,False,1.0,1.0001)
    electron_temperature = ct.c_double(300.)
    electron_density = ct.c_double(1e20)
    effective_charge = ct.c_double(1.)
    electric_field = ct.c_double(3.66)
    magnetic_field = ct.c_double(1.)
    inv_asp_ratio = ct.c_double(0.30303)
    rate_values = (ct.c_double * 4)(0.,0.,0.,0.)
    
    eq = (fp.TransientTerm() == fp.DiffusionTerm(coeff=diffCoeff)
          - fp.ConvectionTerm(coeff=convCoeff))
    for i in range(nt):
        for j in range(nr):
            plasma_local = PLASMA(ct.c_double(mesh.x[j]),
                                  electron_density,
                                  electron_temperature,
                                  effective_charge,
                                  electric_field,
                                  magnetic_field,
                                  ct.c_double(n.value[j]))
            n.value[j] = adv_RE_pop(ct.byref(plasma_local),dt,inv_asp_ratio,ct.c_double(mesh.x[j]),ct.byref(modules),rate_values)
        
        eq.solve(var=n, dt=dt)
        solution[i,0:nr,1]=copy.deepcopy(n.value)

    plot_solution(solution,ticks=ticks,levels=levels,logdiff=logdiff,figsize=figsize,
                  duration=duration, nt=nt, saveplot=saveplot)
    if plotcoeff == True:
        coeff_plot(conv_i=conv_i, diff_i=diff_i)
    else:
        1
    
    if hdf5 == True:
        hdf5_save(fname="uniformIC", solution=solution, conv=conv_i, diff=diff_i, duration=duration)
    else:
        1
    
    return solution
    
def solve_avalanche(saveplot = False, R_from = 0.7, R_to = 1.0, nr = 1000, duration = 0.001, nt = 1000,
                    conv_file = 'convC.txt', diff_file = 'diffC.txt', plotcoeff = False,
                    levels = 300, logdiff = 5, ticks = None, figsize=(10,5), hdf5 = False):
    
    dr = (R_to - R_from) / nr  ## distance between the centers of the mesh cells
    dt = duration / nt  ## length of one timestep
    solution = np.zeros((nt,nr,2))
    for j in range(nr):
        solution[:,j,0] = (j * dr) + (dr / 2) + R_from

    mesh = fp.CylindricalGrid1D(dx=dr, nx=nr)  ## 1D mesh based on the radial coordinates 
    mesh = mesh + (R_from,)  ## translation of the mesh to R_from
    n = fp.CellVariable(mesh=mesh)  ## fipy.CellVariable for the density solution in each timestep
    conv_data = np.genfromtxt(conv_file, delimiter=',')
    diff_data = np.genfromtxt(diff_file, delimiter=',')
    conv_i = np.zeros((nr, 2))
    diff_i = np.zeros((nr, 2))
    for i in range(conv_i.shape[0]):
        conv_i[i, 0] = R_from + (i * dr) + (dr / 2)

    for i in range(diff_i.shape[0]):
        diff_i[i, 0] = R_from + (i * dr) + (dr / 2)

    conv_i[:,1] = np.interp(conv_i[:,0],conv_data[:,0],conv_data[:,1])
    diff_i[:,1] = np.interp(diff_i[:,0],diff_data[:,0],diff_data[:,1])
    dC = diff_i[:,1]
    diffCoeff = fp.CellVariable(mesh=mesh, value=dC)
    cC = conv_i[:,1]
    convCoeff = fp.CellVariable(mesh=mesh, value=[cC])
    n.setValue(1e15)
    gradLeft = (0.,)  ## density gradient (at the "left side of the radius") - must be a vector
    valueRight = 0.  ## density value (at the "right end of the radius")
    n.faceGrad.constrain(gradLeft, where=mesh.facesLeft)  ## applying Neumann boundary condition
    n.constrain(valueRight, mesh.facesRight)  ## applying Dirichlet boundary condition
    convCoeff.setValue(0, where=mesh.x<(R_from + dr))  ## convection coefficient 0 at the inner edge
    diffCoeff.setValue(0.001, where=mesh.x<(R_from + dr))  ## diffusion coefficient almost 0 at inner edge
    
    modules = MODULE(b"",False,b"rosenbluth_putvinski",False,False,1.0,1.0001)
    electron_temperature = ct.c_double(300.)
    electron_density = ct.c_double(1e20)
    effective_charge = ct.c_double(1.)
    electric_field = ct.c_double(3.66)
    magnetic_field = ct.c_double(1.)
    inv_asp_ratio = ct.c_double(0.30303)
    rate_values = (ct.c_double * 4)(0.,0.,0.,0.)
    
    eq = (fp.TransientTerm() == fp.DiffusionTerm(coeff=diffCoeff)
          - fp.ConvectionTerm(coeff=convCoeff))
    for i in range(nt):
        for j in range(nr):
            plasma_local = PLASMA(ct.c_double(mesh.x[j]),
                                  electron_density,
                                  electron_temperature,
                                  effective_charge,
                                  electric_field,
                                  magnetic_field,
                                  ct.c_double(n.value[j]))
            n.value[j] = adv_RE_pop(ct.byref(plasma_local),dt,inv_asp_ratio,ct.c_double(mesh.x[j]),ct.byref(modules),rate_values)
        
        eq.solve(var=n, dt=dt)
        solution[i,0:nr,1]=copy.deepcopy(n.value)

    plot_solution(solution,ticks=ticks,levels=levels,logdiff=logdiff,figsize=figsize,
                  duration=duration, nt=nt, saveplot=saveplot)
    if plotcoeff == True:
        coeff_plot(conv_i=conv_i, diff_i=diff_i)
    else:
        1
    
    if hdf5 == True:
        hdf5_save(fname="uniformIC", solution=solution, conv=conv_i, diff=diff_i, duration=duration)
    else:
        1
    
    return solution
    
def solve_Dreicer_withI(saveplot = False, R_from = 0.7, R_to = 1.0, nr = 1000, duration = 0.001, nt = 1000,
                    conv_file = 'convC.txt', diff_file = 'diffC.txt', plotcoeff = False,
                    levels = 300, logdiff = 5, ticks = None, figsize=(10,5), hdf5 = False):
    
    dr = (R_to - R_from) / nr  ## distance between the centers of the mesh cells
    dt = duration / nt  ## length of one timestep
    solution = np.zeros((nt,nr,2))
    for j in range(nr):
        solution[:,j,0] = (j * dr) + (dr / 2) + R_from

    mesh = fp.CylindricalGrid1D(dx=dr, nx=nr)  ## 1D mesh based on the radial coordinates 
    mesh = mesh + (R_from,)  ## translation of the mesh to R_from
    n = fp.CellVariable(mesh=mesh)  ## fipy.CellVariable for the density solution in each timestep
    conv_data = np.genfromtxt(conv_file, delimiter=',')
    diff_data = np.genfromtxt(diff_file, delimiter=',')
    conv_i = np.zeros((nr, 2))
    diff_i = np.zeros((nr, 2))
    for i in range(conv_i.shape[0]):
        conv_i[i, 0] = R_from + (i * dr) + (dr / 2)

    for i in range(diff_i.shape[0]):
        diff_i[i, 0] = R_from + (i * dr) + (dr / 2)

    conv_i[:,1] = np.interp(conv_i[:,0],conv_data[:,0],conv_data[:,1])
    diff_i[:,1] = np.interp(diff_i[:,0],diff_data[:,0],diff_data[:,1])
    dC = diff_i[:,1]
    diffCoeff = fp.CellVariable(mesh=mesh, value=dC)
    cC = conv_i[:,1]
    convCoeff = fp.CellVariable(mesh=mesh, value=[cC])
    
    idata = np.genfromtxt('island_data.csv', delimiter=',')
    islands_ratio = np.zeros((nr,2))
    for i in range(nr):
        islands_ratio[i, 0] = R_from + (i * dr) + (dr / 2)
    
    islands_ratio[:,1] = np.interp(islands_ratio[:,0],idata[:,0],idata[:,1])
    w_length = math.ceil(nr/20)
    if (w_length % 2) == 0:
        w_length = w_length + 1
    else:
        pass
    
    islands_ratio[:,1] = savgol_filter(islands_ratio[:,1],w_length, 3)
    re_in_islands = islands_ratio[:,1]
    
    gradLeft = (0.,)  ## density gradient (at the "left side of the radius") - must be a vector
    valueRight = 0.  ## density value (at the "right end of the radius")
    n.faceGrad.constrain(gradLeft, where=mesh.facesLeft)  ## applying Neumann boundary condition
    n.constrain(valueRight, mesh.facesRight)  ## applying Dirichlet boundary condition
    convCoeff.setValue(0, where=mesh.x<(R_from + dr))  ## convection coefficient 0 at the inner edge
    diffCoeff.setValue(0.001, where=mesh.x<(R_from + dr))  ## diffusion coefficient almost 0 at inner edge
    
    modules = MODULE(b"hc_formula_63",False,b"",False,False,1.0,1.0001)
    electron_temperature = ct.c_double(300.)
    electron_density = ct.c_double(1e20)
    effective_charge = ct.c_double(1.)
    electric_field = ct.c_double(3.66)
    magnetic_field = ct.c_double(1.)
    inv_asp_ratio = ct.c_double(0.30303)
    rate_values = (ct.c_double * 4)(0.,0.,0.,0.)
    
    eq = (fp.TransientTerm() == fp.DiffusionTerm(coeff=diffCoeff)
          - fp.ConvectionTerm(coeff=convCoeff))
    for i in range(nt):
        for j in range(nr):
            plasma_local = PLASMA(ct.c_double(mesh.x[j]),
                                  electron_density,
                                  electron_temperature,
                                  effective_charge,
                                  electric_field,
                                  magnetic_field,
                                  ct.c_double(n.value[j]))
            n.value[j] = adv_RE_pop(ct.byref(plasma_local),dt,inv_asp_ratio,ct.c_double(mesh.x[j]),ct.byref(modules),rate_values)
        
        eq.solve(var=n, dt=dt)
        if i == 0:
            solution[i,0:nr,1] = copy.deepcopy(n.value) + re_in_islands * copy.deepcopy(n.value)
        else:
            re_in_islands = re_in_islands / (re_in_islands + copy.deepcopy(n.value))
            solution[i,0:nr,1] = copy.deepcopy(n.value) + re_in_islands * copy.deepcopy(n.value)
            re_in_islands = solution[i,0:nr,1] * re_in_islands

    plot_solution(solution,ticks=ticks,levels=levels,logdiff=logdiff,figsize=figsize,
                  duration=duration, nt=nt, saveplot=saveplot)
    if plotcoeff == True:
        coeff_plot(conv_i=conv_i, diff_i=diff_i)
    else:
        1
    
    if hdf5 == True:
        hdf5_save(fname="uniformIC", solution=solution, conv=conv_i, diff=diff_i, duration=duration)
    else:
        1
    
    return solution
    
def solve_avalanche_withI(saveplot = False, R_from = 0.7, R_to = 1.0, nr = 1000, duration = 0.001, nt = 1000,
                    conv_file = 'convC.txt', diff_file = 'diffC.txt', plotcoeff = False, n0 = 1,
                    levels = 300, logdiff = 5, ticks = None, figsize=(10,5), hdf5 = False):
    
    dr = (R_to - R_from) / nr  ## distance between the centers of the mesh cells
    dt = duration / nt  ## length of one timestep
    solution = np.zeros((nt,nr,2))
    for j in range(nr):
        solution[:,j,0] = (j * dr) + (dr / 2) + R_from

    mesh = fp.CylindricalGrid1D(dx=dr, nx=nr)  ## 1D mesh based on the radial coordinates 
    mesh = mesh + (R_from,)  ## translation of the mesh to R_from
    n = fp.CellVariable(mesh=mesh)  ## fipy.CellVariable for the density solution in each timestep
    conv_data = np.genfromtxt(conv_file, delimiter=',')
    diff_data = np.genfromtxt(diff_file, delimiter=',')
    conv_i = np.zeros((nr, 2))
    diff_i = np.zeros((nr, 2))
    for i in range(conv_i.shape[0]):
        conv_i[i, 0] = R_from + (i * dr) + (dr / 2)

    for i in range(diff_i.shape[0]):
        diff_i[i, 0] = R_from + (i * dr) + (dr / 2)

    conv_i[:,1] = np.interp(conv_i[:,0],conv_data[:,0],conv_data[:,1])
    diff_i[:,1] = np.interp(diff_i[:,0],diff_data[:,0],diff_data[:,1])
    dC = diff_i[:,1]
    diffCoeff = fp.CellVariable(mesh=mesh, value=dC)
    cC = conv_i[:,1]
    convCoeff = fp.CellVariable(mesh=mesh, value=[cC])
    
    idata = np.genfromtxt('island_data.csv', delimiter=',')
    islands_ratio = np.zeros((nr,2))
    for i in range(nr):
        islands_ratio[i, 0] = R_from + (i * dr) + (dr / 2)
    
    islands_ratio[:,1] = np.interp(islands_ratio[:,0],idata[:,0],idata[:,1])
    w_length = math.ceil(nr/20)
    if (w_length % 2) == 0:
        w_length = w_length + 1
    else:
        pass
    
    islands_ratio[:,1] = savgol_filter(islands_ratio[:,1],w_length, 3)
    islands_ratio[islands_ratio < 0] = 0
    n.setValue(n0)
    n.setValue(n.value[:] - (islands_ratio[:,1] * n0))
    re_in_islands = islands_ratio[:,1] * n0
    
    gradLeft = (0.,)  ## density gradient (at the "left side of the radius") - must be a vector
    valueRight = 0.  ## density value (at the "right end of the radius")
    n.faceGrad.constrain(gradLeft, where=mesh.facesLeft)  ## applying Neumann boundary condition
    n.constrain(valueRight, mesh.facesRight)  ## applying Dirichlet boundary condition
    convCoeff.setValue(0, where=mesh.x<(R_from + dr))  ## convection coefficient 0 at the inner edge
    diffCoeff.setValue(0.001, where=mesh.x<(R_from + dr))  ## diffusion coefficient almost 0 at inner edge
    
    modules = MODULE(b"",False,b"rosenbluth_putvinski",False,False,1.0,1.0001)
    electron_temperature = ct.c_double(300.)
    electron_density = ct.c_double(1e20)
    effective_charge = ct.c_double(1.)
    electric_field = ct.c_double(3.66)
    magnetic_field = ct.c_double(1.)
    inv_asp_ratio = ct.c_double(0.30303)
    rate_values = (ct.c_double * 4)(0.,0.,0.,0.)
    
    eq = (fp.TransientTerm() == fp.DiffusionTerm(coeff=diffCoeff)
          - fp.ConvectionTerm(coeff=convCoeff))
    for i in range(nt):
        for j in range(nr):
            plasma_local = PLASMA(ct.c_double(mesh.x[j]),
                                  electron_density,
                                  electron_temperature,
                                  effective_charge,
                                  electric_field,
                                  magnetic_field,
                                  ct.c_double(n.value[j]))
            n.value[j] = adv_RE_pop(ct.byref(plasma_local),dt,inv_asp_ratio,ct.c_double(mesh.x[j]),ct.byref(modules),rate_values)
        
        eq.solve(var=n, dt=dt)
        re_in_islands = re_in_islands / (re_in_islands + copy.deepcopy(n.value))
        solution[i,0:nr,1] = copy.deepcopy(n.value) + re_in_islands * copy.deepcopy(n.value)
        re_in_islands = solution[i,0:nr,1] * re_in_islands

    plot_solution(solution,ticks=ticks,levels=levels,logdiff=logdiff,figsize=figsize,
                  duration=duration, nt=nt, saveplot=saveplot)
    if plotcoeff == True:
        coeff_plot(conv_i=conv_i, diff_i=diff_i)
    else:
        1
    
    if hdf5 == True:
        hdf5_save(fname="uniformIC", solution=solution, conv=conv_i, diff=diff_i, duration=duration)
    else:
        1
    
    return solution

def fhelp():
    fname = "help.txt"
    with open(fname, 'r') as fin:
        print(fin.read(), end="")
    
    return

print("For usage help call fhelp()")