import ctypes as ct
from os import path
import platform

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

basepath = path.dirname("python_constructor.py")
if platform.system() == 'Windows':
	library_path = path.abspath(path.join(basepath, "..", "build/src/libRunaphys.dll"))
else:
	library_path = path.abspath(path.join(basepath, "..", "build/src/libRunaphys.so"))
	
lib_Runaphys = ct.CDLL(library_path)
adv_RE_pop = lib_Runaphys.Runaphys_advance_runaway_population
adv_RE_pop.argtypes = ct.POINTER(PLASMA),ct.c_double,ct.c_double,ct.c_double,ct.POINTER(MODULE),ct.POINTER(ct.c_double)
adv_RE_pop.restype = ct.c_double
rate_values_type = (ct.c_double * 4)

print("MODULE and PLASMA structs are ready to use. ")
print("rate_values_type and adv_RE_pop() are ready to use.")
