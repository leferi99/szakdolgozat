import ctypes as ct
from os import path
import unittest

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


modules = MODULE(b"hc_formula_63",False,b"rosenbluth_putvinski",False,False,1.0,0.85)

reference_rho = ct.c_double(0.95)
reference_electron_density = ct.c_double(1e21)
reference_electron_temperature = ct.c_double(1e5)
reference_effective_charge = ct.c_double(1.0)
reference_electric_field = ct.c_double(1.2 * 0.935882)
reference_magnetic_field = ct.c_double(2.0)
reference_runaway_density_before = ct.c_double(3.12292e15)
reference_rate_values = (ct.c_double * 4)(0.,0.,0.,0.)
reference_timestep = ct.c_double(1.82129e-3)
reference_inv_asp_ratio = ct.c_double(0.30303)
reference_rho_tor_norm = ct.c_double(0.65)
reference_runaway_density_after63 = 7.85003e16

reference_plasma_local = PLASMA(reference_rho,
		         	reference_electron_density,
		      		reference_electron_temperature,
		      		reference_effective_charge,
		      		reference_electric_field,
		      		reference_magnetic_field,
		      		reference_runaway_density_before)

basepath = path.dirname("python_help.py")
library_path = path.abspath(path.join(basepath, "..", "build/src/libRunaphys.so"))

lib_Runaphys = ct.CDLL(library_path)
adv_RE_pop = lib_Runaphys.runaphys_advance_runaway_population
adv_RE_pop.argtypes = ct.POINTER(PLASMA), ct.c_double, ct.c_double, ct.c_double, ct.POINTER(MODULE), ct.POINTER(ct.c_double)

adv_RE_pop.restype = ct.c_double

class Python_wrapper_test(unittest.TestCase):
	
		
	def test_python_advance_runaway_population_63(self):
		self.assertAlmostEqual(reference_runaway_density_after63, adv_RE_pop(ct.byref(reference_plasma_local),reference_timestep,reference_inv_asp_ratio,reference_rho_tor_norm,ct.byref(modules),reference_rate_values), delta=reference_runaway_density_after63 * 1e-2)
		

if __name__ == '__main__':
    unittest.main()
    