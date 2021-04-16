from ctypes import *
lib = cdll.LoadLibrary('MyLib.so')

#lib.MyClass_int_set.restype = c_int
#lib.MyClass_int_set.argtypes = [c_double]

class MyClass(object):
	def __init__(self):
		self.obj = lib.newMyClass()
	
	def MyClass_int_set(self, i):
		lib.MyClass_int_set(self.obj, i)
	
	def MyClass_int_get(self):
		return lib.MyClass_int_get(self.obj)
	
	def deleteMyClass(self):
		lib.deleteMyClass(self.obj)

f = MyClass()
f.MyClass_int_set(3)
print(f.MyClass_int_get())
f.deleteMyClass()

print("The code ran until the end")
	
	