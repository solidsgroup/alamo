import os, subprocess

alamo_home = os.getcwd()

result = subprocess.run(['mpicxx', '--showme:compile'],capture_output=True,text=True).stdout.strip()
openmpi_includes = [r.replace('-I','') for r in result.split()]

f = open("setup.py",'w')
f.write("""
from setuptools import setup

setup(
    name='alamo',
    py_modules=['alamo'],
    install_requires=['cppyy'],
)
""")
f.close()

os.makedirs('alamo',exist_ok=True)

f = open("alamo.py",'w')
f.write(f"""
import cppyy
import types
import sys

""")

for openmpi_include in openmpi_includes:
    f.write(f"cppyy.add_include_path('{openmpi_include}')\n")
    #cppyy.add_include_path('/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi')

f.write(f"""
cppyy.add_include_path('{alamo_home}/ext/amrex/2d-debug-g++-25.07/include/')
cppyy.add_include_path('{alamo_home}/src/')
cppyy.load_library("{alamo_home}/ext/amrex/2d-debug-g++-25.07/lib/libamrex.so")
cppyy.load_library("{alamo_home}/lib/libalamo-2d-debug-g++.so")

# Create proxy module
class AlamoModule(types.ModuleType):
    def __getattr__(self, name):
        return getattr(cppyy.gbl, name)

    def include(self, path):
        cppyy.include("{alamo_home}/src/" + path)

# Replace current module with AlamoModule instance
sys.modules[__name__].__class__ = AlamoModule
""")
