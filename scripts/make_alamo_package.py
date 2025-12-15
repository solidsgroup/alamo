import os, subprocess, argparse

parser = argparse.ArgumentParser(description='For use by Makefile generating Alamo python library');
parser.add_argument('--alamo_home',   default=os.getcwd(), help='Alamo home directory')
parser.add_argument('--postfix',      required=True, help="Build postfix")
parser.add_argument('--amrex',required=True, help='AMReX version')
args=parser.parse_args()

result = subprocess.run(['mpicxx', '--showme:compile'],capture_output=True,text=True).stdout.strip()
openmpi_includes = [r.replace('-I','') for r in result.split()]


def get_version():
    """Return a PEP 440–compliant version string based on git tags and commit hash."""
    try:
        # Most recent tag reachable from HEAD
        tag = subprocess.run(
            ["git", "describe", "--tags", "--abbrev=0"],
            capture_output=True, text=True, check=True
        ).stdout.strip()

        # Number of commits since the tag
        commits = int(subprocess.run(
            ["git", "rev-list", f"{tag}..HEAD", "--count"],
            capture_output=True, text=True, check=True
        ).stdout.strip())

        # Current commit hash (short)
        sha = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, check=True
        ).stdout.strip()

        # Check if working tree is dirty
        dirty = bool(subprocess.run(
            ["git", "diff", "--quiet", "HEAD"],
            capture_output=True
        ).returncode)

        version = tag
        if commits > 0:
            version += f".post{commits}"
        if dirty:
            version += ".dev0"
        if commits > 0 or dirty:
            version += f"+{sha}"

        return version
    except subprocess.CalledProcessError:
        # fallback if not in a git repo
        return "0.0.0.dev0"

version = get_version()

f = open("setup.py",'w')
f.write(f"""
from setuptools import setup

setup(
    name='alamo',
    version='{version}',
    py_modules=['alamo'],
    install_requires=['cppyy','sympy'],
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
cppyy.add_include_path('{args.alamo_home}/{args.amrex}/include/')
cppyy.add_include_path('{args.alamo_home}/src/')
cppyy.load_library("{args.alamo_home}/{args.amrex}/lib/libamrex.so")
cppyy.load_library("{args.alamo_home}/lib/libalamo-{args.postfix}.so")

# Create proxy module
class AlamoModule(types.ModuleType):
    def __getattr__(self, name):
        return getattr(cppyy.gbl, name)

    def include(self, path):
        cppyy.include("{args.alamo_home}/src/" + path)

# Replace current module with AlamoModule instance
sys.modules[__name__].__class__ = AlamoModule
""")
