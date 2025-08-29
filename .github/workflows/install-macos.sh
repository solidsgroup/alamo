set -eu -o pipefail 

#
# DEPENDENCIES
# ============
#

#
# Use the Homebrew package manager for macOS to install dependencies
# [ you should only need to do this once ]
#
brew update
brew reinstall --force --binaries gcc openmpi eigen libpng

#
# CONFIGURATION
# =============
#

#
# This command updates the include path to be able to find
# the above libraries
#
# [ you need to do this in every new shell OR add to your shell config file (like .bashrc) ]
#
export CPLUS_INCLUDE_PATH="$(brew --prefix)/include/c++:$(brew --prefix eigen)/include:$(brew --prefix libpng)/include"

#
# Create version-agnostic executables
#
(cd /opt/homebrew/bin && for tool in gcc g++ gfortran; do latest=$(ls -1 ${tool}-[0-9]* 2>/dev/null | sort -V | tail -n 1); [ -n "$latest" ] && ln -sf "$latest" "$tool"; done)

#
# Configure OpenMPI to use the correct compiler executables
#
export OMPI_CC=gcc
export OMPI_CXX=g++ 
export OMPI_FC=gfortran

#
# In the alamo directory, run this command with any additional arguments. 
#
# [ you need to include these arguments every time you configure ]
#
./configure --link "$(brew --prefix libpng)/lib" "$(dirname $(mpifort -print-file-name=libgfortran.dylib))" "$(dirname $(mpicxx -print-libgcc-file-name))"

#
# Compile the code by running make
#
make

#
# Run the unit test suite
#
./bin/test-3d-g++

