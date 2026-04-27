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
brew reinstall --force --binaries llvm gfortran openmpi eigen libpng

#
# [optional] Needed for HDF5 output
#
brew reinstall --force --binaries hdf5-mpi

#
# CONFIGURATION
# =============
#
EIGEN_PREFIX="$(brew --prefix eigen)"
LIBPNG_PREFIX="$(brew --prefix libpng)"
LLVM_PREFIX="$(brew --prefix llvm)"
GCC_PREFIX="$(brew --prefix gcc)"
HDF5_PREFIX="$(brew --prefix hdf5-mpi)"

#
# This command updates the include path to be able to find
# the above libraries
#
# [ you need to do this in every new shell OR add to your shell config file (like .bashrc) ]
#
export CPLUS_INCLUDE_PATH="$EIGEN_PREFIX/include:$LIBPNG_PREFIX/include:$HDF5_PREFIX/include"

#
# Add the LLVM directory first in the path to avoid using Apple Clang
#
export PATH="$LLVM_PREFIX/bin:$PATH"

#
# In the alamo directory, run this command with any additional arguments. 
#
# [ you need to include these arguments every time you configure ]
#
./configure --comp clang++ --link "$LLVM_PREFIX/lib/c++" "$GCC_PREFIX/lib/gcc/current" "$LIBPNG_PREFIX/lib" "$HDF5_PREFIX/lib"

#
# Compile the code by running make
#
make

#
# Run the unit test suite
#
./bin/test-3d-clang++

