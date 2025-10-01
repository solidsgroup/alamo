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
# CONFIGURATION
# =============
#
EIGEN_PREFIX="$(brew --prefix eigen)"
LIBPNG_PREFIX="$(brew --prefix libpng)"
LLVM_PREFIX="$(brew --prefix llvm)"
GCC_PREFIX="$(brew --prefix gcc)"

#
# This command updates the include path to be able to find
# the above libraries
#
# [ you need to do this in every new shell OR add to your shell config file (like .bashrc) ]
#
export CPLUS_INCLUDE_PATH="$EIGEN_PREFIX/include:$LIBPNG_PREFIX/include"

#
# Add the LLVM directory first in the path to avoid using Apple Clang
#
export PATH="$LLVM_PREFIX/bin:$PATH"

#
# In the alamo directory, run this command with any additional arguments. 
#
# [ you need to include these arguments every time you configure ]
#
./configure --comp clang++ --link "$LLVM_PREFIX/lib/c++" "$GCC_PREFIX/lib/gcc/current" "$LIBPNG_PREFIX/lib"

#
# Compile the code by running make
#
make

# TODO: remove this line before merge
otool -l bin/alamo-3d-clang++

#
# Run the unit test suite
#
./bin/test-3d-clang++

