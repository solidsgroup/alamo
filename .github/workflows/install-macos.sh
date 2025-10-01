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

#
# This command updates the include path to be able to find
# the above libraries
#
# [ you need to do this in every new shell OR add to your shell config file (like .bashrc) ]
#
export CPLUS_INCLUDE_PATH="$(brew --prefix eigen)/include:$(brew --prefix libpng)/include"

#
# Create version-agnostic executables
#
ls -alh "$(brew --prefix llvm)/bin"

#
# In the alamo directory, run this command with any additional arguments. 
#
# [ you need to include these arguments every time you configure ]
#
./configure --comp clang++ --link "$(brew --prefix llvm)/lib/c++" "$(brew --prefix gcc)/lib/gcc/current" "$(brew --prefix libpng)/lib"

#
# Compile the code by running make
#
make

#
# Run the unit test suite
#
./bin/test-3d-clang++

