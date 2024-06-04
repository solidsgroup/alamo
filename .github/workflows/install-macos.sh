set -eu -o pipefail 

#
# Use Mac OS HomeBrew system to install mpich and eigen
# [ you should only need to do this once ]
#
brew install mpich
brew install eigen
brew install libpng

#
# This command updates the include path to be able to find
# the above libraries
#
# [ you need to do this in every new shell OR add to your shell config file (like .bashrc) ]
#
export CPLUS_INCLUDE_PATH=${CPLUS_INCLUDE_PATH}:$(brew --prefix)/include:$(brew --prefix libpng)/include

#
# In the alamo directory, run this command with any additional arguments. 
#
# [ you need to include these arguments every time you configure ]
#
./configure --macos --link $(brew --prefix)/lib/gcc/current/ $(brew --prefix libpng)/lib/

#
# Compile the code by running make
#
make

#
# Executables should now be available under ./bin
#
ls ./bin/

#
# Run the unit test suite
#
./bin/test-3d-g++

