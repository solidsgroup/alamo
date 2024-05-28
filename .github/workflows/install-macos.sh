set -eu -o pipefail 

#
# Use Mac OS HomeBrew system to install mpich and eigen
#
brew install mpich
brew install eigen

#
# This command updates the include path to be able to find
# the above libraries
#
export CPLUS_INCLUDE_PATH=$(brew --prefix)/include

#
# In the alamo directory, run this command with any additiona/
# arguments. On MacOS, PNG is not supported.
#
./configure --macos --link=$(brew --prefix)/lib/gcc/current/ 

#
# Compile the code by running make
#
make

#
# Executables should now be available under ./bin
#
ls ./bin/
