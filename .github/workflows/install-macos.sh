set -eu -o pipefail 

#
# DEPENDENCIES
# ============
#

#
# Use Mac OS HomeBrew system to install mpich and eigen
# [ you should only need to do this once ]
#
brew update
brew install gfortran || true
brew install mpich || true
brew install eigen || true
brew install libpng || true


#
# If your system is unable to find gfortran (for instance, if you are
# a github action), you can try creating an alias using the following
# command:
#
# ln -s /opt/homebrew/bin/gfortran-14 /opt/homebrew/bin/gfortran

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
export CPLUS_INCLUDE_PATH=$(brew --prefix)/include:$(brew --prefix libpng)/include

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


#
# PYTHON [OPTIONAL]
# =================
#
# These are not required tu run alamo, but are needed to use alamo scripts such as
# the regression test script.
# The following commands work on the Github VM, but your configuration may vary.
#
# (If you already have another way of installing python packages, use it - the
# "break-system-packages" is kind of dangerous and is only needed for this CI script)
#

#
# Install packages needed for regression test script
#
pip3 install sympy yt matplotlib numpy pandas --break-system-packages --user

