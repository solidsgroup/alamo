set -eu -o pipefail 

#
# DEPENDENCIES
# ============
#

#
# Use apt (or apt-get) to install these packages
# [ you should only need to do this once ]
#
sudo apt install mpich
sudo apt install libeigen3-dev
sudo apt install libpng-dev

#
# CONFIGURATION
# =============
#

#
# In the alamo directory, run this command with any additional arguments. 
#
./configure 

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

#
# Install packages needed for regression test script
#
pip3 install sympy yt matplotlib numpy pandas

