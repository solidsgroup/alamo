set -eu -o pipefail 

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
# Run the unit test suite in serial using the regression test script
#
scripts/runtests.py --dim=3 --serial tests/Unit
