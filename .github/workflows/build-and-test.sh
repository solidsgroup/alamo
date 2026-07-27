set -eu -o pipefail 

#
# In the alamo directory, run this command with any additional arguments. 
#
compiler=clang++
for arg in "$@"; do
	case "$arg" in
		--comp=*) compiler="${arg#--comp=}" ;;
	esac
done
./configure "$@"

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
scripts/runtests.py --dim=3 --serial --comp="$compiler" tests/Unit
