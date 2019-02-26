set -e
bin/test
mpirun -np 2 bin/test
