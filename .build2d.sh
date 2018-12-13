set -e
make clean
make AMREX=./amrex2d MPICHFORT=mpichf90 bin/test

