set -e
make clean
make info
make AMREX=./amrex$1d MPICHFORT=mpichf90 bin/test
