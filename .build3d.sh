set -e
make clean
make AMREX=./amrex3d MPICHFORT=mpichf90 bin/fem

