set -e
make clean
make info
./configure --dim $1 --build-amrex
make
