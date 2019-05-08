set -e
make clean
make info
./configure --dim $1 --amrex ./amrex$1d
make
