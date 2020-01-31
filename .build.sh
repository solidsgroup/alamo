set -e
./configure --dim $1 --build-amrex --no-diff
make
