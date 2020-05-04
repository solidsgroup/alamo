set -e
./configure --dim $1 --no-diff --create-gfortran-link
make
