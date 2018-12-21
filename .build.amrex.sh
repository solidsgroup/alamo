set -e
mkdir amrex$1
git clone https://github.com/AMReX-Codes/amrex.git
cd amrex
git checkout 18.11
./configure --dim=$1 --prefix=../amrex$1d/ --debug=yes
make
make install
