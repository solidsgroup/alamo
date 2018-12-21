set -e
mkdir amrex$1
git clone https://github.com/AMReX-Codes/amrex.git
cd amrex
if [ $# -eq 2 ]
then
    git checkout $2
fi
./configure --dim=$1 --prefix=../amrex$1d/ --debug=yes
make
make install
