set -e
mkdir amrex3d
git clone https://github.com/AMReX-Codes/amrex.git
cd amrex
git checkout development
./configure --dim=3 --prefix=../amrex3d/ --debug=yes
make
make install
cd ..
make AMREX=./amrex3d 
rm -rf amrex3d amrex
