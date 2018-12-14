mkdir amrex2d amrex3d
git clone https://github.com/AMReX-Codes/amrex.git
cd amrex
git checkout 18.11
./configure --dim=2 --prefix=../amrex2d/ --debug=yes
make
make install
make clean
rm -rf tmp_build_dir
./configure --dim=3 --prefix=../amrex3d/ --debug=yes
make
make install
cd ..
