set -eu -o pipefail 

#
# Generate documentation with
#
./configure --dim=2 --comp=clang++
make -j2
./configure --dim=3 --comp=clang++
make -j2
make docs

# Open 
#    docs/build/html/index.html
# in a browser to view the built documentation
