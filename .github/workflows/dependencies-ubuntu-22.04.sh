set -eu -o pipefail 

#
# Make sure your system is up-to-date with standard build software
#
sudo apt update
sudo apt install build-essential g++ gfortran

#
# Use apt (or apt-get) to install these packages
# [ you should only need to do this once ]
#
sudo apt install libmpich-dev libmpich12 mpich
sudo apt install libeigen3-dev
sudo apt install libpng-dev

#
# Install these packages if compiling with clang
#
sudo apt install clang libstdc++-14-dev libgfortran-14-dev

#
# These are needed for regression test scripts and other code infrastructure
#
sudo apt install python3-pip npm lcov doxygen

# Needed for regression testing
pip3 install yt matplotlib numpy pandas

# Needed to build documentation
pip3 install sphinx breathe m2r sphinx_rtd_theme linuxdoc sphinx_design
pip3 install sphinx-copybutton sphinxcontrib-bibtex xmltodict

# Needed for editorconfig linting
npm install -g eclint

