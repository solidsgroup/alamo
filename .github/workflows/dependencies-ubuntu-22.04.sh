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
sudo apt install clang

#
# These are needed for regression test scripts and other code infrastructure
#
sudo apt install python3-pip npm lcov doxygen

# Needed for regression testing
pip3 install yt matplotlib numpy pandas

# Needed to build documentation
pip3 install sphinx sphinx_rtd_theme sphinx_design
pip3 install sphinx-copybutton sphinxcontrib-bibtex xmltodict

# Needed for linting and code scraping
# You need to install python clang with apt, otherwise you must manually
# provide the clang library location in the python script.
sudo apt install python3-clang
sudo apt install npm
npm install -g eclint

