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
sudo apt install libopenmpi-dev 
sudo apt install libeigen3-dev
sudo apt install libpng-dev

#
# Install these packages if compiling with clang
#
sudo apt install clang libgfortran-14-dev libstdc++-14-dev

#
# These are needed for regression test scripts and other code infrastructure
#
sudo apt install lcov doxygen

# Needed for regression testing
sudo apt install python3-yt python3-matplotlib python3-numpy python3-pandas

# Needed to build documentation
sudo apt install python3-sphinx
sudo apt install python3-sphinx-rtd-theme python3-sphinx-design
sudo apt install python3-sphinx-copybutton python3-sphinxcontrib.bibtex
sudo apt install python3-xmltodict

# Needed for editorconfig linting
sudo apt install npm
npm install -g eclint

