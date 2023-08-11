set -eu -o pipefail

# Install basic C++ compilers and build tools
sudo yum install -y gcc gcc-c++ make

# Install mpich
sudo yum install -y mpich mpich-devel mpich-autoload

# Manually add MPICH to paths
# If you are installing Alamo yourself, you will want to add these two 
# lines to your .bashrc file
export PATH=${PATH}:/usr/lib64/mpich/bin/
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/lib64/mpich/lib/

# Install eigen
sudo yum install -y eigen3-devel

# Install python and node.js
sudo yum install -y python3 python3-pip npm

# Requirements for regression test scripts
pip3 install yt matplotlib numpy pandas 

# Requirements for regression test scripts
pip3 install yt matplotlib numpy pandas

# Requirements for building documentation
pip3 install sphinx breathe m2r sphinx_rtd_theme linuxdoc sphinx_design sphinx-copybutton

# Editor config linter
npm install -g eclint
