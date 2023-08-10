set -eu -o pipefail

# Install basic C++ compilers and build tools
sudo yum install -y gcc gcc-c++ make

# Install mpich
sudo yum install -y mpich2 mpich2-devel mpich-autoload

# Install python and node.js
sudo yum install -y python3 npm

# Install pip
sudo python3 -m pip install pip

# Requirements for regression test scripts
pip3 install yt matplotlib numpy pandas 

# Requirements for regression test scripts
pip3 install yt matplotlib numpy pandas

# Requirements for building documentation
pip3 install sphinx breathe m2r sphinx_rtd_theme linuxdoc sphinx_design sphinx-copybutton

# Editor config linter
npm install -g eclint
