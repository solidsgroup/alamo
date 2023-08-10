set -eu -o pipefail

sudo yum install -y gcc gcc-c++ make

sudo yum install mpich2 mpich2-devel mpich-autoload

sudo yum install rh-python36 npm

sudo python3 -m pip install pip

# Requirements for regression test scripts
pip3 install yt matplotlib numpy pandas 

# Requirements for regression test scripts
pip3 install yt matplotlib numpy pandas

# Requirements for building documentation
pip3 install sphinx breathe m2r sphinx_rtd_theme linuxdoc sphinx_design sphinx-copybutton

# Editor config linter
npm install -g eclint
