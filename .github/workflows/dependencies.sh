set -eu -o pipefail

sudo apt-get update

sudo apt-get install -y --no-install-recommends \
  build-essential \
  g++ gfortran \
  libmpich-dev libmpich12 mpich libeigen3-dev libpng12-dev \
  python3-pip npm lcov

# Requirements for regression test scripts
pip3 install yt matplotlib numpy pandas

# Requirements for building documentation
pip3 install sphinx breathe m2r sphinx_rtd_theme linuxdoc sphinx_design sphinx-copybutton sphinxcontrib-bibtex

npm install -g eclint
