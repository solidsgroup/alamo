set -eu -o pipefail

sudo apt-get update

sudo apt-get install -y --no-install-recommends \
  build-essential \
  g++ gfortran \
  libmpich-dev libmpich12 libeigen3-dev \
  python3-pip

pip3 install yt matplotlib numpy
