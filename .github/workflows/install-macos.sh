set -eu -o pipefail

#
# Assumes Xcode Command Line Tools and Homebrew are installed.
#   - Xcode CLT: xcode-select --install
#   - Homebrew:  https://brew.sh
#

#
# DEPENDENCIES
#
brew update
brew install openmpi eigen libpng

#
# CONFIGURE
#
./configure --comp clang++

#
# BUILD
#
make

#
# TEST
#
./bin/test-3d-clang++
