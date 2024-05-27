set -eu -o pipefail 

#
# Use Mac OS HomeBrew system to 
#
brew install mpich

brew install eigen

export CPLUS_INCLUDE_PATH=/opt/homebrew/include

./configure --no-png

make
