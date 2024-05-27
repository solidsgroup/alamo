set -eu -o pipefail 

brew install mpich

brew install eigen

export CPLUS_INCLUDE_PATH=/opt/homebrew/include

./configure --no-png

