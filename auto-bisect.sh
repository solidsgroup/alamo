#!/bin/bash
# ./auto-bisect.sh
# Build and run the hydro test for git bisect automation

# Exit immediately if any command fails
set -e

# 1. Clean and build the code (adjust for your build system)
make clean
make -j8  

# 2. Run the test
# Redirect both stdout and stderr to a log file
./bin/hydro-2d-clang++ tests/FlowDrivenCavity/input > bisect.log 2>&1

# 3. Check for failure
# If the program exits with non-zero, mark as bad
if [ $? -ne 0 ]; then
    echo "Test failed on commit $(git rev-parse HEAD)"
    exit 1   # non-zero exit signals "bad" to git bisect
else
    echo "Test passed on commit $(git rev-parse HEAD)"
    exit 0   # zero exit signals "good" to git bisect
fi
