set -eu -o pipefail 

# Must use debug mode using g++ for python
./configure --comp=g++ --debug

# Build python package
make py

# Set up a virtual environment and activate
mkdir alamovenv 
python -m venv alamovenv
source alamovenv/bin/activate

# Install alamo using pip (-e for editable mode is optional)
pip install -e .

# Run the python regression test
./scripts/runtests.py --python --comp=g++ tests/Py 
