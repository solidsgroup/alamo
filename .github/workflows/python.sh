set -eu -o pipefail 

# Set up a virtual environment and activate
mkdir alamovenv 
python -m venv alamovenv
source alamovenv/bin/activate


# Install alamo using pip (-e for editable mode is optional)
pip install -e .

# Run the python regression test
./scripts/runtests.sh --python tests/Py 
