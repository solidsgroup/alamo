#
# Use strict mode only when this file is executed as a script. When sourced as
# installation documentation, strict mode would leak into the caller's shell;
# then a later failing command such as `make` can close the terminal session.
#
ALAMO_DEPS_SOURCED=0
if [ -n "${BASH_VERSION:-}" ] && [ "${BASH_SOURCE[0]}" != "$0" ]; then
    ALAMO_DEPS_SOURCED=1
fi
if [ "${ALAMO_DEPS_SOURCED}" -eq 0 ]; then
    set -eu -o pipefail
fi

#
# Make sure your system is up-to-date with standard build software
#
sudo apt update
sudo apt install build-essential g++ ca-certificates wget

#
# Use apt (or apt-get) to install these packages
# [ you should only need to do this once ]
#
sudo apt install libopenmpi-dev 
sudo apt install libeigen3-dev
sudo apt install libpng-dev

#
# [optional] Install these packages if compiling with clang
#
sudo apt install clang libstdc++-14-dev

#
# [optional] These are needed for regression test scripts and
#            other code infrastructure
#
sudo apt install lcov doxygen

# [optional] Needed for regression testing
sudo apt install python3-yt python3-matplotlib python3-numpy python3-pandas

# [optional] Needed to build documentation
sudo apt install python3-sphinx
sudo apt install python3-sphinx-rtd-theme python3-sphinx-design
sudo apt install python3-sphinx-copybutton python3-sphinxcontrib.bibtex
sudo apt install python3-xmltodict

# [optional] Needed for linting and code scraping
sudo apt install python3-clang
sudo apt install npm
npm install -g eclint

# [optional] Needed if compiling for NVidia GPUs.
# Use NVIDIA's versioned CUDA repository instead of Ubuntu's nvidia-cuda-toolkit,
# which can lag behind AMReX's minimum supported nvcc version.
CUDA_TOOLKIT_VERSION="${CUDA_TOOLKIT_VERSION:-12-9}"
CUDA_KEYRING_DEB="cuda-keyring_1.1-1_all.deb"
CUDA_REPO_URL="https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2404/x86_64"
if ! dpkg-query -W -f='${Status}' cuda-keyring 2>/dev/null | grep -q "install ok installed"; then
    tmpdir="$(mktemp -d)"
    wget -q -O "${tmpdir}/${CUDA_KEYRING_DEB}" "${CUDA_REPO_URL}/${CUDA_KEYRING_DEB}"
    sudo apt install -y "${tmpdir}/${CUDA_KEYRING_DEB}"
    rm -rf "${tmpdir}"
    sudo apt update
fi
sudo apt install "cuda-toolkit-${CUDA_TOOLKIT_VERSION}"
CUDA_HOME="/usr/local/cuda-${CUDA_TOOLKIT_VERSION/-/.}"
if [ -d "${CUDA_HOME}/bin" ]; then
    export CUDA_HOME
    export PATH="${CUDA_HOME}/bin:${PATH}"
    export LD_LIBRARY_PATH="${CUDA_HOME}/lib64:${LD_LIBRARY_PATH:-}"
    hash -r 2>/dev/null || true

    # GitHub Actions reads these files after this script exits, so later
    # workflow steps will see the same CUDA environment.
    if [ -n "${GITHUB_PATH:-}" ]; then echo "${CUDA_HOME}/bin" >> "${GITHUB_PATH}"; fi
    if [ -n "${GITHUB_ENV:-}" ]; then
        echo "CUDA_HOME=${CUDA_HOME}" >> "${GITHUB_ENV}"
        echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}" >> "${GITHUB_ENV}"
    fi

    cat <<EOF
CUDA toolkit configured for this session:
  CUDA_HOME=${CUDA_HOME}
  nvcc=$(command -v nvcc)

For persistent interactive use, add this to ~/.bashrc:
  export CUDA_HOME=${CUDA_HOME}
  export PATH="\$CUDA_HOME/bin:\$PATH"
  export LD_LIBRARY_PATH="\$CUDA_HOME/lib64:\${LD_LIBRARY_PATH:-}"
EOF
fi

# [optional] Needed for automated profiling
sudo apt install google-perftools libgoogle-perftools-dev
