#!/usr/bin/env bash
# Source this to use the project-local CUDA toolkit.
CUDA_HOME="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/.local/cuda"
export CUDA_HOME
export CUDA_PATH="${CUDA_HOME}"
export PATH="${CUDA_HOME}/bin:${PATH}"
export LD_LIBRARY_PATH="${CUDA_HOME}/lib64:${CUDA_HOME}/lib:${CUDA_HOME}/lib/x86_64-linux-gnu:${CUDA_HOME}/lib/cuda/lib64:${LD_LIBRARY_PATH:-}"
export LIBRARY_PATH="${CUDA_HOME}/lib/stubs:${CUDA_HOME}/lib64/stubs:${LIBRARY_PATH:-}"
