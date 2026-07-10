make clean
make -j8 ./bin/hydro
./scripts/runtests.py ./tests/FlowSource1D
