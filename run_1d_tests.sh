make clean
make -j8 ./bin/hydro
./bin/hydro-2d-g++ tests/FlowSource1D/input m0.ic.expression.region0 = "10.0" u0.ic.expression.region0 = "0.0" q.ic.expression.region0 = "-0.001" plot_file=tests/FlowSource1D/output_2026-07-10_12.50.23_pepe_mass-source-term
./bin/hydro-2d-g++ tests/FlowSource1D/input m0.ic.expression.region0 = "10.0" u0.ic.expression.region0 = "-1.0" q.ic.expression.region0 = "0.0" plot_file=tests/FlowSource1D/output_2026-07-10_12.50.23_pepe_momemtum-source-term
./bin/hydro-2d-g++ tests/FlowSource1D/input m0.ic.expression.region0 = "10.0" u0.ic.expression.region0 = "0.0" q.ic.expression.region0 = "-10.0" plot_file=tests/FlowSource1D/output_2026-07-10_12.50.23_pepe_energy-source-term
