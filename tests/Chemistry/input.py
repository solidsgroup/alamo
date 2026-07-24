import csv
from pathlib import Path
import sys

import numpy


ROOT = Path(__file__).resolve().parents[2]
MECHANISM = Path(__file__).with_name("h2o2.yaml")
CONFIGURATION = ROOT / ".make/Makefile.pre.conf"


def yaml_configuration_error(problem):
    message = f"""
ERROR: tests/Chemistry requires YAML-enabled Alamo Python bindings.

{problem}

From the Alamo repository, activate the Python environment used for testing,
then configure, build, and install alamopy with:

  ./configure --comp=g++ --debug --yaml
  make -j py
  python -m pip install -e .

Then rerun:

  python tests/Chemistry/input.py
"""
    print(message.strip(), file=sys.stderr)
    raise SystemExit(2)


try:
    yaml_enabled = "-DALAMO_YAML" in CONFIGURATION.read_text()
except OSError:
    yaml_enabled = False

if not yaml_enabled:
    yaml_configuration_error("The active Alamo build was not configured with --yaml.")

try:
    import cppyy
    import alamo
except (ImportError, OSError, RuntimeError) as error:
    yaml_configuration_error(f"alamopy is not available: {error}")

# These headers contain YAML-dependent inline functions, so Cling needs the
# same include path and feature definition as the YAML-enabled Alamo library.
cppyy.add_include_path(str(ROOT / "ext/yaml-cpp/include"))
cppyy.cppdef("#define ALAMO_YAML")

alamo.include("Util/Util.H")
alamo.include("IO/ParmParse.H")
alamo.include("Model/Gas/Gas.H")
alamo.include("Model/Chemistry/FiniteRate.H")


MOLECULAR_WEIGHTS = [
    "2.016_g/mol",
    "1.008_g/mol",
    "15.999_g/mol",
    "31.998_g/mol",
    "17.007_g/mol",
    "18.015_g/mol",
    "33.006_g/mol",
    "34.014_g/mol",
    "39.948_g/mol",
    "28.013_g/mol",
]
SPECIES = ["H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "AR", "N2"]
LJ_DIAMETERS = [
    "2.92_ang",
    "2.05_ang",
    "2.75_ang",
    "3.458_ang",
    "2.75_ang",
    "2.605_ang",
    "3.458_ang",
    "3.458_ang",
    "3.33_ang",
    "3.621_ang",
]
LJ_WELL_DEPTHS = [
    "38.0_K",
    "145.0_K",
    "80.0_K",
    "107.4_K",
    "80.0_K",
    "572.4_K",
    "107.4_K",
    "107.4_K",
    "136.5_K",
    "97.53_K",
]
EXPECTED_WDOT = numpy.array(
    [
        -4.025974566628345e-06,
        2.0129872913589804e-06,
        1.107058498505617e-16,
        -6.390036392279833e-05,
        0.0,
        0.0,
        6.591335119795699e-05,
        0.0,
        0.0,
        0.0,
    ]
)
EXPECTED_QDOT = -461.905471792037


def cpp_vector(cpp_type, values):
    result = alamo.std.vector[cpp_type]()
    for value in values:
        result.push_back(value)
    return result


def add_string(pp, name, value):
    pp.add(name, alamo.std.string(str(value)))


def add_strings(pp, name, values):
    pp.addarr(name, cpp_vector("std::string", values))


alamo.Util.Initialize()
try:
    # The mechanism parser and rate coefficients use kmol-based quantities.
    alamo.Unit.setAmountUnit("kmol")
    alamo.Set.Constant.SetGlobalConstants()

    pp = alamo.IO.ParmParse()
    add_strings(pp, "gas.mw", MOLECULAR_WEIGHTS)
    add_string(pp, "gas.thermo.type", "nasa7")
    add_string(pp, "gas.thermo.nasa7.type", "yaml")
    add_string(pp, "gas.thermo.nasa7.yaml", MECHANISM)
    add_string(pp, "gas.transport.type", "mixture_averaged")
    add_string(pp, "gas.transport.mixture_averaged.type", "LJ")
    add_strings(pp, "gas.transport.mixture_averaged.LJdiameter", LJ_DIAMETERS)
    add_strings(pp, "gas.transport.mixture_averaged.LJwelldepth", LJ_WELL_DEPTHS)
    add_string(pp, "gas.eos.type", "tpg")

    gas = alamo.Model.Gas.Gas(pp, "gas")

    chemistry_prefix = "chemistry.model.finite_rate"
    add_string(pp, f"{chemistry_prefix}.yaml", MECHANISM)
    chemistry = alamo.Model.Chemistry.FiniteRate()
    alamo.Model.Chemistry.FiniteRate.Parse(
        chemistry, alamo.IO.ParmParse(chemistry_prefix), len(MOLECULAR_WEIGHTS)
    )

    partial_densities = alamo.Model.Chemistry.SpeciesArray()
    partial_densities[0] = 0.016378785527602464
    partial_densities[3] = 0.12998223693259516

    source = chemistry.ComputeChemistrySources(
        101325.0, 1000.0, partial_densities, 0.0, gas
    )
    wdot = numpy.array([source.first[n] for n in range(len(MOLECULAR_WEIGHTS))])

    numpy.testing.assert_allclose(
        wdot,
        EXPECTED_WDOT,
        rtol=2.0e-8,
        atol=2.0e-12,
        err_msg="finite-rate species sources differ from the Cantera reference",
    )
    numpy.testing.assert_allclose(
        source.second,
        EXPECTED_QDOT,
        rtol=2.0e-8,
        atol=2.0e-8,
        err_msg="finite-rate heat release differs from the Cantera reference",
    )
    numpy.testing.assert_allclose(
        wdot.sum(),
        0.0,
        rtol=0.0,
        atol=2.0e-12,
        err_msg="finite-rate chemistry does not conserve mass",
    )

    empty_state = alamo.Model.Chemistry.SpeciesArray()
    empty_source = chemistry.ComputeChemistrySources(
        101325.0, 1000.0, empty_state, 0.0, gas
    )
    numpy.testing.assert_array_equal(
        [empty_source.first[n] for n in range(len(MOLECULAR_WEIGHTS))],
        numpy.zeros(len(MOLECULAR_WEIGHTS)),
    )
    assert empty_source.second == 0.0

    output = Path(__file__).with_name("output")
    output.mkdir(exist_ok=True)
    rows = []
    for species, reference, actual in zip(SPECIES, EXPECTED_WDOT, wdot):
        rows.append(
            {
                "quantity": f"wdot[{species}]",
                "units": "kg/m^3/s",
                "reference": reference,
                "alamo": actual,
                "relative_error": (
                    abs(actual - reference) / abs(reference)
                    if reference != 0.0
                    else ""
                ),
            }
        )
    rows.append(
        {
            "quantity": "qdot",
            "units": "W/m^3",
            "reference": EXPECTED_QDOT,
            "alamo": source.second,
            "relative_error": (
                abs(source.second - EXPECTED_QDOT) / abs(EXPECTED_QDOT)
            ),
        }
    )

    csv_path = output / "chemistry_sources.csv"
    with csv_path.open("w", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    print(f"CSV: {csv_path}")
    print(f"Mass conservation residual: {wdot.sum():.10e} kg/m^3/s")
    print("PASS: Chemistry sources match the reference data")
finally:
    alamo.Util.Finalize()
