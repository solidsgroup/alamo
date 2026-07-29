import ctypes
import math
import os
from pathlib import Path
import sys
import time

SUITE_START = time.perf_counter()
os.environ.setdefault("MPLCONFIGDIR", "/tmp/alamo-matplotlib")

import numpy
import pylab


ROOT = Path(__file__).resolve().parents[2]
MECHANISM = Path(__file__).with_name("h2o2.yaml")
CONFIGURATION = ROOT / ".make/Makefile.pre.conf"
OUTPUT = Path(__file__).with_name("output")


def yaml_configuration_error(problem):
    message = f"""
ERROR: tests/Chemistry requires YAML-enabled Alamo Python bindings.

{problem}

From the Alamo repository, activate the Python environment used for testing,
then configure, build, and install alamopy with:

  ./configure --comp=g++ --debug --yaml
  make -j
  make -j py
  python -m pip install -e .

Then rerun:

  scripts/runtests.py tests/Chemistry --python
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
    import alamo
except (ImportError, OSError, RuntimeError) as error:
    yaml_configuration_error(f"alamopy is not available: {error}")

for header in [
    "Util/Util.H",
    "IO/ParmParse.H",
    "Model/Gas/Gas.H",
    "Model/Chemistry/Chemistry.H",
    "Model/Chemistry/FiniteRate.H",
    "Model/Chemistry/Frozen.H",
    "Model/Chemistry/GrossModel.H",
    "Model/Chemistry/GrossModel_Aluminized.H",
]:
    alamo.include(header)


FINITE_SPECIES = [
    "H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "AR", "N2"
]
GROSS_MODEL_SPECIES = ["AP", "HTPB", "Mono", "Premixed", "Primary", "Final"]
GROSS_MODEL_ALUMINIZED_SPECIES = GROSS_MODEL_SPECIES + ["Al_gas", "Al2O3_gas"]

FINITE_MOLECULAR_WEIGHTS = [
    "2.016_g/mol", "1.008_g/mol", "15.999_g/mol", "31.998_g/mol",
    "17.007_g/mol", "18.015_g/mol", "33.006_g/mol", "34.014_g/mol",
    "39.948_g/mol", "28.013_g/mol",
]
LJ_DIAMETERS = [
    "2.92_ang", "2.05_ang", "2.75_ang", "3.458_ang", "2.75_ang",
    "2.605_ang", "3.458_ang", "3.458_ang", "3.33_ang", "3.621_ang",
]
LJ_WELL_DEPTHS = [
    "38.0_K", "145.0_K", "80.0_K", "107.4_K", "80.0_K", "572.4_K",
    "107.4_K", "107.4_K", "136.5_K", "97.53_K",
]

FINITE_INITIAL_PARTIAL_DENSITIES = numpy.array(
    [0.016378785527602464, 0.0, 0.0, 0.12998223693259516,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
)
GROSS_MODEL_INITIAL_PARTIAL_DENSITIES = numpy.array(
    [0.28154378279353426, 0.03530772294877530, 0.0, 0.0, 0.0, 0.0]
)
GROSS_MODEL_ALUMINIZED_INITIAL_PARTIAL_DENSITIES = numpy.append(
    GROSS_MODEL_INITIAL_PARTIAL_DENSITIES, [0.05, 0.03]
)

EXPECTED_WDOT = numpy.array(
    [-4.025974566628345e-06, 2.0129872913589804e-06,
     1.107058498505617e-16, -6.390036392279833e-05, 0.0, 0.0,
     6.591335119795699e-05, 0.0, 0.0, 0.0]
)
EXPECTED_QDOT = -461.905471792037
EXPECTED_SUBSTEP_WDOT = numpy.array(
    [-2.6446576525951393e-06, 0.0, 3.037109989146865e-06,
     -5.413475945292845e-05, 3.8658137269011934e-06,
     1.0938310771028209e-05, 3.885837065026453e-05,
     7.981216158242902e-08, 0.0, 0.0]
)
EXPECTED_SUBSTEP_QDOT = 83.43352674578557

# Cantera 3.2 IdealGasConstPressureReactor reference sampled every 10 us.
CANTERA_TIMES = numpy.linspace(0.0, 2.0e-4, 21)
CANTERA_ACTIVE_MASS_FRACTIONS = numpy.array(
    [
        [0.11190674437968359, 0.0, 0.0, 0.8880932556203164, 0.0, 0.0, 0.0, 0.0],
        [0.11190674333629981, 1.6010985051590677e-10, 2.729086132004047e-10, 0.8880932397763611, 9.563189178230012e-11, 4.637249759846221e-09, 1.1697256137378567e-08, 2.4182728938495634e-11],
        [0.11190673964590007, 4.798887039728584e-10, 8.406517459859773e-10, 0.8880931869885857, 2.9176609403271535e-10, 2.5077915760024245e-08, 4.650977632437549e-08, 1.655155081466767e-10],
        [0.11190673061031398, 1.1256139378080005e-09, 1.986933304461735e-09, 0.888093059716687, 6.87854381642264e-10, 7.740667533485318e-08, 1.278630270605499e-07, 6.02894952263017e-10],
        [0.11190671076540318, 2.4318471274544905e-09, 4.30556320872896e-09, 0.8880927818874947, 1.4893633819751628e-09, 1.9421034247265078e-07, 3.0321171591663386e-07, 1.6982698734290899e-09],
        [0.11190666902993597, 5.077478172511156e-09, 9.001423931612401e-09, 0.888092199319878, 3.1139359233843083e-09, 4.4160454132620636e-07, 6.686647476368403e-07, 4.1880594038015955e-09],
        [0.11190658289692659, 1.0443130967724033e-08, 1.8524526333265225e-08, 0.888090999251989, 6.413842194274512e-09, 9.539977633174247e-07, 1.418889217510191e-06, 9.582604431063898e-09],
        [0.11190640650969884, 2.1348980886471888e-08, 3.787824014141533e-08, 0.8880885460757859, 1.314206544900258e-08, 2.0057393083768606e-06, 2.9482821165319805e-06, 2.1023804580792404e-08],
        [0.1119060458462003, 4.3607242425618385e-08, 7.736918179287824e-08, 0.8880835432727396, 2.6961140620166975e-08, 4.161397131736082e-06, 6.056324010013824e-06, 4.5222355098514084e-08],
        [0.11190530567108958, 8.941032607521893e-08, 1.5859700770759399e-07, 0.8880733262668692, 5.575997271767918e-08, 8.601774125823993e-06, 1.2365315096183084e-05, 9.720551565420519e-08],
        [0.11190377026193618, 1.852327628587626e-07, 3.2837693444059957e-07, 0.8880523357374576, 1.1752932449931751e-07, 1.7876404684068712e-05, 2.5173127815805545e-05, 2.1332908892875346e-07],
        [0.11190051004537586, 3.9243064929179983e-07, 6.948434352480188e-07, 0.8880086201656685, 2.576382680233757e-07, 3.783393045114412e-05, 5.119929281583254e-05, 4.916533398071107e-07],
        [0.11189324649944524, 8.707413813750577e-07, 1.5379540635175393e-06, 0.8879149381385929, 6.106038503973578e-07, 8.344619308113996e-05, 0.00010411832411266323, 1.2315454765090418e-06],
        [0.11187540572473335, 2.12410234337126e-06, 3.733821363898278e-06, 0.8877019139211862, 1.6808619772346244e-06, 0.00020078809461258692, 0.0002109057526621935, 3.4477211256843932e-06],
        [0.1118222590004655, 6.275102296641505e-06, 1.0940983757325725e-05, 0.8871533568623987, 6.033563041984108e-06, 0.0005772919860937222, 0.0004132021863830178, 1.0640315567363552e-05],
        [0.11159999792686527, 2.640960581492398e-05, 4.5788359272319566e-05, 0.8853070237829361, 3.111951480145264e-05, 0.0022896243613752944, 0.0006690947294005572, 3.094171953768456e-05],
        [0.1100302486541149, 0.00021223380442126906, 0.00038285460228379783, 0.8738451518097383, 0.00026482370660678645, 0.014513350933452141, 0.0006998812053173501, 5.145528406846902e-05],
        [0.02461867984931559, 0.01338577409063175, 0.07591736541540003, 0.17390756149742154, 0.10993549636646685, 0.6021390972679007, 9.301362411999944e-05, 3.0118887493816884e-06],
        [0.023825096387086166, 0.008295863520261654, 0.055852595570013486, 0.13359573607801514, 0.1390452243857622, 0.6392988740492612, 8.266773377677175e-05, 3.942275829402184e-06],
        [0.023296327517989934, 0.007354856532221269, 0.05097022145591259, 0.1252615333527489, 0.14231099839388528, 0.6506963216683095, 0.00010382549367419575, 5.915585264900898e-06],
        [0.023213203922539963, 0.007230186158089977, 0.05030084299627257, 0.12410934889698887, 0.14267171946957258, 0.6523609243568937, 0.00010748630189849035, 6.287897749864566e-06],
    ]
)
CANTERA_MASS_FRACTIONS = numpy.pad(
    CANTERA_ACTIVE_MASS_FRACTIONS, ((0, 0), (0, 2))
)
CANTERA_TEMPERATURES = numpy.array(
    [1000.0, 1000.0000078961748, 1000.0000806894886,
     1000.0002843225067, 1000.0007524669115, 1000.0017563420429,
     1000.0038473662493, 1000.0081512698921, 1000.0169854342088,
     1000.035199960352, 1000.073280011222, 1000.155329339466,
     1000.3432343961042, 1000.8279595137572, 1002.385282058233,
     1009.3933927541337, 1055.2312085798892, 2376.51545692503,
     2997.3317879691244, 3137.5919460842397, 3156.9585418678357]
)

# Backward Euler with ten substeps provides the gas-only GrossModel reference;
# condensed-phase heat is exercised separately through LowMach PhaseChange.
GROSS_MODEL_FINAL_MASS_FRACTIONS = numpy.array(
    [0.0, 0.0, 0.466528083, 0.0138619938, 0.207818794, 0.311791129]
)
GROSS_MODEL_FINAL_TEMPERATURE = 3313.7700943090126

CHEMISTRY_TYPE = alamo.Model.Chemistry.Chemistry[
    alamo.Model.Chemistry.Frozen,
    alamo.Model.Chemistry.FiniteRate,
    alamo.Model.Chemistry.GrossModel,
    alamo.Model.Chemistry.GrossModel_Aluminized,
]


def cpp_vector(cpp_type, values):
    result = alamo.std.vector[cpp_type]()
    for value in values:
        result.push_back(value)
    return result


def add_string(pp, name, value):
    pp.add(name, alamo.std.string(str(value)))


def add_strings(pp, name, values):
    pp.addarr(name, cpp_vector("std::string", values))


def make_finite_gas(pp):
    prefix = "finite_gas"
    add_strings(pp, f"{prefix}.mw", FINITE_MOLECULAR_WEIGHTS)
    add_string(pp, f"{prefix}.thermo.type", "nasa7")
    add_string(pp, f"{prefix}.thermo.nasa7.type", "yaml")
    add_string(pp, f"{prefix}.thermo.nasa7.yaml", MECHANISM)
    add_string(pp, f"{prefix}.transport.type", "mixture_averaged")
    add_string(pp, f"{prefix}.transport.mixture_averaged.type", "LJ")
    add_strings(
        pp, f"{prefix}.transport.mixture_averaged.LJdiameter", LJ_DIAMETERS
    )
    add_strings(
        pp, f"{prefix}.transport.mixture_averaged.LJwelldepth", LJ_WELL_DEPTHS
    )
    add_string(pp, f"{prefix}.eos.type", "tpg")
    return alamo.Model.Gas.Gas(pp, prefix)


def make_gross_model_gas(pp):
    prefix = "gross_model_gas"
    add_strings(pp, f"{prefix}.mw", ["26.0_g/mol"] * len(GROSS_MODEL_SPECIES))
    add_string(pp, f"{prefix}.thermo.type", "gross_model")
    add_string(pp, f"{prefix}.transport.type", "gross_model")
    add_string(pp, f"{prefix}.eos.type", "gross_model")
    return alamo.Model.Gas.Gas(pp, prefix)


def make_gross_model_aluminized_gas(pp):
    prefix = "gross_model_aluminized_gas"
    add_strings(
        pp, f"{prefix}.mw",
        ["26.0_g/mol"] * 6 + ["26.9815385_g/mol", "101.96_g/mol"],
    )
    add_string(pp, f"{prefix}.thermo.type", "gross_model")
    add_string(pp, f"{prefix}.transport.type", "gross_model")
    add_string(pp, f"{prefix}.eos.type", "gross_model")
    return alamo.Model.Gas.Gas(pp, prefix)


def make_integrator(pp, prefix, model, solver, nsubsteps, nspecies):
    add_string(pp, f"{prefix}.model.type", model)
    if model == "finite_rate":
        add_string(pp, f"{prefix}.model.finite_rate.yaml", MECHANISM)
    add_string(pp, f"{prefix}.solver.type", solver)
    add_string(pp, f"{prefix}.solver.{solver}.nsubsteps", nsubsteps)
    chemistry = CHEMISTRY_TYPE()
    chemistry.Define(nspecies)
    CHEMISTRY_TYPE.Parse(chemistry, alamo.IO.ParmParse(prefix))
    return chemistry


def make_finite_source_model(pp, prefix, substeps=None):
    add_string(pp, f"{prefix}.yaml", MECHANISM)
    if substeps is not None:
        add_string(pp, f"{prefix}.integration", "substep")
        add_string(pp, f"{prefix}.nsubsteps", substeps)
    chemistry = alamo.Model.Chemistry.FiniteRate()
    alamo.Model.Chemistry.FiniteRate.Parse(
        chemistry, alamo.IO.ParmParse(prefix), len(FINITE_SPECIES)
    )
    return chemistry


def species_array(values):
    result = alamo.Model.Chemistry.SpeciesArray()
    for index, value in enumerate(values):
        result[index] = value
    return result


def integrate_composition(
    chemistry, gas, initial_partial_densities, timestep, output_times
):
    nspecies = len(initial_partial_densities)
    mixture_density = float(initial_partial_densities.sum())
    state = species_array(initial_partial_densities)
    temperature = ctypes.c_double(1000.0)
    mass_fractions = [initial_partial_densities / mixture_density]
    temperatures = [temperature.value]

    for time_lo, time_hi in zip(output_times[:-1], output_times[1:]):
        interval = float(time_hi - time_lo)
        nsteps = round(interval / timestep)
        assert numpy.isclose(nsteps * timestep, interval)
        for _ in range(nsteps):
            result = chemistry.Advance(
                timestep, 101325.0, mixture_density, nspecies,
                state, temperature, gas,
            )
            assert result.converged, (
                f"chemistry integration failed: dt={result.failed_dt}, "
                f"residual={result.residual_norm}"
            )

        densities = numpy.array([state[n] for n in range(nspecies)])
        mass_fractions.append(densities / densities.sum())
        temperatures.append(temperature.value)

    return numpy.array(mass_fractions), numpy.array(temperatures)


def validate_mass_fractions(name, mass_fractions):
    numpy.testing.assert_allclose(
        mass_fractions.sum(axis=1), 1.0, rtol=0.0, atol=2.0e-14,
        err_msg=f"{name} mass fractions do not sum to one",
    )


def run_case(name, function):
    start = time.perf_counter()
    try:
        result = function()
    except Exception:
        elapsed = time.perf_counter() - start
        print(f"{name:.<62}[FAIL] {elapsed:8.3f} s", flush=True)
        raise
    elapsed = time.perf_counter() - start
    print(f"{name:.<62}[PASS] {elapsed:8.3f} s", flush=True)
    return result


def plot_histories(
    output_path, title, species, histories, times,
    reference_mass_fractions=None, reference_temperatures=None,
):
    ncols = 4
    nrows = math.ceil((len(species) + 1) / ncols)
    figure, axes = pylab.subplots(
        nrows, ncols, figsize=(18.0, 3.2 * nrows), squeeze=False
    )
    line_styles = ["--", "-.", ":", "--"]
    markers = ["o", "s", "^", "D"]
    scaled_time = times * 1.0e6

    for species_index, (axis, species_name) in enumerate(
        zip(axes.flat, species)
    ):
        if reference_mass_fractions is not None:
            axis.plot(
                scaled_time, reference_mass_fractions[:, species_index],
                color="black", linewidth=2.0, label="Cantera",
            )
        for result_index, (label, mass_fractions, _) in enumerate(histories):
            axis.plot(
                scaled_time, mass_fractions[:, species_index],
                color=f"C{result_index}", linestyle=line_styles[result_index],
                marker=markers[result_index], markevery=2, markersize=3.5,
                markerfacecolor="none", linewidth=1.2, label=label,
            )
        axis.set_title(species_name)
        axis.set_xlabel("Time (microseconds)")
        axis.set_ylabel("Mass fraction")
        axis.grid(True, alpha=0.25)

    temperature_axis = axes.flat[len(species)]
    if reference_temperatures is not None:
        temperature_axis.plot(
            scaled_time, reference_temperatures,
            color="black", linewidth=2.0, label="Cantera",
        )
    for result_index, (label, _, temperatures) in enumerate(histories):
        temperature_axis.plot(
            scaled_time, temperatures,
            color=f"C{result_index}", linestyle=line_styles[result_index],
            marker=markers[result_index], markevery=2, markersize=3.5,
            markerfacecolor="none", linewidth=1.2, label=label,
        )
    temperature_axis.set_title("Temperature")
    temperature_axis.set_xlabel("Time (microseconds)")
    temperature_axis.set_ylabel("Temperature (K)")
    temperature_axis.grid(True, alpha=0.25)

    for axis in axes.flat[len(species) + 1:]:
        axis.axis("off")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    figure.legend(
        handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.975),
        ncols=len(labels), fontsize=9, frameon=False,
    )
    figure.suptitle(title, y=0.995)
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.935))
    figure.savefig(output_path, dpi=160)
    pylab.close(figure)


alamo.Util.Initialize()
try:
    alamo.Unit.setAmountUnit("kmol")
    alamo.Set.Constant.SetGlobalConstants()
    pp = alamo.IO.ParmParse()

    finite_gas = make_finite_gas(pp)
    gross_model_gas = make_gross_model_gas(pp)
    gross_model_aluminized_gas = make_gross_model_aluminized_gas(pp)
    source_model = make_finite_source_model(pp, "source_finite")
    substep_source_model = make_finite_source_model(
        pp, "source_finite_substep", substeps=10
    )

    frozen_integrator = make_integrator(
        pp, "frozen", "frozen", "forward_euler", 1, len(FINITE_SPECIES)
    )
    finite_cases = [
        ("Finite rate: Forward Euler, 1 substep", "finite_forward_1",
         "forward_euler", 1, 1.0e-8),
        ("Finite rate: Forward Euler, 10 substeps", "finite_forward_10",
         "forward_euler", 10, 1.0e-7),
        ("Finite rate: Backward Euler, 1 substep", "finite_backward_1",
         "backward_euler", 1, 1.0e-8),
    ]
    finite_integrators = [
        make_integrator(
            pp, prefix, "finite_rate", solver, nsubsteps,
            len(FINITE_SPECIES),
        )
        for _, prefix, solver, nsubsteps, _ in finite_cases
    ]
    gross_model_cases = [
        ("GrossModel: Forward Euler, 1 substep", "gross_model_forward_1",
         "forward_euler", 1),
        ("GrossModel: Forward Euler, 10 substeps", "gross_model_forward_10",
         "forward_euler", 10),
        ("GrossModel: Backward Euler, 1 substep", "gross_model_backward_1",
         "backward_euler", 1),
        ("GrossModel: Backward Euler, 10 substeps", "gross_model_backward_10",
         "backward_euler", 10),
    ]
    gross_model_integrators = [
        make_integrator(
            pp, prefix, "gross_model", solver, nsubsteps,
            len(GROSS_MODEL_SPECIES),
        )
        for _, prefix, solver, nsubsteps in gross_model_cases
    ]
    gross_model_aluminized_integrator = make_integrator(
        pp, "gross_model_aluminized_backward", "gross_model_aluminized",
        "backward_euler", 1, len(GROSS_MODEL_ALUMINIZED_SPECIES),
    )

    print("\nChemistry cases")

    def check_finite_sources():
        source = source_model.ComputeChemistrySources(
            101325.0, 1000.0,
            species_array(FINITE_INITIAL_PARTIAL_DENSITIES), 0.0, finite_gas,
        )
        wdot = numpy.array(
            [source.first[n] for n in range(len(FINITE_SPECIES))]
        )
        numpy.testing.assert_allclose(wdot, EXPECTED_WDOT, rtol=2.0e-8, atol=2.0e-12)
        numpy.testing.assert_allclose(source.second, EXPECTED_QDOT, rtol=2.0e-8, atol=2.0e-8)
        numpy.testing.assert_allclose(wdot.sum(), 0.0, rtol=0.0, atol=2.0e-12)

    run_case("Finite-rate chemistry sources", check_finite_sources)

    def check_substep_sources():
        source = substep_source_model.ComputeChemistrySources(
            101325.0, 1000.0,
            species_array(FINITE_INITIAL_PARTIAL_DENSITIES), 1.0e-4, finite_gas,
        )
        wdot = numpy.array(
            [source.first[n] for n in range(len(FINITE_SPECIES))]
        )
        numpy.testing.assert_allclose(
            wdot, EXPECTED_SUBSTEP_WDOT, rtol=2.0e-8, atol=2.0e-12
        )
        numpy.testing.assert_allclose(
            source.second, EXPECTED_SUBSTEP_QDOT, rtol=2.0e-8, atol=2.0e-8
        )
        numpy.testing.assert_allclose(wdot.sum(), 0.0, rtol=0.0, atol=2.0e-12)

    run_case("Finite-rate chemistry sources, 10 substeps", check_substep_sources)

    frozen_initial = FINITE_INITIAL_PARTIAL_DENSITIES / FINITE_INITIAL_PARTIAL_DENSITIES.sum()

    def check_frozen():
        mass_fractions, temperatures = integrate_composition(
            frozen_integrator, finite_gas, FINITE_INITIAL_PARTIAL_DENSITIES,
            2.0e-4, numpy.array([0.0, 2.0e-4]),
        )
        numpy.testing.assert_allclose(
            mass_fractions[-1], frozen_initial, rtol=0.0, atol=2.0e-15,
            err_msg="frozen chemistry changed the composition",
        )
        numpy.testing.assert_allclose(
            temperatures[-1], temperatures[0], rtol=0.0, atol=2.0e-12,
            err_msg="frozen chemistry changed the temperature",
        )

    run_case("Frozen chemistry", check_frozen)

    def check_gross_model_aluminized():
        initial = GROSS_MODEL_ALUMINIZED_INITIAL_PARTIAL_DENSITIES
        mass_fractions, temperatures = integrate_composition(
            gross_model_aluminized_integrator, gross_model_aluminized_gas, initial,
            1.0e-7, numpy.array([0.0, 1.0e-6]),
        )
        validate_mass_fractions("GrossModel aluminized", mass_fractions)
        initial_aluminum_fractions = initial[-2:] / initial.sum()
        numpy.testing.assert_allclose(
            mass_fractions[:, -2:],
            numpy.broadcast_to(initial_aluminum_fractions,
                               mass_fractions[:, -2:].shape),
            rtol=0.0, atol=2.0e-14,
            err_msg="GrossModel_Aluminized reacted an inert aluminum gas species",
        )
        if numpy.allclose(mass_fractions[-1, :6], mass_fractions[0, :6]):
            raise RuntimeError(
                "GrossModel_Aluminized did not advance its AP/HTPB chemistry"
            )
        if not numpy.all(numpy.isfinite(temperatures)):
            raise RuntimeError("GrossModel_Aluminized produced non-finite temperature")

    run_case("GrossModel aluminized inert aluminum", check_gross_model_aluminized)

    finite_histories = []
    for case, integrator in zip(finite_cases, finite_integrators):
        name, _, _, _, timestep = case

        def check_finite(integrator=integrator, timestep=timestep, name=name):
            mass_fractions, temperatures = integrate_composition(
                integrator, finite_gas, FINITE_INITIAL_PARTIAL_DENSITIES,
                timestep, CANTERA_TIMES,
            )
            validate_mass_fractions(name, mass_fractions)
            numpy.testing.assert_allclose(
                mass_fractions[-1], CANTERA_MASS_FRACTIONS[-1],
                rtol=5.0e-2, atol=2.0e-6,
                err_msg=f"{name} final composition differs from Cantera",
            )
            numpy.testing.assert_allclose(
                temperatures[-1], CANTERA_TEMPERATURES[-1],
                rtol=5.0e-2, atol=1.0,
                err_msg=f"{name} final temperature differs from Cantera",
            )
            return mass_fractions, temperatures

        history = run_case(name, check_finite)
        finite_histories.append((name.removeprefix("Finite rate: "), *history))

    gross_model_times = numpy.linspace(0.0, 1.0e-5, 21)
    gross_model_histories = []
    for case, integrator in zip(gross_model_cases, gross_model_integrators):
        name, _, _, nsubsteps = case

        def check_gross_model(
            integrator=integrator, name=name, nsubsteps=nsubsteps
        ):
            mass_fractions, temperatures = integrate_composition(
                integrator, gross_model_gas, GROSS_MODEL_INITIAL_PARTIAL_DENSITIES,
                1.0e-7, gross_model_times,
            )
            validate_mass_fractions(name, mass_fractions)
            numpy.testing.assert_allclose(
                mass_fractions[-1], GROSS_MODEL_FINAL_MASS_FRACTIONS,
                rtol=5.0e-2 if nsubsteps == 10 else 1.5e-1,
                atol=2.0e-3,
                err_msg=f"{name} final composition differs from the reference",
            )
            numpy.testing.assert_allclose(
                temperatures[-1], GROSS_MODEL_FINAL_TEMPERATURE,
                rtol=5.0e-2, atol=1.0,
                err_msg=f"{name} final temperature differs from the reference",
            )
            return mass_fractions, temperatures

        history = run_case(name, check_gross_model)
        gross_model_histories.append((name.removeprefix("GrossModel: "), *history))

    OUTPUT.mkdir(exist_ok=True)
    plot_histories(
        OUTPUT / "finiterate.png",
        "Finite-rate chemistry integration compared with Cantera",
        FINITE_SPECIES, finite_histories, CANTERA_TIMES,
        CANTERA_MASS_FRACTIONS, CANTERA_TEMPERATURES,
    )
    plot_histories(
        OUTPUT / "gross_model.png",
        "GrossModel stoichiometric AP/HTPB chemistry integration",
        GROSS_MODEL_SPECIES, gross_model_histories, gross_model_times,
    )
    print(f"\nPlots: {OUTPUT / 'finiterate.png'}")
    print(f"       {OUTPUT / 'gross_model.png'}")
finally:
    alamo.Util.Finalize()

print(f"Total chemistry test runtime: {time.perf_counter() - SUITE_START:.3f} s")
