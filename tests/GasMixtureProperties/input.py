import csv
from pathlib import Path

import cppyy

import alamo
import numpy


alamo.include("Util/Util.H")
alamo.include("IO/ParmParse.H")
alamo.include("Model/Gas/Gas.H")

cppyy.cppdef(
    r"""
    namespace AlamoGasMixturePropertiesTest
    {
    struct Composition
    {
        explicit Composition(const std::vector<Set::Scalar>& a_values)
            : values(a_values)
        {}

        Set::Scalar operator()(int, int, int, int species) const
        {
            return values[species];
        }

        std::vector<Set::Scalar> values;
    };
    }
    """
)


EXPECTED = {
    "viscous-thermal-constant": {
        "molecular_weights": ["44.01_g/mol", "32.00_g/mol", "28.02_g/mol"],
        "partial_densities": [
            0.243454875129296,
            0.0519074926855928,
            0.964970271235302,
        ],
        "temperature": 292.99955533880296,
        "transport": {
            "type": "constant",
            "mu": ["1462E-7_g/cm/s", "2031E-7_g/cm/s", "1754E-7_g/cm/s"],
            "k": ["383e-7_cal/cm/s/K", "612e-7_cal/cm/s/K", "627e-7_cal/cm/s/K"],
        },
        "properties": [
            1.7142991458692116e-05,
            0.024446070520350752,
            1.3957619455684473e-05,
        ],
    },
    "viscous-thermal-LJ": {
        "molecular_weights": ["44.01_g/mol", "32.00_g/mol", "28.02_g/mol"],
        "partial_densities": [
            0.243454875129296,
            0.0519074926855928,
            0.964970271235302,
        ],
        "temperature": 292.99955533880296,
        "transport": {
            "type": "LJ",
            "LJdiameter": ["3.996_ang", "3.433_ang", "3.667_ang"],
            "LJwelldepth": ["190.0_K", "113.0_K", "99.8_K"],
        },
        "properties": [
            1.6884124446477103e-05,
            0.0240145886371462,
            1.4573514272511353e-05,
        ],
    },
    "species-constant": {
        "molecular_weights": ["28.01_g/mol", "44.01_g/mol", "28.02_g/mol"],
        "partial_densities": [0.574470635009588, 0.902622372251766, 0.0],
        "temperature": 297.0966170009438,
        "transport": {
            "type": "constant",
            "mu": ["1740E-7_g/cm/s", "1462E-7_g/cm/s", "1754E-7_g/cm/s"],
            "k": ["554e-7_cal/cm/s/K", "383e-7_cal/cm/s/K", "627e-7_cal/cm/s/K"],
            "d_fact": 1.4,
        },
        "properties": [
            1.5815104346700175e-05,
            0.01907231139644361,
            1.5060017708001669e-05,
        ],
    },
    "species-LJ": {
        "molecular_weights": ["28.01_g/mol", "44.01_g/mol", "28.02_g/mol"],
        "partial_densities": [0.574470635009588, 0.902622372251766, 0.0],
        "temperature": 297.0966170009438,
        "transport": {
            "type": "LJ",
            "LJdiameter": ["3.590_ang", "3.996_ang", "3.667_ang"],
            "LJwelldepth": ["100.0_K", "190.0_K", "99.8_K"],
        },
        "properties": [
            1.620866188510972e-05,
            0.019795463972170718,
            1.5249863744230164e-05,
        ],
    },
}

PROPERTY_NAMES = [
    "viscosity [Pa s]",
    "conductivity [W/(m K)]",
    "diffusion species 1 [m^2/s]",
]


def cpp_vector(cpp_type, values):
    result = alamo.std.vector[cpp_type]()
    for value in values:
        result.push_back(value)
    return result


def add_string(pp, name, value):
    pp.add(name, alamo.std.string(value))


def add_strings(pp, name, values):
    pp.addarr(name, cpp_vector("std::string", values))


def make_gas(pp, prefix, case):
    add_strings(pp, f"{prefix}.mw", case["molecular_weights"])
    add_string(pp, f"{prefix}.thermo.type", "cpconstant")
    add_strings(pp, f"{prefix}.thermo.cpconstant.cp_moles", ["32.0_J/mol/K"] * 3)
    add_strings(pp, f"{prefix}.thermo.cpconstant.h0", ["0.0_J/mol"] * 3)
    add_strings(pp, f"{prefix}.thermo.cpconstant.s0", ["0.0_J/mol/K"] * 3)
    add_strings(pp, f"{prefix}.thermo.cpconstant.Tref", ["0.0_K"] * 3)
    add_string(pp, f"{prefix}.transport.type", "mixture_averaged")

    transport_prefix = f"{prefix}.transport.mixture_averaged"
    transport = case["transport"]
    add_string(pp, f"{transport_prefix}.type", transport["type"])
    for name in ("mu", "k", "LJdiameter", "LJwelldepth"):
        if name in transport:
            add_strings(pp, f"{transport_prefix}.{name}", transport[name])
    if "d_fact" in transport:
        pp.add(f"{transport_prefix}.d_fact", transport["d_fact"])

    add_string(pp, f"{prefix}.eos.type", "cpg")
    return alamo.Model.Gas.Gas(pp, prefix)


alamo.Util.Initialize()
try:
    # The reference data use kg/kmol and J/kmol, matching the original tests.
    alamo.Unit.setAmountUnit("kmol")
    alamo.Set.Constant.SetGlobalConstants()

    pp = alamo.IO.ParmParse()
    Composition = alamo.AlamoGasMixturePropertiesTest.Composition
    results = []

    for case_number, (name, case) in enumerate(EXPECTED.items()):
        gas = make_gas(pp, f"gas_mixture_case{case_number}", case)

        molecular_weights = numpy.array(
            [float(value.split("_")[0]) for value in case["molecular_weights"]]
        )
        mole_amounts = numpy.array(case["partial_densities"]) / molecular_weights
        mole_fractions = mole_amounts / mole_amounts.sum()
        composition = Composition(cpp_vector("double", mole_fractions))

        temperature = case["temperature"]
        pressure = 101325.0
        actual = numpy.array(
            [
                gas.dynamic_viscosity(temperature, composition, 0, 0, 0),
                gas.thermal_conductivity(temperature, composition, 0, 0, 0),
                gas.diffusion_coefficient(
                    temperature, pressure, composition, 0, 0, 0, 0
                ),
            ]
        )
        numpy.testing.assert_allclose(
            actual,
            case["properties"],
            rtol=2.0e-12,
            atol=2.0e-15,
            err_msg=f"{name} gas mixture properties differ from the reference",
        )

        for property_name, reference, value in zip(
            PROPERTY_NAMES, case["properties"], actual
        ):
            results.append((name, property_name, reference, value))

    output = Path(__file__).with_name("output")
    output.mkdir(exist_ok=True)
    rows = [
        {
            "case": name,
            "property": property_name,
            "reference": reference,
            "alamo": actual,
            "relative_error": (
                abs(actual - reference) / abs(reference)
            ),
        }
        for name, property_name, reference, actual in results
    ]
    csv_path = output / "gas_mixture_properties.csv"
    with csv_path.open("w", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    print(f"CSV: {csv_path}")
    print("PASS: Gas mixture properties match the reference data")
finally:
    alamo.Util.Finalize()
