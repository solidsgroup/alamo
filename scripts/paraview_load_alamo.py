from paraview.simple import *
import re
import sys

# ------------------------------------------------------------
# Read Chombo file list from celloutput.visit
# ------------------------------------------------------------

output = "celloutput.visit"

files = []

with open(output, "r") as f:
    for line in f:
        line = line.strip()
        if line and not line.startswith("!"):
            files.append(line)


# ------------------------------------------------------------
# Create Chombo reader
# ------------------------------------------------------------

reader = VisItChomboReader(FileName=files)

reader.UpdatePipelineInformation()


# ------------------------------------------------------------
# Enable all available cell variables
# ------------------------------------------------------------

all_cell_arrays = []

try:
    cell_array_info = reader.GetProperty("CellArrayInfo")

    if cell_array_info is not None:
        all_cell_arrays = [
            cell_array_info.GetElement(i)
            for i in range(cell_array_info.GetNumberOfElements())
        ]

except Exception as e:
    print("CellArrayInfo unavailable:")
    print(e)


# Fallback
if len(all_cell_arrays) == 0:

    print("Using CellArrayStatus fallback")

    status = reader.GetProperty("CellArrayStatus")

    all_cell_arrays = [
        status.GetElement(i)
        for i in range(status.GetNumberOfElements())
    ]


if len(all_cell_arrays) == 0:
    raise RuntimeError("No cell arrays found")

reader.CellArrayStatus = all_cell_arrays

reader.UpdatePipelineInformation()


# ------------------------------------------------------------
# Detect vector components
#
# Supports:
#   velocityx velocityy velocityz
#   solid.momentumx solid.momentumy solid.momentumz
#
# Also supports:
#   velocity.x velocity.y velocity.z
#   solid.momentum.x ...
# ------------------------------------------------------------

vectors = {}

for name in all_cell_arrays:

    # Case 1: suffix x/y/z
    m = re.match(r"^(.*)(x|y|z)$", name)

    if m:
        base, component = m.groups()
        vectors.setdefault(base, {})[component] = name
        continue

    # Case 2: suffix .x/.y/.z
    m = re.match(r"^(.*)\.(x|y|z)$", name)

    if m:
        base, component = m.groups()
        vectors.setdefault(base, {})[component] = name


vector_defs = {}

for base, comps in vectors.items():

    if "x" in comps and "y" in comps:
        vector_defs[base] = comps


# ------------------------------------------------------------
# Chain Calculator filters
# ------------------------------------------------------------

current_input = reader

for vec_name, comps in vector_defs.items():

    calc = Calculator(
        Input=current_input,
        AttributeType="Cell Data",
        ResultArrayName=vec_name
    )

    if "z" in comps:

        calc.Function = (
            f'iHat*"{comps['x']}" + '
            f'jHat*"{comps['y']}" + '
            f'kHat*"{comps['z']}"'
        )

    else:

        calc.Function = (
            f'iHat*"{comps['x']}" + '
            f'jHat*"{comps['y']}"'
        )

    current_input = calc


# ------------------------------------------------------------
# Enable animation
# ------------------------------------------------------------

animation = GetAnimationScene()
animation.UpdateAnimationUsingDataTimeSteps()

timekeeper = GetTimeKeeper()


# ------------------------------------------------------------
# Convert cell data to point data
# ------------------------------------------------------------

current_input = CellDatatoPointData(Input=current_input)
current_input.UpdatePipeline()

display = Show(current_input)
display.Representation = "Surface"

SetActiveSource(current_input)
Render()
