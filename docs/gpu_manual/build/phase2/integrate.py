#!/usr/bin/env python3
import csv
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
BUILD = ROOT / "docs/gpu_manual/build"
PHASE2 = BUILD / "phase2"

PATTERNS = {
    "GPU-001": ("Annotate the full device call chain", "src/Numeric/Stencil.H", "correctness"),
    "GPU-002": ("Replace runtime polymorphism with static value dispatch", "src/Model/Gas/Gas.H", "correctness"),
    "GPU-003": ("Forward constructor arguments through static selection", "src/IO/ParmParse.H", "correctness"),
    "GPU-004": ("Devirtualize boundary conditions inside kernels", "src/BC/Operator/Elastic/Elastic.H", "correctness"),
    "GPU-005": ("Name and expose extended-lambda types", "src/Integrator/Flame.H", "correctness"),
    "GPU-006": ("Hoist kernel captures out of host objects", "src/IC/Expression.H", "correctness"),
    "GPU-007": ("Copy geometry into device-safe values", "src/Operator/Elastic.cpp", "correctness"),
    "GPU-008": ("Stage host buffers and fence local device storage", "src/Util/BMP.H", "correctness"),
    "GPU-009": ("Hold asynchronous FArrayBox temporaries with Elixir", "src/Operator/Operator.cpp", "correctness"),
    "GPU-010": ("Quarantine host-only loops in GPU builds", "src/IC/Laminate.H", "scaffolding"),
    "GPU-011": ("Select paired CPU and GPU loop implementations", "src/Solver/Nonlocal/Newton.H", "correctness"),
    "GPU-012": ("Report kernel failures through device-safe flags", "src/Integrator/Flame.cpp", "correctness"),
    "GPU-013": ("Move and fuse kernel reductions with ReduceOps", "src/Integrator/Base/Mechanics.H", "correctness"),
    "GPU-015": ("Fuse component launches into one multidimensional kernel", "src/Operator/Operator.cpp", "performance"),
    "GPU-016": ("Route spectral transforms through the GPU wrapper", "src/Integrator/CahnHilliard.cpp", "correctness"),
    "GPU-017": ("Materialize chained matrix expressions before device use", "src/Model/Solid/Finite/NeoHookean.H", "correctness"),
    "GPU-018": ("Hoist reused field and tensor loads", "src/Operator/Elastic.cpp", "performance"),
    "GPU-019": ("Defer branch-only expensive work", "src/Operator/Elastic.cpp", "performance"),
    "GPU-020": ("Specialize only the tensor contraction consumed", "src/Set/Matrix4_Major.H", "performance"),
    "GPU-021": ("Use dimension-aware compile-time constants", "src/Numeric/Stencil.H", "correctness"),
    "GPU-022": ("Replace host globals with device-safe sentinel access", "src/Set/Base.H", "correctness"),
    "GPU-023": ("Initialize device-backed MultiFabs on device", "src/Operator/Elastic.cpp", "correctness"),
    "GPU-024": ("Const-qualify captured model accessors", "src/Model/Propellant/Propellant.H", "correctness"),
    "GPU-025": ("Declare and enforce the supported GPU closure", "src/GPU/IntegratorPolicy.mk", "scaffolding"),
    "NUM-001": ("Propagate conservative face-flux mode end to end", "src/Operator/Elastic.cpp", "correctness"),
    "NUM-002": ("Preserve coarse-fine ghost-row policy", "src/Operator/Operator.cpp", "correctness"),
    "NUM-003": ("Carry component layouts through AMR transfers", "src/Operator/Elastic.cpp", "correctness"),
    "NUM-004": ("Dampen Newton updates with bounded line search", "src/Solver/Nonlocal/Newton.H", "correctness"),
    "GPU-030": ("Initialize device-local aggregates before use", "src/Util/PNG.H", "correctness"),
}

# Amendment 22: these device-independent solver/numerics transforms retire to
# FEATURES.md after the full Tier 1 support audit.
RETIRED_PATTERNS = {"NUM-001", "NUM-002", "NUM-003", "NUM-004"}
ACTIVE_PATTERNS = {
    pattern_id: metadata
    for pattern_id, metadata in PATTERNS.items()
    if pattern_id not in RETIRED_PATTERNS
}

KEY_MAP = {
    "gpu_const_qualify_model_methods": "GPU-024",
    "gpu_device_error_flag": "GPU-012",
    "gpu_device_parse_forwarding": "GPU-003",
    "gpu_device_safe_captures": "GPU-006",
    "gpu_device_safe_geometry": "GPU-007",
    "gpu_dimension_aware_constants": "GPU-021",
    "gpu_elixir_capture": "GPU-009",
    "gpu_face_gradient_stencil": "NUM-001",
    "gpu_fft_operator_migration": "GPU-016",
    "gpu_host_device_stencil_helpers": "GPU-001",
    "gpu_ic_static_dispatch_extension": "GPU-025",
    "gpu_integrator_safety_guards": "GPU-025",
    "gpu_integrator_source_closure": "GPU-025",
    "gpu_model_host_device_annotations": "GPU-001",
    "gpu_numeric_host_device_annotations": "GPU-001",
    "gpu_reduction_race_avoidance": "GPU-013",
    "gpu_static_model_dispatch": "GPU-002",
    "gpu_stream_lifetime_synchronization": "GPU-008",
    "BC_DEVICE_CAPTURE": "GPU-006",
    "BC_DEVICE_STATE": "GPU-008",
    "BC_HOST_DEVICE_ANNOTATION": "GPU-001",
    "BC_LOOP_DISPATCH": "GPU-011",
    "BC_STATIC_DISPATCH": "GPU-004",
    "BMP_DEVICE_STAGING": "GPU-008",
    "DEVICE_CAPTURE_HOIST": "GPU-006",
    "DEVICE_ERROR_API": "GPU-012",
    "DEVICE_ERROR_FLAG": "GPU-012",
    "DEVICE_SAFE_VECTOR": "GPU-007",
    "ELASTIC_API_CONTROLS": "NUM-001",
    "ELASTIC_COMPONENT_INTERPOLATION": "NUM-003",
    "ELASTIC_COMPONENT_LAYOUT": "NUM-003",
    "ELASTIC_CONSERVATIVE_FLUX": "NUM-001",
    "ELASTIC_KERNEL_DISPATCH": "GPU-011",
    "ELIXIR_LIFETIME": "GPU-009",
    "GPU_DEVICE_LAUNCH": "GPU-023",
    "HOST_DEVICE_DIAGNOSTIC_API": "GPU-012",
    "HOST_LOOP_FALLBACK": "GPU-010",
    "HOST_LOOP_GUARD": "GPU-010",
    "IC_DEVICE_CAPTURE": "GPU-006",
    "MATRIX4_UNROLL": "GPU-020",
    "NEWTON_CONSERVATIVE_SOLVE": "NUM-001",
    "NEWTON_DAMPED_STEP": "NUM-004",
    "NEWTON_LOOP_DISPATCH": "GPU-011",
    "NEWTON_STATIC_BC": "GPU-004",
    "OPERATOR_GHOST_POLICY": "NUM-002",
    "OPERATOR_INTERPOLATION": "NUM-003",
    "PNG_DEVICE_API": "GPU-001",
    "SET_DEVICE_GARBAGE": "GPU-022",
    "SET_FIELD_DEVICE_COPY": "NUM-003",
    "SET_HOST_DEVICE_ANNOTATION": "GPU-001",
}

FORCE_ONEOFF_KEYS = {
    "gpu_finest_level_spectral_guard",
    "gpu_residual_api_update",
    "gpu_restart_field_copy_api",
    "DIAGNOSTIC_INSTRUMENTATION",
    "OPERATOR_GHOST_SANITIZE",
}

OVERRIDES = {
    "src/Model/Gas/Gas.H:112-118": ["ONEOFF"],
    "src/Model/Propellant/Homogenize.H:8-14": ["ONEOFF"],
    "src/Model/Propellant/PowerLaw.H:31-36": ["ONEOFF"],
    "src/Model/Propellant/Propellant.H:35-40": ["ONEOFF"],
    "src/Model/Solid/Finite/NeoHookean.H:14-34": ["GPU-001", "ONEOFF"],
    "src/Numeric/Function.H:1-6": ["ONEOFF"],
    "src/Numeric/Stencil.H:83-107": ["GPU-001", "ONEOFF"],
    "src/Numeric/Stencil.H:785-821": ["NUM-001", "GPU-001", "ONEOFF"],
    "src/Numeric/Stencil.H:1434-1445": ["GPU-021", "GPU-001"],
    "src/Set/Base.H:16-21": ["ONEOFF"],
    "src/Set/Base.H:188-193": ["ONEOFF"],
    "src/Integrator/Base/Mechanics.H:145-151": ["GPU-005"],
    "src/Integrator/Base/Mechanics.H:198-218": ["GPU-008", "ONEOFF"],
    "src/Integrator/Flame.H:45-51": ["GPU-005"],
    "src/Integrator/Flame.H:86-106": ["GPU-005", "ONEOFF"],
    "src/Integrator/Flame.H:113-144": ["GPU-005", "ONEOFF"],
    "src/Integrator/Flame.cpp:407-418": ["GPU-013"],
    "src/Integrator/Flame.cpp:484-534": ["GPU-012", "ONEOFF"],
    "src/Integrator/Flame.cpp:539-554": ["GPU-012"],
    "src/Integrator/Flame.cpp:582-592": ["GPU-012"],
    "src/Integrator/Flame.cpp:595-693": ["GPU-007", "GPU-006", "ONEOFF"],
    "src/Integrator/Integrator.H:256-262": ["GPU-005"],
    "src/Integrator/Integrator.H:267-273": ["GPU-005"],
    "src/Model/Solid/Finite/NeoHookean.H:36-90": ["GPU-001", "GPU-017", "ONEOFF"],
    "src/BC/Constant.cpp:1-4": ["GPU-008"],
    "src/BC/Operator/Elastic/Expression.H:6-11": ["ONEOFF"],
    "src/IC/Expression.H:34-39": ["ONEOFF"],
    "src/IC/Laminate.H:26-31": ["GPU-010"],
    "src/IC/Laminate.H:5-10": ["ONEOFF"],
    "src/IC/PNG.H:48-53": ["GPU-010"],
    "src/IC/PNG.H:4-9": ["ONEOFF"],
    "src/IC/PSRead.H:6-11": ["ONEOFF"],
    "src/IC/PSRead.H:128-134": ["NOISE"],
    "src/IC/Trig.H:9-14": ["ONEOFF"],
    "src/Operator/Elastic.H:108-115": ["ONEOFF"],
    "src/Operator/Elastic.H:129-167": ["GPU-005", "GPU-006"],
    "src/Operator/Elastic.cpp:173-264": [
        "GPU-004", "GPU-006", "GPU-007", "GPU-012", "GPU-018",
        "GPU-019", "GPU-020", "NUM-001",
    ],
    "src/Operator/Elastic.cpp:148-154": ["GPU-007", "GPU-012"],
    "src/Operator/Elastic.cpp:317-347": ["GPU-006", "GPU-018"],
    "src/Operator/Elastic.cpp:299-315": ["GPU-007", "GPU-012", "NUM-002"],
    "src/Operator/Elastic.cpp:352-377": ["GPU-004", "GPU-012", "GPU-018", "NUM-001"],
    "src/Operator/Elastic.cpp:722-728": ["GPU-010"],
    "src/Operator/Operator.H:57-101": ["GPU-005", "ONEOFF"],
    "src/Operator/Operator.cpp:116-150": ["GPU-015", "NUM-002"],
    "src/Operator/Operator.cpp:368-380": ["ONEOFF", "ONEOFF"],
    "src/Integrator/PFC.cpp:10-15": ["ONEOFF"],
    "src/Set/Base.H:7-12": ["ONEOFF"],
    "src/Set/Matrix4_Major.H:54-93": ["GPU-022", "ONEOFF"],
    "src/Set/Matrix4_Major.H:242-247": ["GPU-001"],
    "src/Set/Matrix4_Major.H:472-483": ["GPU-001"],
    "src/Set/Matrix4_MajorMinor.H:48-54": ["GPU-022"],
    "src/Set/Matrix4_MajorMinor.H:233-246": ["GPU-001"],
    "src/Set/Matrix4_MajorMinor.H:624-637": ["GPU-001"],
    "src/Util/PNG.H:155-161": ["GPU-001", "GPU-024"],
    "src/Util/PNG.H:168-174": ["GPU-030"],
    "src/Util/PNG.H:199-206": ["NOISE"],
    "src/Solver/Nonlocal/Newton.H:60-67": ["GPU-007", "GPU-012"],
    "src/Solver/Nonlocal/Newton.H:177-188": ["GPU-012", "NUM-001", "NUM-003"],
    "src/Util/Util.H:15-20": ["ONEOFF"],
    "src/Util/Util.H:51-62": ["GPU-012", "GPU-001"],
    "src/Util/Util.H:77-88": ["GPU-001"],
    "src/Util/Util.H:101-112": ["GPU-012", "GPU-001"],
    "src/Util/Util.H:123-133": ["GPU-001"],
    "src/Util/Util.H:146-151": ["GPU-001"],
    "src/Util/Util.H:182-189": ["GPU-001"],
    "src/Util/Util.H:196-206": ["GPU-001"],
    "src/Util/Util.H:219-230": ["GPU-001"],
    "src/Util/Util.H:247-260": ["GPU-001"],
    "src/Util/Util.H:276-281": ["GPU-001"],
    "src/Util/Util.cpp:28-33": ["ONEOFF"],
    "src/Util/Util.cpp:138-143": ["ONEOFF"],
    "src/Util/Util.cpp:265-271": ["GPU-001"],
    "src/Operator/Elastic.H:122-127": ["NUM-002"],
}


def parent_id(hunk_id):
    return re.sub(r"\s*#\d+$", "", hunk_id)


def load_rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


index_rows = load_rows(PHASE2 / "HUNK_INDEX.csv")
index_ids = [row["hunk_id"] for row in index_rows]
index_set = set(index_ids)
if len(index_ids) != 432 or len(index_set) != 432:
    raise SystemExit("HUNK_INDEX must contain 432 unique contextual parents")

proposals = defaultdict(list)
for part in ("A", "B"):
    for row in load_rows(PHASE2 / f"{part}_HUNKS.csv"):
        parent = parent_id(row["hunk_id"])
        if parent not in index_set:
            raise SystemExit(f"worker row not in HUNK_INDEX: {row['hunk_id']}")
        classification = row["classification"]
        key = row["candidate_key"]
        if classification == "NOISE":
            value = "NOISE"
        elif classification == "ONEOFF" or key in FORCE_ONEOFF_KEYS:
            value = "ONEOFF"
        else:
            if key not in KEY_MAP:
                raise SystemExit(f"unmapped candidate key: {key} ({row['hunk_id']})")
            value = KEY_MAP[key]
        if value not in proposals[parent]:
            proposals[parent].append(value)

missing = index_set - proposals.keys()
if missing:
    raise SystemExit(f"parents missing worker classifications: {sorted(missing)}")

resolved = {}
for parent in index_ids:
    values = list(OVERRIDES.get(parent, proposals[parent]))
    if not values:
        raise SystemExit(f"empty resolution: {parent}")
    for value in values:
        if value not in PATTERNS and value not in {"ONEOFF", "NOISE"}:
            raise SystemExit(f"bad resolution {value}: {parent}")
    resolved[parent] = values

# Amendments 18-20: independently review every provisional ONEOFF subentry
# against the BASE-port decision test. FEATURE and NUM tagging are orthogonal.
review_rows = []
for part in ("A", "B"):
    review_rows.extend(load_rows(PHASE2 / "feature_review" / f"{part}.csv"))
review_by_id = {}
for row in review_rows:
    hunk_id = row["hunk_id"]
    if hunk_id in review_by_id:
        raise SystemExit(f"duplicate FEATURE-review row: {hunk_id}")
    if row["classification"] not in {"ONEOFF", "FEATURE", "NOISE"}:
        raise SystemExit(f"bad FEATURE-review classification: {hunk_id}")
    if row["num_tag"] not in {"yes", "no"}:
        raise SystemExit(f"bad FEATURE-review NUM tag: {hunk_id}")
    if row["ambiguous"] not in {"yes", "no"}:
        raise SystemExit(f"bad FEATURE-review ambiguity flag: {hunk_id}")
    if not row["summary"] or not row["decision_reason"]:
        raise SystemExit(f"incomplete FEATURE-review rationale: {hunk_id}")
    if row["ambiguous"] == "yes":
        if row["classification"] != "ONEOFF" or "conservative default after two attempts" not in row["decision_reason"]:
            raise SystemExit(f"ambiguous review must use documented conservative ONEOFF default: {hunk_id}")
    review_by_id[hunk_id] = row

provisional_oneoffs = set()
for parent, values in resolved.items():
    for number, value in enumerate(values, 1):
        hunk_id = parent if len(values) == 1 else f"{parent}#{number}"
        if value == "ONEOFF":
            provisional_oneoffs.add(hunk_id)
if set(review_by_id) != provisional_oneoffs:
    missing = sorted(provisional_oneoffs - set(review_by_id))
    extra = sorted(set(review_by_id) - provisional_oneoffs)
    raise SystemExit(f"FEATURE review mismatch: missing={missing}, extra={extra}")

index_files = {row["hunk_id"]: row["file"] for row in index_rows}
entry_num_tags = {}
for parent, values in resolved.items():
    for number, value in enumerate(list(values), 1):
        hunk_id = parent if len(values) == 1 else f"{parent}#{number}"
        if value != "ONEOFF":
            entry_num_tags[(parent, number)] = "no"
            continue
        review = review_by_id[hunk_id]
        if review["file"] != index_files[parent]:
            raise SystemExit(f"FEATURE-review file mismatch: {hunk_id}")
        values[number - 1] = review["classification"]
        entry_num_tags[(parent, number)] = review["num_tag"]

retired_entries = []
for parent, values in resolved.items():
    for number, value in enumerate(list(values), 1):
        if value in RETIRED_PATTERNS:
            hunk_id = parent if len(values) == 1 else f"{parent}#{number}"
            retired_entries.append((hunk_id, index_files[parent], value))
            values[number - 1] = "FEATURE"
            entry_num_tags[(parent, number)] = "yes"

with (BUILD / "HUNK_MAP.csv").open("w", newline="") as handle:
    writer = csv.writer(handle, lineterminator="\n")
    writer.writerow(["hunk_id", "file", "classification", "num_tag"])
    for row in index_rows:
        parent = row["hunk_id"]
        values = resolved[parent]
        for number, value in enumerate(values, 1):
            hunk_id = parent if len(values) == 1 else f"{parent}#{number}"
            writer.writerow([hunk_id, row["file"], value, entry_num_tags[(parent, number)]])

with (PHASE2 / "RETIRED_PATTERNS.csv").open("w", newline="") as handle:
    writer = csv.writer(handle, lineterminator="\n")
    writer.writerow(["hunk_id", "file", "retired_pattern"])
    writer.writerows(retired_entries)

pattern_parents = defaultdict(set)
parent_classes = {}
for parent, values in resolved.items():
    for value in values:
        if value in ACTIVE_PATTERNS:
            pattern_parents[value].add(parent)
    if "ONEOFF" in values:
        parent_classes[parent] = "ONEOFF"
    elif "FEATURE" in values:
        parent_classes[parent] = "FEATURE"
    elif any(value in ACTIVE_PATTERNS for value in values):
        parent_classes[parent] = "CANDIDATE"
    else:
        parent_classes[parent] = "NOISE"

missing_patterns = set(ACTIVE_PATTERNS) - pattern_parents.keys()
if missing_patterns:
    raise SystemExit(f"patterns with no hunk support: {sorted(missing_patterns)}")

with (PHASE2 / "CLUSTER_CHECKPOINT.csv").open("w", newline="") as handle:
    writer = csv.writer(handle, lineterminator="\n")
    writer.writerow(["pattern_id", "candidate_name", "class", "hunk_count", "example_file"])
    for pattern_id, (name, example, pattern_class) in ACTIVE_PATTERNS.items():
        writer.writerow([pattern_id, name, pattern_class, len(pattern_parents[pattern_id]), example])

counts = defaultdict(int)
for classification in parent_classes.values():
    counts[classification] += 1
oneoff_denominator = len(index_ids) - counts["NOISE"] - counts["FEATURE"]
oneoff_pct = 100.0 * counts["ONEOFF"] / oneoff_denominator
with (PHASE2 / "COUNTS.txt").open("w") as handle:
    handle.write(f"contextual_parents={len(index_ids)}\n")
    handle.write(f"ledger_rows={sum(len(v) for v in resolved.values())}\n")
    handle.write(f"candidate_parents={counts['CANDIDATE']}\n")
    handle.write(f"oneoff_parents={counts['ONEOFF']}\n")
    handle.write(f"feature_parents={counts['FEATURE']}\n")
    handle.write(f"noise_parents={counts['NOISE']}\n")
    handle.write(f"oneoff_denominator={oneoff_denominator}\n")
    handle.write(f"oneoff_percent={oneoff_pct:.3f}\n")
    handle.write(f"clusters={len(ACTIVE_PATTERNS)}\n")
    handle.write(f"retired_patterns={len(RETIRED_PATTERNS)}\n")
    handle.write("single_evidence_patterns=GPU-015,GPU-017,GPU-019,GPU-030\n")

print((PHASE2 / "COUNTS.txt").read_text(), end="")
if counts["ONEOFF"] / oneoff_denominator > 0.20:
    raise SystemExit("ONEOFF stop condition fired")
