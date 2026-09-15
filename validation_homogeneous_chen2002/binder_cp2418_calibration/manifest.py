#!/usr/bin/env python3
"""Record frozen inputs before launch; does not execute simulations."""
import json
from audit import HERE, STUDY, digest

parameters_path = HERE / "parameters_frozen.json"
parameters = json.loads(parameters_path.read_text())
freeze = json.loads((HERE / "freeze_record.json").read_text())
assert freeze["parameters_sha256"] == digest(parameters_path)
cases = []
for folder in sorted((STUDY / "runs").glob("*_cp2418_frozen*")):
    case = json.loads((folder / "case.json").read_text())
    assert case["parameters"] == parameters
    assert case["input_sha256"] == digest(folder / "input")
    cases.append(dict(name=folder.name,
        q_cal_cm2_s=case["reference"]["heat_flux_cal_cm2_s"],
        ap_volume_fraction=case["reference"]["ap_volume_fraction"],
        relaxations=case["relaxations"], input_sha256=case["input_sha256"]))
assert len(cases) == 17
expected = {(q, t, 6) for q in (200., 500., 1000.) for t in (0., .2, .4, .6, .8, 1.)
            if not (t == 0 and q in (200., 1000.))} | {(500., .4, 12)}
assert {(c["q_cal_cm2_s"], c["ap_volume_fraction"], c["relaxations"]) for c in cases} == expected
out = HERE / "launch_manifest.json"
assert not out.exists(), "Preserve original launch manifest"
out.write_text(json.dumps(dict(binary_sha256=digest(STUDY.parent / "bin/lowmach-2d-clang++"),
    parameters_sha256=digest(parameters_path), cases=cases), indent=2)+"\n")
print(out, len(cases), "cases")
