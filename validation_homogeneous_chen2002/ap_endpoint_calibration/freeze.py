#!/usr/bin/env python3
"""Freeze an independently verified AP endpoint fit before running the holdout."""
import argparse
import json
from pathlib import Path

from audit import HERE, audit, digest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--parameters", type=Path, required=True)
    p.add_argument("--summary", type=Path, required=True)
    args = p.parse_args()
    parameters = json.loads(args.parameters.read_text())
    rows = json.loads(args.summary.read_text())
    assert len(rows) == 2
    assert {(r["q_cal_cm2_s"], r["ap_volume_fraction"]) for r in rows} == {(200., 1.), (1000., 1.)}
    evidence = audit(rows, parameters)
    assert max(abs(r["error_figure4_solid_percent"]) for r in rows) < 1.
    assert max(abs(r["late_speed_drift_percent"]) for r in rows) < .2
    destination = HERE / "parameters_frozen.json"
    assert not destination.exists(), "Preserve frozen parameters"
    destination.write_bytes(args.parameters.read_bytes())
    (HERE / "freeze_record.json").write_text(json.dumps(dict(
        source_parameters=str(args.parameters), source_summary=str(args.summary),
        parameters_sha256=digest(destination), summary_sha256=digest(args.summary),
        criteria=dict(max_endpoint_error_percent=1., max_late_drift_percent=.2),
        observed_max_endpoint_error_percent=max(abs(r["error_figure4_solid_percent"]) for r in rows),
        heldout_flux_cal_cm2_s=500., execution_audit=evidence), indent=2)+"\n")
    for r in rows:
        print(r["name"], r["r_cm_s"], "cm/s; error", r["error_figure4_solid_percent"], "%")
    print("Frozen:", destination)


if __name__ == "__main__":
    main()
