import csv
import tempfile
import unittest
from pathlib import Path

from scan import FIELDS, TABLE_FIELDS, previous_rows, read_table, scan


def table(path, candidate="BAD\\s+value", converted="GOOD\\s+value", exclude=""):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=TABLE_FIELDS)
        w.writeheader()
        w.writerow({"schema_version": "3", "pattern_id": "GPU-001", "type": "regex",
                     "candidate_expression": candidate, "converted_expression": converted,
                     "exclude_expression": exclude, "confirmation": "inspect", "notes": "fixture"})


class ScannerTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        (self.root / "src").mkdir()
        self.table_path = self.root / "table.csv"

    def tearDown(self):
        self.tmp.cleanup()

    def test_candidate_converted_and_mixed_file(self):
        table(self.table_path)
        (self.root / "src/a.cpp").write_text("BAD value;\nGOOD value;\n")
        rows = scan(self.root / "src", read_table(self.table_path), {}, "fixture", "rev1")
        self.assertEqual({r["state"] for r in rows}, {"candidate", "converted"})
        self.assertEqual([r["line"] for r in rows], ["1", "2"])

    def test_same_site_candidate_wins_over_converted(self):
        table(self.table_path, candidate="VALUE", converted="VALUE")
        (self.root / "src/a.cpp").write_text("VALUE;\n")
        rows = scan(self.root / "src", read_table(self.table_path), {}, "fixture", "rev1")
        self.assertEqual([(r["line"], r["state"]) for r in rows], [("1", "candidate")])

    def test_same_revision_line_movement_preserves_reviewed_disposition(self):
        table(self.table_path)
        source = self.root / "src/a.cpp"
        source.write_text("BAD value;\n")
        first = scan(self.root / "src", read_table(self.table_path), {}, "fixture", "rev1")
        reviewed = dict(first[0], state="false-positive", evidence="reviewed")
        source.write_text("// moved\nBAD value;\n")
        rows = scan(self.root / "src", read_table(self.table_path), {reviewed["site_id"]: reviewed}, "fixture", "rev1")
        self.assertEqual(rows[0]["state"], "false-positive")
        self.assertEqual(rows[0]["line"], "2")

    def test_same_match_in_two_files_has_distinct_identity(self):
        table(self.table_path)
        (self.root / "src/a.cpp").write_text("BAD value;\n")
        (self.root / "src/b.cpp").write_text("BAD value;\n")
        rows = scan(self.root / "src", read_table(self.table_path), {}, "fixture", "rev1")
        self.assertEqual(len({r["site_id"] for r in rows if r["state"] == "candidate"}), 2)

    def test_default_root_skips_generated_components(self):
        table(self.table_path)
        (self.root / "src/a.cpp").write_text("BAD value;\n")
        (self.root / "docs").mkdir()
        (self.root / "docs/d.cpp").write_text("BAD value;\n")
        (self.root / "ext").mkdir()
        (self.root / "ext/e.cpp").write_text("BAD value;\n")
        rows = scan(self.root, read_table(self.table_path), {}, "fixture", "rev1")
        self.assertEqual({r["file"] for r in rows if r["state"] == "candidate"}, {"src/a.cpp"})

    def test_exclusion_emits_not_applicable(self):
        table(self.table_path, exclude="ALLOW")
        (self.root / "src/a.cpp").write_text("BAD value; // ALLOW\n")
        rows = scan(self.root / "src", read_table(self.table_path), {}, "fixture", "rev1")
        self.assertEqual(rows[0]["state"], "not-applicable")

    def test_invalid_schema_rejected(self):
        self.table_path.write_text("pattern_id,type\nGPU-001,regex\n")
        with self.assertRaises(ValueError):
            read_table(self.table_path)

    def test_metadata_and_port_scoped_review_carryover(self):
        table(self.table_path)
        source = self.root / "src/a.cpp"
        source.write_text("BAD value;\n")
        first = scan(self.root / "src", read_table(self.table_path), {}, "hydro", "rev-a")
        self.assertEqual(first[0]["port_id"], "hydro")
        self.assertEqual(first[0]["source_revision"], "rev-a")
        reviewed = self.root / "previous.csv"
        with reviewed.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=FIELDS)
            w.writeheader()
            w.writerow(dict(first[0], state="false-positive", evidence="reviewed"))
        self.assertEqual(previous_rows(reviewed, "hydro", "rev-a")[first[0]["site_id"]]["state"], "false-positive")
        self.assertEqual(previous_rows(reviewed, "fracture", "rev-a"), {})
        self.assertEqual(previous_rows(reviewed, "hydro", "rev-b"), {})
        carried = scan(self.root / "src", read_table(self.table_path),
                       previous_rows(reviewed, "hydro", "rev-b"), "hydro", "rev-b")
        self.assertEqual(carried[0]["state"], "candidate")
        self.assertEqual(carried[0]["source_revision"], "rev-b")

    def test_production_shapes_are_identifier_independent(self):
        with (Path(__file__).parent / "table.csv").open(newline="") as source_table:
            rows = [row for row in csv.DictReader(source_table)
                    if row["pattern_id"] in {"GPU-002", "GPU-004", "GPU-013"}]
        with self.table_path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=TABLE_FIELDS); w.writeheader(); w.writerows(rows)
        (self.root / "src/shapes.cpp").write_text("""std::vector<double> coeffs;\namrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i,int j,int k) {\n  coeffs[i].value();\n  *dt_handle = std::min(*dt_handle, cfl);\n  driving_force_norm += value;\n});\n""")
        hits = {(r["pattern_id"], r["state"]) for r in scan(self.root / "src", read_table(self.table_path), {}, "hydro", "rev")}
        self.assertIn(("GPU-002", "candidate"), hits)
        self.assertIn(("GPU-004", "candidate"), hits)
        self.assertIn(("GPU-013", "candidate"), hits)

    def test_generic_fixture_rows_for_selection_load_branch_sentinel_and_init(self):
        with (Path(__file__).parent / "table.csv").open(newline="") as source_table:
            table_rows = list(csv.DictReader(source_table))
        wanted = {"GPU-003", "GPU-018", "GPU-019", "GPU-022", "GPU-030"}
        with self.table_path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=TABLE_FIELDS); w.writeheader()
            w.writerows(r for r in table_rows if r["pattern_id"] in wanted)
        (self.root / "src/generic.cpp").write_text("""select_default<Model>(value);
amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i,int j,int k) {
  patch.load(i,j); patch.load(i,j);
  Set::Vector tmp = make_vector();
  if (enabled) use(tmp);
  static int sentinel;
  return Set::Scalar();
  Set::Vector local; local(0) = 1.;
});
""")
        ids = {r["pattern_id"] for r in scan(self.root / "src", read_table(self.table_path), {}, "fixture", "rev")}
        self.assertTrue(wanted.issubset(ids))


if __name__ == "__main__":
    unittest.main()
