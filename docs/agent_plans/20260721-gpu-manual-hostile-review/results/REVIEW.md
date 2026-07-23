# Focused hostile review

Reviewer: independent reviewer agent, read-only.

The first pass found two high-severity evidence-scope defects and one medium
carryover defect:

- the permanent manual validator incorrectly depended on the task's full-tree
  audit CSV;
- that dependency could accept stale or unrelated class hits rather than exact
  reproduced aggregate sites;
- prior false-positive/not-applicable dispositions survived source revisions.

Corrections removed task coverage from the permanent validator, added the
task-local `verify_current_coverage.py` oracle with live-HEAD and exact-site
checks, retained regenerate-and-byte-compare in the task oracle, and restricted
carryover to identical port and source revision.

Recheck result: scoped validator, 10 scanner tests, fresh scan, exact-site
verification, and byte comparison all pass. The reviewer found no remaining
high- or medium-severity defect in the reviewed blocker items or corrections.
