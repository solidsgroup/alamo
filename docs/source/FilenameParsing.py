#!/usr/bin/python
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts" / "builder"))

from reference import leading_documentation


source = REPO_ROOT / "src/IO/FileNameParse.H"
documentation = leading_documentation(source)
if not documentation:
    raise RuntimeError(f"No leading documentation found in {source}")

Path(__file__).with_suffix(".rst").write_text(
    documentation + "\n",
    encoding="utf-8",
)
