#!/usr/bin/env python3
import os
import tarfile
import re
from pathlib import Path
import shutil

REPORT_HTML = Path("report.html")
REPORT_DIR = REPORT_HTML.parent
TARBALL = Path("report.tar.gz")

# Regex to find image and link paths
IMG_RE = re.compile(r'<img[^>]+src=["\']([^"\']+)["\']')
LINK_RE = re.compile(r'<a[^>]+href=["\']([^"\']+)["\']')

def extract_paths_from_html(html_file):
    with open(html_file, "r", encoding="utf-8") as f:
        html = f.read()
    paths = set()
    paths.update([(html_file.parent/Path(p)).absolute().resolve() for p in IMG_RE.findall(html)])
    paths.update([(html_file.parent/Path(p)).absolute().resolve() for p in LINK_RE.findall(html)])


    subpaths = set()
    for path in paths:
        if str(path).endswith("html"):
            subpaths.update(extract_paths_from_html(Path(path).absolute()))
    # Always include the HTML file itself
    paths.update(subpaths)
    paths.add(html_file)

    return [p for p in paths if p.exists()]

def create_tarball(files, output_path):
    with tarfile.open(output_path, "w:gz") as tar:
        for file in files:
            arcname = file.relative_to(REPORT_DIR.absolute())
            tar.add(file, arcname=arcname)

def main():
    if not REPORT_HTML.exists():
        print("No report found.")
        return

    print("Scraping report.html...")
    files = extract_paths_from_html(REPORT_HTML.absolute())

    for f in files: print(f)

    print(f"Found {len(files)} files. Creating {TARBALL}...")
    create_tarball(files, TARBALL)

    print("Done.")

if __name__ == "__main__":
    main()
