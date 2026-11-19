#!/usr/bin/env python3
import os
import tarfile
import re
from pathlib import Path
import shutil
import argparse

parser = argparse.ArgumentParser(description='Bundle report into a single directory');
parser.add_argument('--name',default="report.tar.gz",help='Name of output. If ends in .tar.gz, output will be taball.')
args=parser.parse_args()


REPORT_HTML = Path("report.html")
REPORT_DIR = REPORT_HTML.parent


# Regex to find resource paths in HTML attributes
IMG_RE = re.compile(r'<img[^>]+src=["\']([^"\']+)["\']', re.IGNORECASE)
LINK_RE = re.compile(r'<a[^>]+href=["\']([^"\']+)["\']', re.IGNORECASE)
CSS_RE = re.compile(r'<link[^>]+href=["\']([^"\']+)["\']', re.IGNORECASE)
SCRIPT_RE = re.compile(r'<script[^>]+src=["\']([^"\']+)["\']', re.IGNORECASE)
SOURCE_RE = re.compile(r'<source[^>]+src=["\']([^"\']+)["\']', re.IGNORECASE)
IFRAME_RE = re.compile(r'<iframe[^>]+src=["\']([^"\']+)["\']', re.IGNORECASE)

SKIP_PREFIXES = ("http://", "https://", "mailto:", "data:", "javascript:")

def _normalize_ref(ref):
    ref = ref.strip()
    if not ref or ref.startswith(SKIP_PREFIXES):
        return None
    # Drop fragment and query string for filesystem lookup
    ref = ref.split("#", 1)[0].split("?", 1)[0]
    if not ref:
        return None
    return ref

def extract_paths_from_html(html_file, seen=None):
    if seen is None:
        seen = set()
    html_file = Path(html_file).absolute().resolve()
    if html_file in seen or not html_file.exists():
        return set()
    seen.add(html_file)

    with open(html_file, "r", encoding="utf-8", errors="ignore") as f:
        html = f.read()

    refs = []
    for regex in (IMG_RE, LINK_RE, CSS_RE, SCRIPT_RE, SOURCE_RE, IFRAME_RE):
        refs.extend(regex.findall(html))

    paths = set()
    for ref in refs:
        ref = _normalize_ref(ref)
        if not ref:
            continue
        resolved = (html_file.parent / Path(ref)).absolute().resolve()
        if resolved.exists():
            paths.add(resolved)

    # Recurse into linked HTML files
    subpaths = set()
    for path in paths:
        if path.is_file() and path.suffix.lower() in {".html", ".htm"}:
            subpaths.update(extract_paths_from_html(path, seen))

    paths.update(subpaths)
    paths.add(html_file)
    return paths

def create_tarball(files, output_path):
    with tarfile.open(output_path, "w:gz") as tar:
        for file in files:
            arcname = file.relative_to(REPORT_DIR.absolute())
            tar.add(file, arcname=arcname)

def copy_to_directory(files, output_path):
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    base_dir = REPORT_DIR.resolve()
    for file in files:
        # Preserve the relative path (same as your tarball version)
        try:
            rel = file.relative_to(base_dir)
        except ValueError:
            # Fallback: place files outside the base at the top level.
            rel = file.name
        dest = output_path / rel
        # Ensure parent directories exist
        dest.parent.mkdir(parents=True, exist_ok=True)
        # Copy file
        shutil.copy2(file, dest)

def main():
    if not REPORT_HTML.exists():
        print("No report found.")
        return

    print("Scraping report.html...")
    files = set()
    if REPORT_HTML.exists():
        files.update(extract_paths_from_html(REPORT_HTML.absolute()))
    # Support report/report.html as an entry point if present
    report_summary = Path("report") / "report.html"
    if report_summary.exists():
        files.update(extract_paths_from_html(report_summary.absolute()))
    files = [p for p in files if p.exists()]

    for f in files: print(f)

    print(f"Found {len(files)} files. Creating {args.name}...")
    if args.name.endswith(".tar.gz"):
        create_tarball(files, args.name)
    else:
        copy_to_directory(files,args.name)
    

    print("Done.")

if __name__ == "__main__":
    main()
