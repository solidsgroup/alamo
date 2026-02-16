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
    files = extract_paths_from_html(REPORT_HTML.absolute())

    for f in files: print(f)

    print(f"Found {len(files)} files. Creating {args.name}...")
    if args.name.endswith(".tar.gz"):
        create_tarball(files, args.name)
    else:
        copy_to_directory(files,args.name)
    

    print("Done.")

if __name__ == "__main__":
    main()
