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
    paths.update(IMG_RE.findall(html))
    paths.update(LINK_RE.findall(html))
    # Always include the HTML file itself
    paths.add(str(REPORT_HTML.relative_to(REPORT_DIR)))
    return [REPORT_DIR / p for p in paths if (REPORT_DIR / p).exists()]

def create_tarball(files, output_path):
    with tarfile.open(output_path, "w:gz") as tar:
        for file in files:
            arcname = file.relative_to(REPORT_DIR.parent)
            tar.add(file, arcname=arcname)

def copy_to_directory(files, output_path):
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    for file in files:
        # Preserve the relative path (same as your tarball version)
        rel = file.relative_to(REPORT_DIR.parent)
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
    files = extract_paths_from_html(REPORT_HTML)

    print(f"Found {len(files)} files. Creating {args.name}...")
    if args.name.endswith(".tar.gz"):
        create_tarball(files, args.name)
    else:
        copy_to_directory(files,args.name)
    

    print("Done.")

if __name__ == "__main__":
    main()
