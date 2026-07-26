#!/usr/bin/env python3

import argparse
import html
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import stat
import tempfile
import urllib.parse
import urllib.request
import zipfile


class GitHub:
    def __init__(self, repository, token):
        self.repository = repository
        self.headers = {
            "Accept": "application/vnd.github+json",
            "User-Agent": "alamo-pages-assembler",
            "X-GitHub-Api-Version": "2022-11-28",
        }
        if token:
            self.headers["Authorization"] = f"Bearer {token}"
        self.run_cache = {}
        self.run_artifact_cache = {}
        self.named_artifact_cache = {}

    def request(self, endpoint):
        request = urllib.request.Request(
            f"https://api.github.com/repos/{self.repository}/{endpoint}",
            headers=self.headers,
        )
        with urllib.request.urlopen(request) as response:
            return json.load(response)

    def paginate(self, endpoint, key):
        separator = "&" if "?" in endpoint else "?"
        page = 1
        values = []
        while True:
            result = self.request(f"{endpoint}{separator}per_page=100&page={page}")
            batch = result[key] if isinstance(result, dict) else result
            values.extend(batch)
            if len(batch) < 100:
                return values
            page += 1

    def run(self, run_id):
        run_id = int(run_id)
        if run_id not in self.run_cache:
            self.run_cache[run_id] = self.request(f"actions/runs/{run_id}")
        return self.run_cache[run_id]

    def artifacts_named(self, name):
        if name not in self.named_artifact_cache:
            encoded = urllib.parse.quote(name)
            artifacts = self.paginate(f"actions/artifacts?name={encoded}", "artifacts")
            self.named_artifact_cache[name] = sorted(
                (item for item in artifacts if not item["expired"]),
                key=lambda item: item["created_at"],
                reverse=True,
            )
        return self.named_artifact_cache[name]

    def run_artifacts(self, run_id):
        run_id = int(run_id)
        if run_id not in self.run_artifact_cache:
            artifacts = self.paginate(
                f"actions/runs/{run_id}/artifacts", "artifacts"
            )
            self.run_artifact_cache[run_id] = {
                item["name"]: item for item in artifacts if not item["expired"]
            }
        return self.run_artifact_cache[run_id]

    def open_pull_requests(self):
        return self.paginate("pulls?state=open", "")

    def extract_artifact(self, artifact, destination):
        destination = destination.resolve()
        destination.mkdir(parents=True, exist_ok=True)
        request = urllib.request.Request(
            artifact["archive_download_url"], headers=self.headers
        )
        with tempfile.TemporaryFile() as archive:
            with urllib.request.urlopen(request) as response:
                shutil.copyfileobj(response, archive)
            archive.seek(0)
            with zipfile.ZipFile(archive) as zipped:
                total_size = sum(item.file_size for item in zipped.infolist())
                if total_size > 2 * 1024**3:
                    raise RuntimeError(
                        f"Artifact {artifact['name']} exceeds the extraction limit"
                    )
                for item in zipped.infolist():
                    relative = PurePosixPath(item.filename)
                    if relative.is_absolute() or ".." in relative.parts:
                        raise RuntimeError(
                            f"Unsafe path in artifact {artifact['name']}: {item.filename}"
                        )
                    mode = item.external_attr >> 16
                    if stat.S_ISLNK(mode):
                        raise RuntimeError(
                            f"Symlink in artifact {artifact['name']}: {item.filename}"
                        )
                    target = destination.joinpath(*relative.parts)
                    if destination not in target.resolve().parents and target != destination:
                        raise RuntimeError(
                            f"Path escapes artifact destination: {item.filename}"
                        )
                    if item.is_dir():
                        target.mkdir(parents=True, exist_ok=True)
                        continue
                    target.parent.mkdir(parents=True, exist_ok=True)
                    with zipped.open(item) as source, target.open("wb") as output:
                        shutil.copyfileobj(source, output)


def output(name, value):
    output_file = os.environ.get("GITHUB_OUTPUT")
    if output_file:
        with open(output_file, "a", encoding="utf-8") as stream:
            stream.write(f"{name}={value}\n")
    else:
        print(f"{name}={value}")


def valid_push_run(run, branch):
    return (
        run["event"] == "push"
        and run["head_branch"] == branch
        and run["conclusion"] == "success"
    )


def valid_publisher_run(run):
    return (
        run.get("path") == ".github/workflows/docs-pages.yml"
        and run["name"] == "Publish Documentation"
        and run["conclusion"] == "success"
    )


def valid_branch_docs(client, branch, trusted_run_id=None):
    name = f"docs-{branch}"
    if trusted_run_id:
        artifact = client.run_artifacts(trusted_run_id).get(name)
        if artifact:
            return artifact
    for artifact in client.artifacts_named(name):
        run_id = artifact["workflow_run"]["id"]
        run = client.run(run_id)
        if valid_push_run(run, branch) or valid_publisher_run(run):
            return artifact
    return None


def branch_version(client, branch, trusted_run_id=None):
    docs_name = f"docs-{branch}"
    coverage_name = f"coverage-{branch}"

    for coverage in client.artifacts_named(coverage_name):
        run_id = coverage["workflow_run"]["id"]
        run = client.run(run_id)
        if not valid_push_run(run, branch):
            continue
        docs = client.run_artifacts(run_id).get(docs_name)
        if docs:
            return {
                "slug": branch,
                "label": branch,
                "kind": "branch",
                "sha": run["head_sha"],
                "run_id": run_id,
                "docs": docs,
                "coverage": coverage,
            }

    docs = valid_branch_docs(client, branch, trusted_run_id)
    if not docs:
        return None
    run = client.run(docs["workflow_run"]["id"])
    source_sha = run["head_sha"]
    if trusted_run_id and int(run["id"]) == int(trusted_run_id):
        source_sha = client.request(f"commits/{urllib.parse.quote(branch)}")["sha"]
    return {
        "slug": branch,
        "label": branch,
        "kind": "branch",
        "sha": source_sha,
        "run_id": run["id"],
        "docs": docs,
        "coverage": None,
    }


def preview_allowed(pull, repository):
    source = pull.get("head", {}).get("repo")
    if source and source.get("full_name") == repository:
        return True
    return any(label["name"] == "docs-preview" for label in pull["labels"])


def pull_request_version(client, pull, repository):
    if not preview_allowed(pull, repository):
        return None

    number = pull["number"]
    head_sha = pull["head"]["sha"]
    docs_name = f"docs-pr-{number}"
    coverage_name = f"coverage-pr-{number}"
    for docs in client.artifacts_named(docs_name):
        if docs["workflow_run"]["head_sha"] != head_sha:
            continue
        run_id = docs["workflow_run"]["id"]
        run = client.run(run_id)
        run_pull_numbers = {item["number"] for item in run.get("pull_requests", [])}
        if (
            run["event"] != "pull_request"
            or run["conclusion"] != "success"
            or run["head_sha"] != head_sha
            or number not in run_pull_numbers
        ):
            continue
        coverage = client.run_artifacts(run_id).get(coverage_name)
        if not coverage:
            continue
        return {
            "slug": f"pr-{number}",
            "label": f"PR #{number}",
            "kind": "pull_request",
            "number": number,
            "title": pull["title"],
            "html_url": pull["html_url"],
            "sha": head_sha,
            "run_id": run_id,
            "docs": docs,
            "coverage": coverage,
        }
    return None


def write_redirect(path, target, title):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                "<!doctype html>",
                '<meta charset="utf-8">',
                f'<meta http-equiv="refresh" content="0; url={html.escape(target)}">',
                f'<link rel="canonical" href="{html.escape(target)}">',
                f"<title>{html.escape(title)}</title>",
                "",
            ]
        ),
        encoding="utf-8",
    )


def write_versions_page(site, versions):
    rows = []
    for version in versions:
        description = "Published branch"
        if version["kind"] == "pull_request":
            description = (
                f'<a href="{html.escape(version["html_url"])}">'
                f'{html.escape(version["title"])}</a>'
            )
        reports = [f'<a href="docs/{version["slug"]}/Builder.html">Input builder</a>']
        reports.append(
            f'<a href="docs/{version["slug"]}/cov/">Coverage</a>'
            if version["coverage"]
            else "Coverage unavailable"
        )
        rows.append(
            "<tr>"
            f'<td><a href="docs/{version["slug"]}/">{html.escape(version["label"])}</a></td>'
            f"<td>{description}</td>"
            f"<td><code>{html.escape(version['sha'][:12])}</code></td>"
            f"<td>{' &middot; '.join(reports)}</td>"
            "</tr>"
        )

    site.joinpath("versions.html").write_text(
        "\n".join(
            [
                "<!doctype html>",
                '<html lang="en"><head><meta charset="utf-8">',
                '<meta name="viewport" content="width=device-width, initial-scale=1">',
                "<title>Alamo documentation versions</title>",
                "<style>",
                "body{font:15px system-ui,sans-serif;max-width:960px;margin:40px auto;padding:0 20px;color:#263238}",
                "table{width:100%;border-collapse:collapse}th,td{text-align:left;padding:10px;border-bottom:1px solid #d8dde2}",
                "th{background:#f3f6f8}a{color:#2878b5}code{font-size:13px}",
                "</style></head><body>",
                "<h1>Alamo documentation versions</h1>",
                "<table><thead><tr><th>Version</th><th>Source</th><th>Commit</th><th>Report</th></tr></thead><tbody>",
                *rows,
                "</tbody></table></body></html>",
            ]
        ),
        encoding="utf-8",
    )


def install_version(client, site, version):
    destination = site / "docs" / version["slug"]
    client.extract_artifact(version["docs"], destination)
    if not (destination / "index.html").is_file():
        raise RuntimeError(f"{version['docs']['name']} has no index.html")
    manifest = destination / "alamo-docs-manifest.json"
    if version["kind"] == "branch" and manifest.is_file():
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        version["sha"] = metadata["source_sha"]
    if version["coverage"]:
        client.extract_artifact(version["coverage"], destination / "cov")
        if not (destination / "cov" / "index.html").is_file():
            raise RuntimeError(f"{version['coverage']['name']} has no index.html")


def add_compatibility_paths(site):
    builders = site / "docs" / "development" / "_static" / "input-builders"
    schemas = site / "docs" / "development" / "_static" / "input-schemas"
    if not builders.is_dir() or not schemas.is_dir():
        raise RuntimeError("Development documentation has no input builder assets")
    shutil.copytree(builders, site / "inputs", dirs_exist_ok=True)
    shutil.copytree(schemas, site / "input-schemas", dirs_exist_ok=True)
    for builder in (site / "inputs").glob("*.html"):
        text = builder.read_text(encoding="utf-8")
        builder.write_text(
            text.replace("../../doxygen/", "../docs/development/doxygen/"),
            encoding="utf-8",
        )


def discover(args):
    client = GitHub(args.repository, args.token)
    output("master_docs_found", str(bool(valid_branch_docs(client, "master"))).lower())
    output(
        "development_docs_found",
        str(bool(valid_branch_docs(client, "development"))).lower(),
    )


def assemble(args):
    client = GitHub(args.repository, args.token)
    site = Path(args.site)
    site.mkdir(parents=True, exist_ok=True)

    versions = []
    for branch in ("master", "development"):
        version = branch_version(client, branch, args.trusted_run_id)
        if not version:
            raise RuntimeError(f"No valid documentation artifact for {branch}")
        versions.append(version)

    for pull in client.open_pull_requests():
        version = pull_request_version(client, pull, args.repository)
        if version:
            versions.append(version)

    for version in versions:
        install_version(client, site, version)

    public_versions = [
        {key: value for key, value in version.items() if key not in {"docs", "coverage"}}
        | {"coverage": bool(version["coverage"])}
        for version in versions
    ]
    (site / "versions.json").write_text(
        json.dumps(public_versions, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    write_versions_page(site, versions)
    write_redirect(site / "index.html", "docs/development/", "Alamo documentation")
    write_redirect(site / "docs" / "index.html", "development/", "Alamo documentation")
    add_compatibility_paths(site)
    (site / ".nojekyll").touch()

    preview = next(
        (
            version
            for version in versions
            if version["kind"] == "pull_request"
            and args.trigger_run_id
            and int(version["run_id"]) == int(args.trigger_run_id)
        ),
        None,
    )
    trigger_pull_number = ""
    if args.trigger_run_id:
        trigger_run = client.run(args.trigger_run_id)
        trigger_pulls = trigger_run.get("pull_requests", [])
        if trigger_pulls:
            trigger_pull_number = trigger_pulls[0]["number"]

    output("preview_published", str(bool(preview)).lower())
    output(
        "preview_pr_number",
        preview["number"] if preview else trigger_pull_number,
    )
    output(
        "preview_slug",
        preview["slug"]
        if preview
        else (f"pr-{trigger_pull_number}" if trigger_pull_number else ""),
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("discover", "assemble"))
    parser.add_argument("--repository", required=True)
    parser.add_argument("--token", required=True)
    parser.add_argument("--site", default="site")
    parser.add_argument("--trusted-run-id")
    parser.add_argument("--trigger-run-id")
    args = parser.parse_args()
    if args.command == "discover":
        discover(args)
    else:
        assemble(args)


if __name__ == "__main__":
    main()
