#!/usr/bin/env python3
"""Generate the Sphinx input reference from Alamo input schemas."""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import os
import re
import shutil
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable


ICONS = {
    "BC": ":fas:`border-top-left;fa-fw`",
    "IC": ":fas:`circle-right;fa-fw`",
    "IO": ":fas:`print;fa-fw`",
    "Integrator": ":fas:`gear;fa-fw`",
    "Model": ":fas:`panorama;fa-fw`",
    "Numeric": ":fas:`calculator;fa-fw`",
    "Operator": ":far:`map;fa-fw`",
    "Set": ":fas:`braille;fa-fw`",
    "Solver": ":fas:`diamond-turn-right;fa-fw`",
    "Util": ":fas:`sliders;fa-fw`",
}
HTML_ICONS = {
    "BC": "fa-border-top-left",
    "IC": "fa-circle-right",
    "IO": "fa-print",
    "Integrator": "fa-gear",
    "Model": "fa-panorama",
    "Numeric": "fa-calculator",
    "Operator": "fa-map",
    "Set": "fa-braille",
    "Solver": "fa-diamond-turn-right",
    "Util": "fa-sliders",
}
BADGE_DEFAULT = "input-reference-badge input-reference-badge-default"
BADGE_OPTION = "input-reference-badge input-reference-badge-option"
BADGE_META = "input-reference-badge input-reference-badge-meta"
BADGE_REQUIRED = "input-reference-badge input-reference-badge-required"


@dataclass
class Entry:
    name: str
    kind: str
    directive: str
    source_file: str
    source_line: int
    description: str = ""
    required: bool = False
    has_default: bool = False
    default_value: str = ""
    options: set[str] = field(default_factory=set)
    conditional: bool = False
    members: tuple[str, ...] = ()
    count: int = 0
    resolved_names: dict[str, set[str]] = field(
        default_factory=lambda: defaultdict(set)
    )


def parse_args() -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=repo_root)
    parser.add_argument(
        "--schema-dir",
        type=Path,
        default=repo_root / "docs/source/_static/input-schemas",
    )
    parser.add_argument("--schema", type=Path, action="append", default=[])
    parser.add_argument(
        "--output",
        type=Path,
        default=repo_root / "docs/source/Inputs.generated.rst",
        help="Generated toctree fragment included by Inputs.rst.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=repo_root / "docs/source/InputsReference",
        help="Directory for generated namespace and class pages.",
    )
    return parser.parse_args()


def normalized_source_file(value: str) -> str:
    value = value.replace("\\", "/")
    while value.startswith("./"):
        value = value[2:]
    return value


def owner_for_file(source_file: str) -> str:
    path = normalized_source_file(source_file)
    if path.startswith("src/"):
        path = path[4:]
    path = re.sub(r"\.(?:H|HH|HPP|h|hh|hpp|cpp)$", "", path)
    return path.replace("/", "::")


def source_path(repo_root: Path, source_file: str) -> Path | None:
    path = Path(normalized_source_file(source_file))
    candidates = [path] if path.is_absolute() else [repo_root / path]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return None


def load_doxygen_sources(doxygen_dir: Path) -> dict[str, str]:
    sources: dict[str, str] = {}
    title_pattern = re.compile(r"<title>Alamo: (.+) Source File</title>")
    for page in doxygen_dir.glob("*source.html"):
        contents = page.read_text(encoding="utf-8", errors="replace")
        match = title_pattern.search(contents)
        if match:
            source_file = normalized_source_file(html.unescape(match.group(1)))
            sources[source_file] = page.name
    return sources


def doxygen_source_url(
    page: Path,
    docs_source_dir: Path,
    source_file: str,
    source_line: int,
    doxygen_sources: dict[str, str],
) -> str:
    filename = doxygen_sources.get(normalized_source_file(source_file))
    if not filename:
        return ""
    html_page = page.relative_to(docs_source_dir).with_suffix(".html")
    target = Path("doxygen") / filename
    url = Path(os.path.relpath(target, start=html_page.parent)).as_posix()
    if source_line > 0:
        url += f"#l{source_line:05d}"
    return url


def leading_documentation(path: Path) -> str:
    parts: list[str] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("///"):
            parts.append(line[3:])
        elif line.startswith("// "):
            parts.append(line[3:])
        elif line.startswith("//"):
            parts.append(line[2:])
        else:
            break

    normalized: list[str] = []
    math_content = False
    for line in parts:
        stripped = line.strip()
        if stripped.startswith(".. ") and normalized and normalized[-1].strip():
            normalized.append("")
        if stripped == ":list: none":
            line = line.replace(":list: none", ":list: bullet")
        if math_content and stripped and not line.startswith((" ", "\t")):
            line = "   " + line
        normalized.append(line)
        if stripped == ".. math::":
            math_content = True
        elif math_content and stripped:
            math_content = False

    documentation = "\n".join(normalized).strip()
    return re.sub(
        r"(`[^`\n]+ <https?://[^>`]+>`)_",
        r"\1__",
        documentation,
    )


def root_objects(node: Any) -> Iterable[dict[str, Any]]:
    if not isinstance(node, dict):
        return
    yield node
    for child in node.get("children", []):
        yield from root_objects(child)


def is_conditional(value: Any) -> bool:
    return isinstance(value, list) and any(bool(context) for context in value)


def local_member_name(path: str) -> str:
    return path.rsplit(".", 1)[-1]


def declared_name(
    repo_root: Path,
    source_file: str,
    source_line: int,
    fallback: str,
) -> str:
    path = source_path(repo_root, source_file)
    if path is None or source_line < 1:
        return fallback
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if source_line > len(lines):
        return fallback
    call_lines: list[str] = []
    for line in lines[source_line - 1 : source_line + 5]:
        call_lines.append(line)
        if ";" in line:
            break
    call = " ".join(call_lines)
    match = re.search(r'"((?:\\.|[^"\\])*)"', call)
    if match is None:
        return fallback
    return match.group(1)


def merge_entry(entries: dict[tuple[Any, ...], Entry], entry: Entry) -> None:
    key = (
        entry.name,
        entry.kind,
        entry.directive,
        entry.source_file,
        entry.source_line,
        entry.members,
    )
    current = entries.get(key)
    if current is None:
        entries[key] = entry
        return
    if len(entry.description) > len(current.description):
        current.description = entry.description
    current.required = current.required or entry.required
    current.has_default = current.has_default or entry.has_default
    if entry.default_value:
        current.default_value = entry.default_value
    current.options.update(entry.options)
    current.conditional = current.conditional or entry.conditional
    for executable, names in entry.resolved_names.items():
        current.resolved_names[executable].update(names)


def load_entries(
    repo_root: Path,
    schemas: list[Path],
) -> tuple[dict[str, list[Entry]], dict[str, list[str]]]:
    grouped: dict[str, dict[tuple[Any, ...], Entry]] = {}
    ignored: dict[str, list[str]] = {}

    for schema_path in schemas:
        executable = schema_path.name.removesuffix(".schema.json")
        schema = json.loads(schema_path.read_text(encoding="utf-8"))
        for node in root_objects(schema.get("root", {})):
            source = node.get("source")
            directive = str(node.get("directive", ""))
            if (
                isinstance(source, dict)
                and directive
                and node.get("kind") != "scope"
                and directive != "query_exactly"
            ):
                source_file = normalized_source_file(str(source.get("file", "")))
                if source_file.endswith((".cc", ".cxx")):
                    continue
                source_line = int(source.get("line", 0))
                name = declared_name(
                    repo_root,
                    source_file,
                    source_line,
                    str(node.get("name", "")),
                )
                entry = Entry(
                    name=name,
                    kind=str(node.get("kind", "parameter")),
                    directive=directive,
                    source_file=source_file,
                    source_line=source_line,
                    description=str(node.get("doc") or node.get("description") or ""),
                    required=bool(node.get("required")),
                    has_default=bool(node.get("has_default")),
                    default_value=str(node.get("default_value", "")),
                    options={str(option) for option in node.get("options", [])},
                    conditional=is_conditional(node.get("contexts", [])),
                )
                entry.resolved_names[executable].add(str(node.get("path", name)))
                owner = owner_for_file(source_file)
                merge_entry(grouped.setdefault(owner, {}), entry)

            for constraint in node.get("constraints", []):
                if not isinstance(constraint, dict):
                    continue
                constraint_source = constraint.get("source")
                if not isinstance(constraint_source, dict):
                    continue
                source_file = normalized_source_file(
                    str(constraint_source.get("file", ""))
                )
                if source_file.endswith((".cc", ".cxx")):
                    continue
                resolved_members = tuple(
                    str(member) for member in constraint.get("members", [])
                )
                members = tuple(local_member_name(member) for member in resolved_members)
                entry = Entry(
                    name=" / ".join(members),
                    kind="constraint",
                    directive=str(constraint.get("kind", "constraint")),
                    source_file=source_file,
                    source_line=int(constraint_source.get("line", 0)),
                    description=str(
                        constraint.get("doc") or constraint.get("description") or ""
                    ),
                    required=True,
                    conditional=bool(constraint.get("conditions", [])),
                    members=members,
                    count=int(constraint.get("count", 0)),
                )
                entry.resolved_names[executable].add(" / ".join(resolved_members))
                owner = owner_for_file(source_file)
                merge_entry(grouped.setdefault(owner, {}), entry)

        for item in schema.get("traversal_ignored", []):
            source = item.get("source", {}) if isinstance(item, dict) else {}
            source_file = normalized_source_file(str(source.get("file", "")))
            if not source_file or source_file.endswith((".cc", ".cxx")):
                continue
            owner = owner_for_file(source_file)
            note = str(item.get("note", "")).strip()
            if note and note not in ignored.setdefault(owner, []):
                ignored[owner].append(note)

    return {
        owner: sorted(
            values.values(),
            key=lambda item: (item.source_line, item.name.lower(), item.directive),
        )
        for owner, values in grouped.items()
    }, ignored


def inline_doc(value: str) -> str:
    value = re.sub(r"\s+", " ", value).strip()
    math: list[str] = []

    def preserve_math(match: re.Match[str]) -> str:
        math.append(match.group(1))
        return f"ALAMOMATHTOKEN{len(math) - 1}END"

    value = re.sub(r":math:`([^`]+)`", preserve_math, value)
    value = re.sub(r"\\\((.+?)\\\)", preserve_math, value)
    value = html.escape(value)
    value = re.sub(r":code:`([^`]+)`", r"<code>\1</code>", value)
    value = re.sub(r":(?:ref|doc):`([^`]+)`", r"\1", value)
    for index, expression in enumerate(math):
        rendered = (
            '<span class="math notranslate nohighlight">'
            f"\\({html.escape(expression)}\\)</span>"
        )
        value = value.replace(f"ALAMOMATHTOKEN{index}END", rendered)
    return value


def entry_anchor(entry: Entry) -> str:
    identity = "|".join(
        [
            entry.source_file,
            str(entry.source_line),
            entry.directive,
            entry.name,
            *entry.members,
        ]
    )
    digest = hashlib.sha1(identity.encode("utf-8")).hexdigest()[:10]
    name = re.sub(r"[^a-z0-9]+", "-", entry.name.lower()).strip("-")
    return f"input-{name[:40] or 'parameter'}-{digest}"


def badge(css_class: str, text: str) -> str:
    return f'<span class="{css_class}">{html.escape(text)}</span>'


def metadata_html(entry: Entry) -> str:
    values: list[str] = []
    if entry.kind == "constraint":
        values.append(badge(BADGE_REQUIRED, f"specify exactly {entry.count}"))
    elif entry.required:
        values.append(badge(BADGE_REQUIRED, "required"))
    if entry.has_default:
        values.append(badge(BADGE_DEFAULT, entry.default_value or "default"))
    for option in sorted(entry.options):
        if entry.has_default and option == entry.default_value:
            continue
        values.append(badge(BADGE_OPTION, option))
    if entry.kind == "sequence":
        values.append(badge(BADGE_META, "repeatable"))
    if entry.directive == "query_file":
        values.append(badge(BADGE_META, "file path"))
    if entry.conditional:
        values.append(badge(BADGE_META, "conditional"))
    return " ".join(values)


def resolved_names_html(entry: Entry) -> str:
    executable_count = len(entry.resolved_names)
    use_count = sum(len(names) for names in entry.resolved_names.values())
    noun = "executable" if executable_count == 1 else "executables"
    use_noun = "use" if use_count == 1 else "uses"
    lines = [
        '<details class="input-reference-resolved">',
        (
            "  <summary>Names by executable "
            f'<span>{executable_count} {noun}, {use_count} {use_noun}</span></summary>'
        ),
        '  <div class="input-reference-resolved-body">',
    ]
    for executable, names in sorted(entry.resolved_names.items()):
        lines.extend(
            [
                '    <div class="input-reference-executable">',
                f"      <strong>{html.escape(executable)}</strong>",
                '      <div class="input-reference-paths">',
            ]
        )
        for name in sorted(names):
            lines.append(f"        <code>{html.escape(name)}</code>")
        lines.extend(["      </div>", "    </div>"])
    lines.extend(["  </div>", "</details>"])
    return "\n".join(lines)


def input_panels(
    entries: list[Entry],
    page: Path,
    docs_source_dir: Path,
    doxygen_sources: dict[str, str],
) -> list[str]:
    lines = [".. raw:: html", "", '   <div class="input-reference-list">']
    for entry in entries:
        description = inline_doc(entry.description)
        if not description:
            description = "<em>No documentation available.</em>"
        metadata = metadata_html(entry)
        source_url = doxygen_source_url(
            page,
            docs_source_dir,
            entry.source_file,
            entry.source_line,
            doxygen_sources,
        )
        if source_url:
            input_name = (
                f'<a class="input-reference-name" href="{html.escape(source_url)}">'
                f"<code>{html.escape(entry.name)}</code></a>"
            )
        else:
            input_name = (
                f'<code class="input-reference-name">{html.escape(entry.name)}</code>'
            )
        lines.extend(
            [
                (
                    f'     <article class="input-reference-item" '
                    f'id="{entry_anchor(entry)}">'
                ),
                '       <header class="input-reference-item-header">',
                f"         {input_name}",
                f'         <div class="input-reference-badges">{metadata}</div>',
                "       </header>",
                f'       <p class="input-reference-description">{description}</p>',
            ]
        )
        lines.extend(f"       {line}" for line in resolved_names_html(entry).splitlines())
        lines.extend(["     </article>"])
    lines.extend(["   </div>", ""])
    return lines


def title_lines(title: str) -> list[str]:
    return [title, "=" * len(title), ""]


def owner_page_path(output_dir: Path, owner: str, is_namespace: bool) -> Path:
    parts = owner.split("::")
    if is_namespace:
        return output_dir.joinpath(*parts, "index.rst")
    return output_dir.joinpath(*parts[:-1], parts[-1] + ".rst")


def write_owner_content(
    page: Path,
    owner: str,
    entries: list[Entry],
    notes: list[str],
    repo_root: Path,
    docs_source_dir: Path,
    doxygen_sources: dict[str, str],
) -> list[str]:
    lines: list[str] = []
    source_files = sorted({entry.source_file for entry in entries})
    for source_file in source_files:
        source_url = doxygen_source_url(
            page, docs_source_dir, source_file, 0, doxygen_sources
        )
        if source_url:
            lines.extend(
                [f":bdg-link-secondary-line:`{source_file} <{source_url}>`", ""]
            )
        else:
            lines.extend([f":bdg-secondary-line:`{source_file}`", ""])
        path = source_path(repo_root, source_file)
        documentation = leading_documentation(path) if path else ""
        if documentation:
            lines.extend([documentation, ""])
    for note in notes:
        lines.extend([".. warning::", "", f"   {note}", ""])
    if entries:
        if any(
            re.search(r":math:`|\\\(", entry.description)
            for entry in entries
        ):
            lines.extend(
                [
                    ".. container:: input-reference-math-loader",
                    "",
                    r"   :math:`\phantom{0}`",
                    "",
                ]
            )
        lines.extend(
            input_panels(entries, page, docs_source_dir, doxygen_sources)
        )
    return lines


def directory_html(
    page: Path,
    owners: list[str],
    output_dir: Path,
    namespace_owners: set[str],
    grouped: dict[str, list[Entry]],
) -> list[str]:
    lines = [
        ".. raw:: html",
        "",
        '   <nav class="input-reference-directory" aria-label="Input reference sections">',
    ]
    for owner in owners:
        is_namespace = owner in namespace_owners
        target_page = owner_page_path(output_dir, owner, is_namespace)
        target = target_page.with_suffix(".html").relative_to(page.parent).as_posix()
        count = sum(
            len(entries)
            for entry_owner, entries in grouped.items()
            if entry_owner == owner or entry_owner.startswith(owner + "::")
        )
        noun = "input" if count == 1 else "inputs"
        icon = HTML_ICONS.get(owner, "fa-folder-tree" if is_namespace else "fa-file-lines")
        lines.extend(
            [
                f'     <a class="input-reference-directory-item" href="{html.escape(target)}">',
                f'       <span class="fas {icon} fa-fw" aria-hidden="true"></span>',
                "       <span>",
                f"         <strong>{html.escape(owner)}</strong>",
                f"         <small>{count} {noun}</small>",
                "       </span>",
                '       <span class="fas fa-chevron-right" aria-hidden="true"></span>',
                "     </a>",
            ]
        )
    lines.extend(["   </nav>", ""])
    return lines


def generate_pages(
    repo_root: Path,
    output: Path,
    output_dir: Path,
    grouped: dict[str, list[Entry]],
    ignored: dict[str, list[str]],
) -> tuple[int, int]:
    docs_source_dir = output.parent
    doxygen_sources = load_doxygen_sources(
        repo_root / "docs/build/html/doxygen"
    )
    owners = sorted(set(grouped) | set(ignored))
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not owners:
        output.write_text(
            ".. note::\n\n"
            "   No generated input schemas are available. Run "
            "``make docs-input-builders`` before publishing the input reference.\n",
            encoding="utf-8",
        )
        return 0, 0

    namespace_children: dict[str, set[str]] = defaultdict(set)
    namespace_owners: set[str] = set()
    for owner in owners:
        parts = owner.split("::")
        for depth in range(1, len(parts)):
            namespace = "::".join(parts[:depth])
            child = "::".join(parts[: depth + 1])
            namespace_owners.add(namespace)
            namespace_children[namespace].add(child)

    roots = sorted({owner.split("::", 1)[0] for owner in owners})
    root_lines = [
        ".. raw:: html",
        "",
        '   <h2 class="input-reference-directory-title">Input namespaces</h2>',
        "",
    ]
    root_lines.extend(
        directory_html(output, roots, output_dir, namespace_owners, grouped)
    )
    root_lines.extend(
        [
            ".. toctree::",
            "   :maxdepth: 10",
            "   :hidden:",
            "",
        ]
    )
    for root in roots:
        root_is_namespace = root in namespace_owners
        root_page = owner_page_path(output_dir, root, root_is_namespace)
        root_target = root_page.with_suffix("").relative_to(output.parent)
        root_lines.append(f"   {root_target.as_posix()}")
    output.write_text("\n".join(root_lines) + "\n", encoding="utf-8")

    page_count = 0
    for owner in owners:
        is_namespace = owner in namespace_owners
        page = owner_page_path(output_dir, owner, is_namespace)
        page.parent.mkdir(parents=True, exist_ok=True)
        title = owner
        if "::" not in owner:
            title = f"{ICONS.get(owner, '')} {owner}".strip()
        lines = title_lines(title)
        lines.extend(
            write_owner_content(
                page,
                owner,
                grouped.get(owner, []),
                ignored.get(owner, []),
                repo_root,
                docs_source_dir,
                doxygen_sources,
            )
        )
        children = sorted(namespace_children.get(owner, set()))
        if children:
            lines.extend(
                directory_html(page, children, output_dir, namespace_owners, grouped)
            )
            lines.extend(
                [".. toctree::", "   :maxdepth: 10", "   :hidden:", ""]
            )
            for child in children:
                child_is_namespace = child in namespace_owners
                child_page = owner_page_path(output_dir, child, child_is_namespace)
                target = child_page.with_suffix("").relative_to(page.parent)
                lines.append(f"   {target.as_posix()}")
            lines.append("")
        page.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
        page_count += 1

    for namespace in sorted(namespace_owners - set(owners)):
        page = owner_page_path(output_dir, namespace, True)
        page.parent.mkdir(parents=True, exist_ok=True)
        title = namespace
        if "::" not in namespace:
            title = f"{ICONS.get(namespace, '')} {namespace}".strip()
        lines = title_lines(title)
        children = sorted(namespace_children.get(namespace, set()))
        lines.extend(
            directory_html(page, children, output_dir, namespace_owners, grouped)
        )
        lines.extend([".. toctree::", "   :maxdepth: 10", "   :hidden:", ""])
        for child in children:
            child_page = owner_page_path(
                output_dir, child, child in namespace_owners
            )
            target = child_page.with_suffix("").relative_to(page.parent)
            lines.append(f"   {target.as_posix()}")
        lines.append("")
        page.write_text("\n".join(lines), encoding="utf-8")
        page_count += 1

    return page_count, sum(len(entries) for entries in grouped.values())


def main() -> None:
    args = parse_args()
    repo_root = args.repo_root.resolve()
    schema_dir = args.schema_dir
    if not schema_dir.is_absolute():
        schema_dir = repo_root / schema_dir
    schemas = [
        path if path.is_absolute() else repo_root / path for path in args.schema
    ]
    if not schemas:
        schemas = sorted(schema_dir.glob("*.schema.json"))
    else:
        schemas = [path for path in schemas if path.is_file()]

    grouped, ignored = load_entries(repo_root, schemas)
    output = args.output if args.output.is_absolute() else repo_root / args.output
    output_dir = (
        args.output_dir
        if args.output_dir.is_absolute()
        else repo_root / args.output_dir
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    page_count, entry_count = generate_pages(
        repo_root, output, output_dir, grouped, ignored
    )
    print(
        f"Wrote {page_count} input reference pages from {len(schemas)} schemas "
        f"({entry_count} source declarations)"
    )


if __name__ == "__main__":
    main()
