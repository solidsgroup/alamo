#!/usr/bin/python
from __future__ import annotations
import html
import json
import os 
import glob
import re
import sys
from os import listdir
from os.path import isfile, join
import io
import configparser
from collections import OrderedDict
from pathlib import Path
import fnmatch

from pygments import highlight
from pygments.formatters import HtmlFormatter
from pygments.lexers import MakefileLexer

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from make_schema_reference import (
    entry_anchor,
    load_entries,
    normalized_source_file,
    owner_page_path,
    root_objects,
)


def load_input_reference_targets():
    schema_dir = REPO_ROOT / "docs/source/_static/input-schemas"
    schemas = sorted(schema_dir.glob("*.schema.json"))
    grouped, _ = load_entries(REPO_ROOT, schemas)
    owners = set(grouped)
    namespace_owners = set()
    for owner in owners:
        parts = owner.split("::")
        for depth in range(1, len(parts)):
            namespace_owners.add("::".join(parts[:depth]))

    targets = {}
    for owner, entries in grouped.items():
        page = owner_page_path(
            Path("InputsReference"),
            owner,
            owner in namespace_owners,
        ).with_suffix(".html")
        for entry in entries:
            target = f"../{page.as_posix()}#{entry_anchor(entry)}"
            for executable, names in entry.resolved_names.items():
                executable_targets = targets.setdefault(executable, {})
                for name in names:
                    executable_targets.setdefault(name, []).append(
                        {
                            "target": target,
                            "source_file": entry.source_file,
                            "source_line": entry.source_line,
                            "directive": entry.directive,
                            "contexts": [],
                        }
                    )

    defaults = {}
    for schema_path in schemas:
        executable = schema_path.name.removesuffix(".schema.json")
        schema = json.loads(schema_path.read_text(encoding="utf-8"))
        default_candidates = {}
        for node in root_objects(schema.get("root", {})):
            path = str(node.get("path", ""))
            if path and node.get("has_default"):
                default_value = node.get("default_value")
                if default_value is None and node.get("options"):
                    default_value = node["options"][0]
                default_candidates.setdefault(path, set()).add(
                    normalize_input_value(str(default_value or ""))
                )

            source = node.get("source")
            directive = str(node.get("directive", ""))
            if not isinstance(source, dict) or not directive or not path:
                continue
            source_file = normalized_source_file(str(source.get("file", "")))
            source_line = int(source.get("line", 0))
            for record in targets.get(executable, {}).get(path, []):
                if (
                    record["source_file"] == source_file
                    and record["source_line"] == source_line
                    and record["directive"] == directive
                ):
                    record["contexts"].extend(node.get("contexts", []))

        defaults[executable] = {
            path: next(iter(values))
            for path, values in default_candidates.items()
            if len(values) == 1
        }
    return targets, defaults


def normalize_input_value(value):
    value = value.strip().strip("\"'").strip().lower()
    if value in {"true", "yes", "on"}:
        return "1"
    if value in {"false", "no", "off"}:
        return "0"
    return value


INPUT_REFERENCE_TARGETS, INPUT_DEFAULTS = load_input_reference_targets()


def parsed_input_values(path):
    values = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = re.match(r"\s*([^\s=#]+)\s*=\s*(.*)$", line)
        if not match:
            continue
        value = re.split(r"\s+#", match.group(2), maxsplit=1)[0]
        values[match.group(1)] = normalize_input_value(value)
    return values


def test_runs(config, path):
    base_values = parsed_input_values(path)
    runs = []
    for section in config:
        if section == "DEFAULT":
            continue
        executable = config[section].get("exe", "alamo")
        dimension = config[section].get("dim", "3")
        prefix = f"{executable}-{dimension}d-"
        schemas = {
            name for name in INPUT_REFERENCE_TARGETS if name.startswith(prefix)
        }
        for schema in schemas:
            values = dict(INPUT_DEFAULTS.get(schema, {}))
            values.update(base_values)
            for argument in config[section].get("args", "").splitlines():
                match = re.match(r"\s*([^\s=#]+)\s*=\s*(.*)$", argument)
                if match:
                    values[match.group(1)] = normalize_input_value(match.group(2))
            runs.append((schema, values))
    if runs:
        return runs
    return [
        (schema, {**INPUT_DEFAULTS.get(schema, {}), **base_values})
        for schema in INPUT_REFERENCE_TARGETS
    ]


def context_matches(contexts, values):
    if not contexts:
        return True
    return any(
        all(
            condition.get("path") not in values
            or values[condition["path"]]
            == normalize_input_value(str(condition.get("value", "")))
            for condition in context
        )
        for context in contexts
    )


def input_reference_target(name, runs):
    matches = set()
    for schema, values in runs:
        for record in INPUT_REFERENCE_TARGETS.get(schema, {}).get(name, []):
            if context_matches(record["contexts"], values):
                matches.add(record["target"])
    if len(matches) == 1:
        return matches.pop()
    return None


def write_linked_input(testdocfile, path, caption, runs, block_id):
    source = path.read_text(encoding="utf-8", errors="replace")
    rendered = highlight(
        source,
        MakefileLexer(),
        HtmlFormatter(nowrap=True),
    )
    names = {
        match.group(1)
        for line in source.splitlines()
        if (match := re.match(r"\s*([^\s=#]+)\s*=", line))
    }
    for name in names:
        target = input_reference_target(name, runs)
        if not target:
            continue
        token = f'<span class="nv">{html.escape(name)}</span>'
        link = (
            f'<a class="test-input-reference" href="{html.escape(target)}" '
            f'title="View input documentation">{token}</a>'
        )
        rendered = rendered.replace(token, link)

    contents = [
        f'<div class="literal-block-wrapper docutils container" id="{block_id}">',
        '<div class="code-block-caption">',
        f'<span class="caption-text">{html.escape(caption)}</span>',
        (
            f'<a class="headerlink" href="#{block_id}" '
            'title="Link to this code">&#182;</a>'
        ),
        "</div>",
        '<div class="highlight-makefile notranslate">',
        '<div class="highlight"><pre><span></span>',
        *rendered.splitlines(),
        "</pre></div>",
        "</div>",
        "</div>",
    ]
    testdocfile.write(".. raw:: html\n\n")
    for line in contents:
        testdocfile.write(f"   {line}\n")
    testdocfile.write("\n")


#
# Special order from SO - dictionary allows for keys to be specified multiple times and 
# config parser will read it all. (Copy/pasted from runtests.sh)
class MultiOrderedDict(OrderedDict):
    def __setitem__(self, key, value):
        if isinstance(value, list) and key in self:
            self[key].extend(value)
        else:
            super().__setitem__(key, value)


def find_files_ignore_case(target_dir: Path, file_name: str) -> list[Path]:
    """Get paths in `target_dir` that match `file_name`, ignoring case.

    Parameters
    ----------
    target_dir
        The directory to search for matching files.
    file_name
        The file name to search for in `target_dir`.

    Returns
    -------
    list[Path]
        File paths that match `file_name` in `target_dir`.

    """
    try:
        return [f for f in target_dir.iterdir()
            if f.name.lower() == file_name.lower() and Path.is_file(f)]
    except OSError:
        print(f"'{target_dir}' is not a directory or otherwise inaccessible.")
        return []

#def icon(str)

docfile    = open("Tests.rst","w")
docfile.write(r"""

.. _tests:

==============================
:fas:`flask;fa-fw` Tests
==============================

""")

headerchar = ["=","*","-","~","."]
written_headers = []

num_tot = 0
num_doc = 0


docfile.write(r"""


""")

#docfile.write(r"- A regular icon: :material-outlined:`data_exploration;2em`, some more text")
#docfile.write(r""":raw:'<span class="material-symbols-outlined">search</span>'""")
#docfile.write(r""":html:'<span class="material-symbols-outlined">add</span>'""")

docfile.write("\n\n")

docfile.write(".. list-table:: \n")
docfile.write("    :widths: 3 15 10 10 10\n")
docfile.write("    :header-rows: 1\n\n")
docfile.write("    * - Status\n")
docfile.write("      - Name\n")
docfile.write("      - Sections\n")
docfile.write("      - Dimension\n")
docfile.write("      - Validation\n")

if not os.path.isdir("Tests"):
    os.mkdir("Tests")

toctreestr  = ".. toctree::\n"
toctreestr += "   :hidden:\n\n"

for testdirname in sorted(glob.glob("../../tests/*")):
    if not os.path.isdir(testdirname): continue
    
    testname = os.path.basename(testdirname)

    if os.path.isfile(testdirname+"/input.py"):
        docfile.write(f"    * - :fab:`python;sd-text-success fa-fw fa-lg`\n\n")
        docfile.write("      - :ref:`{}`\n".format(testname))
        docfile.write(f"      - \n\n")
        docfile.write(f"      - \n\n")
        docfile.write(f"      - \n\n")
        with open("Tests/{}.rst".format(testname),"w") as testdocfile:
            toctreestr += "   Tests/{}\n".format(testname)
            testdocfile.write(testname + "\n")
            testdocfile.write("="*len(testname) + "\n")

            readmes = find_files_ignore_case(Path(testdirname), "README.rst")
            for readme in readmes: testdocfile.write(f".. include:: ../{testdirname}/{readme.name}\n")
            if readmes: testdocfile.write("\n\n")

            testdocfile.write(".. literalinclude:: ../{}/input.py\n".format(testdirname))
            testdocfile.write("   :caption: Input file ({}/input.py)\n".format(testdirname))
            testdocfile.write("   :language: python\n")
        continue

    if not os.path.isfile(testdirname+"/input"):
        docfile.write(f"    * - :fas:`circle-xmark;sd-text-danger fa-fw fa-lg`\n\n")
        docfile.write(f"      - {testname}\n\n")
        docfile.write(f"      - \n\n")
        docfile.write(f"      - \n\n")
        docfile.write(f"      - \n\n")
        continue




    # Parse the input file ./tests/MyTest/input containing #@ comments.
    # Everything commeneted with #@ will be interpreted as a "config" file
    # The variable "config" is a dict of dicts where each item corresponds to
    # a test configuration.
    cfgfile = io.StringIO()
    input = open(testdirname + "/input")
    for line in input.readlines():
        if line.startswith("#@"):
            cfgfile.write(line.replace("#@",""))
    cfgfile.seek(0)
    config = configparser.ConfigParser(dict_type=MultiOrderedDict,strict=False)
    config.read_file(cfgfile)

    ## this block of code should be the same as in ./scripts/runtests.py
    sections = [s for s in config.sections() if "*" not in s]
    section_wildcards = [s for s in config.sections() if "*" in s]
    for desc in sections:
        desc_wildcards = [s for s in section_wildcards if fnmatch.fnmatch(desc,s)]
        for wc in desc_wildcards:
            new = dict(config[wc])
            newargs = None
            if "args" in config[desc].keys() and "args" in config[wc].keys():
                newargs = new["args"] + "\n" + config[desc]["args"]
            new.update(config[desc])
            if newargs: new["args"] = newargs
            config[desc] = new
    for wc in section_wildcards:
        config.remove_section(wc)
    
    readmes = find_files_ignore_case(Path(testdirname), "README.rst")

    if len(config) <= 1:
        docfile.write("    * - :fas:`triangle-exclamation;sd-text-secondary fa-fw fa-lg`\n")
        if readmes:
            docfile.write("      - :ref:`{}`\n".format(testname))
        else: 
            docfile.write("      - {}\n\n".format(testname))
    else:
        docfile.write("    * - :fas:`circle-check;sd-text-success fa-fw fa-lg`\n")
        docfile.write("      - :ref:`{}`\n".format(testname))
        docfile.write("      - {}\n".format(str(len(config)-1)))
    
        has2D = False
        has3D = False
        for c in config:
            if "dim" in config[c]:
                if config[c]["dim"] == "2": has2D = True
                if config[c]["dim"] == "3": has3D = True
        
        dimstr = ""
        if has2D: dimstr += ":fas:`maximize;fa-fw fa-lg sd-text-secondary` "
        if has3D: dimstr += ":fab:`unity;fa-fw fa-lg sd-text-secondary` "
        docfile.write("      - {}\n".format(dimstr))
        
        if os.path.isfile(testdirname+"/test"):
            docfile.write("      - :fas:`medal;fa-fw fa-lg sd-text-secondary`\n")
        else:
            docfile.write("      - \n")
        docfile.write("\n")

    if len(config) <= 1 and len(readmes) == 0:
        docfile.write("      - \n")
        docfile.write("      - \n")
        docfile.write("      - \n")
        continue
    with open("Tests/{}.rst".format(testname),"w") as testdocfile:
        toctreestr += "   Tests/{}\n".format(testname)

        testdocfile.write(testname + "\n")
        testdocfile.write("="*len(testname) + "\n")

        for readme in readmes:
            testdocfile.write(
                f".. include:: ../{testdirname}/{readme.name}\n"
            )
        if readmes:
            testdocfile.write("\n\n")

        for c in config:
            if c == "DEFAULT": continue
            #testsectionname = "[{}] {}".format(testname,c)
            testsectionname = c
            testdocfile.write(testsectionname+"\n")
            testdocfile.write("-"*len(testsectionname)+"\n")

            testdocfile.write(".. list-table:: \n")
            testdocfile.write("    :widths: 10 90\n")
            testdocfile.write("    :header-rows: 0\n\n")

            #
            # DIMENSION
            #
            if config[c]["dim"] == "2":
                #testdocfile.write("    * - :icon:`2d`\n")
                testdocfile.write("    * - :fas:`maximize;fa-fw fa-lg`\n")
                testdocfile.write("      - Two-dimensional\n")
            else:
                config[c]["dim"] = "3"
                #testdocfile.write("    * - :icon:`3d_rotation`\n")
                testdocfile.write("    * - :fab:`unity;fa-fw fa-lg`\n")
                testdocfile.write("      - Three-dimensional\n")

            #
            # PARALLELISM
            #
            if "nprocs" in config[c] and int(config[c]["nprocs"]) > 1:
                #testdocfile.write("    * - :icon:`grid_view`\n")
                testdocfile.write("    * - :fas:`cubes;fa-fw fa-lg`\n")
                testdocfile.write("      - Parallel ({} procs)\n".format(config[c]["nprocs"]))
            else:
                #testdocfile.write("    * - :icon:`square`\n")
                testdocfile.write("    * - :fas:`cube;fa-fw fa-lg`\n")
                testdocfile.write("      - Serial\n")
            
            #
            # TESTING OR NOT TESTING
            #
            if not os.path.isfile("{}/test".format(testdirname)) or ("check" in config[c] and config[c]["check"] in {"no","No","false","False","0"}):
                #testdocfile.write("    * - :icon:`report_off`\n")
                testdocfile.write("    * - :fas:`question;fa-fw fa-lg`\n")
                testdocfile.write("      - Not validated\n")
            else:
                #testdocfile.write("    * - :icon:`verified`\n")
                testdocfile.write("    * - :fas:`medal;fa-fw fa-lg`\n")
                testdocfile.write("      - Validated using check script\n")

            #
            # BENCHMARK TIME
            #
            if any(["benchmark-" in key for key in config[c]]):
                #testdocfile.write("    * - :icon:`timer`\n")
                testdocfile.write("    * - :fas:`stopwatch;fa-fw fa-lg`\n")
                testdocfile.write("      - ")
                for key in config[c]:
                    if "benchmark-" in key:
                        if "\n" in config[c][key]:
                            raise Exception("Error reading benchmark time for test {} section {}".format(testname,c))
                        testdocfile.write(config[c][key] + "s ({}) ".format(key.replace("benchmark-","")))
                testdocfile.write("\n")

            #testdocfile.write("    * - :icon:`play_circle`\n")
            testdocfile.write("    * - :fas:`circle-play;fa-fw fa-lg`\n")
            cmd = ""
            if "nprocs" in config[c] and int(config[c]["nprocs"]) > 1:
                cmd += "mpiexec -np {} ".format(config[c]["nprocs"])
            exe = "alamo"
            if "exe" in config[c]: exe = config[c]["exe"]
            cmd += "./bin/{}-{}d-g++".format(exe,config[c]["dim"])
            cmd += " {}/input".format(testdirname.replace("../../",""))
            if "args" in config[c]:
                cmd += " "
                cmdargs = [s.replace("= ","=").replace(" =","=") for s in config[c]["args"].split("\n")]
                for s in cmdargs:
                    if len(s.split('=')) == 2:
                        cmd += ' {}="{}"'.format(s.split('=')[0], s.split('=')[1])
            if "ignore" in config[c]:
                cmd += ' ignore="{}"'.format(config[c]["ignore"])
            if "restart" in config[c]:
                cmd += ' restart='+config[c]["restart"]

            testdocfile.write("      - .. code-block:: bash \n\n             {}\n".format(cmd))

            testdocfile.write("\n\n")
        
        #print(os.path.isfile("../../../{}/input".format(testdirname)))
        write_linked_input(
            testdocfile,
            Path(testdirname) / "input",
            "Input file ({}/input)".format(testdirname),
            test_runs(config, Path(testdirname) / "input"),
            "test-input-{}".format(re.sub(r"[^a-z0-9]+", "-", testname.lower())),
        )

        
        
docfile.write("\n\n")
docfile.write(toctreestr)    
docfile.close()
