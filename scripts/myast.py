import clang.cindex
from clang.cindex import CursorKind
import pathlib
from pathlib import Path
from glob import glob


def scan(root,path):
    absroot = Path(root).resolve()

    args = f'-xc++ --std=c++17 -I../src/ -DAMREX_SPACEDIM=2 -fimplicit-modules -fimplicit-module-maps'.split()

    idx = clang.cindex.Index.create()
    tu = idx.parse(root+"/"+path, args=args)

    templateclasses = dict()

    def get_fully_qualified_name(node):
        if node is None or node.kind == clang.cindex.CursorKind.TRANSLATION_UNIT:
            return ""
        parent_name = get_fully_qualified_name(node.semantic_parent)
        if parent_name:
            return f"{parent_name}::{node.spelling}"
        return node.spelling

    def extract_template_arguments(node):
        if not node.kind == clang.cindex.CursorKind.TRANSLATION_UNIT:
            if not node.location.file:
                return

            absfilename = str(Path(node.location.file.name).resolve())

            if str(absroot) not in absfilename:
                return

        if "../src" in str(node.location):
            if node.kind in [clang.cindex.CursorKind.CLASS_TEMPLATE,
                             clang.cindex.CursorKind.CLASS_DECL]:
                fq_name = get_fully_qualified_name(node)

                templateclasses[fq_name]={}
                templateclasses[fq_name]['templates'] = []
                templateclasses[fq_name]['parents'] = []

                ## Iterate over template parameters
                for child in node.get_children():
                    if child.kind == clang.cindex.CursorKind.CXX_BASE_SPECIFIER:
                        templateclasses[fq_name]['parents'].append(get_fully_qualified_name(child))
                    if child.kind == clang.cindex.CursorKind.TEMPLATE_TYPE_PARAMETER:
                        templateclasses[fq_name]['templates'].append(child.spelling)

        # Recurse through all nodes
        for child in node.get_children():
            extract_template_arguments(child)

    extract_template_arguments(tu.cursor)
    return templateclasses
    
