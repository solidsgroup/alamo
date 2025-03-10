import clang.cindex
from clang.cindex import CursorKind

fqs = []

if False:
    idx = clang.cindex.Index.create()
    tu = idx.parse("../src/alamo.cc", args='-xc++ --std=c++17 -I./src/ -DAMREX_SPACEDIM=2'.split())
    def traverse_ast(node, indent=0):
        try:
            if "../src/Integrator/Base/Mechanics" in str(node.location):
                #if node.kind == clang.cindex.CursorKind.NAMESPACE:
                if "Mechanics" in node.spelling:
                    print('  ' * indent + str(node.spelling), node.kind, node.location)
            for child in node.get_children():
                traverse_ast(child, indent + 1)
        except Exception as e:
            print(e)
    traverse_ast(tu.cursor)
    exit(0)

def ast_classes(path):
    
    idx = clang.cindex.Index.create()
    tu = idx.parse(path, args='-xc++ --std=c++17 -I./src/ -DAMREX_SPACEDIM=2'.split())

    templateclasses = dict()

    def get_fully_qualified_name(node):
        if node is None or node.kind == clang.cindex.CursorKind.TRANSLATION_UNIT:
            return ""
        parent_name = get_fully_qualified_name(node.semantic_parent)
        if parent_name:
            return f"{parent_name}::{node.spelling}"
        return node.spelling

    def extract_template_arguments(node):
        """Recursively find class templates and extract their template parameters."""
        try:
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
        except Exception as e:
            print("EXCEPTION", e)


    extract_template_arguments(tu.cursor)

    return templateclasses
    
if False:
    classes = ast_classes("../src/mechanics.cc")
    for cl in classes:
        print(cl)
        print("  ",classes[cl])
