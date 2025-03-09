import clang.cindex
from clang.cindex import CursorKind

#def fully_qualified(c):
#    try:
#        if c is None:
#            return ''
#        elif c.kind == CursorKind.TRANSLATION_UNIT:
#            return ''
#        else:
#            res = fully_qualified(c.semantic_parent)
#            if res != '':
#                return res + '::' + c.spelling
#        return c.spelling
#    except Exception as e:
#        print(e)

idx = clang.cindex.Index.create()
tu = idx.parse('../src/sfi.cc', args='-xc++ --std=c++17 -I./src/'.split())
fqs = []

#def get_calling_function(node):
#    """
#    Given a node with kind MEMBER_REF_EXPR, find and return the cursor for the calling function.
#    """
#    # Traverse upwards to find the calling function
#    parent_node = node
#    while parent_node:
#        # Look for FUNCTION_DECL (regular functions) or CXX_METHOD (methods)
#        if parent_node.kind == clang.cindex.CursorKind.FUNCTION_DECL:
#            print(f"Found function declaration: {parent_node.spelling} at {parent_node.location}")
#            return parent_node
#        elif parent_node.kind == clang.cindex.CursorKind.CXX_METHOD:
#            print(f"Found class method: {parent_node.spelling} at {parent_node.location}")
#            return parent_node
#        
#        # Move to the parent of the current node (this gives the enclosing function or method)
#        parent_node = parent_node.semantic_parent
#    
#    # If we didn't find anything, return None
#    print("No calling function found")
#    return None



def walk(node,stack=[]):
    #print(len(stack))
    if "query" in node.spelling:
        print(node.kind,node.spelling,"({}:{})".format(node.location.file,node.location.line))
        #for thing in reversed(stack):
        #    print(len(stack),thing.kind, thing.spelling, "({}:{})".format(thing.location.file,thing.location.line))
    for child in node.get_children():
        try:
            child.kind
            stack.append(node)
            walk(child,stack)
            stack.pop()
        except ValueError:
            True

walk(tu.cursor)



exit(0)




def getFirst(first=True,verbose=False):
    for c in tu.cursor.walk_preorder():
        try:
            #if c.kind in [
            #        CursorKind.CALL_EXPR,CursorKind.CXX_MEMBER_CALL]:
            
            if "src/IC/BMP.H" in str(c.location.file):
                if not c.referenced: continue
                if c.location.line == 206:
                    print(c.referenced.spelling)
                    print("\t",c.kind)
                    print("\t",c.semantic_parent)
                    print("\t",c.lexical_parent)
                    print("\t",c.parent)
                    for c2 in tu.cursor.walk_preorder():
                        if c in c2.get_children():
                            print("\t\t==>", c2)
                            print("\t\t==>", c2.kind)
                            print("\t\t==>", c2.location)
                            print("\t\t==>", [str(myc2.spelling) for myc2 in c2.get_children()])
                            continue

            continue

            if c.kind not in [CursorKind.OVERLOADED_DECL_REF]:
                continue
            #fq = fully_qualified(c.referenced)
            if not c.referenced:
                continue
            fq = c.referenced.spelling
            #fq = c.referenced.spelling
            #if fq: fqs.append(fq)
            if "ParmParse" in str(c.location.file):
                continue
            if "query_default" in fq.lower():
                get_calling_function(c)
                if verbose:
                    print(fq)
                    print("\tparent:   ",c.semantic_parent)
                    print("\tkind:     ",c.kind)
                    print("\tchildren: ",c.location)
                    for child in c.get_children():
                        print("\t\t",child.kind)
                        print("\t\t",child.spelling)
                if first:
                    return c
        except ValueError as e:
            print(e)

getFirst(False,True)


#exit(0)





## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## import clang.cindex
## 
## # #print(str(clang.cindex.CursorKind))
## # print(clang.cindex.CursorKind.get_all_kinds())
## # CursorKinds = 
## # 
## # # Print all available CursorKinds
## # for kind in clang.cindex.CursorKind.get_all_kinds():
## #     print(f"{kind} : {clang.cindex.CursorKind[kind]}")
## # 
## # exit(0)
## # 
## # Initialize libclang
## clang.cindex.Config.set_library_file('/usr/lib/llvm-18/lib/libclang.so.1')  # Update this with your actual libclang path
## 
## # Function to track function calls and handle recursion
## # Function to track function calls and handle recursion
## class FunctionCallTracker:
##     def __init__(self, target_function):
##         self.target_function = target_function  # The function we're tracking
##         self.call_stack = []  # Stack to track recursion
##         self.call_sites = []  # List to store call sites
## 
##     def visit_node(self, node):
##         print(f"Visiting node: {node.kind} -> {node.spelling if node.spelling else 'No Spelling'}")
##         
##         # Skip unsupported node kinds by checking the valid node types
## #        if node.kind not in [
## #            clang.cindex.CursorKind.CALL_EXPR,
## #            clang.cindex.CursorKind.FUNCTION_DECL,
## #            clang.cindex.CursorKind.VAR_DECL,  # Example of other valid nodes
## #            clang.cindex.CursorKind.TEMPLATE_TYPE_PARAMETER  # Add more if necessary
## #        ]: return
## 
##         # If the node is a function call
##         if node.kind == clang.cindex.CursorKind.CALL_EXPR:
##             print(f"Checking function call: {node.spelling}")  # Debugging output
## 
##             # Check for function name match (accounting for class methods)
##             if node.spelling == self.target_function or \
##                (self.target_function in node.displayname):  # Check for method name with class name
##                 # Check for recursion (if we're currently inside the target function)
##                 if self.target_function in self.call_stack:
##                     print(f"Recursive call to {self.target_function} at {node.location}")
##                 else:
##                     print(f"Found call to {self.target_function} at {node.location}")
##                 self.call_sites.append(node.location)
## 
##         # Recursively visit the children of the node
##         for child in node.get_children():
##             print(f"Visiting child node: {child.kind} -> {child.spelling if child.spelling else 'No Spelling'}")
##             
##             # Check if the child node is a valid CursorKind
##             if child.kind in [clang.cindex.CursorKind.CALL_EXPR,
##                                clang.cindex.CursorKind.FUNCTION_DECL,
##                                clang.cindex.CursorKind.METHOD,
##                                clang.cindex.CursorKind.VAR_DECL]:  # Extend as needed
##                 # Push the current function onto the stack before visiting it
##                 if node.kind == clang.cindex.CursorKind.FUNCTION_DECL or \
##                    node.kind == clang.cindex.CursorKind.METHOD:
##                     self.call_stack.append(node.spelling)
##                 try:
##                     child.kind
##                     self.visit_node(child)
##                 except ValueError as e:
##                     print(e)
##                 # Pop the function from the stack after visiting it
##                 if node.kind == clang.cindex.CursorKind.FUNCTION_DECL or \
##                    node.kind == clang.cindex.CursorKind.METHOD:
##                     self.call_stack.pop()
## 
## # Function to parse the C++ code and analyze the AST
## def analyze_code(file_path, target_function):
##     # Initialize the Clang index
##     index = clang.cindex.Index.create()
## 
##     # Parse the file and get the translation unit (AST)
##     translation_unit = index.parse(file_path)
## 
##     # Create a tracker for the target function
##     tracker = FunctionCallTracker(target_function)
## 
##     # Start traversing the AST from the root
##     tracker.visit_node(translation_unit.cursor)
## 
##     # Print out all recorded call sites
##     print(f"\nTotal occurrences of {target_function}: {len(tracker.call_sites)}")
##     for site in tracker.call_sites:
##         print(f" - {site}")
## 
## # Main execution
## if __name__ == '__main__':
##     file_path = 'src/alamo.cc'  # Path to your C++ file
##     target_function = 'main'  # Name of the function to track
## 
##     analyze_code(file_path, target_function)
## 
