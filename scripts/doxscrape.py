import os
import xmltodict

# Function to parse a single XML file
def parse_single_doxygen_file(file_path):
    with open(file_path, "r", encoding="utf-8") as file:
        return xmltodict.parse(file.read())

# Function to extract data from all Doxygen XML files
def parse_all_doxygen_files(doxygen_dir):
    # Parse index.xml first
    index_path = os.path.join(doxygen_dir, "index.xml")
    index_data = parse_single_doxygen_file(index_path)

    # Get all referenced XML files from index.xml
    compounds = index_data.get("doxygenindex", {}).get("compound", [])
    if isinstance(compounds, dict):  # Handle case where there's only one compound
        compounds = [compounds]

    # Collect all data
    all_inner_files = []
    all_classes = []
    all_methods = []  # For methods outside of classes

    for compound in compounds:
        refid = compound.get("@refid")
        kind = compound.get("@kind")
        if not refid or not kind:  # Skip invalid entries
            continue

        # Each compound is stored in a separate XML file
        compound_file_path = os.path.join(doxygen_dir, f"{refid}.xml")
        if not os.path.exists(compound_file_path):
            print(f"Warning: File {compound_file_path} not found!")
            continue

        compound_data = parse_single_doxygen_file(compound_file_path)

        # Navigate to <compounddef>
        compounddefs = compound_data.get("doxygen", {}).get("compounddef", [])
        if isinstance(compounddefs, dict):  # Handle case with only one compounddef
            compounddefs = [compounddefs]

        for compounddef in compounddefs:
            # Extract file location from <location> tag in <compounddef>
            compound_location = compounddef.get("location", {}).get("@file", "Unknown")

            # Extract inner files
            inner_files = compounddef.get("innerfile", [])
            if isinstance(inner_files, dict):
                inner_files = [inner_files]
            all_inner_files.extend(
                [{"name": f.get("#text"), "refid": f.get("@refid")} for f in inner_files]
            )

            # Extract class documentation (brief and detailed)
            class_docs = {}
            if compounddef.get("briefdescription", {}):
                class_docs["brief"] = compounddef.get("briefdescription", {}).get("#text", None)
            else:
                class_docs["brief"] = None
            class_docs["detailed"] = None
                        
            #{
            #    "brief":  compounddef.get("briefdescription", {}).get("#text", None),
            #    "detailed": None, #compounddef.get("detaileddescription", {}).get("#text", None),
            #}

            # Store classes with their methods and variables
            if compounddef.get("@kind") in {"class", "struct"}:
                class_data = {
                    "name": compounddef.get("compoundname"),
                    "id": compounddef.get("@id"),
                    "type": compounddef.get("@kind"),
                    "file": compound_location,
                    "brief_description": compounddef["briefdescription"], 
                    "detailed_description": compounddef["detaileddescription"], 
                    "methods": [],
                    "variables": [],
                }

                # Iterate through sectiondefs for methods and variables
                sectiondefs = compounddef.get("sectiondef", [])
                if isinstance(sectiondefs, dict):
                    sectiondefs = [sectiondefs]

                for sectiondef in sectiondefs:
                    # Methods
                    if sectiondef.get("@kind") in {"function", "public-func", "private-func", "protected-func"}:
                        methods = sectiondef.get("memberdef", [])
                        if isinstance(methods, dict):
                            methods = [methods]
                        for method in methods:
                            method_location = method.get("location", {})
                            method_brief_description = None
                            if method.get("briefdescription", {}):
                                method_brief_description = method.get("briefdescription", {}).get("#text", 'none provided'),
                            method_detailed_description = None
                            if method.get("detaileddescription", {}):
                                method_detailed_description = method.get("detaileddescription", {}).get("#text", 'none provided'),
                            #print(method_brief_description,method_detailed_description)
                            method_data = {
                                "name": method.get("name"),
                                "id": method.get("@id"),
                                "type": method.get("type"),
                                "argsstring": method.get("argsstring"),
                                "protection": method.get("prot"),
                                "file": method_location.get("@file", compound_location),
                                "line": method_location.get("@line", "Unknown"),
                                "brief_description": method["briefdescription"],
                                "detailed_description": method["detaileddescription"], 
                            }
                            class_data["methods"].append(method_data)

                    # Variables
                    if sectiondef.get("@kind") in {"variable", "public-attrib", "private-attrib", "protected-attrib"}:
                        variables = sectiondef.get("memberdef", [])
                        if isinstance(variables, dict):
                            variables = [variables]
                        for var in variables:
                            var_location = var.get("location", {})
                            
                            variable_data = {
                                "name": var.get("name"),
                                "id": var.get("@id"),
                                "type": var.get("type"),
                                "protection": var.get("prot"),
                                "definition": var["definition"],
                                "file": var_location.get("@file", compound_location),
                                "line": var_location.get("@line", "Unknown"),
                                "brief_description": var["briefdescription"],
                                "detailed_description": var["detaileddescription"],
                            }
                            class_data["variables"].append(variable_data)

                all_classes.append(class_data)

            # Methods outside of classes (global functions)
            sectiondefs = compounddef.get("sectiondef", [])
            if isinstance(sectiondefs, dict):
                sectiondefs = [sectiondefs]

            for sectiondef in sectiondefs:
                if sectiondef.get("@kind") == "function":
                    methods = sectiondef.get("memberdef", [])
                    if isinstance(methods, dict):
                        methods = [methods]
                    for method in methods:
                        method_location = method.get("location", {})
                        method_data = {
                            "name": method.get("name"),
                            "id": method.get("@id"),
                            "type": method.get("type"),
                            "argsstring": method.get("argsstring"),
                            "protection": method.get("prot"),
                            "file": method_location.get("@file", compound_location),
                            "line": method_location.get("@line", "Unknown"),
                            "brief_description": method["briefdescription"],
                            "detailed_description": method["detaileddescription"],
                        }
                        all_methods.append(method_data)

    return {"files": all_inner_files, "classes": all_classes, "methods": all_methods}

## Directory containing the Doxygen XML files
#doxygen_dir = doxygen_xml #"path/to/doxygen/output/xml"
#
#myfile = "src/Integrator/Integrator.H"
#
## Parse all files and print results
#doxygen_data = parse_all_doxygen_files(doxygen_dir)
#
#
#
#
##print("Inner Files:")
##for file in doxygen_data["files"]:
##    print(f"- File Name: {file['name']}, RefID: {file['refid']}")
#
#print("\nClasses:")
#for cls in doxygen_data["classes"]:
#    if cls["file"]  != myfile: continue
#
#    print(f"- Class Name: {cls['name']}, ID: {cls['id']}, Type: {cls['type']}, File: {cls['file']}")
#    print(f"  Brief Description: {cls['brief_description']}")
#    print(f"  Detailed Description: {cls['detailed_description']}")
#
#    
#    # Methods
#    print("  Methods:")
#    for method in cls["methods"]:
#        #print(f"    - Method Name: {method['name']}, ID: {method['id']}, Return Type: {method['type']}, "
#        #      f"Args: {method['argsstring']}, Protection: {method['protection']}, File: {method['file']}, "
#        #      f"Line: {method['line']}")
#        print(f"      Name: {method['name']}")
#        if method['brief_description']:
#            print(f"      --> Brief Description: {method['brief_description']['para']}")
#        if method['detailed_description']:
#            print(f"      --> Detailed Description: {method['detailed_description']['para']}")
#
#        #print(f"      Detailed Description: {method['detailed_description']}")
#    
##    # Variables
##    print("  Variables:")
##    for var in cls["variables"]:
##        print(f"    - Variable Name: {var['name']}, ID: {var['id']}, Type: {var['type']}, "
##              f"Protection: {var['protection']}, File: {var['file']}, Line: {var['line']}")
##        print(f"      Brief Description: {var['brief_description']}")
##        print(f"      Detailed Description: {var['detailed_description']}")
#
### Global Methods (outside classes)
##print("\nGlobal Methods (Outside Classes):")
##for method in doxygen_data["methods"]:
##
##    if method["file"]  != myfile: continue
##
##    print(f"- Method Name: {method['name']}, ID: {method['id']}, Return Type: {method['type']}, "
##          f"Args: {method['argsstring']}, Protection: {method['protection']}, File: {method['file']}, "
##          f"Line: {method['line']}")
##    print(f"  Brief Description: {method['brief_description']}")
##    print(f"  Detailed Description: {method['detailed_description']}")


