import yaml
import os
import sys


def write_header_by_level(output_file, string, level=0):
    """
    Writes an appropriate header out
    """
    if level is 0:
        output_file.write("{}\n====================\n".format(string))
    elif level is 1:
        output_file.write("{}\n--------------------\n".format(string))
    else:
        output_file.write("{}\n^^^^^^^^^^^^^^^^^^^^\n".format(string))
    return


def write_str_indent(output_file, string, indent="  "):
    """
    Writes a passed string to the output file, with an indent for multiple lines if necessary
    """
    lines = string.splitlines()
    output_file.write("{}\n".format(lines[0]))
    for line in lines[1:]:
        output_file.write("{}{}\n".format(indent, line))
    output_file.write("\n")
    return


def output_parameter(output_file, parameter, level=0):
    """
    Output a parameter to file
    """
    write_header_by_level(output_file, parameter['name'], level)

    type_file.write("{}\n".format(parameter['description']))

    if parameter.get('type'):
        type_file.write("**Type:** {}\n\n".format(parameter['type']))
    if parameter.get('unit'):
        type_file.write("**Unit:** {}\n\n".format(parameter['unit']))

    if parameter.get('values'):
        if isinstance(parameter['values'], dict):
            type_file.write("**Values:**\n\n")
            for key, value in parameter['values'].items():
                if isinstance(value, str):
                    write_str_indent(type_file, "{}. {}".format(key, value.strip()), indent="   ")
                elif isinstance(value, list):
                    list_text = ', '.join([str(x) for x in value])
                    write_str_indent(type_file, "{}. {}".format(key, list_text, indent="   "))
                else:
                    type_file.write("{}\n. {}\n".format(key, value))

        elif isinstance(parameter['values'], list):
            # If this is a list of values, write each out as a bullet-point
            type_file.write("**Values:**\n\n")
            for value in parameter['values'].items():
                write_str_indent(type_file, "* {}".format(value), indent="  ")

        else:
            type_file.write("**Value:** {}\n".format(parameter['values']))
        type_file.write("\n")

    if parameter.get('parent'):
        if isinstance(parameter['parent'], dict):
            type_file.write("**Parent(s):**\n")
            for key, value in parameter['parent'].items():
                if isinstance(value, str):
                    write_str_indent(type_file, "  *{}:* {}".format(key, value), indent="    ")
                elif isinstance(value, list):
                    list_text = ', '.join([str(x) for x in value])
                    write_str_indent(type_file, "  *{}:* {}".format(key, list_text, indent="    "))
                else:
                    type_file.write("  *{}:* {}\n\n".format(key, value))
            type_file.write("\n")

        elif isinstance(parameter['parent'], list):
            # If this is a list of parents, write each out as a bullet-point
            type_file.write("**Parent(s):**\n")
            for value in parameter['parent'].items():
                write_str_indent(type_file, "* {}".format(value), indent="  ")
        else:
            type_file.write("**Parent(s):** {}\n\n".format(parameter['parent']))

    type_file.write("**File:** {}\n\n\n".format(parameter['file']))

    # Go through all the children and output them
    if parameter.get('children'):
        for child_name, child in parameter['children'].items():
            output_parameter(output_file, child, level+1)
    return


input_folder = os.path.join(os.environ["PYTHON"], "docs", "parameters")
output_folder = os.path.join(os.environ["PYTHON"], "docs")
dox_root = {}
dox_all = {}

# First, read the directory to get all the input files
input_files = []
for filename in os.listdir(input_folder):
    # For each file in the input folder
    if filename.endswith(".yaml") and 'reverb' in filename:
        # This is a valid input file
        input_files.append(filename)
    else:
        continue

# Second, read in all parameters from the input files
for input_file in input_files:
    with open(os.path.join(input_folder, input_file), "r") as file_object:
        # For each key, open a yaml file and load from it.
        try:
            parameter = yaml.load(file_object)
            parameter_name = parameter['name']

            # Get the type from the name (e.g. reverb.mode and set no children)
            parameter_type = parameter['name'].split('.')[0]
            parameter['root_type'] = parameter_type.strip()
            parameter['children'] = {}

            # Register this on the full list of all parameters
            dox_all[parameter_name] = parameter

        except yaml.YAMLError as error:
            print(error)

dox_structured = {}

# For each possible root type, set up a root dictionary for it
for parameter_name, parameter in dox_all.items():
    dox_structured[parameter['root_type']] = {}

# We want to build a tree with children etc. to respect parentage
# Now, for the unstructured list of all parameters
for parameter_name, parameter in dox_all.items():
    # For each parameter in the list of all parameters, check its parents
    parameter_type = parameter['root_type']
    if parameter.get('parent'):
        # If it has parents
        found_parent = False
        for parent_name, value in parameter['parent'].items():
            # Check each parent.
            if parent_name in dox_all.keys():
                # If it's parent is another parameter (i.e. in .keys())
                if parameter_type == dox_all[parent_name]['root_type']:
                    # If both the parameter and its parents are in the same class
                    # We assign this parameter to the children's list of its parent
                    dox_all[parent_name]["children"][parameter_name] = parameter
                    # And then 'continue' on, skipping other parents
                    found_parent = True
                    continue

        if not found_parent:
            # If we didn't find the parameter's parent elsewhere,
            # then add it to the root.
            dox_structured[parameter_type][parameter_name] = parameter

    else:
        # If it has no parents, add it to the top layer
        dox_structured[parameter_type][parameter_name] = parameter

# Go through the parameter tree and write it out
for parameter_type, type_parameters in dox_structured.items():
    with open(os.path.join(output_folder, "{}.rst".format(parameter_type)), 'w') as type_file:
        type_file.write("\n====================\n{}\n====================\n\n".format(parameter_type.title()))
        for parameter_name, parameter in type_parameters.items():
            output_parameter(type_file, parameter, level=0)
