import yaml
import os
import sys
from collections import OrderedDict


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


def output_parameter(output_file, parameter):
    return


input_folder = os.path.join(os.environ["PYTHON"], "docs", "parameters")
output_folder = os.path.join(os.environ["PYTHON"], "docs")
dox_root = {}

input_files = []
for filename in os.listdir(input_folder):
    # For each file in the input folder
    if filename.endswith(".yaml"):
        # This is a valid input file
        input_files.append(filename)
    else:
        continue

for input_file in input_files:
    with open(os.path.join(input_folder, input_file), "r") as file_object:
        # For each key, open a yaml file and load from it.
        try:
            parameter = yaml.load(file_object)
            parameter_name = parameter['name']
            parameter_type = parameter['name'].split('.')[0]

            if parameter_type not in dox_root:
                dox_root[parameter_type] = {}
            dox_root[parameter_type][parameter_name] = parameter

        except yaml.YAMLError as error:
            print(error)

dox_unstructured = dox_root.copy()
dox_structured = {}
for key in dox_unstructured.keys():
    dox_structured[key] = {}

for parameter_type, parameters in dox_unstructured.items():
    for parameter_name in list(parameters.keys()):
        parameter = parameters[parameter_name]
        if parameter.get('parent'):
            for parent_name, value2 in parameter['parent'].items():
                # If this parameter has a parent parameter that is *also* in this file
                parent_type = parent_name.split('.')[0]
                if parent_type in parameter_type:
                    continue
                else:
                    dox_structured[parameter_type][parent_name] = parameter
                    del dox_unstructured[parameter_type][parameter_name]

for parameter_type, type_parameters in dox_root.items():
    with open(os.path.join(output_folder, "{}.rst".format(parameter_type)), 'w') as type_file:
        type_file.write("\n=====\n{}\n=====\n\n".format(parameter_type.title()))
        for parameter_name, parameter in type_parameters.items():
            type_file.write("{}\n==============================\n\n".format(parameter['name']))
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
                            write_str_indent(type_file, "{}. {}".format(key, value), indent="   ")
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
                            type_file.write("  *{}:* {}\n".format(key, value))

                elif isinstance(parameter['parent'], list):
                    # If this is a list of parents, write each out as a bullet-point
                    type_file.write("**Parent(s):**\n")
                    for value in parameter['parent'].items():
                        write_str_indent(type_file, "* {}".format(value), indent="  ")
                else:
                    type_file.write("**Parent(s):** {}\n".format(parameter['parent']))
                type_file.write("\n")

            type_file.write("**File:** {}\n\n".format(parameter.get('file')))
