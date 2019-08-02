"""
This script goes through the created yaml documentation and converts them to
.rst files, then adds the new RST files to the git repository. It then calls
Sphinx to build the docs, in the documentation/html folder, and opens them in
the browser to test.

Arguments:
    Any
        Prints this documentation
"""
#! /usr/bin/env python
from typing import TextIO
import os
import sys
import webbrowser
from subprocess import call
import yaml


def sanitise(text: str) -> str:
    """

    :param text:
    :return:
    """
    while '$' in text:
        equation_start = text.find('$')
        equation_end = text.find('$', equation_start + 1)
        text = text.replace(
            text[equation_start:equation_end + 1],
            image_link_from_latex(text[equation_start + 1:equation_end])
        )

    while ":ref:" in text:
        ref_start = text.find(":ref:")
        name_start = text.find('`', ref_start)
        name_end = text.find('`', name_start+1)

        text = text.replace(
            text[ref_start:name_end+1],
            link_from_name(text[name_start+1:name_end])
        )

    return text


def write_definition_list(output_file: TextIO, definitions: dict):
    """
    Converts a list of key-value pairs to a definition list

    :param output_file:
    :param definitions:
    :return:
    """
    for key, value in definitions.items():
        write_str_indent(
            output_file, "* `{}`: {}".format(key, value),
            indent="  "
        )

def image_link_from_latex(equation: str) -> str:
    """
    Converts a LaTeX equation into an image

    :param equation: The string to convert into an image link
    :return: A markdown image link to codecogs
    """
    return '![{}](https://latex.codecogs.com/gif.latex?{})'.format(equation, equation)


def link_from_name(name: str) -> str:
    """
    Converts a name into a mediawiki link

    :param name: The name to link to
    :return: A mediawiki format link [[display name|link path]]
    """
    return '[[{}|{}#{}]]'.format(
        name,
        name.split('.')[0]+' Parameters',
        name.lower().replace('.', '')
    )


def link_from_file(name: str) -> str:
    """
    Converts a filename to a github link

    :param name:
    :return:
    """
    if isinstance(name, list):
        return ', '.join(['[{}]({})'.format(
            name, 'https://github.com/agnwinds/python/blob/master/source/' + one_name) for one_name in name])
    else:
        return '[{}]({})'.format(
            name, 'https://github.com/agnwinds/python/blob/master/source/'+name
        )


def write_header_by_level(output_file: TextIO, string: str, level: int = 0):
    """
    Writes a passed string as a header, with an appropriate type of underscore depending on
    nesting level (e.g. heading -> subheading)

    :param output_file: The file to write to
    :param string: The text to write
    :param level: The level of heading, 0=heading, 1=subheading, 2=subsubheading, 3=paragraph heading
    """
    output_file.write('#{} {}\n'.format('#'*level, string))


def write_str_indent(output_file: TextIO, string: str, indent: str = "  ", all: bool = False):
    """
    Writes a passed string to the output file. If it's split over multiple lines,
    indent the subsequent lines by the specified amount.

    :param output_file: The file to write to
    :param string: The text to write
    :param indent: The indent to apply if the text splits over multiple lines
    :param all: Whether or not all lines should be indented
    """
    lines = sanitise(string).splitlines()
    if all:
        output_file.write("{}{}\n".format(indent, lines[0]))
    else:
        output_file.write("{}\n".format(lines[0]))

    for line in lines[1:]:
        output_file.write("{}{}\n".format(indent, line))
    output_file.write("\n")
    return


def output_parameter(output_file: TextIO, parameter: dict, level: int = 0):
    """
    Output a parameter to file

    :param output_file: The file to write to
    :param parameter: The parameter to write to file
    :param level: The level of the parameter, e.g. heading, subheading
    """
    write_header_by_level(output_file, parameter['name'], level)
    output_file.write("{}\n".format(sanitise(parameter['description'])))

    if parameter.get('type'):
        output_file.write("##### Type\n {}\n\n".format(sanitise(parameter['type'])))
    if parameter.get('unit'):
        output_file.write("##### Unit\n {}\n\n".format(sanitise(parameter['unit'])))

    if parameter.get('values'):
        output_file.write("##### Values\n\n")
        if isinstance(parameter['values'], dict):
            write_definition_list(output_file, parameter['values'])
        elif isinstance(parameter['values'], list):
            # If this is a list of values, write each out as a bullet-point
            for value in parameter['values'].items():
                write_str_indent(output_file, "* {}".format(value), indent="  ")
        else:
            output_file.write("{}\n".format(sanitise(parameter['values'])))

        output_file.write("\n")

    if parameter.get('parent'):
        if isinstance(parameter['parent'], dict):
            output_file.write("##### Parent(s)\n\n")
            for key, value in parameter['parent'].items():
                if isinstance(value, str):
                    if ' ' not in value:
                        output_file.write("* {}: ``{}``\n".format(link_from_name(key), value))
                    else:
                        write_str_indent(
                            output_file, "* {}: {}".format(link_from_name(key), value),
                            indent="  "
                        )
                elif isinstance(value, list):
                    list_text = ', '.join(['``'+str(x)+'``' for x in value])
                    write_str_indent(output_file, "* {}: {}".format(link_from_name(key), list_text, indent="  "))
                else:
                    output_file.write("* {}: ``{}``\n\n".format(link_from_name(key), value))
            output_file.write("\n")

        elif isinstance(parameter['parent'], list):
            # If this is a list of parents, write each out as a bullet-point
            output_file.write("##### Parent(s)\n")
            for value in parameter['parent'].items():
                write_str_indent(output_file, "* {}".format(value), indent="  ")
        else:
            output_file.write("##### Parent(s) {}\n\n".format(sanitise(parameter['parent'])))

    output_file.write("**File:** {}\n\n\n".format(link_from_file(parameter['file'])))

    # Go through all the children and output them
    if parameter.get('children'):
        for child in parameter['children'].values():
            output_parameter(output_file, child, level+1)
    return


def read_yaml(input_folder: str) -> dict:
    """
    Read all the .yaml files in to dicts. Categorise by type, e.g. 'reverb.mode' is type 'reverb'

    :param input_folder: Path to the folder of input yaml files.
    :return: Dictionary of parameter dicts
    """
    dox_all = {}

    # Read in all parameters from the input files
    for input_file in [filename for filename in os.listdir(input_folder) if filename.endswith(".yaml")]:
        with open(os.path.join(input_folder, input_file), "r") as file_object:
            # Open each yaml file and load from it.
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

    return dox_all


def write_md(output_folder: str, dox_structured: dict):
    """
    Write the structured documentation out to md files

    :param output_folder: The folder to write to
    :param dox_structured: The tree structured documentation
    """
    # Go through the parameter tree and write it out
    for parameter_type, type_parameters in dox_structured.items():
        with open(os.path.join(output_folder, "{} Parameters.md".format(parameter_type)), 'w') as type_file:
            for parameter in type_parameters.values():
                output_parameter(type_file, parameter, level=0)


def autogenerate_wiki_pages():
    """
    Write the wiki files to disk and add them to git.
    """
    output_folder = os.path.join(os.environ["PYTHON"], "wiki", "params")

    # If we can't create the folder, clear out the old autogenerated files
    try:
        for filename in os.listdir(output_folder):
            if " Parameters.md" in filename:
                os.remove(os.path.join(output_folder, filename))
    except:
        os.mkdir(output_folder)

    call(["git", "rm", "-f", os.path.join(output_folder, "* Parameters.md")],
         cwd=os.path.join(os.environ["PYTHON"], "wiki"))

    try:
        os.mkdir(output_folder)
    except:
        pass

    dox_all = read_yaml(
        os.path.join(os.environ["PYTHON"], "docs", "parameters")
    )

    # This is a tree-structured dictionary e.g.
    # reverb: {
    #    children: {
    #        mode: {},
    #        visualisation: {
    #            children: {
    #                angular_bins: {}
    #            }
    #        }
    #    }
    # }
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
            if isinstance(parameter['parent'], dict):
                for parent_name in parameter['parent'].keys():
                    # Check each parent.
                    if parent_name in dox_all.keys():
                        # If it's parent is another parameter (i.e. in .keys())
                        if parameter_type == dox_all[parent_name]['root_type']:
                            # If both the parameter and its parents are in the same class
                            # We assign this parameter to the children's list of its parent
                            dox_all[parent_name]["children"][parameter_name] = parameter
                            # And then 'continue' on, skipping other parents
                            found_parent = True
                            break

            if not found_parent:
                # If we didn't find the parameter's parent elsewhere,
                # then add it to the root.
                dox_structured[parameter_type][parameter_name] = parameter

        else:
            # If it has no parents, add it to the top layer
            dox_structured[parameter_type][parameter_name] = parameter

    write_md(output_folder, dox_structured)

    call(["git", "add", os.path.join(output_folder, "* Parameters.md")],
         cwd=os.path.join(os.environ["PYTHON"], "wiki"))


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    if len(sys.argv) > 1:
        print(__doc__)
    else:
        autogenerate_wiki_pages()
