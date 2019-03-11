"""
Script for automatically generating documentation from Python *.c files
"""
# -*- coding: <utf-8> -*-
import yaml
import os
import sys
import typing
from collections import OrderedDict


def parse_param_to_dict(found_parameters: typing.List[dict], text: str, input_file: str, param_type: str) -> bool:
    """
    Given a string that's the text of a c function call for a parameter, decide
    if this is a parameter and return so.

    Arguments:
        found_parameters: Array of found parameter entries so far.
        text: Parameter text e.g. "foo(bar, baz)".
        input_file: File the parameter is in e.g. "setup.c".
        param_type: Parameter type e.g. Integer, Enumerator.

    Returns:
        False: I'm not sure why.
    """
    param_dict = OrderedDict([
        ('name', 'Name'),
        ('description', 'Multi-line description, must keep indentation.\n'),
        ('type', param_type),
        ('unit', 'None'),
        ('values', 'Condition e.g. greater than 0 or list e.g. [1, 2, 5]'),
        ('parent', {"parameter": "Condition e.g. greater than 0 or list e.g. [1, 2, 5]"}),
        ('file', input_file)
    ])

    # If this is a string parameter, it has no units or enumerator values
    if param_type in ["rdstr", "rdline"]:
        param_dict.pop('unit', None)
        param_dict.pop('values', None)

    # If this is an 'advanced' parameter, record so
    if '@' in text[0]:
        param_dict["advanced"] = True
        text = text[1:]

    text_finished = False

    # If this parameter doesn't have a description, we have full name!
    if '(' not in text:
        param_dict["name"] = text.strip()
        text_finished = True
    # If this parameter *does* have a description. record the name
    else:
        if ')' in text:
            # If this is a properly formatted description with a )
            # We want to cut out the whole bracketed section e.g. "foo(bar,baz)_blah" -> "foo_blah"
            description = text[text.find('('):text.rfind(')')+1].strip()
        else:
            # Otherwise, this is a mangled name e.g. "foo(bar,baz" -> "foo"
            description = text[text.find('('):].strip()
        param_dict["name"] = text.replace(description, '').strip()
        text = description.replace('(', '').replace(')', '').strip()

    # Record we've found it
    found_parameters.append(param_dict["name"])

    # If this parameter is called 'exponent', it has no units
    if param_dict["name"].endswith('exponent'):
        param_dict.pop("unit", None)

    # If this parameter isn't already documented, record!
    if param_dict["name"] not in existing_documentation:
        output_dict[param_dict["name"]] = param_dict
    # If this parameter is in the existing documentation, stop and find another
    else:
        return False

    # If this parameter had a name only, no description, stop here
    if text_finished:
        return False

    # Otherwise, we need to parse the remaining descriptive text e.g
    # the (bar,baz) from foo(bar,baz)
    for key, value in bool_types.items():
        # Check to see if this is a de-facto boolean
        if key in text.lower():
            param_dict["type"] = value.strip()
            param_dict.pop("values", None)
            param_dict.pop("unit", None)
            text_finished = True
            break
    if text_finished:
        return False

    # So: Is this string (blahblahblah) or is it (foo,bar,baz)?
    if ',' not in text:
        # If this is not an enumerated list, it's a unit or description.
        for key, value in unit_types.items():
            # Is the description a unit?
            if key in text.lower():
                param_dict["unit"] = value.strip()
                return False

        # If it's not a unit, and not an enumerator, it's just a description
        param_dict["description"] = text.strip()+'\n'

    else:
        # If it is enumerated, it has no units and a list of options
        del param_dict["unit"]
        param_dict['values'] = {
            option: 'Multi-line description, must keep indentation.\n' for option in text.split(',')
        }

    return False


# ==============================================================================
# SETUP
# ==============================================================================

# Set up input and output folders. Requires PYTHON environment variable.
input_folder = os.path.join(os.environ["PYTHON"], "source")
output_folder = os.path.join(os.environ["PYTHON"], "docs", "parameters")
output_old_folder = os.path.join(os.environ["PYTHON"], "docs", "parameters", "old")

# Do not generate documentation for inputs in files containing these substrings
blacklist = ["py_wind", "plot_roche", "py_grid", "synonyms", "t_bilinear", "rdpar"]

# Types of parameter read functions and the type they correspond to.
# Must include 'wrappers' e.g. get_spectype.
param_types = {"get_spectype": "Int",
               "rdint": "Int",
               "rdflo": "Float",
               "rddoub": "Double",
               "rdstr": "String",
               "rdchoice": "Enumerator",
               "rdline": "Line"}

# Types of unit in the description, and what to write in the output.
# The first matching entry will be applied (e.g. msol/yr before msol).
unit_types = OrderedDict([("cm", "cm"),
                          ("ergs/s", "ergs/s"),
                          ("ev", "eV"),
                          ("deg", 'Degrees'),
                          ("msol/yr", u"M☉/year"),
                          ("msol", u"M☉"),
                          ("hr", "Hours"),
                          ("r_star", "co.radius"),
                          ("rstar", "co.radius"),
                          ("wd_rad", "co.radius"),
                          ("r_g", "co.gravitational_radius"),
                          ("cgs", "cgs"),
                          ("vescap", "Escape velocity")])

# Ways boolean variables are shown in the description and their types.
bool_types = {"1=yes": "Boolean (1/0)",
              "y=1": "Boolean (1/0)",
              'anything_else': "Boolean (1/0)",
              'y/n': 'Boolean (Y/N)',
              'yes,no': 'Boolean (yes/no)'}


# Make a list of output files. First, try making output folder
try:
    os.makedirs(os.fsencode(output_folder))
except OSError:
    pass
try:
    os.makedirs(os.fsencode(output_old_folder))
except OSError:
    pass


# ==============================================================================
# FIRST STEP: Read existing documentation and source files.
# ==============================================================================

# Look through the output folder, and record every .yaml file in it.
existing_documentation = []
for filename in os.listdir(output_folder):
    parameter = os.path.splitext(filename)
    if 'yaml' in parameter[1]:
        existing_documentation.append(parameter[0])

# Look through the source folder, and track each file that isn't in the blacklist.
input_files = []
for filename in os.listdir(input_folder):
    # For each file in the input folder
    if filename.endswith(".c") and not any(black in filename for black in blacklist):
        # This is a valid input file
        input_files.append(filename)
    else:
        continue

# existing_documentation = []
# input_files = ['autogen_test.c']

# ==============================================================================
# SECOND STEP: Process the source files and pull the parameters from them.
# ==============================================================================

# For each input file, process and add to found parameters list
# Found parameters is just a list of names. This is for if a parameter is in multiple files.
found_parameters = []
# Output dict is the actual structured dict of documentation dicts to output.
output_dict = OrderedDict()

for input_file in input_files:
    # For each input file with parameters in
    file_object = open(os.path.join(input_folder, input_file), "r")

    # file_text = file_object.read()
    #
    # # Remove all the block comments
    # while "/*" in file_text:
    #     start = file_text.find('/*')
    #     end = file_text.find('*/', start)
    #     comment = file_text[start:end+3]
    #     file_text = file_text.replace(comment, '')
    #
    # # Remove all the in-line comments
    # while '//' in file_text:
    #     start = file_text.find('//')
    #     end = file_text.find('\n', start)
    #     comment = file_text[start:end+2]
    #     file_text = file_text.replace(comment, '')
    #
    # # Remove the linebreaks, then split based on command terminators.
    # file_text = file_text.replace('\n', '')
    # file_lines = file_text.split(';')

    while not file_object.closed:
        # Whilst this file is open, read a line in to find a new parameter.
        line = file_object.readline()
        if not line:
            # If it's "", this is the end of the file. Close it.
            file_object.close()
            break
        elif "//" in line.strip()[:2]:
            # Comment line, ignore
            pass
        elif "/*" in line.strip()[:2]:
            # Multi-line comment, keep reading extra lines.
            # Until we get to the other end of the comment block.
            while "*/" not in line.strip()[-2:]:
                line = file_object.readline()
        else:
            # This is a line of actual c code.
            # It may be split over multiple lines, contain comments, and be messy.
            # We need to clean this line up and send it on to be parsed.

            # First, we check the dict of parameter read functions to see if this line
            # contains one of them.
            for rd_function, param_type in param_types.items():
                if rd_function in line.lower():
                    # We've found a rdpar line with a readable entry,
                    # now let's build the text to parse
                    while ";" not in line:
                        # If this line is spread across multiple lines:
                        # Staple lines together until we have a full "foo (bar \n baz);"
                        line += file_object.readline()

                    # The line may be structured as "foo (bar /* blah */ \n baz);"
                    # Remove any block comments from the line
                    while '/*' in line:
                        comment_start = line.find('/*')
                        comment_finish = line.find('*/')+1
                        comment = line[comment_start:comment_finish]
                        line = line.replace(comment, '')

                    # Remove any single-line comments
                    while '//' in line:
                        comment_start = line.find('//')
                        comment_finish = line.find('\n', comment_start)+1
                        comment = line[comment_start:comment_finish]
                        line = line.replace(comment, '')

                    # Remove any linebreaks
                    line = line.replace('\n', '')

                    if '(' not in line or '"' not in line:
                        # This isn't anything useful
                        continue

                    # Now we've got the full line, extract the parameter text
                    # e.g. rdchoice("var(foo,bar,baz)"", value, choices)
                    # we want to pass "var(foo,bar,baz)" for analysis.
                    start = line.find('"')
                    end = line.find('"', start+1)
                    text = line[start+1:end]
                    parse_param_to_dict(found_parameters=found_parameters, text=text,
                                        input_file=input_file, param_type=param_type)

# Intersections of documentation
deprecated_documentation = sorted(list(set(existing_documentation) - set(found_parameters)),
                                  key=lambda s: s.lower())
new_documentation = sorted(list(set(found_parameters) - set(existing_documentation)),
                           key=lambda s: s.lower())

# If we're not running in write mode
if len(sys.argv) is 1:
    print("Documentation for parameters that no longer exist:")
    for param in deprecated_documentation:
        print(" - {}".format(param))
    print("New, previously undocumented parameters:")
    for param in new_documentation:
            print(" - {}".format(param))
    print("Rerun with:\n" +
          "'-w' to move old parameters to 'parameters/old' and write new files.\n" +
          "'-p' to print new files to screen.")
    sys.exit(1)


# From https://stackoverflow.com/questions/8640959/how-can-i-control-what-scalar-form-pyyaml-uses-for-my-data
# Use literal outputs for multiline
def should_use_block(value):
    """
    Whether or not this string is multi-line and should use the '|' literal string format

    Args:
        value (string): String to test
    Returns:
        bool: Whether or not this is a multi-line string
    """
    for character in u"\u000a\u000d\u001c\u001d\u001e\u0085\u2028\u2029":
        if character in value:
            return True
    return False


def my_represent_scalar(self, tag, value, style=None):
    """
    Replacement 'representer' for strings that intelligently switches between single- and
    multi-line output

    Args:
        tag (string): YAML tag
        value (string): YAML value to assess
    Returns:
        function: The representation function
    """
    if style is None:
        if should_use_block(value):
            style = '|'
        else:
            style = self.default_style

    node = yaml.representer.ScalarNode(tag, value, style=style)
    if self.alias_key is not None:
        self.represented_objects[self.alias_key] = node
    return node


# https://stackoverflow.com/questions/5121931/in-python-how-can-you-load-yaml-mappings-as-ordereddicts
# Dump ordered dicts properly
def dict_representer(dumper, data):
    """
    Function for representing ordered dictionaries like regular ones

    Args:
        dumper (class): Type of object to be dumped
        data(iterator): Output representation of class
    Returns:
        function: Way a dictionary should be represented
    """
    return dumper.represent_dict(data.items())


yaml.representer.BaseRepresenter.represent_scalar = my_represent_scalar
yaml.add_representer(OrderedDict, dict_representer)

# Output to file and print to screen
if '-p' in sys.argv[1:]:
    print("Printing new documentation to screen...\n")
if '-w' in sys.argv[1:]:
    print("Writing out...")
    for filename in deprecated_documentation:
        original_path = os.path.join(output_folder, '{}.yaml'.format(filename))
        moved_path = os.path.join(output_old_folder, '{}.yaml'.format(filename))
        os.rename(original_path, moved_path)


for key, value in output_dict.items():
    print("File: {}/{}.yaml".format(output_folder, key))
    if '-p' in sys.argv[1:] and key in new_documentation:
        print(yaml.dump(value, default_flow_style=False, allow_unicode=True))

    if '-w' in sys.argv[1:] and key in new_documentation:
        with open(os.path.join(output_folder, "{}.yaml".format(key)), "w") as file_object:
            # For each key, open a yaml file and write to it.
            yaml.dump(value, file_object, default_flow_style=False, allow_unicode=True)
