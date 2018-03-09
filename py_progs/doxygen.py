#!/usr/bin/env python
'''
                    Space Telescope Science Institute

Synopsis:

Replace current informational headers with one that can
be parsed with doxygen


Command line usage (if any):

    usage: doxygen.py whatever.c

Description:

    This is an attempt to write a parser which will read one of the
    .c files used in Python and locate the headers that are supposed
    to be at the top of every function call.  Parsing these headers
    and also using information from cproto the routine attemps to
    create a partially populated doxygen block for each routine.

    The routine writes out a new file, new_whatever.c containging
    doxygen blocks just in fromt of every function call.

    It is up to the user to copy new_whatever.c back to whatever.c once
    he/she has thoroughly looked at the file.

Primary routines:

    doit

Notes:

    The routine has to be run in a directory that contains the header
    files, e.g atomic.h, so that cproto will perform properly.

    The routine will not work with 'non-standard' old headers, so
    one should carefully look at all of the created blocks to see that
    they make sense.

    The old headers are commented out with //OLD and at once we are fairly
    sure they are no longer needed.

    If there was not old_header, there will still be a partially filled out
    doxygen block in front of each routine in the file

    Right now there is a lot of repetive code for populating the various
    blocks.  ksl suspects this should be systemetized


History:

180123 ksl Coding begun
180127 ksl Added the file blocks, which are needed for doxygen to actually
           work

'''

import sys
import subprocess


def is_installed(program):
    """
    Tests to see if a program has been installed. Code taken from:
    https://stackoverflow.com/questions/377017/test-if-executable-exists-in-python

    Args:
        program (str): Program name to test
    Returns:
        Either returns True or raises an OSError
    """
    import os

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return True
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return True

    raise(OSError("Executable '{}' could not be found!".format(program)))


def read_file(filename):
    '''
    Read a file

    Arguments:
        filename (str): Name of the file to read
    Returns:
        list: Lines within the file
    '''

    try:
        f = open(filename, 'r')
        xlines = f.readlines()
        f.close()
    except IOError as e:
        print("The file %s does not exist" % filename)
        return []

    return xlines


def write_header(function):
    """
    Given a function dictionary, writes it to an array of strings

    Arguments:
        function: Dictionary describing the function header

    Returns:
        list (string): The header to be written to file
    """
    header = module_string_start.format(
        function['name'],
        function['synopsis'].replace('\n', '\n * ')
    )

    for name, argument in function['arguments'].items():
        header += module_string_param.format(
            argument['type'], argument['name'],
            argument['description'].replace('\n', '\n * ')
        )

    header += module_string_end.replace(r'%s', r'{}').format(
        function['returns'].replace('\n', '\n * '),
        function['description'].replace('\n', '\n * '),
        function['notes'].replace('\n', '\n * ')
    )
    return header.splitlines(True)


def parse_header_block(header, title, strip=True, default=None):
    """
    Given a list of strings describing a function header, extracts the specified block.

    Does so by iterating over lines in the array, stapling them together into a single string.

    Arguments:
        header (list): List of strings making up the header
        title (string): Block title

    Returns:
        list: List of string making up the block
    """
    titles = ['Arguments:', 'Returns:', 'Notes:', 'Description:',
              'History:', 'Synopsis:', '**********']
    block = None
    for line in header:
        # For each line in the header
        if not block and title in line:
            # We've found the start of the arguments block
            block = line[line.find(title)+len(title):]
        elif block:
            if any(title in line for title in titles):
                # We've reached the start of another block
                if not block.strip():
                    # If this is just a block of whitespace with no text
                    return default
                elif strip:
                    # If we're stripping whitespace out
                    return block.strip()
                else:
                    # If we're not, well, leave it be
                    return block
            else:
                # We're still going through the block, so append this
                block += line
    # We never found the title, so give up and return the default
    return default


def parse_header(header, func_dict):
    """
    Given a list of strings describing a function header, parses them into info.

    Arguments:
        header (list): List of strings making up this header
        func_dict (dict): Dictionary describing this functions assembled by cproto

    Returns:
        dict: Updated func_dict. Could just edit in place, but it's more pythonic?
    """
    # Seek through to find the start of the arguments block
    block = parse_header_block(header, 'Arguments:', strip=False)

    print('Parsing header for {}'.format(func_dict['name']))

    if block:
        # If there is an argument block
        block_test = block.replace('\t', ' ').lower()

        for name, parameter in func_dict['arguments'].items():
            # We want to find this name in the line. BUT! It could be, like 'a'. That's no good.
            # So we need to find a string like ' a; blah' or ' a: blah' or ' a blah'
            # Made a bit more messy given there's some listed as 'a,b,c'

            name_test = ' {};'.format(name).lower()
            name_index = block_test.find(name_test)
            if name_index == -1:
                name_test = ' {}:'.format(name).lower()
                name_index = block_test.find(name_test)
            if name_index == -1:
                name_test = ' {} '.format(name).lower()
                name_index = block_test.find(name_test)
            if name_index == -1:
                name_test = ',{};'.format(name).lower()
                name_index = block_test.find(name_test)
            if name_index == -1:
                name_test = ',{}:'.format(name).lower()
                name_index = block_test.find(name_test)
            if name_index == -1:
                name_test = ',{} '.format(name).lower()
                name_index = block_test.find(name_test)

            if name_index > -1:
                # If we *did* find the substring, add the length of the variable name to it
                name_index += len(name_test)
                comma_index = block.find(',', name_index)
                break_index = block.find('\n', name_index)
                if comma_index > -1 and break_index > -1:
                    if comma_index < break_index:
                        description = block[name_index:comma_index].strip()
                    else:
                        description = block[name_index:break_index].strip()
                elif comma_index > -1:
                    description = block[name_index:comma_index].strip()
                elif break_index > -1:
                    description = block[name_index:break_index].strip()
                else:
                    description = block[name_index:].strip()

                # Did we ever actually find a description?
                if description:
                    parameter['description'] = description

    func_dict['synopsis'] = parse_header_block(header, 'Synopsis:', default='??? SYNOPSIS ???')\
        .replace(func_dict['name'], '').strip()
    func_dict['description'] = parse_header_block(header, 'Description:',
                                                  default='??? DESCRIPTION ???')
    func_dict['returns'] = parse_header_block(header, 'Returns:', default='??? RETURNS ???')
    func_dict['notes'] = parse_header_block(header, 'Notes:', default='??? NOTES ???')
    func_dict['history'] = parse_header_block(header, 'History:', default='??? HISTORY ???')
    return func_dict


def split_definition(definition):
    """
    Taking a definition as 'float* blah', splits it into the name and type

    Args:
        definition (str): String containing definition (function or variable)

    Returns:
        tuple: Name and type of the definition
    """
    last_space = definition.rfind(' ')
    last_pointer = definition.rfind('*')
    # We need to see if this variable is type name or type* name or type *name
    if last_pointer > last_space:
        # If there's a pointer assignment, use that in the type
        definition_name = definition[last_pointer+1:].strip()
        definition_type = definition[:last_pointer+1].strip()
    else:
        definition_name = definition[last_space+1:].strip()
        definition_type = definition[:last_space].strip()
    return (definition_name, definition_type)


def get_modules(filename='emission.c'):
    '''
    use cproto to capture the functions etc that
    are contained in the routine.  Split what is
    returned into something which can be incorporated
    into a search

    Args:
        filename (str): Name of the file to be read

    Returns:
        List of records and dict of records
    '''

    from collections import OrderedDict

    proc = subprocess.Popen('cproto %s ' % filename, shell=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    stdout = stdout.decode().split('\n')

    records_dict = OrderedDict()
    for line in stdout[1:]:
        if line and not ('bad character' in line or 'syntax error' in line):
            # Find the substring of arguments in 'type func_name(type arg1, type arg2);'
            function = {'arguments': OrderedDict()}
            function['name'], function['type'] = split_definition(line[:line.find('(')])
            argument_lines = line[line.find('(')+1:line.find(')')].split(',')
            if 'void' not in argument_lines:
                for argument_line in argument_lines:
                    argument = {'description': '???'}
                    argument['name'], argument['type'] = split_definition(argument_line)
                    function['arguments'][argument['name']] = argument

            records_dict[function['name']] = function

    records = []
    i = 1
    while i < len(stdout)-1:
        z = stdout[i].replace('(', ' ')
        z = z.replace(')', ' ')
        z = z.replace(',', ' ')
        z = z.replace(';', ' ')
        words = z.split()
        # print(words)
        one_record = [stdout[i]] + words
        records.append(one_record)
        # records.append([stdout[i],words[0],words[1]])

        i += 1

    return (records, records_dict)


file_string = '''
/***********************************************************/
/** @file  %s
 * @Author ksl
 * @date   January, 2018
 *
 * @brief  ???
 *
 * ???
 ***********************************************************/
'''

module_string_start = '''
/**********************************************************/
/** @name      {}
 * @brief      {}
 *
 * <NOTE: The [in out] tag describes if the value of a parameter is used or altered. If it is used but not altered, delete 'OUT'. If the original value is not used and it is written to, delete 'IN'.>
'''
module_string_param = ''' * @param [in out] {}  {}   {}
'''
module_string_end = ''' * @return     {}
 *
 * @details
 * {}
 *
 * ### Notes ###
 * {}
 *
 **********************************************************/

'''


def doit(filename='emission.c', outputfile=None):
    '''
    Do something magnificent

    Description:

    Notes:

    For cproto to work properly, files like atomic.h and python.h
    must be in the same directory

    History:


    '''

    if not outputfile:
        outputfile = 'new_'+filename

    lines = read_file(filename)

    if len(lines) == 0:
        raise(EOFError('File {} has no lines!').format(filename))

    modules, mod_dict = get_modules(filename)

    print('These are the functions found in the file\n')
    for one in modules:
        print(one[0])
    print('\n')

    i = 0  # Index for iterating over the array of lines
    j = 0  # Index for the position in the header array
    header_start = []  # Position of the header start lines
    header_end = []    # Position of the header end lines
    header_assigned = []  # Has this header been assigned to a function?
    module_start = []  # Position of the module start lines
    module_end = []    # Position of the module end lines
    module_name = []   # Name of the module for each start-end

    # This section scans through the file, looking for header starts and ends.
    found_header = False
    found_module = False
    while i < len(lines)-1:
        line = lines[i].strip()

        if line[0:24] == '/***********************':
            # print('This is the beginning of a header')
            header_start.append(i)
            found_header = True

        elif found_header:
            if line[0:20] == '********************':
                # print('Found the end of a header')
                header_end.append(i)
                found_header = False

        elif found_module and line:
            # If we've found the start of a module, we want to know where it ends.
            # If this line isn't just whitespace...
            if lines[i][0] == '}':
                # If we find a line that starts with '}', we've found the end (almost certainly)
                module_end.append(i)
                found_module = False

        elif not found_header and not found_module and j < len(modules):
            # To try to match this header to a module, we peek at pairs of lines
            # As most of our functions are defined as 'type \n name (arguments)'
            x = line.strip()
            y = lines[i+1].strip()
            # print('test',x,y)
            x = x.split()
            y = y.split()
            # print('test2',x,'y',y)
            if len(x) > 0 and len(y) > 0 and x[0] == modules[j][1] and y[0] == modules[j][2]:
                # If we've found the first two lines of the module, we append the start line of it
                # and its name to the lists
                module_start.append(i)
                module_name.append(y[0])
                j += 1
                found_module = True

        i += 1

    if len(module_end) != len(module_start):
        # Since the loop only goes to the line *before* the end of the file,
        # if we never caught the end of the last function it must be the last line.
        module_end.append(i)

    print('line where header starts   :', header_start)
    print('line_where header_ends     :', header_end)
    print('lines where function_starts:', module_start)
    print('lines where function_ends:', module_end)

    # Now we need to try to match the current headers with the modules because
    # some may be missing

    # xmatch is the array used to cross-match the headers in the code and their associated modules
    xmatch = [-1] * len(module_start)
    header_assigned = [False] * len(header_start)

    for header_index, header_line in enumerate(header_end):
        # We iterate through the array of headers
        for module_index, module_line in enumerate(module_start):
            # If the start line of this module is after the end of the current header, *and*
            # there's no other module between it and the header.
            if not header_assigned[header_index]:
                # If this header hasn't already been assigned to another module
                if header_line < module_line:
                    # Is it before this module?
                    header_assigned[header_index] = True
                    xmatch[module_index] = header_index
                    break

    print('xmatch betwen headers and functions', xmatch)

    # Now we have the start and end points of each header, and the crossmatch between headers and
    # functions, we can parse the input file.
    for index, module in enumerate(mod_dict.values()):
        # We iterate over the dictionary of modules, looking up what header range corresponds to
        # them using the xmatch array, and passing that to parse_header.
        if xmatch[index] > -1:
            # If there *is* a header for this module, pass it the text
            module = parse_header(lines[header_start[xmatch[index]]:header_end[xmatch[index]]], module)
        else:
            # Else, pass it a blank string and let it set up a default one
            module = parse_header('', module)

    x = open(outputfile, 'w')

    x.write(file_string % (outputfile))

    # Prepend '//OLD' to all the lines between all the identified header start-end pairs
    for start, end in zip(header_start, header_end):
        for i in range(start, end+1):
            lines[i] = '//OLD '+lines[i]

    # Now go through the list and stick our new headers into place
    for index, start in reversed(list(enumerate(module_start))):
        # We go through in reverse order to avoid having to deal with changing line numbers
        lines = lines[:start]+write_header(list(mod_dict.values())[index])+lines[start:]

    x.write(''.join(lines))
    return


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys

    is_installed('cproto')

    if len(sys.argv) == 1 or sys.argv[1] == '-h':
        print(__doc__)
    elif len(sys.argv) > 1:
        doit(sys.argv[1])
    else:
        print('usage: doxygen.py [-h] filename')
