#!/usr/bin/env python
#
# PRECONFIG, a versatile configuration file generator
#
# Copyright Francois J. Nedelec, EMBL 2010--2018
# This is PRECONFIG version 1.1, last modified on 22.3.2019

"""
# SYNOPSIS

   Preconfig generates files from a template by evaluating doubly-bracketed Python code.
   
# DESCRIPTION
   
   Preconfig reads the template file from top to bottom, identifying snippets
   of code which are surrounded by double square brackets. It then executes this
   code using the python interpreter, proceeding recursively whenever multiple
   values are specified. Values are eventually converted to their string
   representation, and substituted in place of the code snippet. In this way,
   Preconfig will generate all the possible combinations following the order in
   which these combinations were specified in the file. Importantly, any ac-
   -companying text in the template file is copied verbatim to the output file,
   such that any syntax present in the configuration file can be maintained
   during the process.

   At least one template file should be specified, and other arguments are optional.
   If several template files are specified, they will be processed sequentially.
   The names of the produced files are built from the name of the template
   by removing any second extension, and inserting an integer of constant width.
   
   Examples:
   
   - config.cym.tpl --> config0000.cym, config0001.cym, config0002.cym, etc.
   - config.txt.tpl --> config0000.txt, config0001.txt, config0002.txt, etc.
   - model.xml.tpl --> model0000.xml, model0001.xml, model0002.xml, etc.
   
   The width of the variable part (default=4) can be changed on the command line.
   For instance, to specify a width of 2 characters, invoke "preconfig -2 ...".
   
# SYNTAX
   
   preconfig [OPTIONS] TEMPLATE_FILE [ADDITIONAL_TEMPLATE_FILES]

# OPTIONS
   
   - if a positive integer REPEAT is specified, each template file will be
   processed REPEAT times, for example: `preconfig 3 config.cym.tpl` will parse
   the template three times and generate three times more files.
   
   - if the path to an existing directory is specified, files will be created
   in this directory, for example: `preconfig dir config.cym.tpl`
   
   - DEFINITIONS can be specified on the command line as 'name=value' or 
   'name=sequence', with no space around the '='. They are added to the 
   dictionary used to evaluate the code snippets found inside the template file,
   for example: `preconfig n_molecules=100 config.cym.tpl`
   
   - if a negative integer is specified, this will set the width of the integer
   that is used to build the file namess.
   For example: `preconfig -2 config.cym.tpl` will create 'config00.cym', etc.
   
   - if a '-' is specified, all accessory output is suppressed
   
   - if a '+' is specified, more detailed information on the parsing is provided.
   
   - if '++' or 'log' is specified, a file 'log.csv' will be created containing one
   line for each file created, containing the substitutions operated for this file.
   
   - if '--help' is specified, this documentation will be printed.
   
# CODE SNIPPETS
   
   Any plain python code can be embedded in the file, and functions from the
   [Random Module](https://docs.python.org/library/random.html) can be used.
   It is possible to use multiple bracketed expressions in the same file, and
   to define variables in the python environment. An integer 'n', starting at
   zero and corresponding to the file being generated is automatically defined.
   
## Example 1
   
   Generate all combinations with multiple values for 2 parameters:
   
    rate = [[ [1, 10, 100] ]]
    speed = [[ [-1, 0, 1] ]]
   
   Command: `preconfig TEMPLATE_FILE` with the appropriate file name.
   In this case, Preconfig will generate 9 files.

## Example 2
   
   Regularly scan 2 parameters with 10 values each,
   one according to a linear scale, and the other with a geometric scale:
   
    [[ x = range(10) ]]
    [[ y = range(10) ]]
    reaction_rate = [[ 1 + 0.5 * x ]]
    diffusion_rate = [[ 1.0 / 2**y ]]
   
   Command: `preconfig TEMPLATE_FILE`
   In this case, Preconfig will generate 100 files, one for each combination.
   
## Example 3

   Scan multiple parameters values randomly:
   
    diffusion_rate = [[ random.uniform(0,1) ]]
    binding_rate = [[ round(random.uniform(0,1), 3) ]]
    reaction_rate = [[ random.choice([1, 10, 100]) ]]
    abundance = [[ random.randint(0, 1000) ]]
   
   Command: `preconfig 256 TEMPLATE_FILE`
   In this case, Preconfig is instructed to generate 256 files.
   
## Example 4
   
   Randomize two parameters while keeping their ratio constant:
   
    [[ x = random.uniform(0,1) ]]
    binding_rate = [[ 10.0 * x ]]
    unbinding_rate = [[ x ]]
   
   Command: `preconfig 256 TEMPLATE_FILE` to make 256 files.

## Example 5
   
   Randomize one parameter, using 256 values in ascending order:
   
    [[ x = sorted([random.uniform(0.10, 0.25) for i in range(256)]) ]]
    binding_rate = [[ x ]]
   
   Command: `preconfig TEMPLATE_FILE`
   In this case, the number of files (256) is specified in the template

## Example 6
   
   Boolean variables can be used to introduce qualitative differences:

    [[ enable = random.choice([0, 1]) ]]
    feeback = [[ random.uniform(0, 1) if (enable) else 0  ]]
   
   Command: `preconfig 256 TEMPLATE_FILE` to make 256 files.
   
## Example 7
   
   Randomize a value, and print this value as a comment in the file.
   The second line below the [[...]] prints a comment, from which the value
   of 'x' can be read. This can be useful to process the results later.
   
    [[ x = random.uniform(0,1) ]]
    % preconfig.x = [[ x ]];
    binding_rate = [[ 10*x ]]
    unbinding_rate = [[ 2*x ]]

   Command: `preconfig 256 TEMPLATE_FILE` to make 256 files.

## Acknowledgments:

We wish to thank the members of the Nedelec group, and all users of 
Cytosim for their feedback which has contributed greatly to this development.
We thanks Shaun Jackman and Steven Andrews for valuable feedback!

Copyright Francois J. Nedelec, 2010--2018
This is Free Software with absolutely no WARANTY.
Preconfig is distributed under GPL3.0 Licence (see LICENCE)
"""


import sys

try:
    import os, io
    globals = { 'random': __import__('random'), 'math': __import__('math') }
except ImportError:
    sys.stderr.write("Error: Preconfig could not load necessary python modules\n")
    sys.exit()

#-------------------------------------------------------------------------------

# code snippets are surrounded by double square brackets:
SNIPPET_OPEN = '['
SNIPPET_CLOSE = ']'

# streams for output (all output is hidden by default):
out = open(os.devnull, 'w')
log = []

# motif used to compose file names
pattern = 'config%04i.txt'

# number of digits used to compose `pattern`
nb_digits = 4

# index of file being generated
file_index = 0

# list of files generated
files_made = []

# name of current input file being processed (used for error reporting)
template = ''

#-------------------------------------------------------------------------------

def set_pattern(name, destination):
    """
    Extract the root and the extension of the file
    """
    global pattern, nb_digit, file_index
    [main, ext] = os.path.splitext(os.path.basename(name))
    if '.' in main:
        [main, ext] = os.path.splitext(main)
    pattern = main + '%0' + str(nb_digits) + 'i' + ext
    if destination:
        pattern = os.path.join(destination, pattern)
    file_index = 0


def next_file_name():
    """
    Generate the name of the next output file
    """
    global file_index, pattern
    n = pattern % file_index
    file_index += 1
    return n


def make_file(text, values):
    """
    Create a file with the specified text
    """
    global files_made, file_index
    name = next_file_name()
    with open(name, 'w') as f:
        f.write(text)
        files_made.extend([name])
    # fancy ouput:
    out.write("\\"+repr(values)+'\n')
    out.write(" \\"+('> '+name).rjust(78, '-')+'\n')
    # write log:
    if log:
        keys = sorted(values.keys())
        if file_index == 1:
            log.write('%20s' % 'file')
            for k in keys:
                log.write(', %10s' % k)
            log.write('\n')
        log.write('%20s' % name)
        for k in keys:
            log.write(', %10s' % repr(values[k]))
        log.write('\n')
    values['n'] = file_index


#-------------------------------------------------------------------------------

def pop_sequence(dic):
    """
        Remove an entry in the dictionary that has multiple values
    """
    for k in dic:
        v = dic[k]
        try:
            len(v)
            if not isinstance(v, str):
                dic.pop(k)
                return (k, v)
        except:
            pass
    return ('', [])


def try_assignment(cmd):
    """
        Check if `cmd` follows the format of a variable assignent,
        and if that succeeds, return the key and value strings in a tuple.
    """
    try:
        k, v = cmd.split("=")
        if k and v:
            k = k.strip()
            v = v.strip()
            return (k, v)
    except ValueError:
        pass
    return ('', cmd)


def evaluate(cmd, locals):
    """
        Evaluate `cmd` and return the result
    """
    res = cmd
    try:
        res = eval(cmd, globals, locals)
    except Exception as e:
        sys.stderr.write("\033[95m")
        sys.stderr.write("Error in `%s`:\n" % template)
        sys.stderr.write("\033[0m")
        sys.stderr.write("    Could not evaluate [[%s]]\n" % cmd)
        sys.stderr.write("    "+str(e)+'\n')
        sys.exit(1)
    try:
        res = list(res)
    except Exception:
        pass
    return res


def get_block(file, s, e):
    """
    Extract the next block starting with DOUBLE delimiters 'ss' and ending with 'ee'
    Returns a set with 3 values:
        - the text found before the block start
        - the block with its delimiters
        - a boolean EOF indicator
    """
    ch = file.read(1)
    pre = ''
    blk = ''
    dep = 0
    while file and ch:
        pc = ch
        ch = file.read(1)
        if ch == s:
            if pc == s and dep == 0:
                dep = 1
            if dep > 0:
                dep += 1
        if dep > 0:
            blk += pc
        else:
            pre += pc
        #print("%c%c dep %i" %(pc, ch, dep))
        if ch == e:
            if pc == e and dep == 1:
                return (pre, blk+ch, False)
            if dep > 0:
                dep -= 1
    if blk:
        out.write("Error: unclosed bracketted block in:%s\n" % pre);
    return (pre, '', True)


#-------------------------------------------------------------------------------

def process(file, locals, text):
    """
        `process()` will identify and substitute bracketed code blocks
        embedded in the input file, and generate a file at EOF.
    """
    output = text

    while file:
        (pre, code, eof) = get_block(file, SNIPPET_OPEN, SNIPPET_CLOSE)
        #print("text `", pre[0:32], "' of size ", len(pre))
        #print("code [["+pre+"]] EOF=%i" % eof)
        output += pre
        if eof:
            # having exhausted the input, we generate a file:
            make_file(output, locals)
            return
        # remove outer brackets:
        cmd = code[2:-2]
        #print("embedded code '%s'" % code)
        # interpret command:
        (key, vals) = try_assignment(cmd)
        vals = evaluate(vals, locals)
        #print("`"+key+"' is "+repr(vals))
        try:
            # use 'pop()' to test if multiple values were specified...
            # keep last value aside for later:
            val = vals.pop()
            ipos = file.tell()
            for v in vals:
                # fork recursively for all subsequent values:
                #print("forking", v)
                if key:
                    locals[key] = v
                    out.write("|%50s <-- %s\n" % (key, str(v)) )
                    process(file, locals, output)
                else:
                    out.write("|%50s --> %s\n" % (code, str(v)) )
                    process(file, locals, output+str(v))
                file.seek(ipos)
        except (AttributeError, IndexError):
            # a single value was specified:
            val = vals
        # handle remaining value:
        # print("handling", key, val)
        if key:
            locals[key] = val
            out.write("|%50s <-- %s\n" % (key, str(val)) )
        else:
            output += str(val)
            out.write("|%50s --> %s\n" % (code, str(val)) )


def expand_values(file, values, text):
    """
        Call self recursively to remove all entries of the 
        dictionary 'values' that are associated with multiple keys.
    """
    (key, vals) = pop_sequence(values)
    if key:
        ipos = file.tell()
        for v in vals:
            values[key] = v
            #out.write("|%50s <-- %s\n" % (key, str(v)) )
            expand_values(file, values, text)
            file.seek(ipos)
        # restore multiple values on upward recursion
        values[key] = vals
    else:
        process(file, values, text)


def parse(name, values={}, repeat=1, destination=''):
    """
        process one file, and return the list of files generated
    """
    values['n'] = 0
    global files_made, template
    template = name
    set_pattern(name, destination)
    files_made = []
    for x in range(repeat):
        with open(name, 'r') as f:
            expand_values(f, values, '')
    return files_made


#-------------------------------------------------------------------------------


def main(args):
    global out, log, nb_digits
    inputs = []
    locals = {}
    repeat = 1
    verbose = 1
    destination = ''
    
    for arg in args:
        if os.path.isdir(arg):
            destination = arg
        elif os.path.isfile(arg):
            inputs.append(arg)
        elif arg.isdigit():
            repeat = int(arg)
        elif arg == '-':
            verbose = 0
        elif arg == '+':
            out = sys.stderr
            verbose = 0
        elif arg == '++' or arg == 'log':
            log = open('log.csv', 'w')
        elif arg[0] == '-' and arg[1:].isdigit():
            nb_digits = int(arg[1:])
        else:
            (k,v) = try_assignment(arg)
            if k:
                locals[k] = evaluate(v, locals)
            else:
                sys.stderr.write("  Error: unexpected argument `%s'\n" % arg)
                sys.exit()
    
    if not inputs:
        sys.stderr.write("  Error: you must specify an input template file\n")
        sys.exit()

    for i in inputs:
        #out.write("Reading %s\n" % i)
        res = parse(i, locals, repeat, destination)
        if verbose == 1:
            #print("%i files generated from %s:" % (len(res), i))
            for f in res:
                print(f)


#-------------------------------------------------------------------------------


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("You must specify a template file (for instructions, invoke with option '--help')")
    elif sys.argv[1].endswith("help"):
        print(__doc__)
    elif sys.argv[1]=='--version':
        print("This is PRECONFIG version 1.1 (22.3.2019)")
    else:
        main(sys.argv[1:])


