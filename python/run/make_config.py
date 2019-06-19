#!/usr/bin/env python
#
# A configuration file generator
#
# Copyright Francois J. Nedelec, Porquerolles October 2015

"""
    Generates a set of config files from one input file,
    by varying some of the parameters.
    
Syntax:
   
    make_config.py [REPEAT] INPUT_FILE

Description:
    
    Generate files 'config????.cym' by varying parameters.
    The parameter variations are defined within this python script, and the
    template file should contain corresponding placeholders identified with 
    a dollar sign, for example:
    
    set simul system
    {
      time_step = $time_step
      viscosity = 0.02
    }

    F. Nedelec, October 2015

"""

import sys, os, string

# output stream for error messages:
out = sys.stderr

# elements of the file name
file_parts = ['config', '%04i', '.cym']

# index of output file
index = 0

def next_file_name():
    """
        Generate the name for the next output file
    """
    global index, file_parts
    n = file_parts[0] + ( file_parts[1] % index ) + file_parts[2]
    index += 1
    return n

def write_file(name, content):
    """
        Write content to a new file, called `name`
    """
    file = open(name, 'w')
    file.write(content)
    file.close()

#------------------------------------------------------------------------

def process(name, repeat=1):
    """
        process one file, repeating every output `repeat` times
    """
    # read template file:
    input = open(name, 'r').read()
    template = string.Template(input)
    # vary parameters:
    for a in range(10):
        # set value of parameter
        time_step = 2 ** -a
        # substitute parameters by value-strings:
        output = template.substitute(time_step=str(time_step))
        for r in range(repeat):
            write_file(next_file_name(), output)

#------------------------------------------------------------------------

def main(args):
    global index;
    input = []
    repeat = 1
 
    for arg in args:
        if os.path.isfile(arg):
            input = arg
        elif arg.isdigit():
            repeat = int(arg)
        else:
            out.write("  Error: unexpected argument `%s'\n" % arg)
            sys.exit()

    if not input:
        out.write("  Error: you must specify an input file\n")
        sys.exit()

    process(input, repeat)
    out.write("  %i files generated from %s\n" % (index, input))
 
#------------------------------------------------------------------------

if len(sys.argv) < 2 or sys.argv[1].endswith("help"):
    print(__doc__)
else:
    main(sys.argv[1:])

