#!/usr/bin/env python
#
# A configuration file generator
#
# Copyright Francois J. Nedelec, 2013
#

"""
    Generates a set of config files from one input file,
    by applying all patches found in the current directory.
    
Syntax:
   
    patch_config.py [-NB_DIGITS] [REPEAT] [DESTINATION] TEMPLATE_FILE [MORE_TEMPLATE_FILES]

Description:
    
    To generate its file, this finds all `*.patch` in the current directory,
    and apply them to the template_file by invoking the unix command 'patch'
    
    F. Nedelec, July 2013

"""

import sys, shutil, subprocess

# output stream for error messages:
out = sys.stderr

try:
    import os
except ImportError:
    out.write("Error: could not load necessary python modules\n")
    sys.exit()


#------------------------------------------------------------------------

index = 0
file_parts = ['config', '%04i', '.txt']
files_made = []


def set_file_name(name, dest):
    """
        extract the root and the extension of the file
    """
    global file_parts, files_made
    files_made = []
    file = os.path.basename(name)
    [main, ext] = os.path.splitext(file)
    if '.' in main:
        [main, ext] = os.path.splitext(main)
    if dest:
        main = os.path.join(dest, main)
    file_parts[0] = main;
    file_parts[2] = ext;


def next_file_name():
    """
        Generate the name for the next output
    """
    global index, file_parts
    n = file_parts[0] + ( file_parts[1] % index ) + file_parts[2]
    index += 1
    return n

#------------------------------------------------------------------------

def process(iname, patch, dest):
    oname = next_file_name()
    code = subprocess.call(['patch', '-o', oname, iname, patch])
    #shutil.copyfile(iname, oname)
    return oname


def parse(iname, values={}, repeat=1, dest=''):
    """
        process one file, and return the list of files generated
    """
    res = []
    import glob
    patches = glob.glob('*.patch')
    set_file_name(iname, dest)
    for x in range(repeat):
        for patch in patches:
            s = process(iname, patch, dest)
            res.append(s)
    return res


#------------------------------------------------------------------------


def main(args):
    inputs = []
    values = {}
    repeat = 1
    dest = ''

    for arg in args:
        if os.path.isdir(arg):
            dest = arg
        elif os.path.isfile(arg):
            inputs.append(arg)
        elif arg.isdigit():
            repeat = int(arg)
        elif arg[0] == '-' and arg[1:].isdigit():
            global file_parts
            nd = int(arg[1:])
            file_parts[1] = '%%%03ii' % nd
        else:
            out.write("  Error: unexpected argument `%s'\n" % arg)
            sys.exit()
    
    if not inputs:
        out.write("  Error: you must specify an input template file\n")
        sys.exit()

    for i in inputs:
        out.write("Reading %s\n" % i)
        res = parse(i, values, repeat, dest)
        out.write("%i files generated from %s\n" % (len(res), i))
        #out.write(repr(res))



#------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])


