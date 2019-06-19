#!/usr/bin/env python
#
# read_frame.py
#
# read from Cytosim's packed binary files
#
# F. Nedelec, 27.08.2017


"""
Syntax:
    
    read_frame.py FILES
    
    
Description:
    
    This script reads Cytosim files.
    
    
F. Nedelec, 27.08.2017
"""

import sys, os, struct


def read_file(path):
    """
        Read binary file
    """
    format = 'iiifff'
    nbytes = struct.calcsize(format)
    with open(path, "rb") as f:
        bytes = f.read(nbytes)
        while bytes:
            data = struct.unpack(format, bytes)
            # do something here with the data:
            # sphere( type, id, color, x, y, z )
            # radius should be .003575 micrometers
            print(data)
            bytes = f.read(nbytes)
        f.close()


def main(args):
    paths = []
    for arg in args:
        if os.path.isfile(arg):
            paths.append(arg)
        else:
            sys.stderr.write("  Error: unexpected argument `%s'\n" % arg)
            sys.exit()
    for p in paths:
        read_file(p)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])


