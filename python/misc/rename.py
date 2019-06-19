#!/usr/bin/env python
#
# copyright F. Nedelec, 2012--2014

"""
Synopsis:

     Rename files or folders using consecutive numbers
    
Syntax:

     rename.py PATTERN FILE_NAMES
      
Example:

     rename.py image%04i.png image*.png
     
     will rename images as: image0000.png, image0001.png, etc.

F. Nedelec, 2012--2014.
"""

try:
    import sys, os
except ImportError:
    print("  Error: could not load necessary python modules\n")
    sys.exit()

#------------------------------------------------------------------------

def rename(files, pattern):
    """rename files using consecutive numbers"""
    res = []
    cnt = 0
    for file in files:
        while cnt < 10000:
            name = pattern % cnt
            cnt += 1
            if name == file:
                res.append(file)
                break
            if not os.path.exists(name):
                os.rename(file, name)
                res.append(name)
                print("%s -> %s" % ( file, name ))
                break
    return res

#------------------------------------------------------------------------


def main(args):
    """rename files"""
    files = []
    pattern = args.pop(0);
    
    for arg in args:
        if os.path.isfile(arg) or os.path.isdir(arg):
            files.append(arg)
        else:
            print("ignored '%s' on command line" % arg)
        
    try:
        rename(files, pattern)
    except IOError as e:
        print("Error: "+repr(e))


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])


