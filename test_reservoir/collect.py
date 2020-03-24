#!/usr/bin/env python
#
# collect.py
#
# Copyright F. Nedelec, 2007--2017


"""
Synopsis:
    
    Rename files or folders using a pattern containing an index as in 'image0001.png'.
    This number in the file name will be incremented automatically for each file,
    and no existing file will not overwritten, such that this can be used to pool
    together many similar files in a common directory.
    
Syntax:
    
    collect.py PATTERN [INTEGER] [--copy] PATHS

Arguments:
    
    PATTERN is a string containing '%i' or '%0Xi', for example 'image%04i.png'
    if specified, --copy will copy the files/directory instead of moving them
    if specified, INTEGER is the first index to be used (default=0)
    PATHS is a list of file or directories

Examples:
    
    collect.py image%04i.png  *.png
       will rename image files to: image0000.png, image0001.png, etc.
    
    collect.py --copy image%04i.png 1 run*/image.png"
       will copy the image files, starting at index 1
    
F. Nedelec, 2012--2017. Last modified 1.09.2017
"""


import sys, shutil, os, curses.ascii


#------------------------------------------------------------------------

def copy_recursive(src, dst):
    """Copy directory recursively"""
    if os.path.isfile(src):
        shutil.copy2(src, dst)
    elif os.path.isdir(src):
        try:
            os.mkdir(dst)
        except OSError:
            pass
        files = os.listdir(src)
        for f in files:
            s = os.path.join(src, f)
            d = os.path.join(dst, f)
            copy_recursive(s, d)


def main(args):
    """rename files"""
    do_copy = False
    arg = args.pop(0);
    # check if 'copy' specified before pattern
    if arg=='-c' or arg=='--copy' or arg=='copy=1':
        do_copy = True
        pattern = args.pop(0);
    else:
        pattern = arg
    # check validity of the pattern
    if os.path.isfile(pattern):
        sys.stderr.write("Error: first argument should be the pattern used to build output file name")
        return 1
    pattern % 0
    try:
        str = ( pattern % 0 )
    except:
        sys.stderr.write("Error: the pattern should accept an integer: eg. '%i' or '%04i'\n")
        return 1
    for c in str:
        if curses.ascii.isspace(c):
            sys.stderr.write("Error: the pattern includes or generates white space character\n")
            return 1
    paths = []
    idx = 0
    # parse arguments:
    for arg in args:
        if arg=='-c' or arg=='--copy' or arg=='copy=1':
            do_copy = True
        elif args[0].isdigit():
            idx = int(args[0])
        elif os.path.isfile(arg) or os.path.isdir(arg):
            paths.append(arg)
        else:
            sys.stderr.write("Error: '%s' is not a file or directory" % arg)
            return 1
    # process all files
    res = []
    for src in paths:
        while idx < 1000000:
            dst = pattern % idx
            idx += 1
            if dst == src:
                res.append(dst)
                break
            if not os.path.exists(dst):
                #make directory if name include a directory that does not exist:
                dir = os.path.dirname(dst)
                if dir and not os.path.isdir(dir):
                    os.mkdir(dir)
                # process file:
                if do_copy:
                    copy_recursive(src, dst)
                else:
                    os.rename(src, dst)
                res.append(dst)
                print("%s -> %s" % (src, dst))
                break
    return res


#------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1]=='help':
        print(__doc__)
    else:
        main(sys.argv[1:])


