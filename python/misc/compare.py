#!/usr/bin/env python
#
# compare.py
#
# copyright F. Nedelec, December 14th 2007, 14.03.2018

"""
compare.py
    Compare files from two root directories
Usage:
    compare.py root1 root2 [opendiff]
"""

from __future__ import print_function
import sys, os, shutil, time, subprocess

exe = "diff"
diff="diff --side-by-side -W200 -p --suppress-common-lines"


def spacer(info):
    """print a line of width size, with 'info' in the middle"""
    import os
    rows, cols = os.popen('stty size', 'r').read().split()
    sys.stdout.write(chr(27)+"[36;2m"); sys.stdout.flush()
    print(info.center(int(cols), '-'))
    sys.stdout.write(chr(27)+"[0m"); sys.stdout.flush()


def compareFiles(fileL, fileR):
    comp  = os.popen(diff+" -q "+fileL+" "+fileR)
    empty = ( len(comp.read()) == 0 )
    comp.close()
    
    if not empty:
        if exe == 'diff':
            spacer('%s %s'% ( fileL, fileR ))
            comp = os.popen(diff+" "+fileL+" "+fileR)
            for line in comp:
                print(line, end='')
            comp.close()
            
            sys.stdout.write(chr(27)+"[32;2m")
            print("This was %40s" % fileL)
            ans = raw_input('Action? return/left/right/open/q >')
            sys.stdout.write(chr(27)+"[0m")
            
            if ans == "left" or ans == "l":
                shutil.copyfile(fileL, fileR)
            elif ans == "right" or ans == "r":
                shutil.copyfile(fileR, fileL)
            elif ans == "open":
                os.system("opendiff "+fileL+" "+fileR+"&")
            elif ans == "q":
                sys.exit()
        else:
            subprocess.call(["opendiff", fileL, fileR])
            #we wait a bit for the application to start
            time.sleep(0.5)


def interesting(file):
    return ( file.endswith('.py') or file.endswith('.m')
        or file.endswith('.h') or file.endswith('.cc')
        or file.endswith('.md') or file.endswith('.txt')
        or file.startswith('makefile') or file.startswith('.cym') )


def process_dir(roots, dirnameL, files):
    """compare files in the current directory"""
    #print("dirname=%s  args=%s" % (dirnameL,args))
    if 0 <= dirnameL.find('.svn'):
        return
    if 0 <= dirnameL.find('.git'):
        return
    if 0 <= dirnameL.find('DerivedData'):
        return
    if 0 <= dirnameL.find('bin'):
        return
    if 0 <= dirnameL.find('build'):
        return
    dirnameR = dirnameL.replace(roots[0], roots[1])
    spacer(dirnameL)
    for file in files:
        if interesting(file):
            compareFiles( dirnameL+"/"+file, dirnameR+"/"+file)

#------------------------------------------------------------------------

def main(args):
    global exe
    rootL=args[0].rstrip('/')
    rootR=args[1].rstrip('/')

    #parse command-line arguments:    
    for arg in args[2:]:
        if arg == 'opendiff':
            exe=arg
        else:
            print("unknown argument '%s'" % arg)
            sys.exit()

    if rootL == '' or rootR == '':
        print('Error: you must specify root directories!')

    for path, dirs, files in os.walk(rootL, topdown=False):
        process_dir([rootL, rootR], path, files)

if __name__ == "__main__":
    if len(sys.argv) < 3 or sys.argv[1]=='help':
        print(__doc__)
    else:
        main(sys.argv[1:])


