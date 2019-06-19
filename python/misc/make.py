#!/usr/bin/env python
#
# copyright F. Nedelec, July 2011
#
# Unfinished code

"""
Usage:

    make.py command [options] [directories]

Options are specified as 'option=value', without space around the '=' sign.
Existing options and their values:

     walk    0 or 1             default = 0
      
Example:

     make.py 

F. Nedelec, July 2011,
"""

try:
    import sys, os
except ImportError:
    sys.stderr.write("  Error: could not load necessary python modules\n")
    sys.exit()



#------------------------------------------------------------------------

action='none'

dictionnary= { }


#------------------------------------------------------------------------

def process_dir(dirpath, directories, filenames):
    """make something in directory dirpath"""
    global action, dictionnary
    cdir = os.getcwd()
    os.chdir(dirpath)
    print('- '*32+" DIR "+dirpath)
    if action=='none':
        pass
    else:
        res = eval(action, globals(), dictionnary)
        print res
    os.chdir(cdir)



def process(dirpath):
    """make something in directory dirpath"""
    global action, dictionnary
    if os.path.isdir(dirpath):
        files = os.listdir(dirpath)
        for f in files:
            if os.path.isdir(f):
                files.remove(f)
        process_dir(dirpath, [], files)
    else:
        print('- '*32+" FILE "+dirpath)
        if action=='none':
            pass
        else:
            res = eval(action, globals(), dictionnary)
            print res



#------------------------------------------------------------------------

def main(args):
    """process command line arguments"""
    global action, dictionnary
    walk = 0
    places = []
    
    try:
        action = args[1]
    except IndexError:
        print("Error: an action must be specified as first argument")
        sys.exit()

    for arg in args[1:]:
        [key, equal, value] = arg.partition('=')
        
        if key=='' or equal!='=' or value=='':
            if os.path.isfile(arg) or os.path.isdir(arg):
                places.append(arg)
            else:
                print("ignored '%s' on command line" % arg)
        else:
            if key=='walk':
                walk = int(value)
            else:
                dictionnary[key]=value
    
    print('dictionnary = ', dictionnary)
    
    try:
        if walk:
            for root, dirs, files in os.walk('.', topdown=False):
                process_dir(root, dirs, files)
        else:
            if len(places) == 0:
                process(os.getcwd())
            else:
                for place in places:
                    process(place)
    except IOError as e:
        print("Error: "+repr(e))


#------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])


