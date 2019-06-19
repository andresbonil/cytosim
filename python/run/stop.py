#!/usr/bin/env python
#
#
#
# Copyright F. Nedelec, March 16th, 2008 - 2013


import os, sys, time

pid    = ''
action = ''


# parse command-line arguments:
for arg in sys.argv[1:]:
    if os.path.isdir(arg):
        os.chdir(arg)
    else:
        action = arg


# read 'log.txt' to extract information on the process
with open('log.txt') as f:
    for line in f:
        arg = line.split()
        if arg[0] == 'pid':
            pid = arg[1]
        elif arg[0] == 'ppid':
            pid = arg[1]
        elif arg[0] == 'exec':
            exe = arg[1]
        elif arg[0] == 'host':
            host = arg[1]
            if not host == os.getenv('HOSTNAME', 'unknown'):
                print("running on different host '%s'" % host)
                sys.exit()

if not pid:
    print('could not get process id')
    sys.exit()

print('> process %s %s' % (pid, exe))

if action == 'kill':
    os.system('kill -2 '+pid)
    with open('run.txt','a') as f:
        f.write("stopped   %s\n" % time.asctime())
else:
    print('no action taken (to interrupt the simulation, specify "kill" as argument)')


