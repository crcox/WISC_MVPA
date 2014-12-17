#!/usr/bin/env python

import csv
import json
import os
import sys
from datetime import datetime
import tarfile
import shutil

try:
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('master')
    p.add_argument('-r','--run',action='store_true',
        help="In addition to setting up the directory structure, run all analyses. Useful when the analysis is being done locally.")
    args = p.parse_args()
except ImportError:
    class Args:
        def __init__(self):
            self.master = ''
            self.run = False
    args = Args()

    # Log required positional args
    args.master = sys.argv[1]

    # Specify translation from valid flags to valid arguments.
    flagDict  = {'-r':'run','--r':'run'}

    # Specify type of input associated with argument:
    #   If you specify an integer, that declares the expected
    #   number of inputs.
    #   If you define an empty list, the number of inputs is taken
    #   to be > 0 but bounded.
    #   If you specify a logical value, this means that the flag
    #   takes no inputs, and it's presence will flip this default
    #   value.
    tmp = {'run':False}

    # Parse inputs and arguments
    lst = sys.argv[2:]
    while len(lst) > 0:
        flag = lst.pop(0)
        try:
            arg = flagDict[flag]
        except KeyError:
            print '\nERROR: {f} is not a known flag\n'.format(f=flag)
            raise KeyError

        if isinstance(tmp[arg],list):
            while lst[0] not in tmp.keys():
                val = lst.pop(0)
                tmp[arg].append(val)
        elif isinstance(tmp[arg],int):
            n = tmp[arg]
            if n == 1:
                val = lst.pop(0)
                tmp[arg] = val
            elif n > 1:
                tmp[arg]=[0]*n
                for i in range(n):
                    tmp[arg][i] = lst.pop(0)
        elif isinstance(tmp[arg],bool):
            tmp[arg] = not temp[arg]

    # Assign inputs into an args object as would be produced by argparse.
    args.run = tmp['run']

#############################################################
#   Load data and parameters from the "master" json file    #
#############################################################
jsonfile = args.master
with open(jsonfile,'rb') as f:
    jdat = json.load(f)

#############################################################
#     Set a current config that inherits defaults to be     #
#        potentially overwritten by individual configs      #
#############################################################
try:
    allConfigs = jdat['config']
except KeyError:
    # This just ensures there is a list to loop over.
    # The idea is to initialize currentConfig with defaults and then update
    # those with the paramaters from each config dict.
    allConfigs = [jdat['config']]

#############################################################
#  Define a root folder (current directory if not a condor  #
#                           run)                            #
#############################################################
rootdir = ''
try:
    if jdat['condor']:
        rootdir = jdat['version']
except KeyError:
    pass

#############################################################
#                 Setup directory structure                 #
#############################################################
sharedir = os.path.join(rootdir,'shared')
if not os.path.isdir(sharedir):
    os.makedirs(sharedir)

archivedir = os.path.join(rootdir,'archive')
if not os.path.isdir(archivedir):
    os.makedirs(archivedir)

#############################################################
#           Copy source and binary files into place         #
#############################################################
shutil.copy(jdat['binary'],sharedir)
shutil.copytree(jdat['srcdir'],sharedir)

#############################################################
#                     Write URLS file                       #
#############################################################
URLS = os.path.join(sharedir,'URLS')
with open(URLS,'w') as f:
    for x in jdata['data']:
        p = os.path.join(jdat['squid'],x)
        f.write(p+'\n')
    p = os.path.join(jdat['squid'],jdat['metadata'])
    f.write(p+'\n')

#############################################################
#        Loop over configs (if condor, do setup only)       #
#############################################################
for i, cfg in enumerate(allConfigs):
    # Reset to defaults on each loop
    currentConfig = jdat.copy()
    try:
        # Strip out the config list if it exists
        del currentConfig['config']
    except KeyError:
        pass

    # Update defaults with the new config data
    currentConfig.update(cfg)

    jobdir = os.path.join(rootdir,'{cfgnum:03d}'.format(cfgnum=i))
    os.makedirs(jobdir)

    cfgfile = os.path.join(jobdir,'params.json')
    with open(cfgfile,'wb') as f:
        json.dump(currentConfig, f, sort_keys=True, indent=2, separators=(',', ': '))
