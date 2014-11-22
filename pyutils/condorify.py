#!/usr/bin/env  python

import os
import json
import sys
import itertools
try:
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('jfilename')
    p.add_argument('-k','--keys',nargs='+',
        help="List of keys that contain lists that need to be ``condorized''.")
    args = p.parse_args()
except ImportError:
    class Args:
        def __init__(self):
            self.jfilename = ''
            self.keys = []
    args = Args()
    args.jfilename = sys.argv[1]
    args.keys = []
    for arg in sys.argv[2:]:
        if arg in ('-k','--keys'):
            continue
        args.keys.append(arg)

with open(args.jfilename, 'rb') as f:
    jdat = json.load(f)

x = []
for k in args.keys:
    x.append(list(jdat[k]))
    jdat.pop(k,0)

params = jdat.copy()

jdat['config'] = []
for i,config in enumerate(itertools.product(*x)):
    d = dict(zip(args.keys,config))
    jdat['config'].append(d)
    params.update(d)
    try:
        jobdir = '{:03d}'.format(i)
    except:
        jobdir = '%03d' % i

    try:
        os.makedirs(jobdir)
    except OSError:
        pass
    paramsPath = os.path.join(jobdir,'params.json')
    with open(paramsPath,'wb') as f:
        json.dump(params,f,sort_keys=True,indent=2,separators=(',',': '))

with open('MASTER.json','wb') as f:
    json.dump(jdat,f,sort_keys=True,indent=2,separators=(',',': '))
