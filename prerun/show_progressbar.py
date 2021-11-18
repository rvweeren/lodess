#! /usr/bin/python3

import numpy as np
import subprocess as sp
import time
import sys
import glob

data = []
with open('html.txt','r') as handle:
    for line in handle:
        data.append(line)

try:
    deltatime = float(sys.argv[1])
except:
    deltatime = 1

htmllength = len(data)
prevsize = int(sp.check_output('du -s .',shell=True).decode('utf-8')[:-2])
time.sleep(deltatime)
while True:
    size = int(sp.check_output('du -s .',shell=True).decode('utf-8')[:-2])
    delta = (size - prevsize)/1000 # in MB/s
    delta /= deltatime
    deltsteps = np.linspace(0,65,200)
    prstr = '[ '
    for i in deltsteps:
        if delta > i:
            prstr += '#'
        else:
            prstr += '-'
    dls = len(glob.glob('L*/*'))
    prstr += f']  {delta:.02f} MB/s  |||  {dls:03}/{htmllength:03}'
    print(prstr,end='\r')

    time.sleep(deltatime)
    prevsize = size
