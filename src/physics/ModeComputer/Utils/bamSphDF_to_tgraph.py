#!/usr/bin/env python

# bamSphDF_to_tgraph.py - add theta- and phi-coordinates to bam sphere
# output files.

from __future__ import print_function

import sys

######################################################################

# convert a string to float similar to C's atof
def WT_atof(str):
  list = str.split()
  str1 = list[0]
  try:
    fl = float(str1)
  except:
    fl = 0.0
  return fl

# find and get value of a parameter
def getparameter(line, par):
  val = ''
  ok = 0
  EQsign = 0
  p = line.find(par)
  if p > 0:
    before = line[p-1:p]
    if before.isalpha():
      p=-1
  if p >= 0:
    ok = 1
    afterpar = line[p+len(par):]
    # strip white space at end and beginning
    afterpar  = afterpar.rstrip()
    afterpar2 = afterpar.lstrip()
    if len(afterpar2) > 0:
      # if '=' is there
      if afterpar2[0] == '=':
        afterpar2 = afterpar2[1:]
        val = afterpar2.lstrip()
        EQsign = 1
      # if we have a space instead
      elif afterpar[0].isspace():
        val = afterpar.lstrip()
        EQsign = 0

  return (val , ok, EQsign)

######################################################################

if len(sys.argv) != 2:
    print('bamSphDF_to_tgraph.py:')
    print('add theta- and phi-coordinates to bam sphere output files.\n')
    print('Usage:')
    print('bamSphDF_to_tgraph.py filename')
    exit(0)

argvs = sys.argv[1:]
rfile = argvs[0]

with open(rfile, 'rb') as f:
    while True:
        line = f.readline()
        if not line:
            break

        (val,ok,EQsign) = getparameter(line.lower().decode('ascii'), 'gridfunction')
        if ok == 1:
            gridfunction = val

        (val,ok,EQsign) = getparameter(line.lower().decode('ascii'), 'r')
        if ok == 1:
            radius = val

        (val,ok,EQsign) = getparameter(line.lower().decode('ascii'), 'ntheta')
        if ok == 1:
            ntheta = int(val)

        (val,ok,EQsign) = getparameter(line.lower().decode('ascii'), 'nphi')
        if ok == 1:
            nphi = int(val)

        (val,ok,EQsign) = getparameter(line.lower().decode('ascii'), 'theta')
        if ok == 1:
            li = val.split()
            theta0 = WT_atof(li[1])

        (val,ok,EQsign) = getparameter(line.lower().decode('ascii'), 'dtheta')
        if ok == 1:
            dtheta = float(val)

        (val,ok,EQsign) = getparameter(line.lower().decode('ascii'), 'dphi')
        if ok == 1:
            dphi = float(val)

        (val,ok,EQsign) = getparameter(line.lower().decode('ascii'), 'iteration')
        if ok == 1:
            break;

    npoints = ntheta*nphi
    ntimes = 0
    readdata = 0
    tim = []
    dat = []
    for line in f:
        (val,ok,EQsign) = getparameter(line.lower().decode('ascii'), 'time')
        if ok == 1:
            tim.append(val)
            ntimes +=1
            readdata = 1
            ndat = 0
            continue

        if readdata == 1:
            li = line.split()
            for l in li:
                dat.append(float(l))
            ndat = ndat + len(li)
            if ndat >= npoints:
                readdata = 0
            

# ok now we have read all data, now write it back with point coords
print('# gridfunction = ', gridfunction)
print('# radius =', radius)
print('# 3 columns with: theta, phi,', gridfunction)
print()
di = 0
for ti in range(0,ntimes):
    print('# time =', tim[ti])
    for k in range(0,nphi):
        for j in range(0,ntheta):
            theta = theta0 + j*dtheta
            phi   = k*dphi
            print(theta, phi, dat[di])
            di += 1
        print()
    print()
