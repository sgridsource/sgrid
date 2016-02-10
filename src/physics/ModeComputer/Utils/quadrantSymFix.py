#!/usr/bin/env python

# quadrantSymFix.py - fix output bug in old bam waves, assuming quadrant
# symmetry in region 0<theta<pi/4, phi>pi .

from __future__ import print_function

import sys

######################################################################

# convert a string to float similar to C's atof
def WT_atof(str):
  if len(str) == 0:
    return 0.0
  list = str.split()
  str1 = list[0]
  fl = 0.0
  while len(str1)>0:
    try:
      fl = float(str1)
      break
    except:
      str1 = str1[:-1]
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

if len(sys.argv) < 2:
    print('quadrantSymFix.py:')
    print('fix sphere output bug in old bam waves by assuming quadrant symmetry\n')
    print('Usage:')
    print('quadrantSymFix.py file1 file2 ...')
    exit(0)

#PI = math.pi

argvs = sys.argv[1:]
rfile = argvs[0]

for rfile in argvs:
    print('Reading', rfile)
    with open(rfile, 'rb') as f:
        ntheta = 0
        nphi = 0
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

            (val,ok,EQsign) = getparameter(line.lower().decode('ascii'), 'n2')
            if ok == 1:
                n2 = int(WT_atof(val))

            (val,ok,EQsign) = getparameter(line.lower().decode('ascii'), 'n3')
            if ok == 1:
                n3 = int(val)

            (val,ok,EQsign) = getparameter(line.lower().decode('ascii'), 'gridtype')
            if ok == 1:
                gridtype = val

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
                break;

        npoints = ntheta*nphi
        ntimes = 0
        readdata = 0
        tim = []
        dat = []
        itr = []
        for line in f:
            (val,ok,EQsign) = getparameter(line.lower().decode('ascii'), 'iteration')
            if ok == 1:
                iteration = int(val)
            (val,ok,EQsign) = getparameter(line.lower().decode('ascii'), 'time')
            if ok == 1:
                tim.append(val)
                itr.append(iteration)
                ntimes += 1
                readdata = 1
                ndat = 0
                continue

            if readdata == 1:
                li = line.split()
                for l in li:
                    dat.append(l)
                ndat = ndat + len(li)
                if ndat >= npoints:
                    readdata = 0

    # check if we read data
    if npoints == 0:
        print('  Error no data found!')
        continue

    # now that we have read all data, give some info
    print('# gridfunction =', gridfunction)
    print('# r =', radius)
    print('# SphericalDF:   n2 =', n2, ' n3 =', n3)
    print('# gridtype =', gridtype)
    print('# ntheta =', ntheta)
    print('# nphi   =', nphi)
    print('# theta = { ', theta0, ', ', theta0 + (ntheta-1)*dtheta, ' }', sep='')
    print('# phi   = { 0.0,', (nphi-1)*dphi, '}')
    print('# dtheta =', dtheta)
    print('# dphi   =', dphi)

    # ok now we have read all data, now write it back imposing quadrant Sym
    print('Writing', rfile)
    with open(rfile, 'w') as f:
        print('# gridfunction =', gridfunction, file=f)
        print('# r =', radius, file=f)
        print('# SphericalDF:   n2 =', n2, ' n3 =', n3, file=f)
        print('# gridtype =', gridtype, file=f)
        print('# ntheta =', ntheta, file=f)
        print('# nphi   =', nphi, file=f)
        print('# theta = { ', theta0, ', ', theta0 + (ntheta-1)*dtheta, ' }', sep='', file=f)
        print('# phi   = { 0.0,', (nphi-1)*dphi, '}', file=f)
        print('# dtheta =', dtheta, file=f)
        print('# dphi   =', dphi, file=f)
        di = 0
        for ti in range(0,ntimes):
            print('# iteration =', itr[ti], file=f)
            print('# time      =', tim[ti], file=f)
            for k in range(0,nphi):
                for j in range(0,ntheta):
                    # theta = theta0 + j*dtheta
                    # phi   = k*dphi
                    if k>=n3/2 and j<=((n2/2 + 1)/4 + 2): # j<=((n2/2 + 1)/4 + 1): should work too
                        # find index mdi of point at phi-PI
                        mdi = di - ntheta*nphi/2
                        print(dat[mdi], end=' ', file=f)        
                    else:
                        # print(theta, phi, dat[di])
                        print(dat[di], end=' ', file=f)
                    di += 1
                print(file=f)
            print(file=f)
