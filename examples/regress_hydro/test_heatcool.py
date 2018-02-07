#!/usr/bin/env python3




import csv, sys, os, array, warnings
import numpy as np
import copy
from astropy.io import ascii
import pickle

from astropy import constants as consts


serious_err=0
minor_err=0

minmax=pickle.load(open("minmax.dat","rb"))	


nlim=int(sys.argv[2])
fraclim=float(sys.argv[3])


py_data=ascii.read(sys.argv[1])


fixed=['i','j','rcen','thetacen','vol','temp','rho','n_h']



for key in fixed:
	nmin=0
	nmax=0
	for j in range(len(py_data[key])-1):
		one=minmax['min'][key][j]
		two=py_data[key][j]
		three=minmax['max'][key][j]
		if (two<one):
			nmin=nmin+1
		if (two>three):
			nmax=nmax+1
	if (nmax!=0 or nmin!=0):
		print("Serious error "+key+" has changed")
		serious_err=serious_err+1
	
	
for key in py_data.keys():
	nmin=0
	nmax=0
	fracmin=0.0
	fracmax=0.0
	for j in range(len(py_data[key])-1):
		one=minmax['min'][key][j]
		two=py_data[key][j]
		three=minmax['max'][key][j]
		if (two<one):
			test=(one-two)/one
			if test>fracmin:
				fracmin=test
			nmin=nmin+1
		if (two>three):
			test=(two-three)/three
			if test>fracmax:
				fracmax=test
			nmax=nmax+1
	if np.max([nmax,nmin])>nlim:
		print("Minor error "+key+" is outside limits in "+str(np.max([nmax,nmin]))+" cells")
		minor_err=minor_err+1
	elif np.max([fracmax,fracmin])>fraclim:
		minor_err=minor_err+1
		print("Minor error a cell has changed "+key+" by "+str(np.max([fracmax,fracmin])*100.)+"%")

if serious_err==0 and minor_err==0:
	sys.stderr.write("0\n")
else:
	sys.stderr.write("1\n")
	