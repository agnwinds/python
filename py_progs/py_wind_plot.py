#!/usr/bin/env python 
'''
University of Southampton
James Matthews
09 October 2014

py_plot_spectrum.py 

py_plot_spectrum.py creates wind plots 
from a file root.pf.
'''

#import pylab as p 
from pylab import *
import py_read_output as r 
import numpy as np 
import os, sys

standard_cmds = np.array(["1", "n","t", "r","v","1","2","3","-1",\
						 "I", "i", "1","1","1","2","0","i","0", \
						 "1","1","1","2","2","2","2","1","2","3",
						 "6", "3","6","4","6","5","0","q"])


ion_standard_variables = ["ionh1", "ionhe1", "ionhe2", "ionc4", "ionc5"]

ion_standard_variables


def run_py_wind (vers, fname, cmds=None, ilv=None):
	'''
	run version vers of py_wind on file fname.wind_save
	'''

	if cmds == None:
		cmds = standard_cmds

	x = cmds
	np.savetxt("_tempcmd.txt", x, fmt = "%s")


	isys = os.system('py_wind'+vers+' '+fname+' < _tempcmd.txt > tempfile &')
	time.sleep(3)

	# remove temporary file
	#os.system("rm -f _tempcmd.txt")
	return isys

def read_pywind_smart(filename, return_inwind=False):
	'''
	read a py_wind file using np array reshaping and manipulation
	'''

	# first, simply load the filename 
	d = np.loadtxt(filename, comments="#", dtype = "float", unpack = True)

	# our indicies are already stored in the file- we will reshape them in a sec
	zindices = d[-1]
	xindices = d[-2]

	# we get the grid size by finding the maximum in the indicies list 99 => 100 size grid
	zshape = int(np.max(zindices) + 1)
	xshape = int(np.max(zindices) + 1)


	# reshape our indices arrays
	xindices = xindices.reshape(xshape, zshape)
	zindices = zindices.reshape(xshape, zshape)

	# now reshape our x,z and value arrays
	x = d[0].reshape(xshape, zshape)
	z = d[1].reshape(xshape, zshape)

	values = d[2].reshape(xshape, zshape)

	# these are the values of inwind PYTHON spits out
	inwind = d[3].reshape(xshape, zshape)

	# create an inwind boolean to use to create mask
	inwind_bool = (inwind >= 0)
	mask = (inwind < 0)

	# finally we have our mask, so create the masked array
	masked_values = np.ma.masked_where ( mask, values )

	#print xshape, zshape, masked_values.shape

	#return the transpose for contour plots.
	if return_inwind:
		return x, z, masked_values.T, inwind_bool.T
	else:
		return x, z, masked_values.T


def make_windplot(filename, cmds = standard_cmds, variables = standard_variables, run=True):

	'''
	make a four by two wind plot of eight variables 
	'''





