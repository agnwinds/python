#!/usr/bin/env python
'''
	balmer_decrement.py 

runs tests of the Balmer decrement for a one zone 
thin shell Python model. Involves running py_wind on
a wind_save file and reading some output files. Compares
to Osterbrock values.

Usage:
	python BalmerDecrement.py root_filename
	python BalmerDecrement.py -h for help 

Requirements:
	py_wind
	numpy 
	matplotlib 
        py_plot_util, py_read_output from $PYTHON/py_progs in the python path 
Notes:
    This routine is a routine to check the reulsts of running a one zone
    model intended to produce the Balmer decrements.  It does not run the
    model but meery checks it
'''
import numpy as np 
import py_plot_util as util 
import sys, os
import py_read_output as rd 


# some global parameters to test against
osterbrock_ratios = np.array([2.86,1.0,0.470,0.262,0.159,0.107,0.0748,0.0544])
osterbrock_intensity = 8.30e-26
TOLERANCE = 0.3

def BalmerTest(root, plotit=True):
	'''
	runs tests of of the Balmer decrement for a one zone 
	thin shell Python model. Involves running py_wind on
	a wind_save file and reading some output files. Compares
	to Osterbrock values.
	'''

	print ("Running Balmer Test for run {}...".format(root))

	# create the list of commands to run in py wind
	nlevels = 8
	cmds = ["1","s","n","i","1","1","1","1","2","0","M","2"]
	for i in range(nlevels):
		cmds.append("{:d}".format(i+3))
	cmds.append("-1")
	cmds.append("q")

	# run py wind 
	util.run_py_wind(root, cmds=cmds)

	# these could be in principle be used to check absolute emissivity values 
	# ne = rd.read_pywind("{}.ne.dat".format(root), mode="1d")[2][1]
	# vol = rd.read_pywind("{}.vol.dat".format(root), mode="1d")[2][1]
	# nh1 = rd.read_pywind("{}.ioncH1.dat".format(root), mode="1d")[2][1]
	# nh2 = rd.read_pywind("{}.ioncH2.dat".format(root), mode="1d")[2][1]
	# nprot = nh1 + nh2

	# read the emissivities
	ratios = np.zeros(nlevels)
	for i in range(nlevels):
		ratios[i] = rd.read_pywind("{}.lev{}_emiss.dat".format(root,i+3), mode="1d")[2][1]



	if plotit:
		# make a scatter plot of the ratios to Hbeta (Balmer decrement)
		n = np.arange(nlevels)+3
		plt.scatter(n,ratios/ratios[1], c="k")
		plt.scatter(n,osterbrock_ratios, facecolors="None", edgecolors="k", s=80)
		plt.xlabel("$n_u$", fontsize=16)
		plt.ylabel("Balmer Decrement", fontsize=16)
		plt.savefig("BalmerDecrement_{}.png".format(root))

	# define a 0 or 1 pass or fail 
	pass_fail = (np.fabs(ratios/ratios[1] - osterbrock_ratios)/ osterbrock_ratios) 
	print ("\n----------------------------------")
	print ("\n\nArray of line ratio relative errors for Balmer series:\n", pass_fail)

	return (np.all(pass_fail < TOLERANCE))
	


if __name__ == "__main__":

	if len (sys.argv) > 1:
		if sys.argv[1] == "-h":
			print (__doc__)
		else:
			root = sys.argv[1]

			plotit=True
			# only plot if matplotlib installed 
			try:
				import matplotlib.pyplot as plt 
			except:
				print ("No matplotlib, so not making plot.")
				plotit = False

			# run the test 
			ifail = BalmerTest(root, plotit=plotit)

			# Tell the user whether the test is passed. 
			print ("\nTest passed?:", ifail) 
			if ifail == False: # possible this should be an exception instead?
				print ("ERROR: Balmer emissivities did not match those expected\n")
				sys.exit(-1)
	else: 
		print (__doc__)
