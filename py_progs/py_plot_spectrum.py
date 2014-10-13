#!/usr/bin/env python 
'''
University of Southampton
James Matthews
09 October 2014

py_plot_spectrum.py 

py_plot_spectrum.py creates spectrum plots 
from a file root.pf.
'''

#import pylab as p 
from pylab import *
import py_read_output as r 
import numpy as np 
import os, sys

def smooth(x,window_len=20,window='hanning'):

        if x.ndim != 1:
                raise ValueError, "smooth only accepts 1 dimension arrays."

        if x.size < window_len:
                raise ValueError, "Input vector needs to be bigger than window size."

        if window_len<3:
                return x

        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
                raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

        s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]


        if window == 'flat': #moving average
                w = np.ones(window_len,'d')
        else:  
                w = eval('np.'+window+'(window_len)')

        y=np.convolve(w/w.sum(),s,mode='same')

        return y[window_len:-window_len+1]


def make_plot(s, smooth_factor = 10, angles = True, components = False):


	if angles:

		# first make viewing angle plot
		figure(figsize=(8,12))
		nspecs = len(s.spec)

		nx = 1
		ny = nspecs
		if nspecs > 4:
			nx = 2
			ny = (1 + nspecs) / nx 


		print "Making a %i by %i plot, %i spectra" % (nx, ny, nspecs)

		for i in range(nspecs):


			subplot(ny, nx, i+1)

			plot(s.wavelength, smooth(s.spec[i], window_len = smooth_factor))

			xlabel("Wavelength")
			ylabel("Flux")

		savefig("spectrum_%s.png" % (fname))
		clf()

	if components:

		figure(figsize=(8,12))
		subplot(211)
		plot(s.wavelength, smooth(s.created, window_len = smooth_factor), label="Created")
		plot(s.wavelength, smooth(s.emitted, window_len = smooth_factor), label="Emitted")

		subplot(212)
		plot(s.wavelength, smooth(s.censrc, window_len = smooth_factor), label="CenSrc")
		plot(s.wavelength, smooth(s.disk, window_len = smooth_factor), label="Disk")
		plot(s.wavelength, smooth(s.wind, window_len = smooth_factor), label="Wind")
		plot(s.wavelength, smooth(s.hitsurf, window_len = smooth_factor), label="HitSurf")
		plot(s.wavelength, smooth(s.scattered, window_len = smooth_factor), label="Scattered")
		xlabel("Wavelength")
		ylabel("Flux")
		legend()

		savefig("spec_components_%s.png" % (fname))
		clf()

	return 0



try:
	fname = sys.argv[1]
except IndexError:
	print '''No filename provided: 
usage: py_plot_spectrum filename'''
	sys.exit()

s = r.read_spec_file(fname)

make_plot(s, components = True)







