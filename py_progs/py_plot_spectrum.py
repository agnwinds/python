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
import pylab as p
import py_read_output as r 
import numpy as np 
import os, sys

has_astropy = True 
try:
    import astropy
except ImportError:
    has_astropy = False

def smooth(x,window_len=20,window='hanning'):

	'''smooth data x by a factor with window of length window_len'''

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




def make_spec_plot(s, fname, smooth_factor = 10, angles = True, components = False):

	'''
    make a spectrum plot from astropy.table.table.Table object 

    Parameters
    ----------
    s: astropy.table.table.Table
    	table containing spectrum data outputted from Python 

    fname: str 
    	filename to save as e.g. sv

    smooth_factor: int 
    	factor you would like to smooth by, default 10

    angles: Bool 
    	Would you like to plot the viewing angle spectra? 

    components: Bool 
    	would you like to plot the individual components e.g. Disk Wind 
    
    Returns
    ----------
    Success returns 0
    Failure returns 1

    Saves output as "spectrum_%s.png" % (fname)
    '''

	if type(s) != astropy.table.table.Table:
		raise TypeError("make_spec_plot takes astropy.table.table.Table object as first arg")
		return 1

	if s.colnames[8] != "Scattered":
		print "Warning- colnames are not in expected order! %s != Scattered" % (s.colnames[8])

	ncomponents = 9

	if angles:

		# first make viewing angle plot
		p.figure(figsize=(8,12))

		nspecs = len(s.dtype.names) - ncomponents

		nx = 1
		ny = nspecs
		if nspecs > 4:
			nx = 2
			ny = (1 + nspecs) / nx 


		print "Making a %i by %i plot, %i spectra" % (nx, ny, nspecs)

		for i in range(nspecs):

			p.subplot(ny, nx, i+1)

			p.plot(s["Lambda"], smooth(s[s.dtype.names[ncomponents + i]], window_len = smooth_factor))

			p.title(s.dtype.names[ncomponents + i])
			p.xlabel("Wavelength")
			p.ylabel("Flux")

		p.savefig("spectrum_%s.png" % (fname))
		p.clf()

	if components:

		p.figure(figsize=(8,12))
		p.subplot(211)
		p.plot(s["Lambda"], smooth(s["Created"], window_len = smooth_factor), label="Created")
		p.plot(s["Lambda"], smooth(s["Emitted"], window_len = smooth_factor), label="Emitted")

		p.subplot(212)

		for i in range(4,9):
			p.plot(s["Lambda"], smooth(s[s.dtype.names[i]], window_len = smooth_factor), label=s.dtype.names[i])

		p.xlabel("Wavelength")
		p.ylabel("Flux")
		p.legend()

		p.savefig("spec_components_%s.png" % (fname))
		p.clf()
		
	return 0


def make_spec_plot_from_class(s, fname, smooth_factor = 10, angles = True, components = False):

	'''
    make a spectrum plot from py_classes.specclass object

    Parameters
    ----------
    s: specclass object
    	table containing spectrum data outputted from Python 

    fname: str 
    	filename to save as e.g. sv

    smooth_factor: int 
    	factor you would like to smooth by, default 10

    angles: Bool 
    	Would you like to plot the viewing angle spectra? 

    components: Bool 
    	would you like to plot the individual components e.g. Disk Wind 
    
    Returns
    ----------
    Success returns 0
    Failure returns 1

    Saves output as "spectrum_%s.png" % (fname)
    '''


	if angles:

		# first make viewing angle plot
		p.figure(figsize=(8,12))
		nspecs = len(s.spec)

		nx = 1
		ny = nspecs
		if nspecs > 4:
			nx = 2
			ny = (1 + nspecs) / nx 


		print "Making a %i by %i plot, %i spectra" % (nx, ny, nspecs)

		for i in range(nspecs):


			p.subplot(ny, nx, i+1)

			p.plot(s.wavelength, smooth(s.spec[i], window_len = smooth_factor))

			p.xlabel("Wavelength")
			p.ylabel("Flux")

		p.savefig("spectrum_%s.png" % (fname))
		p.clf()

	if components:

		p.figure(figsize=(8,12))
		p.subplot(211)
		p.plot(s.wavelength, smooth(s.created, window_len = smooth_factor), label="Created")
		p.plot(s.wavelength, smooth(s.emitted, window_len = smooth_factor), label="Emitted")

		p.subplot(212)
		p.plot(s.wavelength, smooth(s.censrc, window_len = smooth_factor), label="CenSrc")
		p.plot(s.wavelength, smooth(s.disk, window_len = smooth_factor), label="Disk")
		p.plot(s.wavelength, smooth(s.wind, window_len = smooth_factor), label="Wind")
		p.plot(s.wavelength, smooth(s.hitsurf, window_len = smooth_factor), label="HitSurf")
		p.plot(s.wavelength, smooth(s.scattered, window_len = smooth_factor), label="Scattered")
		p.xlabel("Wavelength")
		p.ylabel("Flux")
		p.legend()

		p.savefig("spec_components_%s.png" % (fname))
		p.clf()

	return 0



# try:
# 	fname = sys.argv[1]
# except IndexError:
# 	print '''No filename provided: 
# usage: py_plot_spectrum filename'''
# 	sys.exit()

# # s = r.read_spec_file(fname)

# # make_plot(s, components = True)







