#!/usr/bin/env python 
'''
	University of Southampton -- JM -- 30 September 2013

				read_output.py

Synopsis:
	this enables one to read outputs from the Python radiative transfer code

Usage:
	

Arguments:
'''

# we need the classes and numpy modules 
import py_classes as cls
import numpy as np
import matplotlib.pyplot as plt
import subprocess


def read_spec_file (filename, new=True):
    
    '''reads a Python .spec file and places in specclass array,
       which is returned'''
    
    if not '.spec' in filename: 
        filename = filename + '.spec'

    if new:
        add = 0
    else:
        add = 1
        
    
    # initialise the spectrum array with blank arrays
    spectrum = cls.specclass ([],[],[],[],[],[],[], [], [], []) 
    
    # first read the file into a temporary storage array
    spectrum_temp = np.loadtxt (filename, comments ='#', unpack=True)
    
    # now set the correct elements of the class to be the temporary array
    #spectrum.tot = spectrum_temp[0]
    spectrum.freq = spectrum_temp[0]
    spectrum.wavelength = spectrum_temp[1]

    if new:
        spectrum.created = spectrum_temp[2]

    spectrum.emitted = spectrum_temp[3-add]
    spectrum.censrc = spectrum_temp[4-add]
    spectrum.disk = spectrum_temp[5-add]
    spectrum.wind = spectrum_temp[6-add] 
    spectrum.scattered = spectrum_temp[8-add]
    spectrum.hitsurf = spectrum_temp[7-add]
    spectrum.spec = spectrum_temp[9-add:]
        
     #finally, return the spectrum class which is a series of named arrays      
    return spectrum



def read_spectot_file (filename, new=True):
    
    '''reads a Python .spec_tot file and places in spectotclass array,
       which is returned'''
    
    if not 'spec_tot' in filename: 
        filename = filename + '.log_spec_tot'

    if new:
        add = 0
    else:
        add = 1
        
    
    # initialise the spectrum array with blank arrays
    spectrum = cls.spectotclass ([],[],[],[],[],[],[], [], []) 
    
    # first read the file into a temporary storage array
    spectrum_temp = np.loadtxt (filename, comments ='#', unpack=True)
    
    # now set the correct elements of the class to be the temporary array
    spectrum.freq = spectrum_temp[0]
    spectrum.wavelength = spectrum_temp[1]

    if new:     # this means the version was post 78a
        spectrum.created = spectrum_temp[2]
    spectrum.emitted = spectrum_temp[3-add]
    spectrum.censrc = spectrum_temp[4-add]
    spectrum.disk = spectrum_temp[5-add]
    spectrum.wind = spectrum_temp[6-add] 
    spectrum.hitsurf = spectrum_temp[7-add]
    
    #finally, return the spectrum class which is a series of named arrays    
    return spectrum






# set some standard parameters for plotting
def setpars():
    
	print 'Setting plot parameters for matplotlib.'
	plt.rcParams['lines.linewidth'] = 1.0
	plt.rcParams['axes.linewidth'] = 1.3
	plt.rcParams['font.family'] = 'serif'
	plt.rcParams['font.serif'] = 'Times New Roman'
	plt.rcParams['text.usetex']='True'
    
	return 0



def read_emissivity ( root ):
    
    '''Read macro atom emissivities from a root diag file. 
       Returns two arrays, kpkt_emiss and matom_emiss.'''
    
    if not "_0.diag" in root:
        root = root + "_0.diag"
    
    # search file for smieeisivities string    
    matom_emiss, matom_abs = [], []
    kpkt_emiss = 0
    with open(root, 'r') as searchfile:
        for line in searchfile:
            
            # check if we have a matom_diagnostics line reporting level emissivities
            if 'emissivities' in line:
                
                data = line.split()
                for i in range(len(data)):
                                       
                    # we now put the appropriate values in corresponding arrays 
                    if data[i] == 'matom_abs':
                        matom_abs.append(float(data[i+1]))
                        
                    if data[i] == 'kpkt_abs':
                        kpkt_abs = float(data[i+1])
                        
                    if data[i] == 'matom_emiss':
                        matom_emiss.append(float(data[i+1]))
                        
                    if data[i] == 'kpkt_emiss':
                        kpkt_emiss = float(data[i+1])
                        

    # convert to numpy arrays
    matom_emiss = np.array(matom_emiss)
    matom_abs = np.array(matom_abs)
    
    return matom_emiss, kpkt_emiss



def thinshell_read ( root ):
    
    '''Read py_wind output filename for thin shell models with one cell'''
    
    inp = open(root, 'r')
    
    for line in inp:
        data = line.split()
        if data[0]!="#":
            if data[2]=="0":
                value = float(data[1])

    return value


def read_convergence (root ):

	''' check convergence in a file '''

	if not "_0.diag" in root:
        	root = root + "_0.diag"
	conv_fraction = []
	with open(root, 'r') as searchfile:
        	for line in searchfile:

           	 	# check if we have a matom_diagnostics line reporting level emissivities
            		if 'Summary  convergence' in line:
                
                		data = line.split()
				conv_fraction.append(float (data[3]))
                print conv_fraction

	final_conv = conv_fraction [-1]
	return final_conv
		



def read_pf(root):

    '''
    reads in a pf file and places in a dictionary
    access the number of observers with 

    dict = read_pf(root)
    nobs = dict["no_observers"]
    '''

    if not ".pf" in root:
        root = root + ".pf"

    params, vals = np.loadtxt(root, dtype="string", unpack=True)

    pf_dict = dict()

    old_param = None 
    old_val = None

    for i in range(len(params)):


        # convert if it is a float
        try:
            val = float(vals[i])

        except ValueError:
            val = vals[i]

        if params[i] == old_param:

            if isinstance(pf_dict[params[i]], list):
                pf_dict[params[i]].append(val)

            else:
                pf_dict[params[i]] = [old_val, val]

        else:
            pf_dict[params[i]] = val

        old_param = params[i]
        old_val = val


    return pf_dict






