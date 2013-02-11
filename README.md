README
***
=========
Precursor python_68h
Changes - This is a major revision to the code, including the ability to define a power law source in the location of the central star, a new ionization mode using the sim factor, a new wind mode to generate a thin shell, and also associated changes to balance.
Several additions to the structures in python.h
New spectype (power law)
new system_type (AGN - just a disk and a power law source)
additions to geo structure for agn power law and luminosity, distance to AGN (for balance), sim weight and alpha for sim code
additions to allow for a thin shell test wind.
Additions to allow for the creation and writing out of spectrum files in logarithmic units
atomic.h
density_min parameter changed to 1e-100 - this is to avoid checks which force the sim code to ignore low density ions which may actually be important once we compute correction factors.
bands.c
New bands to allow arbitrary frequency bands to take in x-rays
extract.c
Small modification to take account of the new system types.
ionization.c
Some external variables (sim_numin, sim_numax and sim_meanfreq) set up to communicate with zbrent to allow sim estimateors to be calculated
zbrent added to the end of the file - this is a duplicate, and perhaps should be dome more cleverly
the function sim_alpha_funch is added to the end of the file. This is the equation 18 in thr 2008 paper and links mean frequency with alpha. It is the equation solved in zbrent.
in ionization, mode 4 is added for the initial run thruogh. This is LTE with sim correction, but weight and alpha set to initial values everywhere in the wind.
Mode 5 added, this calls oneshot then sim correction.
sim_numin and sim_numax are hardwired here to 1.25e15Hz and 1.21e19Hz. This might need correcting in the future.
It is in this file that the calculations are made to compute sim_alpha and sim_w. These are used in the sim correction code.
In one_shot, slight modification to take account of the new ionization modes. Also some printf lines
levels.c
Simple modification to take account of new modes. in if is called with mode 3, it does a non lte level calculation with weight=1 and temperature=t_e
partition.c
Modification to deal with new ionization mode. If called in mode 3, it computes partition functions with w=1 and t=t_e
pdf.c
Lots of modifications 'changed the way in which pdf_gen_from_function works so it will change the number of points that are used in the array pdf_array in situlations where the number of points that are used in the array where the binning is too coarse. October 2010 KSL
photon_gen.c
lines added to allow photons to be generated from AGN. AGN photons added to sum of total nmber of photons, and a call is made to the new function photo_gen_agn (located in agn.c)
python.c
n_ioniz and lum_ioniz are now declared externally so that the total ionizing photon number and luminosity (computed in photon_check) can be reported.
modifications to wind_ionization line to allow sim mode to be selected
modifications to allow the agn system to be defined, along with spectrum type for agn. Some calculations are done to inform the user of sensible radius for the AGN given the schwartzchild radius
A line is here after the call to bands_init to make sure that freqmin and fregmax are updated to take account of any user modified frequency bands.
calls to spectrum_summary added to allow logarithmic spectra to be generated and output.
*

saha.c
A new mode is defined in saha - mode 5. LTE followed by sim correction.
In concentrations, a new mode is also defined, mode 3 which is currently identical to mode 1.
In order to improve the agreement with SS's FE ionisation fractions in his 2010 paper, density(first) is set to 1e-250 instead of 1. This allows a greater dynamic range of SAHA predicted ions, which then allows huge sim correction factors to give more accurate results.
The big parameter which is designed to limit the ratios between adjacent ions in also increased by a factor of 1 million in the line big=big*1e6
There are also a few printf statements to help understand what is