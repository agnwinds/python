README
***
=========
Note: error in Gaunt factor!!

precusror - python_73e
This is intended to be the stable version which is to be used for generating data for the 2012 paper.
The main changes are:
atomic.h - A structure gaunt_total put in to hold Sutherlands gaunt factors
bands.c - the usual changes here, since the power law bands are still hardcoded. There are also new modes to try and match the PL bands.
bb.c - emittance_bb had been changed in 73e to integrate the bb function directly, rather than use tabulated data. This was causing crashes, so changed back in this version.
emission.c - several changes
put in a new counter for lum_adiabatic - there were problems with adiabatic cooling that were not picked up bacause it was not being reported
code put into total_ff and ff to use the sutheralnd gaunt factor data
a new routine gaunt_ff written to compute the gaunt factor for free free from sutherlands data - it interpolates on gsquared.
get_atomicdata.c - lines to read in gaunt factor
levels.c - a new mode to but all levels into ground state - this is an improvement for PL dominated cases but we should really do better
lines.c - some lines added, and then commented out to implement the approximate gaunt factor to line emission as in hazy 2.This is under development
partition.c - made it possible to call partition_functions_2 with a weight, this is to allow more flexibility. For the PL_case, we cann with w=0 to force everything into the ground state
python.c - added ionization options into the question line - set the spectype for the power law to SPECTYPE_POW, rather than defaulting to models, set the innermost stable orbint to 6Rg insead of 12Rg
python.h - added lum_adiabatic into the geo array to store the total adiabatic cooling of the wind, and also kappa_ff_factor into the plasma structure to calculate the gaunt factor for free free heating, mean_ds also put in to see the mean disatnce travlled by a photon packet in a cell. This will allow one to see if ta cell is effectively optically thin...
radiation.c - kappa_ff changed to take account of a gaunt factor. There is also a new subroutine put in to calculate the kappa_ff_factor.


