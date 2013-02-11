README
***
=========
precursor - python_70g (70h got skipped)
emission.c - statement added to allow wind emission from torus to work - probably - not really very well tested.
ionization.c
- disk and agn protons turned into floats to allow actual photon numbers rather than packet numbers
compton and DR cooling subtracted here from the heating of the cell to avoid them going into lum_tot
counters of the three separate convergence criteria zeroed - they get automatically zeroed on the mac but I saw errors when running under linux
photon_gen.c 
- added the weight of photons being created into the reporting - lets you see wether banding is working OK
py_wind_sub -
took account of compton and DR cooling not being in the lum_tot
put in lines to write out the new IPs ?
python.c
lum_ioniz and n_ioniz calls now go to geo structure
put in a call to a new function wind_ip to calculate a simple ionisation parameter
python.h
put lum_ioniz and n_ioniz into the geo structure
changed the photon counters to floats
added ferland_ip and ip to the plasma structure to compute IPs ?
radiation.c
added lines to compute the IP from the photons
added lines to increment photon counters by nphot/energy to give actual numbers rather than packets.
resonate.c - took a print line out
stuart_sim.c - made it so all the diagnostic files only get generated in diag is set to 1


