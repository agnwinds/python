README
***
=========
precursor python_72b
This version contains the experimental code implementing the cardona partition functions and also tracking of photons passing through a given cell and contributing to its spectrum
The main changes are:
atomic.h - a new structure cardon_partition which contains the parameters to calculate the partition functions. There is also an extra member of the nion structure (nxcpart)which points to the correct cardona data. There is also a flag in the ion structure (cpartflag) which says wether or not there are cardona partition functions for this ion.
bands.c - the bands for the photon tracking are set up here for full runs. They need to be changed to agree with the table command in cloudy!
diag.c - lines to check for a file "diag_cells.dat" if we are in diagnostics mode. If the file exists, read in a list of cells for which we are going to write out detailed photon stats.
get_atomicdata.c - lines to read in cardona partition functions. The pointer to this type of data is CPART
partition.c - two new functions, cardona_part_func and cardona_part_func_2 which calculate the partition functions for the given cell (first case) and a pair of ions (second case - for variable temperature). At the moment, this function is separate from the normal partition functions. If we decide it works. it would make sense to incorporate it into the normal partition functions function, choosing this is data is available, or defaulting to the weighted BB form if not. At the moment, if you call cardona_part_func and have no data for a given ion - it simply uses the ground state.
python.h 
- a new file pointer (pstatptr) to hold photon data, a flag (cell_phot_stats) to say wether we are logging photons, an array (ncell_stats[NCSTAT] ? ) which holds the cells we are examining and ncstat - the number of cells we are looking at, up to a maximum of 10.
radiation.c - lines to write out photon stats


