README
***
=========
precursor - python_70i
The main changes are as follows:
elvis.c, photond2d.c and wind.c - Changes were made to the way in which the elvis model is implemented so the offset and the windcone are properly captured. There should no longer be cells that are said to be in the wind with extremely low densities.
py_wind - Changes were made so that the program can be run from the command line using a rdpar file, as well as in the traditional fashion. The goto statement which is quite bad practive has been eliminated. Some new options intended to help with understanding convergence have been added.
py_wind_write.c - The output file format has been changed to include lines for whether a cell is in the wind or not, and to give the i j cell number. This goes with a new version of a plotting routine which Knox has written. Improvements have been made to the way 1-d, that is spherical models are print out.
ionization.c - Code for the power law model has been moved into its own subroutine power.c. Also stuart_sim.c has been renamed power_sub.c. Ultimately we should rename things like sim_alpha to something that doesnot includes Stuart's name, but his has not been done and is not necessary until we decide where we are going with all of that.
Dielectronic recombination has been commented out
bands.c - The code which establishes the bands for calculating the power law factors have been moved here from python.c, and some changes have been made to make sure that the bands for which the power law is calculated do not extend beyond the bands in which photons are generated.
python.h - the parameters that are associated with setting up the bands for the power law spectrum, e.g nxbands and xfreqs have been moved into the geo structure (so they are accessible from py_wind)
Many comments have been added or modified. We need to be more consistent on how we do these.
Note that some routines have been processed with indent in order


