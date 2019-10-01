Model
#####

As Python is run, it repeatedly writes out two binary files that contain essentially all information about the wind as calculated in the ionization phase of the program,
along with status of the program at the last point where the file was written.
These files along with the parameter file are sufficient to restart the program,
if for example, one wants to check point the program after a certain time, and restart where one left off,
or to add spectral cycles to get better spectra.

.wind_save
  A binary file that contains essentially all information about the wind including ion densities,
  temperatures, and velocities in each cell, along with status of the program at the last point where the file was written.

.spec_save
  A binary file that contains all of the information about the spectra that have created.  This file is not of interest to users directly.  It is used when restarting

Two routines exist as part of the Python distribution allow the user to gain insight into the actual model

windsave2table
  Executed from the command line with :code:`windsave2table rootname`.

  Produces a set of standard set ascii tables that that show for each grid cell quantities such as wind velocity,
  :math:`n_e`, temperatures, and densities of prominent ions.

py_wind
  Executed from the command line with :code:`py_wind rootname`

  Allows the user to query for information about the model interactively.  The results can be written to ascii files for future reference
