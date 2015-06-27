140303 - ksl
This is a grid of models calculated for use tiwht Python.  The models are the
same modles as in the disk06 grid. They use Kurucz models for the basis at 
low temperature and Hubeny models at higher temperatures (as soon as it became
possible to calculate Hubeny models).  All of the atmospheres have spectra calculated
with SysnSpec.  The atomic data used to create the spectra was provided to me by
Ivan includes more lines, lower wavelenghts than gfall.dat, the kurucz list used 
to calculate the previous grid.  The final models have been convlved to 
an instrumental profile with a resultion of 0.5 A.  If one is concerned with spectra 
with a resulution less than this, then rotin needs to be rerun.

To use this grid NWAVES needs to be set to 28000
