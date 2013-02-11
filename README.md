README
***
=========
Precursor is python_74a_ksl 
Changes this version
python.c - changed call to pop_kappa_ff_array (); so it only occurs if there is data to allow it to work!
emission.c - changed ff and total_free to allow proper operation if gaunt factor data has not been read in
radiation.c - changed kappa_ff to allow proper operation if no gaunt factors.
variable_temperature.c - changed code so density is computed in a copy of the actual array, and changes to the plasma structure only made on proper completion. If the loop searching for n_e fails, then the old densities are retained.


