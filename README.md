README
***
=========
precursor - python-70c
This version fixes a basic problem with zero emit needed to incorporate adiabatic cooling into the code. Previously the adiabatic cooling did not change as the a one searched for a temperature that would balance heating and cooling. But adiabatic cooling is proportional to t_e and so it should change as you deal with this. The main changes are in the routine ionization.c
Additional minor changes were made to other routine, especially py_wind_subs.c. Various comments were also added to portions of the code.
This needs to be integrated into Nick's version of the code. Note that I did not work on the question of whether the compton torus is actually doing what we want, and so problems with this may remain.


