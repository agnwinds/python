README
***
=========
beta

precursor - python_71b
This is very much a work in progress, but represents the first version of the code to contain a multiple temperature saha equation approach to finding ionisation state. The idea is you pick a good temperature to calculate LTE abundances at - i.e. one that is guaranteed to give you abundances in each of two adjacent states, then you apply corrections to take account of the fact you are not in LTE. This is simply the general case of the corrections we have been doing all along. This version does not work properly, and is put here to make sure we all have a common version to debug
The main changes are as follows
ionization.c - added support for two new ionisation modes 6 - correct for dilute blackbody (strictly redundant since LM should do exactly this, but in at the moment for comparisons) and mode 7 - correct for a broken power law.
partition.c - added a new function partition_functions_2 which calculates the partition functions for a pair of ions at a given temperature
power.c - added support for the new mode 7 - it just calls the new code then oneshot.
python.c - support for the new modes at input
saha.c - support for the new modes
zeta.c - changed the zeta code to do a bit more in the routine. It now also works out where in the ground state table you need to be. It is currently only called by the new code.
There is also an entirely new file - variable_temperature.c that contains all the rest of the new code.
There is also also a file test_saha.c - very rough and ready that lets one explore the saha abundances etc for input parameters. THere ? is support in the makefile to allow make test_saha to build it.

