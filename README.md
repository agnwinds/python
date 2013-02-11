README
***
=========
precursor - python70d (0917 version)
This version fixes the bug that was causing crazy stripes in the data. wind_updates was calling adiabatic_lum with the wrong pointer, and so was getting nonsense results. The previous version fixed the behaviour, but not the root cause.
This version also includes an addition to the convergence valuation which now subtracts the adiabatic luminosity from the heat tot before deciding if the heating/cooling is converging.

There is also some minor tidying up in py_wind and py_wind_sub

