README
***
=========
Precursor - python_68d
Changes
Simply makes modifications that were required for portability to 64bit linux platform. These all have to do with vfprintf. The files that were changed are mostly in the kpar routines but signal.c was also affected.
Note that there may be issues when if one tries to run py_wind on a wind_save file that has been created with a different architecture. This has not yet been fully explored.


