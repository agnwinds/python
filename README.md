README
***
=========
precursor - python_71a
The main changes are as follows
python.c - Modification to make sure the dfudge parameter is set correctly with the new shell wind description. There is still a longer term issue about wether this could be done better and more generally
shell_wind.c - various modifications made so that the shell wind uses stellar wind routines wherever possible to reduce the amount of duplicated code.
spherical.c - slight changes to log some data rather than printf it
stellar_wind.c - removed references to shell code that had now been removed.
wind.c - changed calls that used to go to shell wind routines to call stellar wind routines.


