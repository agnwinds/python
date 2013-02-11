README
***
=========
precursor python_71c (71d's were various test builds)
This version seems stable, and the blackbody variable temperature mode (6) agrees with the lucy mazzali mode (3) although it is much slower due to many integrations. The power law mode still needs more testing.
The main changes are
atomic-h - a new variable in the ion structure called topbase_ground - the ground state level for the ion
get_atomicdata.c - some new lines to work out which level is the ground state, and assign it correctly to topbase_ground.
variable_temperature.c - several changes, all commented to make it work better
wind_updates2d.c - minor changes to write out some more details.


