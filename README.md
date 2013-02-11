README
***
=========
precursor python_71e - 72a's were lots of test versions.
I'm posting this version because I'm running some test cases on Iridis using it, and I want to make sure it remains a stable version. It will be changed almost immediately into 72b2.
The main changes are:
agn.c 
- a new subroutine emittance_bpow - it produces photons from a broken power law - this was implemented to allow python to have precisely the same illuminating spectrum as the cloudy command

table
atomic.h 
- added a new constant - planks constant in Rydbergs - we are now using cloudy lots, and its default energy uint is rydbergs.
bands.c 
- a new mode - mode 5 has been included which attempts to ensure band of the spectrum estimators are alignedd with the locations of the broken power law breaks.
compton.c - implemented the induced compton scattering
photon_gen.c - some changes to support the broken power law mode
power.c - a few changes to what is written out - includes energies with band limits, and tells you what band number you are looking at
python.h - some new elements in the plasma structure to hole induced compton heating
radiation.c - some lines to zero induced compton heating, and write it out
test_saha.c lots of changes, but this is a test code...


