README
***
=========
precursor - python_70a
This version implements dielectronic recombination in the ionization balance calculations but not as a cooling mechanism. There are two new files
zeta.c - this takes the routines for computing zeta, the correction to the saha equation for recombinations not returning to the ground state, combines it with a correction for DR with a switchable mode. This is an attempt to reduce reused code. The sim part of the code now calls this rather than calculating zeta internally. It should also be changed in lucy_mazzali, but not yet done.
dielectronic.c - this is a repository for all the code that will deal with DR. At the moment there is just one routine which takes the coefficients read in and computes a set of rates for a given temperature. * *
also required to get this working is a change to Python/data66/atomic/standard39 : standard39 it just adds a line to refer to the DR file:
the DR file is clist_K.txt
The other main changes in this version are to stuart_sim.c. There was an error in the previous implementaion of the piecewise power law implementation where the only band used was one past the last band! This is fixed, but that in turn showed up problems with doing the frequency check (to see what value of w and alpha should be used) in the integrated function. This made the function discontinuous, and caused major issues with qromb. This check is now outside in xinteg_sim. The code runs.
There are also some more comments here and there.


