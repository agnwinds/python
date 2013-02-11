README
***
=========
precursor - python_71
One change:
in power_sub.c, the four arrays used to report how the code is producing the sim correction factor (fudge_store, num_store, denom_store and ion_store) are set to have dimension NIONS rather than 300. This is to prevent overflows that were probably happening for ages, but not seen because they were not overflowing into anything important.

