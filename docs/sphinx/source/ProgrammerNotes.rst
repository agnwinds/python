Programming Notes
#################

Python is written in C and is normally tested on linux machines and on Macs, where the
compiler usually turns out to be clang.  Certin portions of the code are parellized 
using mpi. 

Version control is (obviously) managed through git.  The stable version is on the master 
branch; the main development is carried out on the dev branch. This is generally the 
branch to start with in developing new code.

If one modifies the code, a developer needs to be sure to have $PYTHON/py_progs both in PYTHONPATH and PATH.  One should also have a version of indent installed, preferably but not 
necessarily gnu_indent installed.  This is because, the Makefile will call a script 
run_indent.py on files that the developer has changed, which enforces a specific indent 
style on the code.

In addition to indent, one should have cproto or something equivalent installed. cproto is used to prototypes for all of the subroutines (through the make command 

* make prototypes

(The many warnings that appear when cproto is run on OSX can so far be ignored. cproto for 
macs is available with brew)

All new routines should have Doxygen headers.

printf statements should be avoided, particular in the main code.  Python has 
its own replacements for these commands, e.g Log and Error which standardize the outputs
and allow for managing what is printed to the screen in multprocessor mode.

Structures
==========

In order to understand the code, one needs to understand the data structures.  

The main header files  are:

* atomic.h - This contains all of the structures that hold atomic data, e.g oscillator 
  strengths, photoionization cross-sections, elemental abundances, etc.  These data are 
  read in at the beginning of the program (see atomicdata.c and other similarly named 
  routines)
* python.h - This contains the structures and other data that comprise the wind as well 
  as the parameters of the model.  (This is fairly well-documented, or should be)


Program Flow
============

Basically, (as decribed from a more scientific perspective elsewhere), the program consists
of a number stages

* Data gathering and initialization: This consiss of reading in all of the parameters 
  for the model to be calculated, reading in the associated atomic data, and setting up 
  program to run.  This procuess involves allocating space for many of the data structures.
* Ionization cycles: During this portion of the program fleets of photons are generated 
  and propogated through the wind, interacting with it in various ways. These photons are
  generated over a large range of frequencies, because their purpose is to allow the program
  to determine the ionization state of the wind.  During this 
  process various estimators are accumulated that describe the interaction of the photons
  with the wind.  Once all of the photons have propagated through the wind the various 
  estimators are used to calculate a new estimate of the ionization state of the various
  wind cells that constitute the wind.  This process is repeated for a number of cycles, 
  by which time, hopefully, the wind will have reached a "steady state".
* Spectral cycles: Once the ionization cycles have been completed, the ionization state 
  of the wind is fixed, and more detailed spectra are calculated. For this, photons are generated
  in a limited spectral range, depending on the interests of the user.  In contrast to
  the ionization state, where "cycles" are  crucial concept, the only reason to have spectrl
  cycles in the "Spectral cycle" phase is to allow one to gradually improve the statistics 
  of the spectrum.  At the end of each spectral cycle, the detailed spectra are written out, 
  so one can see the spectra building up.


Parallel Operation
==================

Python uses MPI to parallelize the most compute intensive portions of the routine.  It has
been run on large machines with 100s of cores without problem.

The portions of the routine that are parallized are:

* Photon generation and transfer: When run in multiprocesser mode, each thread creates only a 
  fraction of the total number of photons.  The weight of the photons in each thread is such
  that the sum of the weights is the total energy expected to be produced in one observer frame second.
  These photons are propagated through the wind, and estimators based on these photons are accumulated.
  At the end of photon transfer by all threads, the various quantities, including the spectra,  that 
  have been accumulated in the separate threads are gathered together and averaged or summed as 
  appropriate.  For ionization cycles, this means that all of the data needed to calculate the
  ionization in any cell is available on each of the threads.
* Ionization calculation:  Although all of the threads have all of the data needed to calculate
  the ionization in any cell, in practice what happens is that the program assigns a different set of
  cells to each thread to calculate the ionization.  After the thread calculates the new ionization 
  state for its assigned cells, the ionization states are then gathered back and broadcast to all
  of the threads, in preparation for the next cycle.  
* Preparation for detailed radiative transfer in the macro-atom mode.  When photons go through the
  grid in the simple-atom mode, photon frequencies do not change a great deal, however in macro-atom 
  mode the frequencies can shift by large amounts. This creates a problem during the detailed spectral 
  generation stage, because one does not know before hand how many photons that started out of the 
  desired band end up in the desired band.  (In macro-atom mode, a photon bundle that starts out at
  8 keV photoionizes an Fe ion can turn (for example) into an Hbeta photon).  To handle this, one 
  needs to estimate how often this happens and include this (effectively as a source function) in 
  radiative transfer involving macro-atoms. This is parallelized, in the same manner as the ionization
  calculation by assigning various cells to various threads and gathering the results back before 
  the radiative transfer step in the detailed spectrum phase.
