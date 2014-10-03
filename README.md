# Python

Python is a (confusingly named) Monte Carlo radiative transfer code which uses the Sobolev approximation.
It has been developed by Knox Long, Christian Knigge, Stuart Sim, Nick Higginbottom and James Matthews. 

The code is not yet ready for use and should not be installed yet. If you are interested in using
Python please contact long@stsci.edu.

## Travis Build Status

[![Build Status](https://travis-ci.org/agnwinds/python.png?branch=dev)](https://travis-ci.org/agnwinds/python)


## Getting the radiative transfer code 'Python'

You can download the required structure under the structure branch. e.g.
git clone https://github.com/agnwinds/python.git -b structure
or simply click on the 'zip' button!


Releases of progs can be found under [tags](https://github.com/agnwinds/python/tags "Wiki").

Consult the [wiki](https://github.com/agnwinds/python/wiki/_pages "Wiki") for how to install Python.


## Installation

Python and the various routines associated are set up in a self-contained directory structure. The basic directory structure and the data files that one needs to run Python need to be retrieved and compiled. 

If you want to obtain a stable (!) release, go to the [Releases](https://github.com/agnwinds/python/releases) page.

If you want to download the latest dev version, you can zip up the git repository by clicking on the zip icon to the right of the GitHub page. Aternatively, you can clone the repository using 

    $ git clone https://github.com/agnwinds/python.git 

If you anticipate contributing to development we suggest Forking the repository and submitting pull requests with any proposed changes.

Once you have the files, you need to cd to the new directory and set your environment variables
    
    $ export PYTHON = /path/to/python/
    $ cd $PYTHON 
    $ ./configure
    $ make install
    $ make clean

note that export syntax is for bash- for csh use 
  
    $ setenv PYTHON /path/to/python/

Atomic data is stored in our [data repository](https://github.com/agnwinds/data) with it's own releases page. one should unzip these files and place them in the $PYTHON/data folder.

## Running python

To run python you need to add the following to your $PATH variable:

    $PYTHON/bin

You can then setup your symbolic links by running 

    $ Setup_Py_Dir

and run the code by typing, e.g.

    $ py root.pf


Please see the (wiki)[https://github.com/agnwinds/python/wiki/Installing-and-Running-Python] and docs folder for how to use the code.

Any comments, email jm8g08@soton.ac.uk or long@stsci.edu.