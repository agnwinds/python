
=============
Input_spectra
=============

----------------------------------------

Input_spectra.model_file
========================
In addition to being able to generate several types of spectra, such
as blackbodies and power laws, Python can read in a series of spectra
which are tabulated and are describable in terms of (usually) temperature
and gravity). This parameter defines the name of the file which gives the
location of the individual spectra and the temperate and gravity associated
with each spectrum. (One may wish to use the same files for several radiation sources, viz the disk and the star)
Python actually only reads in the data the first time.

**Type:** String

**Parent(s):**
  :ref:`parameter`: One is asked this question whenever one chooses this option, regardless of the radiation source.


**File:** setup.c


