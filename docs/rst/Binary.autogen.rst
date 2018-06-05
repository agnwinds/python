
======
Binary
======

Binary.period
=============
The perids of a binary system. Along with a mass, the binary period is 
used to define the Roche lobe of the system, which in turn can be used
to see the effect of eclipses on the spectrum.  Defining the system as
a secondary also initializes the outer radius of the disk.

**Type:** Double

**Unit:** Hours

**Value:** greater than 0

**Parent(s):**
  parameter_: Required whenever the system is defined to be a binary system


**File:** setup_star_bh.c


Binary.mass_sec
===============
In binary systems the mass of the secondary. This is used along
with the period to establish the Roche lobes, so that one can
see the effects of eclipses on the system

**Type:** Double

**Unit:** msol

**Value:** greater than 0

**Parent(s):**
  parameter_: Required whenever the system is a binary


**File:** setup_star_bh.c


