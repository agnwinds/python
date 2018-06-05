
=====
Hydro
=====

Hydro.thetamax
==============
The maximum theta value to be read in from a hydrodynamic snapshot.
This is typically used to excise a dense disk from the midplane of
such a snapshot. Use a negative value to tell the code to use all
the data.

**Type:** Double

**Unit:** Degrees

**Values:**

X. use up to that angle (typically less than 90)

-1. use all data


**Parent(s):**
  parameter_: None


**File:** hydro_import.c


Hydro.file
==========
Relative path to a hydrodynamic snapshot file to be imported.

**Type:** String

**Parent(s):**
  parameter_: None


**File:** hydro_import.c


