
====
Band
====

Band.boundary
=============
The boundaries between bands of energies at which photons will be generated

**Type:** Double

**Unit:** eV

**Value:** Greater than 0

**Parent(s):**
  parameter_: None


**File:** bands.c


Band.minimum_fraction)
======================
The minimium fraction of photons that should be made per band

**Type:** Double

**Unit:** None

**Value:** Greater than 0 and should sum to less than 1 over all bands

**Parent(s):**
  parameter_: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** bands.c


