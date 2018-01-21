
========
Spectrum
========

Spectrum.angle
==============
The inclination angle with respect to the polar axis for
obtaining a spectrum.  This question will be repeated once
for each desired incliniation

**Type:** Double

**Unit:** Degrees

**Value:** Normally betwween 0 and 99 degrees

**Parent(s):**
  parameter_: Spectrum.no_observers


**File:** setup2.c


Spectrum.live.or.die
====================
Normally in creating detailed spectrum Python "extracts" photons in a certain
direction rewithing them to account for the fact that they have been extracted
in a certain direction.  It is possible to just count the photons that are emitted
in a single angle range. The two methods should yield the same or very similar results 
but the extraction method is much more efficient and live or die is basically a 
diagnostic mode.  For historical reaaons the live or die method is called with 0
and anything else rsults in the standard extract method being used.

**Type:** Boolean (1/0)

**Parent(s):**
  parameter_: Called anytime detialed spectra are two be created.


**File:** setup2.c


Spectrum.orbit_phase
====================
For binary systems, the orbital phase at which the spectrum
is to be extracted (so the effects of an eclipse can be taken
into account in creating the spectrum. Phase 0 corresponds to
inferior conjunciton, that is with the secondary in front (or 
depending on inclination angle, partially in front of) the
primary

**Type:** Double

**Unit:** None

**Value:** Normally between 0 and 1

**Parent(s):**
  parameter_: Only required when the system is a described as a binary


**File:** setup2.c


Spectrum.type
=============
The type of spectra that are produced in the final spectra. The current choices are flambda, fnu, or basic,
where basic implies simply summmung up the energy packets that escape within a particularly wavelength/
frequency bin..

**Type:** Enum (Int)

**Values:**

1. flambda

2. fnu

other. basic


**Parent(s):**
  parameter_: Called whenever detailed spectra are generated.


**File:** setup2.c


Spectrum.wavemax
================
The maximum waveleenght of the detailed spectra that are to be produced

**Type:** Double

**Unit:** Angstroms

**Value:** Greater than 0 and greater than Spectrum.wavemin

**Parent(s):**
  parameter_: Condition e.g. >0 or list e.g. [1, 2, 5]


**File:** setup2.c


Spectrum.wavemin
================
The minimum wavelength of the final spectra in Angstroms

**Type:** Double

**Unit:** Angstroms

**Value:** Greater than 0

**Parent(s):**
  parameter_: Condition e.g. >0 or list e.g. [1, 2, 5]


**File:** setup2.c


