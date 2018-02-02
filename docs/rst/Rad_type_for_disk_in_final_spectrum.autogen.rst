
===================================
Rad_type_for_disk_in_final_spectrum
===================================

Rad_type_for_disk_in_final_spectrum
===================================
If the system contains a disk then the spectrum of the disk is generally simulated either a collection of blackbodies
based on the temperature, as an appropriately weighted set of stellar atomosphers.  One can also simply generate
a flat spectrum, in which case the relative contribution of different annulae is based on the wavelength limited 
luminostiy in each annulus, something that is usually done only for diagnostic purposes.

**Type:** Enum (Int)

**Values:**

0. bb

1. models

2. uniform


**Parent(s):**
  parameter_: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** python.c


