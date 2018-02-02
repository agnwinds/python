
=======
Surface
=======

Surface.reflection.or.absorption
================================
When photons hit the disk, there are several options 

**Type:** Enum (Int)

**Values:**

0. no.rerad - the photons are simply lost from the calculation

1. high.albedo - the photons are scatttered back into the wind

2. thermalized.rerad - The photons are absorbed, in the next ionization cycle energy lost is treated as extra heat, and the effective temperature of the ring in the disk will be increased accordingly


**Parent(s):**
  parameter_: This question is asked whenever there is a disk


**File:** setup2.c


