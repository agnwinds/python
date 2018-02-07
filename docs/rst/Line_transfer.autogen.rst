
=============
Line_transfer
=============

Line_transfer
=============
The way in which line transfer and scattering is dealt with 
in the code. Governs whether we adopt any approximations
for radiative transfer, whether to use the indivisible packet 
and macro-atom machinery, and whether to use isotropic or 
anisotropic scattering.

Mode 5 is recommended for non-macro atom runs,
while mode 7 is recommended for macro-atom runs.

**Type:** Enum (Int)

**Values:**

0. *pure.abs*
   
   The spure absortion approximation.

1. *pure.scat*
   
   The spure scattering approximation.

2. *sing.scat*
   
   The single scattering approximation.

3. *escape.prob*
   
   Resonance scattering and electron scattering is dealt with isotropically.
   free-free, compton and bound-free opacities attenuate the weight of the photon
   wind emission produces additional photons, which have their directions chosen isotropically. 
   The amount of radiation produced is attenuated by the escape probability.

4. *escape.prob, approximate anisotropic scattering*
   
   Resonance scattering and electron scattering is dealt with using the semi-analytic isotropic scattering routine.
   free-free, compton and bound-free opacities attenuate the weight of the photon
   wind emission produces additional photons, which have their directions chosen isotropically. 
   The amount of radiation produced is attenuated by the escape probability.

5. *escape.prob, anisotropic scattering*
   
   as mode 4, but we use 
   the 'thermal trapping method' to choose an 
   anistropic direction when an r-packet deactivation 
   or scatter occurs.

6. *macro_atoms*
   
   use macro-atom line transfer.
   Packets are indivisible and thus all opacities are dealt with by activate a macro-atom, scattering, 
   or creating a k-packet.
   the new direction following electron scattering or deactivation of 
   a macro atom is chosen isotropically.

7. *macro_atoms+aniso.scattering*
   
   as mode 6, but we use the 'thermal trapping method' to choose an anistropic direction 
   when an r-packet deactivation or scatter occurs.

8. *macro_atoms, anisotropic Scattering, All simple ions*
   
   as mode 7, but we treat all ions as simple

9. *macro_atoms, anisotropic Scattering, All simple ions*
   
   as mode 8, but we use semi-analytic scattering mode


**Parent(s):**
  parameter_: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** setup_domains.c


