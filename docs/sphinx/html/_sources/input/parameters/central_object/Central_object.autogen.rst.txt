##############
Central_object
##############

Central_object.radiation
========================
A booliean variable stating whether of not the central object should radiate.
This will enable different follow-up questions depending on the system type.

**Type:** Boolean (yes/no)

**File:** setup_star_bh.c


Central_object.geometry_for_source
----------------------------------
If the central source in an AGN/BH system is radiating, what geometry should it radiate from?
 This is applicable even for black-body sources, where the *luminosity* depends on :ref:`Central_object.radius`.

**Type:** Enumerator

**Values:**

lamp_post
  The central source radiates from two point sources
  located on the system axis above and below the disk plane.
  Emission is completely isotropic.

sphere
  The central source radiates from a spherical surface with radius :ref:`Central_object.radius`.
  Emission is cosine-weighted in the direction of the sphere normal at the point of emission.
  Photons that would be spawned in an extended disk (if :ref:`Disk.type` is `vertically.extended`)
  are re-generated.


**Parent(s):**

* :ref:`System_type`: ``agn``, ``bh``

* :ref:`Central_object.radiation`: ``True``


**File:** setup_star_bh.c


Central_object.lamp_post_height
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The distance above and below the disk plane of the two point sources used in the lamp-post model.

**Type:** Double

**Unit:** :ref:`Central_object.radius`

**Values:** Greater than 0

**Parent(s):**

* :ref:`Central_object.geometry_for_source`: ``lamp_post``


**File:** setup_star_bh.c


Central_object.rad_type_in_final_spectrum
-----------------------------------------
Determines the SED of the central object in the spectral cycles. The luminosity is set by the options for the
ionisation cycles, however.

**Type:** Enumerator

**Values:**

bb
  Available for :ref:`System_type` of `star` or `cv`.
  Black-body radiation. The boundary layer radiates as a black-body source with surface luminosity set by its
  effective temperature (:ref:`Central_object.temp`) and resulting in a total luminosity
  proportional to its surface area.

models
  Available for :ref:`System_type` of `star` or `cv`.
  Radiate according to a model. Python can support tabulated models that output with a binned luminosity distribution
  depending on system properties like temperature and gravity. See :ref:`Input_spectra.model_file`. The total
  luminosity will be set by :ref:`Central_object.luminosity`.

uniform
  Available for :ref:`System_type` of `star` or `cv`.
  Photons are generated with a random, uniformly-distributed wavelength between
  :ref:`Spectrum.wavemin` and :ref:`Spectrum.wavemax`. Can in some cases substitute for a Kurcuz spectrum.
  This mode is only available when generating final spectra.

brems
  Available for :ref:`System_type` of `agn` or `bh`.
  Central object radiates with SED of a brehmsstralung spectrum as $L_\nu=\nu^{\alpha}e^{-h\nu/kT}$.
  This was originally developed to allow comparison to spectra generated
  according to Blondin heating and cooling rates.

cloudy
  Available for :ref:`System_type` of `agn` or `bh`.
  Central object radiates with a 'broken' power law, intended largely for testing purposes against Cloudy.
  The SED form is $L_\nu=K\nu^\alpha$, but beyond the provided high and low energy
  breakpoints the luminosity falls off sharply.

power
  Available for :ref:`System_type` of `agn` or `bh`.
  Radiate following a power-law model as $L_\nu=K\nu^\alpha$.
  The total luminosity will be set by :ref:`Boundary_layer.luminosity`.


**Parent(s):**

* :ref:`Central_object.radiation`: ``True``


**File:** python.c


Central_object.temp
===================
Temperature of the central star. Physically, this is used in blackbody radiation, shock heating and disk heating in
YSO models. It is also used to help determine the frequency bands in which photons are emitted.

**Type:** Double

**Unit:** Kelvin

**Values:** Greater than zero

**Parent(s):**

* :ref:`System_type`: ``star``, ``cv``


**File:** setup_star_bh.c


Central_object.rad_type_to_make_wind
====================================
Multi-line description, must keep indentation.

**Type:** Enumerator

**Values:**

bb
  Black-body radiation. The boundary layer radiates as a black-body source with surface luminosity set by its
  effective temperature (:ref:`Central_object.temp`) and resulting in a total luminosity
  proportional to its surface area.

models
  Radiate according to a model. Python can support tabulated models that output with a binned luminosity distribution
  depending on system properties like temperature and gravity. See :ref:`Input_spectra.model_file`. The total
  luminosity will be set by :ref:`Central_object.luminosity`.

brems
  Available for :ref:`System_type` of `agn` or `bh`.
  Central object radiates with SED of a brehmsstralung spectrum as $L_\nu=\nu^{\alpha}e^{-h\nu/kT}$.
  This was originally developed to allow comparison to spectra generated
  according to Blondin heating and cooling rates.

cloudy
  Available for :ref:`System_type` of `agn` or `bh`.
  Central object radiates with a 'broken' power law, intended largely for testing purposes against Cloudy.
  The SED form is $L_\nu=K\nu^\alpha$, but beyond the provided high and low energy
  breakpoints the luminosity falls off sharply.

power
  Available for :ref:`System_type` of `agn` or `bh`.
  Radiate following a power-law model as $L_\nu=K\nu^\alpha$.
  The total luminosity will be set by :ref:`Boundary_layer.luminosity`.


**File:** setup_star_bh.c


Central_object.power_law_cutoff
-------------------------------
Adds a low-frequency cutoff to the power law spectrum.
Whilst this is required for power-law emission modes,
it's set globally and also used in `cloudy` broken-power-law emission modes!

**Type:** Double

**Unit:** Hz

**Values:** Greater than or equal to 0

**Parent(s):**

* :ref:`Central_object.rad_type_to_make_wind`: ``power``


**File:** setup_star_bh.c


Central_object.bremsstrahlung_alpha
-----------------------------------
The frequency exponent 𝛼 in bremstrahlung SED of the form
$L_\nu=\nu^{}\alpha}e^{-h\nu/kT}$

**Type:** Double

**Values:** Any - sign is not assumed so use negative if you want negative

**Parent(s):**

* :ref:`Central_object.rad_type_to_make_wind`: ``brems``


**File:** setup_star_bh.c


Central_object.cloudy.low_energy_break
--------------------------------------
This is a command to define a cloudy type broken power
law SED - mainly used for testing the code against cloudy.
This SED has hardwired frequency exponents of 2.5 below the
low energy break and -2.0 above the high energy break. This
parameter defines the energy of the low energy break.

**Type:** Double

**Unit:** eV

**Values:** Greater than 0

**Parent(s):**

* :ref:`Central_object.rad_type_to_make_wind`: ``cloudy``


**File:** setup_star_bh.c


Central_object.bremsstrahlung_temp
----------------------------------
The temperature T in bremstrahlung SED of the form
$L_\nu=\nu^{\alpha}e^{-h\nu/kT}$

**Type:** Double

**Unit:** K

**Values:** Greater than 0

**Parent(s):**

* :ref:`Central_object.rad_type_to_make_wind`: ``brems``


**File:** setup_star_bh.c


Central_object.blackbody_temp
-----------------------------
If the AGN/BH is radiating as a black body, what temperature should it radiate at?

**Type:** Double

**Unit:** Kelvin

**Values:** Greater than 0

**Parent(s):**

* :ref:`System_type`: ``agn``, ``bh``

* :ref:`Central_object.rad_type_to_make_wind`: ``bb``


**File:** setup_star_bh.c


Central_object.cloudy.high_energy_break
---------------------------------------
This is a command to define a cloudy type broken power
law SED - mainly used for testing the code against cloudy.
This SED has hardwired frequency exponents of 2.5 below the
low energy break and -2.0 above the high energy break. This
parameter defines the energy of the high energy break.

**Type:** Double

**Unit:** eV

**Values:** Greater than :ref:`Central_object.cloudy.low_energy_break`

**Parent(s):**

* :ref:`Central_object.rad_type_to_make_wind`: ``cloudy``


**File:** setup_star_bh.c


Central_object.luminosity
-------------------------
The luminosity of a non-blackbody AGN central source.
This is defined as the luminosity from 2-10keV.

**Type:** Double

**Unit:** ergs/s

**Values:** Greater than 0.

**Parent(s):**

* :ref:`System_type`: ``agn``, ``bh``

* :ref:`Central_object.rad_type_to_make_wind`: ``brems``, ``cloudy``, ``model``, ``power``


**File:** setup_star_bh.c


Central_object.power_law_index
------------------------------
The exponent 𝛼 in a power law SED applied to a power law source of the form $L_\nu=K\nu^\alpha$.

See :ref:`Radiation_types` and :ref:`Boundary_layer.power_law_index`.

**Type:** Double

**Values:** Greater than 0

**Parent(s):**

* :ref:`Central_object.rad_type_to_make_wind`: ``cloudy``, ``power``


**File:** setup_star_bh.c


Central_object.mass
===================
Mass of the central object. This is very important, affecting wind speeds, gravitational heating and such.

**Type:** Double

**Unit:** M☉

**Values:** Greater than 0

**File:** setup_star_bh.c


Central_object.radius
=====================
Radius of the central object in the system, e.g the white dwarf or black hole

**Type:** Double

**Unit:** cm

**Values:** Greater than 0

**File:** setup_star_bh.c


