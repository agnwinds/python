##############
Boundary_layer
##############

Boundary_layer.radiation
========================
Says whether the boundary layer will radiate.

Type
  Boolean (yes/no)

Parent(s)
  :ref:`System_type`: ``star``, ``cv``


File
  `setup_star_bh.c <https://github.com/agnwinds/python/blob/dev/source/setup_star_bh.c>`_


----------------------------------------

Boundary_layer.rad_type_in_final_spectrum
-----------------------------------------
Determines the luminosity and SED of the boundary layer.
The code can cause a source to radiate differently in the ionisation and spectral cycles.
This variable allows a boundary layer source to radiate differently from :ref:`Boundary_layer.rad_type_to_make_wind`
during the cycles used to calculate the output spectra. This can be

Type
  Enumerator

Values
  One of the following:

  bb
    Black-body radiation. The boundary layer radiates as a black-body source with surface luminosity set by its
    effective temperature (:ref:`Boundary_layer.temp`) and resulting in a total luminosity
    proportional to its surface area.

  models
    Radiate according to a model. Python can support tabulated models that output with a binned luminosity distribution
    depending on system properties like temperature and gravity. See :ref:`Input_spectra.model_file`. The total
    luminosity will be set by :ref:`Boundary_layer.luminosity`.

  uniform
    Available for :ref:`System_type` of ``star`` or ``cv``.
    Photons are generated with a random, uniformly-distributed wavelength between
    :ref:`Spectrum.wavemin` and :ref:`Spectrum.wavemax`. Can in some cases substitute for a Kurcuz spectrum.
    This mode is only available when generating final spectra.

Parent(s)
  :ref:`Boundary_layer.radiation`: ``True``

File
  `python.c <https://github.com/agnwinds/python/blob/dev/source/python.c>`_


----------------------------------------

Boundary_layer.rad_type_to_make_wind
------------------------------------
Determines the luminosity and SED of the boundary layer.
The code can cause a source to radiate differently in the ionisation and spectral cycles.
This variable allows a boundary layer source to radiate differently from :ref:`Boundary_layer.rad_type_in_final_spectrum`
during the cycles used to calculate the wind ionisation state and temperature.

Type
  Enumerator

Values
  One of the following:

  bb
    Black-body radiation. The boundary layer radiates as a black-body source with surface luminosity set by its
    effective temperature (:ref:`Boundary_layer.temp`) and resulting in a total luminosity
    proportional to its surface area.

  models
    Radiate according to a model. Python can support tabulated models that output with a binned luminosity distribution
    depending on system properties like temperature and gravity. See :ref:`Input_spectra.model_file`. The total
    luminosity will be set by :ref:`Boundary_layer.luminosity`.

  power
    Radiate following a power-law model as :math:`L_\nu=K\nu^\alpha`. The total luminosity will be set by :ref:`Boundary_layer.luminosity`.

Parent(s)
  :ref:`Boundary_layer.radiation`: ``True``

File
  `setup_star_bh.c <https://github.com/agnwinds/python/blob/dev/source/setup_star_bh.c>`_


----------------------------------------

Boundary_layer.luminosity
^^^^^^^^^^^^^^^^^^^^^^^^^
The luminosity of the boundary layer.

Type
  Double

Unit
  ergs/s

Values
  Greater than 0

Parent(s)
  :ref:`Boundary_layer.rad_type_to_make_wind`: ``models``, ``power``

  :ref:`Boundary_layer.rad_type_in_final_spectrum`: ``models``, ``uniform``


File
  `setup_star_bh.c <https://github.com/agnwinds/python/blob/dev/source/setup_star_bh.c>`_


----------------------------------------

Boundary_layer.power_law_cutoff
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This is a low frequency cutoff for an AGN-style power law spectrum
of a form :math:`L_\nu=K\nu^\alpha`, as applied to the boundary layer of a star.
It prevents the power-law being applied to low frequencies and giving an odd SED.
See :ref:`Radiation_types` and :ref:`Boundary_layer.power_law_cutoff`.

Type
  Double

Unit
  Hz

Values
  Greater than 0

Parent(s)
  :ref:`Boundary_layer.rad_type_to_make_wind`: ``power_law``


File
  `setup_star_bh.c <https://github.com/agnwinds/python/blob/dev/source/setup_star_bh.c>`_


----------------------------------------

Boundary_layer.power_law_index
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The exponent 𝛼 in a power law SED applied to an AGN-style power law source for a non-AGN system.
central source of the form :math:`L_\nu=K\nu^\alpha`.

See :ref:`Radiation_types` and :ref:`Central_object.power_law_index`.

Type
  Double

Values
  Any - but sign is not assumed, so for negative index use a negative value

Parent(s)
  :ref:`Boundary_layer.rad_type_to_make_wind`: ``power_law``


File
  `setup_star_bh.c <https://github.com/agnwinds/python/blob/dev/source/setup_star_bh.c>`_


----------------------------------------

Boundary_layer.temp
^^^^^^^^^^^^^^^^^^^
The temperature of the boundary layer when radiating as a black body.

Type
  Double

Unit
  Kelvin

Values
  Greater than 0

Parent(s)
  :ref:`Boundary_layer.rad_type_to_make_wind`: ``bb``

  :ref:`Boundary_layer.rad_type_in_final_spectrum`: ``bb``


File
  `setup.c <https://github.com/agnwinds/python/blob/dev/source/setup.c>`_


