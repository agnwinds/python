Central_object.rad_type_to_make_wind
====================================
Multi-line description, must keep indentation.

Type
  Enumerator

Values
  bb
    Black-body radiation. The boundary layer radiates as a black-body source with surface luminosity set by its
    effective temperature (:ref:`Central_object.temp`) and resulting in a total luminosity
    proportional to its surface area.

  models
    Radiate according to a model. Python can support tabulated models that output with a binned luminosity distribution
    depending on system properties like temperature and gravity. See :ref:`Input_spectra.model_file`. The total
    luminosity will be set by :ref:`Central_object.luminosity`.

  brems
    Available for :ref:`System_type` of ``agn`` or ``bh``.
    Central object radiates with SED of a brehmsstralung spectrum as $L_\nu=\nu^{\alpha}e^{-h\nu/kT}$.
    This was originally developed to allow comparison to spectra generated
    according to Blondin heating and cooling rates.

  cloudy
    Available for :ref:`System_type` of ``agn`` or ``bh``.
    Central object radiates with a 'broken' power law, intended largely for testing purposes against Cloudy.
    The SED form is $L_\nu=K\nu^\alpha$, but beyond the provided high and low energy
    breakpoints the luminosity falls off sharply.

  power
    Available for :ref:`System_type` of ``agn`` or ``bh``.
    Radiate following a power-law model as $L_\nu=K\nu^\alpha$.
    The total luminosity will be set by :ref:`Boundary_layer.luminosity`.


File
  `setup_star_bh.c <https://github.com/agnwinds/python/blob/master/source/setup_star_bh.c>`_


Parent(s)
  * :ref:`Central_object.radiation`: ``True``


Child(ren)
  * :ref:`Central_object.power_law_cutoff`

  * :ref:`Central_object.bremsstrahlung_alpha`

  * :ref:`Central_object.cloudy.low_energy_break`

  * :ref:`Central_object.bremsstrahlung_temp`

  * :ref:`Central_object.blackbody_temp`

  * :ref:`Input_spectra.model_file`

  * :ref:`Central_object.cloudy.high_energy_break`

  * :ref:`Central_object.luminosity`

  * :ref:`Central_object.power_law_index`

