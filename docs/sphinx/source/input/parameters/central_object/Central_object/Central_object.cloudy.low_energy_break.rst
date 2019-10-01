Central_object.cloudy.low_energy_break
======================================
This is a command to define a cloudy type broken power
law SED - mainly used for testing the code against cloudy.
This SED has hardwired frequency exponents of 2.5 below the
low energy break and -2.0 above the high energy break. This
parameter defines the energy of the low energy break.

Type
  Double

Unit
  eV

Values
  Greater than 0

File
  `setup_star_bh.c <https://github.com/agnwinds/python/blob/master/source/setup_star_bh.c>`_


Parent(s)
  * :ref:`Central_object.rad_type_to_make_wind`: ``cloudy``


