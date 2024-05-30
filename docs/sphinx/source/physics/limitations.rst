Limitations and Caveats
-------------------------

.. todo:: Write these descriptions and complete list.

There are a number of limitations of Python as a code. This page is a non-exhaustive list.

**Extreme optical depths** 

**Pecular emission features in high optical depth models:** As described in `this issue <https://github.com/agnwinds/python/issues/659>`_ 

**Free-free Gaunt factors:** As described in `this issue <https://github.com/agnwinds/python/issues/33>`_, we use a free-free integrated Gaunt factor for both free-free opacity and cooling calculations, but one should really use a frequency-specific (photon by photon) Gaunt factor for the opacity calculation. 

**Compton heating energy leak:** As described in `this issue <https://github.com/agnwinds/python/issues/295>`_, 
there is the potential for a slight energy leak (or poor statistics) for Compton heating in non-macro atom mode.

**Iron Fluorescence:** As described in `this issue <https://github.com/agnwinds/python/issues/499>`_,
Fluorescence data is not dealt with properly in non-macro-atom mode. The user should use a Fe macro-atom data set for
situations where this is important. 

**The Sobolev and related approximations**

**The Size of the Line List**

**General Relativity**