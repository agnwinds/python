Radiation sources
#################

.. todo::
  Fill in. Add description of how to use your own input spectrum. Finish links to keywords. 

External radiation Sources
==========================

In generic terms, there are two main external radiation sources for any Python calculation: a **Central Source** which can be a normal star, a WD, or a BH, and a **disk**.  Even though Python supports the existence of a secondary star for the purposes of calculating when light from a disk system is occulted, the secondary star does not radiate.


Photons for radiation from the central object emerge uniformly over its surface, except when a lamp-post geometry is specified for the ``bh`` or ``agn`` system types. In this lamp-post case, radiation originates from a point source above and below the central object, with a specified height. 

Emission from a boundary layer can also be defined when this is relevant, from which radiation also emerges uniformly over the surface of the central object.

The Wind as a radiation source
==============================

In macro-atom calculations the wind is NOT a radition source.  All of the photons in a macro-atom
calculation are generated externally, and with minor exceptions photons preserve their weight 
throughout their passage through the wind.  (The minor exceptions have to do with processes like
adiabiatic cooling, which result in the loss of photons).

In the simple-atom approach, various processes cause photons passing through the wind to lose energy
as they pass through the wind. This energy heats the plasma.  To account for this, photons are
generated from the wind at the beginning of each cycle.  Processes include, free-free emission, free-bound emission and line emission.  

In non-macro-atom calculations wind radiation can be turned on and off using the :doc:`/input/parameters/wind/Wind/Wind.radiation` keyword. 



Spectra of the external radiation sources
=========================================

For the most part, the various radiation sources can radiate by any of the following process, as appropriate)

1. Blackbody radiation, specified in terms of a temperature.  Depending on the nature of the source, the luminosity
   specified either by the size of the object, or directly as the total luminosity.

2. Bremsstrahlung radiation, specified in terms of a temperature and a luminosity between 2 and 10 keV

3. Power law radiation, specified in terms of a spectral index, and a luminosity between 2 and 10 keV

4. One or more spectral models read from a series of files.  The models must specified in terms of two 
   parameters, usually T and g, each model consists of an ascii file containing a set of wavelengths 
   and a quantity that is understood to be proportional to :math:`F_{\lambda}`.  An example of the ascii files 
   that can be read in is contained in the xdata folder that is part of the distribution.  


In the ionization cycles, the spectra of the central source, boundary layer (if present) and disk are determined by these three keywords:

* :doc:`/input/parameters/central_object/Central_object/Central_object.rad_type_to_make_wind`
* :doc:`/input/parameters/central_object/Boundary_layer/Boundary_layer.rad_type_to_make_wind`
* :doc:`/input/parameters/Disk/Disk.rad_type_to_make_wind`

It is possible to choose different input spectra for the ionization and spectral cycles, so a corresponding keyword of the form :doc:`/input/parameters/Disk/Disk.rad_type_in_final_spectrum` is also needed. 


