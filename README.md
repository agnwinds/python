README
***
=========
precursor - python_70f
This version incorporates a very approximate treatment of dielectric recombination cooling.
dielectronic.c function total_dr included to compute dielectronic recombination cooling
python.c the bands for the power law ionisation have been changed to exactly match the banding using in photon generation in the agn project. It will probably be best to actually force there two arrays to match if this is what we want
emission.c removed compton luminosity from here. It was being added to the quantity total_lum 
- which is then used to generate wind photons. This seems to be incorrect. This change is a quick fix - more thought is needed long term
ionization.c some changes to the values reported in the convergence - the main change is to add compton cooling and DR cooling to the zero_emit function alongside adiabatic cooling. This means that they contribute to the cooling but are not used to calculate wind luminosity.
wind_updates2d.c changes made to how heating and cooling mechanisms are reported to reflect the way that there are now three cooling mechanisms that are not luminosities...
resonate.c - removed a load of reporting lines put in by NSH for an interesting photon ages ago.


