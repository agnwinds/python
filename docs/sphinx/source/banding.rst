Photon Banding
##############

Photon packets are emitted from a number of different radiation sources, such as
the accretion disk or from the wind itself. When a photon is created, it
is defined by its frequency :math:`\nu` and weight :math:`w`. Photons are 
generated at the beginning of each cycle and can either be generated uniformly
over the entire frequency range, or can be generated in pre-defined frequency
bands, where certain frequency bands are biased to have more photons.

Uniform Samlping
================

In the most simple case photon packets are generated uniformly over the entire
frequency range, with equal photon having equal weight. In this uniform scheme, 
the total weight of all photon packets is equal to the luminosity of the system,
where each photon packet has weight,

.. math ::
    w = \frac{\sum_{\text{sources}} \int_{\nu_{\text{min}}}^{\nu_{\text{max}}} L_{\nu}}{N},

where :math:`N` is the total number of photons, :math:`\nu_{\text{min}}` and 
:math:`\nu_{\text{max}}` define the frequency range and :math:`L_{\nu}` is the
luminosity for a radiation source. Note that a summation is used to find the
total luminosity from all radiation sources and that :math:`w` has units
:math:`\text{ergs s}^{-1}`.

Banded Sampling
===============

In practise, uniform sampling has a few limiations.
