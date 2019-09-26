Photon Banding Strategies
#########################

Photon packets are emitted from a number of different radiation sources, such as
the accretion disk or from the wind itself. When a photon is created, it
is defined by its frequency :math:`\nu` and weight :math:`w`. Photons are 
generated at the beginning of each cycle and can either be generated uniformly
over the entire frequency range, or can be generated in pre-defined frequency
bands, where certain frequency bands are biased to have more photons.

Uniform Sampling
================

In the most simple case, the frequency of a photon is sampled uniformly over the
entire frequency range. The total weight of all photons is equal to the 
luminosity of the system and each photon has weight has equal weight given by,

.. math ::
    w_{i} = \frac{\sum_{i = \text{sources}} ~ \int_{\nu_{\text{min}}}^{\nu_{\text{max}}} L_{\nu, i} d\nu}{N},

where :math:`N` is the total number of photons, :math:`\nu_{\text{min}}` and 
:math:`\nu_{\text{max}}` define the frequency range and :math:`L_{\nu, i}` is the
luminosity for a radiation source. Note that a summation is used to find the
luminosity for each radiation source and that :math:`w` has units of 
:math:`\text{ergs s}^{-1}`.

Banded Sampling
===============

In practice, uniform sampling is generally an inefficient approach to generating
photon packets. For example, it is often desirable to produce a sufficient 
number of photons within a specific frequency range, i.e. around a 
photoionisation edge. However, if these frequencies lie on the Wien tail of a
blackbody distribution, it is unlikely that a sufficient number of photons will
be generated as most of the luminosity is generated at lower frequencies.
It is possible to get around this problem by generating an increasingly large
number of photons. But, this is computationally expensive and inefficient.

In order to cope with cases this, Python implements *importance sampling* which
effectively increases the number of photons which are sampled from specific 
frequency bands considered important. Photons are now generated with the weight,

.. math ::
    w_{j} = \frac{\sum_{j = \text{sources}} ~ \int_{\nu_{i}}^{\nu_{i + 1}} L_{\nu, j} ~ d\nu}{f_{i} N},

where, again, this is a summation over all radiation sources. :math:`N` is the
total number of photons, :math:`f_{i}` is the fraction of photons emerging from
frequency band :math:`i`, :math:`\nu_{i}` and :math:`nu_{i+1}` are the lower
and upper frequency boundaries for each frequency band and :math:`L_{\nu, j}` is
the luminosity of the radiation source. Hence, more photons from frequency bands
with a larger fraction :math:`f_{i}` will be generated. However, photons from
*important* bands (where more photons are sampled from) will have reduced 
weight, whilst photons from the *less important* frequency bands will have
increased weight. 

This scheme has the benefit of allowing the generation of a lower number of
photons whilst still sufficiently sampling important frequency ranges, 
decreasing the computational expense of a simulation.

Implemented Banding Schemes
===========================
