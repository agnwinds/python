Macro Atoms
-------------------------

.. todo:: This page shoud contain a description of how macro atoms work. The below is copied from JM's thesis.

The macro-atom scheme was created by Leon Lucy and is outlined in 
his 2002/03 papers. It was implemented in Python by Stuart Sim, initially
for the study of recombination lines in YSOs (Sim et al. 2005).

Lucy (2002,2004) hereafter L02, L03) has shown that it is possible to calculate the emissivity of a gas in statistical equilibrium without approximation for problems with large departures from LTE. His `macro-atom` scheme allows for all possible transition paths from a given level, dispensing with the two-level approximation, and provides a full non-LTE solution for the level populations based on Monte Carlo estimators. The macro-atom technique has already been used to model Wolf-Rayet star winds (Sim 2004), AGN disc winds (Sim et al. 2008), supernovae (Kromer and Sim 2009, Kerzendorf and Sim 2014) and YSOs (Sim et al. 2005). A full description of the approach can be found in L02 and L03. 

Following L02, let us consider an atomic species interacting with a radiation field.
If the quantity :math:`\epsilon_j` represents the ionization plus excitation energy of 
a level :math:`j` then the rates at which the level absorbs and emits radiant energy 
are given by

.. math::

    \dot{A}_{j}^{R} = R_{l j} \epsilon_{j l} 

and

.. math::

    \dot{E}_{j}^{R} = R_{j l} \epsilon_{j l} 


where :math:`\epsilon_{j l} = \epsilon_j - \epsilon_l`.
Here, we have adopted Lucy's convention, in which the subscript 
:math:`l` denotes a summation over all lower states (:math:`<j`), and
(:math:`u`) a summation over all upper states (:math:`>j`).
Similarly, the rates corresponding to _kinetic_ (collisional)
energy transport can then be written as

.. math::

    \dot{A}_{j}^{C} = C_{l j} \epsilon_{j l}

and


.. math::

    \dot{E}_{j}^{C} = C_{j l} \epsilon_{j},

Let us define :math:`{\cal R}` as a total rate, such that
:math:`{\cal R}_{l j}  = R_{l j} + C_{l j}`.
If we now impose statistical equilibrium

.. math::

    ({\cal R}_{l j}-{\cal R}_{j l})+({\cal R}_{u j}-{\cal R}_{ju})=0 \;\;\;,

we obtain 

.. math::

    \dot{E}_{j}^{R}+\dot{E}_{j}^{C}+{\cal R}_{ju}\epsilon_{j}+ {\cal R}_{j \ell}\epsilon_{l}  \nonumber \\  = \dot{A}_{j}^{R}+\dot{A}_{j}^{C}+{\cal R}_{u j} \epsilon_{j} +{\cal R}_{l j} \epsilon_{l}.  

This equation is the starting point for the macro-atom scheme. It shows that, when assuming radiative equilibrium, the energy flows through a system depend only on the transition probabilities and atomic physics associated with the levels the energy flow interacts with. By quantising this energy flow into radiant (:math:`r-`) and kinetic (:math:`k-`) packets, we can simulate the energy transport through a plasma discretised into volume elements (macro-atoms),whose associated transition probabilities govern the interaction of radiant and kinetic energy with the ionization and excitation energy associated with the ions of the plasma.

Although the equation above assumes strict radiative equilbrium,it is trivial to adjust it to include non-radiative source and sink terms. For example, in an expanding parcel of plasma, adiabatic cooling may be included with a simple modification to the RHS. Currently, we include adiabatic cooling by destroying packets with a probability 
:math:`p_{i,\mathrm{destruct}} = {\cal C}_{\rm adiabatic} / {\cal C}_{\rm tot}`.


A Hybrid Scheme
=============================

.. todo:: Describe.


Bound-free Continua of Simple Atoms
=============================================

Historically, when using the indivisible packet form of radiative transfer (`macro_atoms_thermal_trapping`, for example), the bound-free continua of simple atoms were treated in a simplified two-level framework. In this case, simple atoms are those `without` a full macro-atom model atom, usually the metals. In this two-level scheme, whenever a simple atom undergoes a bound-free interaction, it is excited into the continuum state, and this is immediately followed by recombination, and an :math:`r`-packet or :math:`k`-packet is created immediately. As a result, the scheme does not capture the physical situation whereby a recombination cascade can occur from an initial recombination to excited levels, leading to a gradual reddening of the photon if there are many interactions. This situation **is** modelled well by a full macro-atom treatment. 

To try and slightly improve this scheme, we implemented a "total emissivity" upweighting scheme around 2018. The basic idea is that we pay attention to only the heating and cooling. In particular, the rates of all simple atom bound-free emission are governed by the `emissivity` of the bound-free process. 

This result in two changes to the code for ionization cycles: 
   * whenever a k-packet is eliminated via a bound-free channel of a simple macro atom (simulating energy flow from the :math:`k`-packet pool to the radiation pool, :math:`k \to r`), we have that packet carry additional energy corresponding to the required ionization energy for that particular bf process. This means we upweight the energy of the packet by a factor :math:`f_{\rm up} = \nu / (\nu - \nu_0)`, where :math:`\nu` is the frequency of the new bound-free photon and :math:`\nu_0` is the threshold frequency. This quantity is the ratio of the total energy carried by photons in the packet to the energy supplied to photons in the packet from the thermal pool. 
   * whenever an r-packet is “absorbed” by a simple macro atom bound-free process we track explicitly only the flow of energy to the thermal pool. This means we force the creation of a :math:`k`-packet, whereas before there woud be a choice, but we only take the contribution of the absorption to heating only: i.e. we downweight the packet energy by a factor :math:`f_{\rm down} = (\nu - \nu_0) / \nu`.

In the spectral cycles, interactions with simple bound-free continua now kill the photon, and :math:`k \to r` follow the same behaviour as above, because in these cycles we introduce a precalculated band-limited :math:`k`-packet emissivity. 

**It is possible for some numerical problems to occur.** For example, there is nothing to stop the value of :math:`f_{\rm up}` being quite large, if the photon is being emitted close to the edge. This is most likely to happen when the electron temperature $T_e$ is quite low, but there is nothing to stop it happening anywhere. This is most likely to lead to problems when the factor :math:`f_{\rm up}` is comparable to the typical number of photon passages per cell, since then a single photon can dominate the heating or ionization estimators in a given cell and lead to convergence problems by dramatically exacerbating shot noise. 

.. todo:: Finish documentation of upweighting scheme with some basic explanation.
