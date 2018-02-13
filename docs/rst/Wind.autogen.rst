
====
wind
====

wind.radmax
===========
The maximum radial distance to follow photons in a wind.  Beyond
this point photons do not interact with the wind at all. Note that
this is a global parameter and refers to al of the wind domains

**Type:** Double

**Unit:** cm

**Value:** Greater than 0

**Parent(s):**
  parameter_: Whenver there is a wind one will be asked to provide this parameter.


**File:** setup_domains.c


wind.mdot
=========
The mass lost rate of the wind in several of the kinematic
wind models, e.g SV or KWD.

**Type:** Double

**Unit:** Msol/year

**Value:** Greater than 0

**Parent(s):**
  wind_type_: Various of the kinematic mdoes for the wind


**File:** sv.c, knigge.c


wind.t.init
===========
This parameter sets the initial temperature of the wind.  (It currently
refers to all domains).

**Type:** Double

**Unit:** None

**Value:** Greater than 0

**Parent(s):**
  parameter_: This question is asked whenever there is a wind, that is to say essentially always


**File:** setup_domains.c


