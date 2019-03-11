
====
wind
====

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


wind.fixed_concentrations_file
==============================
The filename for the fixed ion concentrations if you have
set Wind_ionization to 2 (fixed). This file has format
[atomic_number  ionizationstage   ion fraction]. 

**Type:** String

**Parent(s):**
  parameter_: Whenever the wind ionization choice is fixed


**File:** setup2.c


