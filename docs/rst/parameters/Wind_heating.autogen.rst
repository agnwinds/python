
============
Wind_heating
============

Wind_heating.extra_processes
============================
Multi-line description, must keep indentation.

**Type:** Enumerator

**Values:**

  ``adiabatic``
    Multi-line description, must keep indentation.

  ``both``
    Multi-line description, must keep indentation.

  ``none``
    Multi-line description, must keep indentation.

  ``nonthermal``
    Multi-line description, must keep indentation.


**File:** setup.c


----------------------------------------

Wind_heating.kpacket_frac
-------------------------
Multi-line description, must keep indentation.

**Type:** Double

**Unit:** None

**Value:** Condition e.g. greater than 0 or list e.g. [1, 2, 5]

**Parent(s):**
  Wind_heating.extra_processes_: nonthermal, both

  Line_transfer_: macro_atoms, macro_atoms_thermal_trapping


**File:** setup.c


----------------------------------------

Wind_heating.extra_luminosity
-----------------------------
This is a very special option put in place for modelling FU Ori stars, and should be used with extreme caution. Determines the shock factor.

**Type:** Double

**Value:** Condition e.g. greater than 0 or list e.g. [1, 2, 5]

**Parent(s):**
  Wind_heating.extra_processes_: nonthermal, both


**File:** setup.c


