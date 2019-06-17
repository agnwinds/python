
==
BH
==

----------------------------------------

BH.lum
======
Multi-line description, must keep indentation.

**Type:** Double

**Unit:** ergs/s

**Values:** Condition e.g. greater than 0 or list e.g. [1, 2, 5]

**Parent(s):**
  :ref:`bh.rad_type_to_make_wind`: Not ``bb``


**File:** setup_star_bh.c


----------------------------------------

BH.power_law_index
==================
Multi-line description, must keep indentation.

**Type:** Double

**Values:** Condition e.g. greater than 0 or list e.g. [1, 2, 5]

**Parent(s):**
  :ref:`parameter`: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** setup_star_bh.c


----------------------------------------

BH.rad_type_to_make_wind
========================
Multi-line description, must keep indentation.

**Type:** Enumerator

**Values:**

  ``bb``
    Multi-line description, must keep indentation.

  ``brems``
    Multi-line description, must keep indentation.

  ``cloudy``
    Multi-line description, must keep indentation.

  ``models``
    Multi-line description, must keep indentation.

  ``power``
    Multi-line description, must keep indentation.


**Parent(s):**
  :ref:`parameter`: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** setup_star_bh.c


BH.blackbody_temp
-----------------
Multi-line description, must keep indentation.

**Type:** Double

**Unit:** Kelvin

**Values:** Greater than 0

**Parent(s):**
  :ref:`BH.rad_type_to_make_wind`: bb


**File:** setup_star_bh.c


----------------------------------------

BH.rad_type_in_final_spectrum
=============================
Multi-line description, must keep indentation.

**Type:** Enumerator

**Values:**

  ``brems``
    Multi-line description, must keep indentation.

  ``cloudy``
    Multi-line description, must keep indentation.

  ``power``
    Multi-line description, must keep indentation.


**Parent(s):**
  :ref:`parameter`: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** python.c


----------------------------------------

BH.power_law_cutoff
===================
Multi-line description, must keep indentation.

**Type:** Double

**Values:** Condition e.g. greater than 0 or list e.g. [1, 2, 5]

**Parent(s):**
  :ref:`parameter`: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** setup_star_bh.c


----------------------------------------

BH.radiation
============
Multi-line description, must keep indentation.

**Type:** Boolean (yes/no)

**Parent(s):**
  :ref:`parameter`: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** setup_star_bh.c


----------------------------------------

BH.geometry_for_pl_source
=========================
Multi-line description, must keep indentation.

**Type:** Enumerator

**Values:**

  ``lamp_post``
    Multi-line description, must keep indentation.

  ``sphere``
    Multi-line description, must keep indentation.


**Parent(s):**
  :ref:`bh.radiation`: True


**File:** setup_star_bh.c


