
==============
Boundary_layer
==============

----------------------------------------

Boundary_layer.radiation
========================
Says whether the boundary layer will radiate 0=no, 1=yes

**Type:** Boolean (yes/no)

**Parent(s):**
  :ref:`parameter`: None


**File:** setup_star_bh.c


----------------------------------------

Boundary_layer.rad_type_to_make_wind
====================================
Multi-line description, must keep indentation.

**Type:** Enumerator

**Values:**

  ``bb``
    Multi-line description, must keep indentation.

  ``models``
    Multi-line description, must keep indentation.

  ``power``
    Multi-line description, must keep indentation.


**Parent(s):**
  :ref:`parameter`: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** setup_star_bh.c


----------------------------------------

Boundary_layer.luminosity
=========================
The luminosity of the boundary layer          

**Type:** Double

**Unit:** ergs/s

**Values:** Condition e.g. greater than 0 or list e.g. [1, 2, 5]

**Parent(s):**
  :ref:`parameter`: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** setup.c


----------------------------------------

Boundary_layer.temp
===================
The temperature of the boundary layer in situations where temperature
is meaningful in generating the spectrum

**Type:** Double

**Unit:** Degrees Kelvin

**Values:** Greater than 0

**Parent(s):**
  :ref:`parameter`: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** setup.c


----------------------------------------

Boundary_layer.rad_type_in_final_spectrum
=========================================
Multi-line description, must keep indentation.

**Type:** Enumerator

**Values:**

  ``bb``
    Multi-line description, must keep indentation.

  ``models``
    Multi-line description, must keep indentation.

  ``uniform``
    Multi-line description, must keep indentation.


**Parent(s):**
  :ref:`parameter`: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** python.c


