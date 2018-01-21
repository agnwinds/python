
==============
Boundary_layer
==============

Boundary_layer.luminosity
=========================
The luminosity of the boundary layer          

**Type:** Double

**Unit:** ergs/s

**Value:** Condition e.g. >0 or list e.g. [1, 2, 5]

**Parent(s):**
  parameter_: Condition e.g. >0 or list e.g. [1, 2, 5]


**File:** setup.c


Boundary_layer.rad_type_in_final_spectrum
=========================================
Multi-line description, must keep indentation.

**Type:** Enum (Int)

**Values:**

0. bb

1. models

2. uniform


**Parent(s):**
  parameter_: Condition e.g. >0 or list e.g. [1, 2, 5]


**File:** python.c


Boundary_layer.rad_type_to_make_wind
====================================
When the system contains a boundary layer, the spectrum of the boundary layer can be simulated
as a blackbody, from a model, or as a power law.

**Type:** Enum (Int)

**Values:**

0. bb

1. models

3. pow


**Parent(s):**
  parameter_: Condition e.g. >0 or list e.g. [1, 2, 5]


**File:** setup.c


Boundary_layer.radiation
========================
Multi-line description, must keep indentation.

**Type:** Boolean (1/0)

**Parent(s):**
  parameter_: Condition e.g. >0 or list e.g. [1, 2, 5]


**File:** setup.c


Boundary_layer.temp
===================
The temperature of the boundary layer in situations where temperature
is meaningful in generating the spectrum

**Type:** Double

**Unit:** Degrees Kelvin

**Value:** Greater than 0

**Parent(s):**
  parameter_: Condition e.g. >0 or list e.g. [1, 2, 5]


**File:** setup.c


