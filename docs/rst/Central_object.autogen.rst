
==============
Central_object
==============

Central_object.rad_type_for_star_to_make_wind
=============================================
The way in which radiation from the central object will be simulated, either as a bb or from models which
have been read in separately.

**Type:** Enum (Int)

**Values:**

0. bb

1. models


**Parent(s):**
  parameter_: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** setup.c


Central_object.radius
=====================
Radius of the central object in the system, e.g the white dwarf or black hole

**Type:** Double

**Unit:** cm

**Value:** Greater than 0

**Parent(s):**
  parameter_: None


**File:** setup.c


Central_object.rad_type_in_final_spectrum
=========================================
The type of spectral models used to simulate radiation from the central object, where bb inplies
bb radiation, models implies spectra from models that are read in, and uniform means a flat spectrum
regardless of parameters like temperature.

**Type:** Enum (Int)

**Values:**

0. bb

1. models

2. uniform


**Parent(s):**
  parameter_: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** python.c


Central_object.radiation
========================
A booliean variable stating whether of not the central object should radiate from its
survace as a star would. 

**Type:** Boolean (1/0)

**Parent(s):**
  parameter_: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** setup.c


Central_object.mass
===================
Mass of the central object

**Type:** Double

**Unit:** Solar masses

**Value:** Greater than 0

**Parent(s):**
  parameter_: None


**File:** setup.c


Central_object.temp
===================
Multi-line description, must keep indentation.

**Type:** rddoub

**Unit:** None

**Value:** Condition e.g. greater than 0 or list e.g. [1, 2, 5]

**Parent(s):**
  parameter_: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** setup_star_bh.c


