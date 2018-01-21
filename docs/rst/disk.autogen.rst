
====
Disk
====

Disk.rad_type_for_disk_to_make_wind
===================================
The disk is generally described in terms of a run of temperature and possibly gravity with radius.  The spectrum
of the disk can be simulated either as a collection of apppriately weighted blackbodies or from stellar
models which are read in and sampled.

**Type:** Enum (Int)

**Values:**

0. bb

1. models


**Parent(s):**
  parameter_: Condition e.g. >0 or list e.g. [1, 2, 5]


**File:** python.c


Disk.radiation
==============
Multi-line description, must keep indentation.

**Type:** Boolean (1/0)

**Parent(s):**
  parameter_: Condition e.g. >0 or list e.g. [1, 2, 5]


**File:** setup_disk.c


Disk.T_profile_file
===================
Multi-line description, must keep indentation.

**Type:** rdstr

**Parent(s):**
  parameter_: Condition e.g. >0 or list e.g. [1, 2, 5]


**File:** setup_disk.c


Disk.temperature.profile
========================
Multi-line description, must keep indentation.

**Type:** Enum (Int)

**Values:**

0. standard

1. readin

2. analytic


**Parent(s):**
  parameter_: Condition e.g. >0 or list e.g. [1, 2, 5]


**File:** setup.c


