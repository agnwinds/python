
====
Disk
====

Disk.radiation
==============
Multi-line description, must keep indentation.

**Type:** Boolean (1/0)

**Parent(s):**
  parameter_: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** setup_disk.c


----------------------------------------

Disk.temperature.profile
------------------------
The choice of disk temperature profile

**Type:** Enum (Int)

**Values:**

0. standard - A Shakura - Sunyaev  disk, with a hard inner boundar

1. readin - read the profile in from a file; the user will be queried for the name of the file

2. analytic - A profile designed for the situation where the disk is being illuminated by star


**Parent(s):**
  Disk.radiation_: This input is requested for all disks that radiate


**File:** setup_disk.c


----------------------------------------

Disk.T_profile_file
^^^^^^^^^^^^^^^^^^^
When the user chooses to read in the temperature profile as a
function of radius, the user is asked the name of the file that
contains the desired profile.

**Type:** String

**Parent(s):**
  Disk.temperature.profile_: Unspecified


**File:** setup_disk.c


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
  parameter_: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** python.c


