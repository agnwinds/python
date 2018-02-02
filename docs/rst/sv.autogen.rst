
==
sv
==

sv.acceleration_exponent
========================
Multi-line description, must keep indentation.

**Type:** Double

**Value:** Condition e.g. greater than 0 or list e.g. [1, 2, 5]

**Parent(s):**
  parameter_: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** sv.c


sv.acceleration_length
======================
Multi-line description, must keep indentation.

**Type:** Double

**Unit:** cm

**Value:** Condition e.g. greater than 0 or list e.g. [1, 2, 5]

**Parent(s):**
  parameter_: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** sv.c


sv.diskmax
==========
Multi-line description, must keep indentation.

**Type:** Double

**Unit:** co.radius

**Value:** Condition e.g. greater than 0 or list e.g. [1, 2, 5]

**Parent(s):**
  parameter_: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** sv.c


sv.diskmin
==========
Multi-line description, must keep indentation.

**Type:** Double

**Unit:** co.radius

**Value:** Condition e.g. greater than 0 or list e.g. [1, 2, 5]

**Parent(s):**
  parameter_: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** sv.c


sv.mdot_r_exponent
==================
Multi-line description, must keep indentation.

**Type:** Double

**Value:** Condition e.g. greater than 0 or list e.g. [1, 2, 5]

**Parent(s):**
  parameter_: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** sv.c


sv.thetamax
===========
The angle at which the wind rises from the outermost launching radius in a Shlossman-Vitello type disk wind.
This angle is measured with respect to the vertical (z) direction i.e. zero describes a vertical wind.
See figure 1 of Shlossman & Vitello 1993, ApJ 409,372.

**Type:** Double

**Unit:** Degrees

**Value:** Greater than sv.thetamin

**Parent(s):**
  wind_type_: 0


**File:** sv.c


sv.thetamin
===========
The angle at which the wind rises from the innermost launching radius in a Shlossman-Vitello type disk wind.
This angle is measured with respect to the vertical (z) direction. I.e. zero descirbes a vertical wind.
See figure 1 of Shlossman & Vitello 1993, ApJ, 409, 372.

**Type:** Double

**Unit:** Degrees

**Value:** Greater than= 0

**Parent(s):**
  parameter_: required when the wind_type is set to 0, a SV wind.


**File:** sv.c


sv.v_infinity
=============
Multi-line description, must keep indentation.

**Type:** Double

**Unit:** Escape velocity

**Value:** Condition e.g. greater than 0 or list e.g. [1, 2, 5]

**Parent(s):**
  parameter_: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** sv.c


