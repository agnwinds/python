
=====
shell
=====

shell.wind_v_at_rmin
====================
The velocity of a shell wind at the inner edge of the 
shell - the variation of the velocity in the shell is
set by the velocity law exponent. It allows a gradient 
to be enforced.

**Type:** Double

**Unit:** cm/s

**Value:** Greater than or equal to zero.

**Parent(s):**
  Wind_type_: 9


**File:** shell_wind.c


shell.wind.radmin
=================
The innermost edge of a diagnostic type of wind made up of a single
(ideally thin) shell.

**Type:** Double

**Unit:** cm

**Value:** greater than zero, less than wind.radmax

**Parent(s):**
  Wind_type_: 9


**File:** shell_wind.c


shell.wind.acceleration_exponent
================================
Exponent beta for the Caster and Lamers description of a stellar wind
v(r)=v_o + (v_inf - v_o) (1+R_s/r)**beta for a shell wind.

**Type:** Double

**Value:** Greater than or equal to 0

**Parent(s):**
  Wind_type_: 9


**File:** shell_wind.c


shell.wind.v_at_rmax
====================
The velocity of a shell wind at the outer edge of the 
shell - the variation of the velocity in the shell is
set by the velocity law exponent. It allows a gradient 
to be enforced.

**Type:** Double

**Unit:** cm/s

**Value:** Greater than or equal to zero.

**Parent(s):**
  Wind_type_: 9


**File:** shell_wind.c


