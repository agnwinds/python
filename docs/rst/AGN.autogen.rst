
===
AGN
===

AGN.bremsstrahlung_temp
=======================
The temperature T in bremstrahlung SED of the form
L_nu=nu**alpha exp(-hnu/kT)

**Type:** Double

**Unit:** K

**Value:** Greater than 0

**Parent(s):**
  parameter_: None


**File:** setup_star_bh.c


AGN.blackbody_temp
==================
The temperature of a blackbody SED to be used for the central AGN source

**Type:** Double

**Unit:** K

**Value:** Greater than 0

**Parent(s):**
  parameter_: None


**File:** setup_star_bh.c


AGN.bremsstrahlung_alpha
========================
The frequency exponent alpha in bremstrahlung SED of the form
L_nu=nu**alpha exp(-hnu/kT)

**Type:** Double

**Unit:** None

**Value:** Any - sign is not asssumed so use negative if you want negative

**Parent(s):**
  parameter_: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** setup_star_bh.c


AGN.power_law_cutoff
====================
An advanced option - this is a low frequency cutoff for an 
AGN power law spectrum. It prevents the powerlaw being 
applied to low frequencies and giving an odd SED.

**Type:** Double

**Unit:** Hz

**Value:** Greater than 0

**Parent(s):**
  parameter_: None


**File:** setup_star_bh.c


AGN.power_law_index
===================
The exponent alpha in a power las SED applied to an AGN
central source of the form L_nu=K nu**alpha

**Type:** Double

**Unit:** None

**Value:** Any - but sign is not assumed, so for negative index use a negative value

**Parent(s):**
  parameter_: None


**File:** setup_star_bh.c


AGN.geometry_for_pl_source
==========================
Choose the geometry for the power law source.

**Type:** Enum (Int)

**Values:**

0. *sphere*
   
   The power law source is a sphere with radius equal to 
   the radius of the central object. In a BH system this is 
   often set to the ISCO.

1. *lamp_post* 
   
   The power law source is a single point at some height above the origin. 
   Photons radiate isotropically from this point. The height is specified in 
   a subsequent parameter, lamp_post.height.


**Parent(s):**
  parameter_: QSO_BH_radiation


**File:** setup_star_bh.c


