
===============
Wind_ionization
===============

Wind_ionization
===============
The method used to calculate the ionization state of the wind. Does
not apply to macro-atoms. Some descriptions of the assumptions and derivations
can be found in the papers referenced below, as well as via this link: 
http://jhmatthews.github.io/pydocs/docs/ionization-schemes-rate.pdf. 

**Type:** Enum (Int)

**Values:**

0. *on.the.spot*
   
   Mazzali & Lucy (1993) modified Saha method but using existing t_e. 
   No heating and cooling balance.

1. *LTE(tr)*
   
   LTE abundances (the Saha equation) using the radiation temperature.

2. *fixed*
   
   Ion fractions are fixed and read from a file.

3. *recalc_bb*
   
   The modified Saha approach, based on dilute blacbody ionization 
   originally described by Abbott & Lucy (1985), modified by Mazzali & Lucy (1993) 
   and also used in the original Long & Knigge (2002) study.

4. *LTE(t_e)*
   
   LTE abundances (the Saha equation) using the electron temperature.

6. *pairwise_bb*
   
   A dilute blackbody intensity is used to calculate the photoionization
   rate, and variable temperature modified Saha equation described by Higginbottom
   et al. (2013) is used.

7. *pairwise_pow*
   
   The modelled mean intensity is used to calculate the photoionization
   rate, and variable temperature modified Saha equation described by Higginbottom
   et al. (2013) is used.

8. *matrix_bb*
   
   A dilute blackbody mean intensity is used to calculate the photoionization
   rate, and the full rate matrix for the ions is solved to give ionization 
   fraction.

9. *matrix_pow*
   
   The modelled mean intensity is used to calculate the photoionization
   rate, and the full rate matrix for the ions is solved to give ionization 
   fraction.


**Parent(s):**
  parameter_: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** setup2.c


