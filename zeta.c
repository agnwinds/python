/**************************************************************************
                    Southampton University
                                                                                                   
                                                                                                   
  Synopsis:

The routines in this file are to do with computing the correction to the saha equation.



  Description:

Zeta is the term in lucy and mazzali which corrects for the proportion of recombinations
going directly to the ground state. I'm experimenting with moving it out of saha and stuart_sim
since the same factor will be used in both, and I'm going to add in some more options to
incorporate a correction factor for dielectronic recombination.

                                                                                                   
  Arguments:  
                                                                                                   
                                                                                                   
  Returns:
                                                                                                   
  Notes:
                                                                                                   
                                                                                                   
                                                                                                   
  History:
	11aug	nsh	Began work
	111211	ksl	Began to comment on it

                                                                                                   
 ************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "atomic.h"
#include "python.h"
#include "recipes.h"




/**************************************************************************
                    Southampton University
                                                                                                   
                                                                                                   
  Synopsis: compute_zeta calculates the correction term in the L+M and sim equations
        for corrected saha.
                                                                                                   
  Description:
	mode 1 is the original version, it takes ilow, ihi and interpfrac and 
        interpolates in the ground state fraction table.
        mode 2 attempts to modify zeta by application of the dielectronic recombination
 	rate
                                                                                                   
  Arguments:  
	temperature - the temperature of the cell
	nion - the ion we are recombining from
  	ilow,ihi - the bracketing elements in the interpfrac table (can be the same in at extremes)
        interpfrac - the distnce between the two bracketing temps
	f1,f2 - the frequency band over which we are interested - NB this doesnt really do that much, 
	since integ_fb ends up resetting the limits.
	mode - this will allow control here, so we can have the 
		1: interpfrac type, 
		2: interpfrac but then incorporating DR as an additional correction.
                                                                                                   
                                                                                                   
  Returns:
	zeta, the correction factor to the recombination rates. 
                                                                                                   
  Notes: When compute_zeta is called, nion is the upper ion in the pair whose abundances are being calculated.
This means, we need to use nion for the DR coefficient, because nion in this case refers to the rate of recombination of the ionised ion into the less ionised state. This equals zero for the totally ionised ion, since there is no electron to get dielectronic with.
We need to use nion-1 for integ_fb, since this is computed via the photoionisation rate from the lower ion into the upper ion via the milne relation. 
 
                                                                                                   
                                                                                                   
                                                                                                   
  History:

                                                                                                   
 ************************************************************************/

double
compute_zeta (temp, nion, ilow, ihi, interpfrac, f1, f2, mode)
     double temp, interpfrac, f1, f2;
     int ilow, ihi, mode, nion;
{
  double zeta, alpha_all, alpha_dr;

  if (mode == 1)
    {
      //This is the old method of getting zeta - just from the ground state tables computed by Christian
      zeta =
	ground_frac[nion - 1].frac[ilow] +
	interpfrac * (ground_frac[nion - 1].frac[ihi] -
		      ground_frac[nion - 1].frac[ilow]);
    }
  else if (mode == 2)
    {
      /* To kick off with, we'll use the old zeta as a basis. 
       * We can tehn multiply it by (recomb to all) / (recomb to all + DR). 
       * Next job will be to compute recomb to ground for all TB ions, and then we can calculate zeta properly for all temperatures
       * , at least for ions where we have TB x sections - i.e. for each level.
       */
      zeta =
	ground_frac[nion - 1].frac[ilow] +
	interpfrac * (ground_frac[nion - 1].frac[ihi] -
		      ground_frac[nion - 1].frac[ilow]);

      /* This initiallizes the FB arrays, from temperatures of 1000K to 1e6, from frequencies between 0 and 1e50.  
       * When these frequency limits are given, they are revised to be more
       * reasonable int he free_bound routines.
       */

      init_freebound (1.e3, 1.e6, 0, 1e50);	//This initialises the FB arrays, we will need them to calculate the rates for radiative recombination.

      alpha_all = integ_fb (temp, 0, 1e50, nion - 1, 2);	//This is the rate for recombinations to all states.
      alpha_dr = dr_coeffs[nion];	//Get the dielectronic recombination coefficient for this ion

      if (alpha_all < 1e-99)
	{
	  Error
	    ("compute_zeta: alpha_all very small value %e (%e) for t of %.1f in element %i in the model?\n", alpha_all, alpha_dr,temp, ion[nion].z);
	  return (zeta);
	}
      zeta *= (alpha_all / (alpha_all + alpha_dr));	//Compute the new value of zeta, corrected for DR.

    }
  return (zeta);
}
