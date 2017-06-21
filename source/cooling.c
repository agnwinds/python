


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

  Description:	

  Arguments:  (Input via .pf file)		


  Returns:

  Notes:


 

  History:
2004	ksl	Coded as better.c

 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


double
cooling (xxxplasma, t)
     PlasmaPtr xxxplasma;
     double t;
{

  /*Original method */
  xxxplasma->t_e = t;



  if (geo.adiabatic)
    {
      if (wmain[xxxplasma->nwind].div_v >= 0.0)
	{
	  /* This is the case where we have adiabatic cooling - we want to retain the old behaviour, 
	     so we use the 'test' temperature to compute it. If div_v is less than zero, we don't do
	     anything here, and so the existing value of adiabatic cooling is used - this was computed 
	     in wind_updates2d before the call to ion_abundances. */
	  xxxplasma->cool_adiabatic =
	    adiabatic_cooling (&wmain[xxxplasma->nwind], t);
	}
    }

  else
    {
      xxxplasma->cool_adiabatic = 0.0;
    }


  /*81c - nsh - we now treat DR cooling as a recombinational process - still unsure as to how to treat emission, so at the moment
     it remains here */

  xxxplasma->cool_dr = total_fb (&wmain[xxxplasma->nwind], t, 0, VERY_BIG, FB_REDUCED, 2);

  /* 78b - nsh adding this line in next to calculate direct ionization cooling without generating photons */

  xxxplasma->cool_di = total_di (&wmain[xxxplasma->nwind], t);

  /* 70g compton cooling calculated here to avoid generating photons */

  xxxplasma->cool_comp = total_comp (&wmain[xxxplasma->nwind], t);

  xxxplasma->cool_tot =
    xxxplasma->cool_adiabatic + xxxplasma->cool_dr + xxxplasma->cool_di +
    xxxplasma->cool_comp + xtotal_emission (&wmain[xxxplasma->nwind], 0.,
					  VERY_BIG);


  return (xxxplasma->cool_tot);
}




/***********************************************************
             Space Telescope Science Institute

Synopsis:  total_emission (one, f1, f2) Calculate the total emission of a single cell 

Arguments:		


Returns:
 
Description:	
	

Notes:
	Total emission gives the total enery loss due to photons.  It does
	not include other coooling sources, e. g. adiabatic expansion.

	It returns the total luminosity, but also stores the luminosity due
	to various types of emssion, e.g ff, fb, lines, compton into the
	Plasms cells

	Comment: The call to this routine was changed when PlasmaPtrs
	were introduced, but it appears that the various routines 
	that were called were not changed.

History:
	97	ksl	Coded
	98oct	ksl     Removed temperature limits within total_emission in
			and attempt to concentrate the frequency limits elsewhere.
	01nov	ksl	Removed superflous code
	01nov	ksl	Changed call so t_e passed through wind structure.  Other
			routines should be changed ALSO, BUT NOT NOW
	04Mar   SS      Modified to divert if Macro Atoms being used.
	04Apr   SS      Modified to include a collisional cooling term
	                - this is getting very messy and confusing - I think
                        I should probably put together a completely separte
			set of routines for doing the heating and cooling in the
			macro atom cases. (SS) since this is now not really a
			total emission but total cooling.
	05may	ksl	57+ -- To use plasma structure for most things.  If vol
			is added to plasma then one can change the call
	11aug	nsh	70 changes made to allow compton heating to be calculated.
 
 
**************************************************************/


double
xtotal_emission (one, f1, f2)
     WindPtr one;               /* WindPtr to a specific cell in the wind */
     double f1, f2;             /* The minimum and maximum frequency over which the emission is
                                   integrated */
{
  double t_e;
  int nplasma;
  PlasmaPtr xplasma;

  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];

  t_e = xplasma->t_e;           // Change so calls to total emission are simpler

  if (f2 < f1)
  {
    xplasma->lum_tot = xplasma->lum_lines = xplasma->lum_ff = xplasma->cool_rr = 0;      //NSH 1108 Zero the new cool_comp variable NSH 1101 - removed
  }
  else
  {
    if (geo.rt_mode == RT_MODE_MACRO)       //Switch for macro atoms (SS)
    {
      xplasma->cool_rr = total_fb_matoms (xplasma, t_e, f1, f2) + total_fb (one, t_e, f1, f2, FB_REDUCED, 1);        //outer shellrecombinations
      //The first term here is the fb cooling due to macro ions and the second gives
      //the fb cooling due to simple ions.
      //total_fb has been modified to exclude recombinations treated using macro atoms.
      xplasma->lum_tot = xplasma->cool_rr;
      //Note: This the fb_matom call makes no use of f1 or f2. They are passed for
      //now in case they should be used in the future. But they could
      //also be removed.
      // (SS)
      xplasma->lum_lines = total_bb_cooling (xplasma, t_e);
      xplasma->lum_tot += xplasma->lum_lines;
      /* total_bb_cooling gives the total cooling rate due to bb transisions whether they
         are macro atoms or simple ions. */
      xplasma->lum_ff = total_free (one, t_e, f1, f2);
      xplasma->lum_tot += xplasma->lum_ff;


    }
    else                        //default (non-macro atoms) (SS)
    {
      xplasma->lum_tot = xplasma->lum_lines = total_line_emission (one, f1, f2);
      xplasma->lum_tot += xplasma->lum_ff = total_free (one, t_e, f1, f2);
      xplasma->lum_tot += xplasma->cool_rr = total_fb (one, t_e, f1, f2, FB_REDUCED, 1);     //outer shell recombinations


    }
  }


  return (xplasma->lum_tot);


}



/***********************************************************
                                       Space Telescope Science Institute

Synopsis:  
	adiabatic_cooling (one, t) determines the amount of 
	adiabatic cooling in a cell, in units of luminosity.

Arguments:		
	WindPtr one;	pointer to wind cell
	double t;		electron temperature

Returns:
 
Description:	
   Adiabatic cooling is clearly the amount of PdV work done by
   a fluid element, per unit time dt. Thus it is equal to P dV/dt.  
   The only real question is whether dV/dt is given by the volume * div v.
   div v here is the divergence of the velocity field.
	
Notes:
  JM 1401 -- I've rederived the expression for dv/dT. It follows
        directly from the continuity equation and is indeed equal 
        to volume * div_v. 

        Note also that this function should only be called
        if geo.adiabatic == 1, in which case it populates
        xplasma->cool_adiabatic. This is used in heating and cooling
        balance. We also use it as a potential destruction choice for 
        kpkts in which case the kpkt is thrown away by setting its istat 
        to P_ADIABATIC.


History:
	04nov	ksl	Stuart had switched adiabatic cooling off
  			by commenting out the entire program.  I have
  			reinstated the calculation, but created a new
  			variable geo.adiabatic to determine whether this
  			routine is called in python.  At present, even
  			if the routine is called it has no effect.  The
  			or more properly, one issue is how to meld adiabatic 
  			cooling with the macro atom approach.
	06may	ksl	57+ -- Adapted to include plasma structure
	11aug	ksl	70 - Adiabatic cooling returns the cooling but
			does not store it.
 
 
**************************************************************/


double
adiabatic_cooling (one, t)
     WindPtr one;
     double t;
{
  double cooling;
  int nplasma, nion;
  double nparticles;
  PlasmaPtr xplasma;

  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];

  //JM 1401 -- here was an old factor of 3/2 which KSL and JM believe to be incorrect. 
  //JM 1601 -- this should include the pressure from all particles, rather than just ne
  //cooling = xplasma->ne * BOLTZMANN * t * xplasma->vol * one->div_v;
  nparticles = xplasma->ne;

  for (nion = 0; nion < nions; nion++)
  {
    /* loop over all ions as they all contribute to the pressure */
    nparticles += xplasma->density[nion];
  }

  cooling = nparticles * BOLTZMANN * t * xplasma->vol * one->div_v;

  return (cooling);
}

