/**************************************************************************
                    Space Telescope Science Institute
                                                                                                   
                                                                                                   
  Synopsis:

The routines in this file all have to do with photoionization and 
recombination rates and emissivities.  

(Note -- These descriptions need to be cleaned up and checked at some point. 
Today 06jul27 -- I just reorganized things more sensibly, but this section
should describe all of the routines generically)

Calculate the free-bound emissivity of a cell.   This is the emissivity,
eg ne * ni *enu.  It includes the electron density

  Description:

This is the version of the program that used detailed balance to
calculate the emissivity.  The specific formulation implemented
is that described most clearly in Hazy, but that comes from
Brown and Matthews 1970.

ion_choice      0- (nions-1) causes fb to return the emissivity of
		a specific ion number, i.e. 0 will be h1, 2 will be
		he1, etc.   
		>=nions returns the total emissivity
		<0  returns the emissivity of all metals, under the
		assumption that H and He comprise the first five
		ions in the input array.
fb_choice	When an ion recombines there are three possible quantities
		of interest.
			0=the emissivity/per unit frequency of the plasma
			1=the heat lost from the electrons (that is the
				fraction of the emissivity that does not
				go into the binding energy of the level
				to which the electron+ion are recombining)
			2=the emissivity in photons/unit frequency
		0 should be used in calculations of the emission spectrum but
		1 and should be used in energy loss and gain  calculations), while
	        2 should be used in ion densities and levels
                                                                                                   
  Arguments:  
                                                                                                   
                                                                                                   
  Returns:
                                                                                                   
  Notes:
                                                                                                   

For reference here is  the freebound structures that can be used for both the
specific emissivity of a free-bound transition.  But the real structue is
in python.h

	#define NTEMPS        100             // The number of temperatures which are stored in each fbstruct
	#define NFB   10              // The maximum number of frequency intervals for which the fb emission is calculated

	struct fbstruc
	{
	  double f1, f2;
	  double emiss[NIONS][NTEMPS];
	}
	freebound[NFB];

	double xnrecomb[NIONS][NTEMPS];       // There is only one set of recombination coefficients
	double fb_t[NTEMPS];

                                                                                                   
  History:
	01oct	ksl	Began work
	01nov	ksl	Adapt to include Verland cross-sections as well.
	02jun	ksl	Modified so could be used to calculate emissivities in
			photons as well as energy.
	06jul	ksl	57h-Standardized the headers and notes.  I have not
			really checked that everything is up to date, just
			put it into a more readable format.
	12jul	nsh	73-Subroutine bad_t_rr coded to generate a total
			recombination rate from badnell type parameters
	14jan	nsh	77a - Added some checks into the integrals, to ensure
			we do not attempt to integrate over a range of 
			frequencies so large that the integrand is zero over
			an excessive range - hence causing QROMB to return
			an answer of zero.
                                                                                                   
 ************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "atomic.h"
#include "python.h"
#include "recipes.h"



// Next line is left in to be sure that nfb is initialized properly at
// the beginning of the program

/* Next line is required for proper initialization */
int nfb = 0;                    // Actual number of freqency intervals calculated

/* FBEMISS was calculated as follows:
x= 2. * PI * MELEC * BOLTZMANN / (H*H);
x=pow(x,-1.5);
x*=8. * PI / (C*C);
x*= H;
*/
#define FBEMISS   7.67413e-62   // Calculated with constants.c



/* These are external structures used primarily because we need to call 
Numerical Recipes routines from fb_verner and fb_topbase */

struct topbase_phot *fb_xtop;   //Topbase description of a photoionization x-section
double fbt;                     // Temperature at which thee emissivity is calculated
int fbfr;                       // fb_choice (see above)




/**************************************************************************
                    Space Telescope Science Institute
                                                                                                   
                                                                                                   
  Synopsis: fb_topbase_partial returns the partial (for a specific ion) emissivity or 
  recombination rate for ions described in terms of Topbase photoionization x-sections.
                                                                                                   
  Description:
                                                                                                   
  Arguments:  

  Some arguments are externally passed, including fbfr which determines whether
  one is computing the total emission (0), the reduced emission (1), or the rate
                                                                                                   
                                                                                                   
  Returns:
                                                                                                   
  Notes:
	The routines fb_verner_partial and fb_topbase_partial return the emissivity 
	at a specific freqency.  Because they are called by some NR recipes routines 
	that integrates over frequency, most of the information for these routines
	 has to be and is passed by the external structures. 
                                                                                                   
  History:
	02jul	ksl	Removed all references to the wind cell.
                                                                                                   
 ************************************************************************/


double
fb_topbase_partial (freq)
     double freq;
{
  int nion;
  double partial;
  double x;
  double gn, gion;
  double fthresh;

  fthresh = fb_xtop->freq[0];
  if (freq < fthresh)
    return (0.0);               // No recombination at frequencies lower than the threshold freq occur

  nion = fb_xtop->nion;

  /* JM -- below lines to address bug #195 */
  gn = 1;
  if (ion[nion].phot_info > 0)  // it's a topbase record
    gn = config[fb_xtop->nlev].g;
  else if (ion[nion].phot_info == 0)    // it's a VFKY record, so shouldn't really use levels
    gn = ion[nion].g;
  else
  {
    Error ("fb_topbase_partial: Did not understand cross-section type %i for ion %i. Setting multiplicity to zero!\n",
           ion[nion].phot_info, nion);
    gn = 0.0;
  }



  gion = ion[nion + 1].g;       // Want the g factor of the next ion up
  x = sigma_phot (fb_xtop, freq);
  // Now calculate emission using Ferland's expression

  partial = FBEMISS * gn / (2. * gion) * pow (freq * freq / fbt, 1.5) * exp (H_OVER_K * (fthresh - freq) / fbt) * x;

  // 0=emissivity, 1=heat loss from electrons, 2=photons emissivity

  if (fbfr == FB_REDUCED)
    partial *= (freq - fthresh) / freq;
  else if (fbfr == FB_RATE)
    partial /= (H * freq);
  
  

  return (partial);
}

/**************************************************************************
                    Space Telescope Science Institute
                                                                                                   
                                                                                                   
  Synopsis: integ_fb calculates the integrated emissivity of the plasma, or the number of 
recombinations per second of a particular ion.  
                                                                                                   
  Description:
                                                                                                   
  Arguments:  
	t		    The temperature at which the emissivity 
			    or recombination rate is calculated
	f1,f2		The frequency limits on the calculation 
			    of the emissivity (ignored if the number 
			    of recombinations is desired.
	nion		The ion for which the emissivity is returned
	fb_choice	A switch which determines exactly what is to
			    be returned: 
			    0- the full emissivity including
			    the energy associated associated with the
			    threshold
			    1- the (reduced) emissivity, e.g. excluding
			    the threshold energy.  This is the energy
			    associated with kinetic energy loss
			    2- the specific recombination rate.
mode	inner or outer shell

                                                                                                   
                                                                                                   
  Returns:
	The routine returns the specific emissivity, e.g. the emissivity and
	or recombination rate per electron and per ion.
                                                                                                   
  Notes:
	As written, in July02, the idea is that if the recombination coefficients
	have been calculated the program is going to return an interpolated
	coefficient, if not, it will calculate it from scratch.  (The later is 
	much slower if one has to do it a lot of times.

	???? Error -- There is definitely an error because the program does
	not support option 0, and this is needed for calculation of the
	relative numbers of photons by fb vs free free photons.  ksl 02jul
 
                                                                                                   
                                                                                                   
  History:
	02jul	ksl	Modified to elimate need to include information about the cell.
			At this point integ_fb yields answers per electron and per ion.
			Moved DENSITY_PHOT_MIN test out of integ_fb
	02jul	ksl	Original routines has been pushed down to xinteg_fb so
			that integ_fb can be modified to use stored values when
			desired.
                                                                                                   
 ************************************************************************/


double
integ_fb (t, f1, f2, nion, fb_choice, mode)
     double t;                  // The temperature at which to calculate the emissivity
     double f1, f2;             // The frequencies overwhich to integrate the emissivity
     int nion;                  // The ion for which the "specific emissivity is calculateed
     int fb_choice;             // 0=full, 1=reduced, 2= rate
     int mode;                  // 1- outer shell 2-inner shell
{
  double fnu;
  int n;

  if (mode == OUTER_SHELL)
  {

    if (fb_choice == FB_FULL)
    {
      for (n = 0; n < nfb; n++)
      {
        /* See if the frequencies correspond to one previously calculated */
        if (f1 == freebound[n].f1 && f2 == freebound[n].f2)
        {
          fnu = get_fb (t, nion, n, fb_choice, mode);
          return (fnu);
        }
      }
      /* If not calculate it here */
      fnu = xinteg_fb (t, f1, f2, nion, fb_choice);
      return (fnu);
    }
    else if (fb_choice == FB_REDUCED)
    {
      for (n = 0; n < nfb; n++)
      {
        /* See if the frequencies correspond to one previously calculated */
        if (f1 == freebound[n].f1 && f2 == freebound[n].f2)
        {
          fnu = get_fb (t, nion, n, fb_choice, mode);
          return (fnu);
        }
      }
      /* If not calculate it here */
      fnu = xinteg_fb (t, f1, f2, nion, fb_choice);
      return (fnu);
    }
    else if (fb_choice == FB_RATE)
    {
      /* See if the frequencies correspond to one previously calculated */
      if (nfb > 0)
      {
        fnu = get_nrecomb (t, nion, mode);
        return (fnu);
      }
      /* If not calculate it here */
      fnu = xinteg_fb (t, f1, f2, nion, fb_choice);
      return (fnu);
    }
    Error ("integ_fb: Unknown fb_choice(%d)\n", fb_choice);
    exit (0);
  }

  else if (mode == INNER_SHELL)           // inner shell
  {
    if (fb_choice == FB_FULL)
    {
      for (n = 0; n < nfb; n++)
      {
        /* See if the frequencies correspond to one previously calculated */
        if (f1 == freebound[n].f1 && f2 == freebound[n].f2)
        {
          fnu = get_fb (t, nion, n, fb_choice, mode);
          return (fnu);
        }
      }
      fnu = xinteg_inner_fb (t, f1, f2, nion, fb_choice);
      return (fnu);
    }
    else if (fb_choice == FB_REDUCED)
    {
      for (n = 0; n < nfb; n++)
      {
        /* See if the frequencies correspond to one previously calculated */
        if (f1 == freebound[n].f1 && f2 == freebound[n].f2)
        {
          fnu = get_fb (t, nion, n, fb_choice, mode);
          return (fnu);
        }
      }
      fnu = xinteg_inner_fb (t, f1, f2, nion, fb_choice);
      return (fnu);
    }
    else if (fb_choice == FB_RATE)
    {
      if (nfb > 0)
      {
        fnu = get_nrecomb (t, nion, mode);
        return (fnu);
      }
      fnu = xinteg_inner_fb (t, f1, f2, nion, fb_choice);
      return (fnu);
    }
    Error ("integ_fb: Unknown fb_choice(%d)\n", fb_choice);
    exit (0);
  }

  Error ("integ_fb: Unknown mode(%d)\n", mode);
  exit (0);

}




/**************************************************************************
                    Space Telescope Science Institute
                                                                                                   
                                                                                                   
  Synopsis: total_fb returns the energy lost from the plasma due to fb emission in a
	single wind cell at a temperature t between the frequncy limits f1 and f2.  
	The energy lost is just the kinetic energy lost from the plasma 
	because it does not include the ionization potential
	associated with each recombination.  Python tracks effectively the kinetic energy
	of the plasma (not the potential energy available if everything recombined. 
                                                                                                   
  Description:
                                                                                                   
  Arguments:  
                                                                                                   
                                                                                                   
  Returns:
                                                                                                   
  Notes:
                                                                                                   
                                                                                                   
  History:
  02jul	ksl	Modified so computes fb contributions of individual ions, instead
 		of just h, he, and z, and so that one does not need to pass
 		the entire wind cell to integ_fb 
  02jul	ksl	Added call to init_freebound so that precalculated values would
  		be used where possible.  This seemed the most conservative place
		to put these calls since total_fb is always called whenever the
		band-limited luminosities are needed.
                                                                                                   
 ************************************************************************/


double
total_fb (one, t, f1, f2, fb_choice, mode)
     WindPtr one;
     double t, f1, f2;
     int fb_choice;
     int mode;                  //inner=2 outer=1
{
  double total;
  int nion;
  int nplasma;
  PlasmaPtr xplasma;

  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];

  if (t < 100. || f2 < f1)
      t=100.;   /* Set the temperature to 100 K so that if there are free electrons emission by this process continues */ 

// Initialize the free_bound structures if that is necessary
  if (mode == OUTER_SHELL)
    init_freebound (100., 1.e9, f1, f2);        //NSH 140121 increased limit to take account of hot plasmas NSH 1706 -


// Calculate the number of recombinations whenever calculating the fb_luminosities
  num_recomb (xplasma, t, mode);

  total = 0;
  xplasma->cool_rr_metals = 0.0;


  for (nion = 0; nion < nions; nion++)
  {
    if (xplasma->density[nion] > DENSITY_PHOT_MIN)
    {
      if (mode == OUTER_SHELL)
      {
		  if (fb_choice == FB_FULL) // we are calculating a luminosity
		  {
             total += xplasma->lum_rr_ion[nion] = xplasma->vol * xplasma->ne * xplasma->density[nion + 1] * integ_fb (t, f1, f2, nion, fb_choice, mode);
	      }
	      else  // we are calculating a cooling rate
	      {
             total += xplasma->cool_rr_ion[nion] = xplasma->vol * xplasma->ne * xplasma->density[nion + 1] * integ_fb (t, f1, f2, nion, fb_choice, mode);
		  	 if (ion[nion].z > 3)
                  xplasma->cool_rr_metals += xplasma->cool_rr_ion[nion];
		  
	      } 
      }
      else if (mode == INNER_SHELL)  // at present we do not compute a luminosity from DR
        total += xplasma->cool_dr_ion[nion] =
          xplasma->vol * xplasma->ne * xplasma->density[nion + 1] * integ_fb (t, f1, f2, nion, fb_choice, mode);

    }

  }
  return (total);
}



/**************************************************************************
                    Space Telescope Science Institute
                                                                                                   
                                                                                                   
  Synopsis: one_fb generates one free bound photon with specific frequency limits
                                                                                                   
  Description:
                                                                                                   
  Arguments:  
	one 	The wind cell in which the photon is being 
		generated
	f1,f2	The frequency limits
                                                                                                   
                                                                                                   
  Returns:
	The frequency of the fb photon that was generated.
                                                                                                   
  Notes:
	57h -- This routine was a major time sync in python57g.  Most of the problem
	was ascociated with the generation of pdfs.  So the program now has a 
	rather multilayered approach to reducing the number of pdf generation
	steps that has to be done.  For a long time, we simply stored a pdf in
 	the pdf_fb array.  For python_57h, I created a new structure, photstoremain
	which parallels plasmamain and with allows one to store photons for future
	use.  
	
	It's possible the time for this routine to take could be reduced by storing
	more than one set of extra photons for different frequncy intervals, since
	our "banding" approach causes the frequencies to shift during the photon
	generation cycle.  Alternatively, it is possible that changing the order
	of photon generation could help.  But the approach adopted eliminate
	most of the problems with the routine, and so I have not pursued that.

	060802 -- ksl
                                                                                                   
                                                                                                   
  History:
	98	ksl	Coded as part of python effort
	98oct	ksl	Removed upper limits on freqency to attempt to resolve problems with
			different frequency limits in total_fb and one_fb
	01oct	ksl	Completely rewritten for topbase x-sections
	06jul	ksl	57h -- Modified to speed the program up by checking whether
			there are free_bound frequencies that have been previsulsy
			calculated for this cell and this condition.  Also, added
			a section which creates and stores multiple photons for
			the same conditions.  This reduces very significantly
			the number of times one has to construct a pdf, which is
			the main time sink for the program
	17jul	NSH - changed references from PDF to CDF
                                                                                                   
 ************************************************************************/


double fb_x[3000], fb_y[3000];
double fb_jumps[NLEVELS];       // There is at most one jump per level
double xfb_jumps[NLEVELS];     // This is just a dummy array that parallels fb_jumpts
int fb_njumps = (-1);

WindPtr ww_fb;
//struct Cdf cdf_fb;
double one_fb_f1, one_fb_f2, one_fb_te; /* Old values */

double
one_fb (one, f1, f2)
     WindPtr one;               /* a single cell */
     double f1, f2;             /* freqmin and freqmax */
{
  double freq, tt, delta;
  int n,nn,nnn;
  double fthresh, dfreq;
  int nplasma;
  PlasmaPtr xplasma;
  PhotStorePtr xphot;

  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];
  xphot = &photstoremain[nplasma];
  

  if (f2 < f1)
  {
    Error ("one_fb: f2 %g < f1 %g Something is rotten  t %g\n", f2, f1, xplasma->t_e);
    exit (0);
  }

/* Check if an apprpriate photon frequency has already been generated, and 
use that instead if possible --  57h */
  tt = xplasma->t_e;
  if (xphot->n < NSTORE && xphot->f1 == f1 && xphot->f2 == f2 && xphot->t == tt)
  {
    freq = xphot->freq[xphot->n];
    (xphot->n)++;
    return (freq);
  }

  delta = 500;                  // Fudge factor to prevent generation of a CDF if t has changed only slightly
  /* Check to see if we have already generated a cdf */
//  if (tt > (one_fb_te + delta) || tt < (one_fb_te - delta) || f1 != one_fb_f1 || f2 != one_fb_f2)
  if (tt!=one_fb_te || f1 != one_fb_f1 || f2 != one_fb_f2)
  
  {

/* Then need to generate a new cdf */

    ww_fb = one;

    /* Create the fb_array */

    /* Determine how many intervals are between f1 and f2.  These need to be
       put in increasing frequency order */

    if (f1 != one_fb_f1 || f2 != one_fb_f2)
    {                           // Regenerate the jumps 
      fb_njumps = 0;
      for (n = 0; n < nphot_total; n++)
      {                         //IS THIS ADDED BRACKET CORRECT? (SS, MAY04)
        fthresh = phot_top_ptr[n]->freq[0];
        if (f1 < fthresh && fthresh < f2)
        {
          fb_jumps[fb_njumps] = fthresh;
          fb_njumps++;
        }
      }                         //IS THIS CORRECT? (SS, MAY04)



      /* The next line sorts the fb_jumps by frequency and eliminates
       * duplicate frequencies which is what was causing the error in
       * cdf.c when more than one jump was intended
       */

	  if (fb_njumps > 1) //We only need to sort and compress if we have more than one jump
	  {
      fb_njumps=sort_and_compress(fb_jumps,xfb_jumps,fb_njumps);
      for (n=0;n<fb_njumps;n++){
          fb_jumps[n]=xfb_jumps[n];
	  }
      }


    }
	

    //!BUG SSMay04
    //It doesn't seem to work unless this is zero? (SS May04)
    // fb_njumps = 0;              // FUDGE (SS, May04)

    /* Note -- Need to fix this to get jumps properly, that is the
       frequencies need to allow for the jumps !! ??? */
	
	/*NSH 1707 - modified the loop below to ensure we have points just below and above any jumps */
	
  nnn=0;   //Zero the index for elements in the flux array
	nn=0;  //Zero the index for elements in the jump array
	  n=0;  //Zero the counting element for equally spaced frequencies
    dfreq = (f2 - f1) / 1999; //This is the frequency spacing for the equally spaced elements
    while (n < (2000))   //We keep going until n=199, which will give the maximum required frequency
    {
		freq=f1 + dfreq * n;  //The frequency of the arrayelement we would make in the normal rin of things
		if (freq > fb_jumps[nn] && nn<fb_njumps) //The element we were going to make has a frequency abouve the jump
		{
			fb_x[nnn]=fb_jumps[nn]*(1.-DELTA_V/(2.*C));  //We make one frequency point 1km/s below the jump
			fb_y[nnn]=fb (xplasma, xplasma->t_e, fb_x[nnn], nions, FB_FULL); //And the flux for that point
			nnn=nnn+1;			//increase the index of the created array
			fb_x[nnn]=fb_jumps[nn]*(1.+DELTA_V/(2*C));  //And one frequency point just above the jump
			fb_y[nnn]=fb (xplasma, xplasma->t_e, fb_x[nnn], nions, FB_FULL); //And the flux for that point
			nn=nn+1;    //We heave dealt with this jump - on to the next one
			nnn=nnn+1;  //And we will be filling the next array element next time
		}
		else  //We haven't hit a jump
		{
			if (freq > fb_x[nnn-1])  //Deal with the unusual case where the upper point in our 'jump' pair is above the next ragular point
				{
      			fb_x[nnn] = freq;   //Set the next array element frequency
      			fb_y[nnn] = fb (xplasma, xplasma->t_e, fb_x[nnn], nions, FB_FULL); //And the flux
				n=n+1;  //Increment the regular grid counter
	  	  		nnn=nnn+1; //Increment the generated array counter
  				}
  	  		else //We dont need to make a new point, the upper frequency pair of the last jump did the trick
  		  		{
	  			n=n+1;  //We only need to increment our regular grid counter
  		  		}
 	 	}
    }
	

	
	
	
	/* At this point, the variable nnn stores the number of points */
	
	
	/* Check to see if the new number of points will exceed the number of points allocated in the cdf - if so, extend */
	



    if (cdf_gen_from_array (&cdf_fb, fb_x, fb_y, nnn, f1, f2) != 0)
    {
      Error ("one_fb after error: f1 %g f2 %g te %g ne %g nh %g vol %g\n",
             f1, f2, xplasma->t_e, xplasma->ne, xplasma->density[1], one->vol);
      Error ("Giving up\n");
      exit (0);
    }
    one_fb_te = xplasma->t_e;
    one_fb_f1 = f1;
    one_fb_f2 = f2;             /* Note that this may not be the best way to check for a previous cdf */
  }

/* OK, we have not created a new cdf actually.  We are in a position to
generate photons */

  //Debug ("one_fb, got here 2\n");

/* First generate the photon we need */
  freq = cdf_get_rand (&cdf_fb);
  if (freq<f1 || freq > f2) {
      Error("one_fb:  freq %e  freqmin %e freqmax %e out of range\n",freq,f1,f2);
  }

/* Now create and store for future use a set of additonal photons */

  for (n = 0; n < NSTORE; n++)
  {
    xphot->freq[n] = cdf_get_rand (&cdf_fb);
  if (xphot->freq[n]<f1 || xphot->freq[n] > f2) {
      Error("one_fb:  freq %e  freqmin %e freqmax %e out of range\n",xphot->freq[n],f1,f2);
  }

  }
  xphot->n = 0;
  xphot->t = tt;
  xphot->f1 = f1;
  xphot->f2 = f2;
  return (freq);
}





/**************************************************************************
                    Space Telescope Science Institute
                                                                                                   
                                                                                                   
  Synopsis: num_recomb calculates the total number of recombinations in (units of #/cm**2/s) 
   in the cell per second for the all ions.    
                                                                                                   
  Description:
                                                                                                   
  Arguments:  

		mode - wether we are coputing inner shell or outer shell rates
                                                                                                   
                                                                                                   
  Returns:
                                                                                                   
  Notes:
   The calculation is made purely direct recombination using the photoionization
   x-sections.  
   
   ?? It is possible that there should be a stimulated emission correction.
   ?? It is possible that a check to see that nfb should have been made
                                                                                                   
                                                                                                   
  History:
	98jun	ck      Error in one rate corrected
	02jun	ksl	This is a completely new version of the routine and
   			uses detailed balance. (python_43.5)
	02jul	ksl	Modified so that integ_fb is the number of 
			recombinations per ne and per ion
	06may	ksl	57+ -- Modified to use plasma structure since on volume
    17jan	nsh 81	Added a mode parameter to allow the same code to work for both inner shell and outer shell recomb
                                                                                                   
 ************************************************************************/

int
num_recomb (xplasma, t_e, mode)
     PlasmaPtr xplasma;
     double t_e;
     int mode;
{
  int nelem;
  int i, imin, imax;
  for (nelem = 0; nelem < nelements; nelem++)
  {
    imin = ele[nelem].firstion;
    imax = imin + ele[nelem].nions;
    for (i = imin; i < imax; i++)
    {
      if (xplasma->density[i] > DENSITY_PHOT_MIN)
      {
        if (mode == OUTER_SHELL)          //outer shell
          xplasma->recomb[i] = xplasma->ne * xplasma->density[i + 1] * integ_fb (t_e, 0.0, VERY_BIG, i, FB_RATE, mode);
        else if (mode == INNER_SHELL)     //innershell
          xplasma->inner_recomb[i] = xplasma->ne * xplasma->density[i + 1] * integ_fb (t_e, 0.0, VERY_BIG, i, FB_RATE, mode);

      }
    }
    xplasma->recomb[imax] = 0.0;        // Can't recombine to highest i-state
    xplasma->inner_recomb[imax] = 0.0;  // Can't recombine to highest i-state

  }

  return (0);
}


/**************************************************************************
                    Space Telescope Science Institute
                                                                                                   
                                                                                                   
  Synopsis: fb calculates the free_bound emissivity of the plasma at a specific frequency
                                                                                                   
  Description:
                                                                                                   
  Arguments:  

	ion_choice	Either the total emissivity or the emissivity for a specific
			ion is caculated depending on whether ion_choice=nions, or 
			a value less than the total number of ions
    fb_choice determines whether what is returned is the emissivity a specific frecuency 0
            the emissivity reduced by the ionization potential 1, or the number of photons per
            unit Hz
                                                                                                   
  Returns:
                                                                                                   
  Notes:
                                                                                                   
                                                                                                   
  History:
	02jul	ksl	Modified to reflect desire to make fb_xtopbase and fb_verner
			independent of the wind cell.
	06may	ksl	57+ -- Switched to plasma structure since no volume
	06jul	ksl	57h -- Cleaned this routine up a bit, in part to avoid
			calling fb_verner_partial when it should not be called.
                                                                                                   
 ************************************************************************/


double
fb (xplasma, t, freq, ion_choice, fb_choice)
     PlasmaPtr xplasma;         // A cell with all its associated density data
     double t;                  // The temperature at which to calculate the emissivity
     double freq;               // The frequency at which to calculate the emissivity
     int ion_choice;            // Selects which ions the emissivity is to be calculated for (see above)
     int fb_choice;             // 0=emissivity in the standard sense, 1 heat loss from electons, 2 number of photons
{
  int n;
  double fnu, x;
  int nmin, nmax;               // These are the photo-ionization xsections that are used
  int nion, nion_min, nion_max;



  if (ion_choice < nions)       //Get emissivity for this specific ion_number
  {
    nion_min = ion_choice;
    nion_max = ion_choice + 1;
  }
  else if (ion_choice == nions) // Get the total emissivity
  {
    nion_min = 0;
    nion_max = nions;
  }
  else
  {
    Error ("fb: This choice %d for ion_choice is not supported\n", ion_choice);
    exit (0);
  }


  fbt = t;                      /* Externally transmitted variable */
  fbfr = fb_choice;             /* Externally transmitted variable */

  fnu = 0.0;                    /* Initially set the emissivity to zero */

  //Debug("in fb for ion_choice %i\n", ion_choice);

  for (nion = nion_min; nion < nion_max; nion++)
  {
    if (ion[nion].phot_info > 0)        // topbase or VFKY+topbase
    {
      nmin = ion[nion].ntop_first;
      nmax = nmin + ion[nion].ntop;
    }
    else if (ion[nion].phot_info == 0)  // VFKY 
    {
      nmin = ion[nion].nxphot;
      nmax = nmin + 1;
    }
    else
      nmin = nmax = 0;          // no XS / ionized - don't do anything 

    //Debug("in fb for ion %i info %i, nmin nmax %i, %i\n", nion, ion[nion].phot_info, nmin, nmax);

    x = 0.0;

    /* Loop over relevent Topbase photoionization x-sections.  If 
       an ion does not have Topbase photoionization x-sections then
       ntmin and ntmax are the same and the loop will be skipped. */

    for (n = nmin; n < nmax; n++)
    {
      fb_xtop = &phot_top[n];   /*Externally transmited to fb_topbase_partial */
      /* We don't want to include fb transitions associated with macro atoms here
         - they are separated out for now. (SS, Apr 04). "If" statement added. */
      if (fb_xtop->macro_info == 0 || geo.macro_simple == 1 || geo.rt_mode == RT_MODE_2LEVEL)
      {
        x += fb_topbase_partial (freq);
      }

//      fnu += xplasma->density[nion] * x;
    }


    /* x is the emissivity from this ion. Add it to the total */
    fnu += xplasma->density[nion+1] * x;
  }

  fnu *= xplasma->ne;           // Correct from specific emissivity to the total fb emissivity

  return (fnu);

}


/**************************************************************************
                    Space Telescope Science Institute
                                                                                                   
                                                                                                   
  Synopsis: init_freebound initializes the structure fb_struc as well as some
associated arrays and variables (found in python.h) that describe
recombination rates and band-limited luminosities.
                                                                                                   
  Description:
                                                                                                   
  Arguments:  
                                                                                                   
	t1, t2	The lower and upper limits for which the fb 
		information is calculated.
	f1, f2	The frequency interval in which the band-limited
		fb information is calculated.
                                                                                                   
  Returns:
                                                                                                   
  Notes:
	The first time the routine is called, both recombination
	rates and band-limited luminosities are calculated.  On
	subsequent calls the routine checks to see whether it has
	already calculated the band-limited freebound emissivities, 
	and if so returns without redoing the calculation.  However, 
	if a new frequency interval is provided, the new luminosities
	are added to the free-bound structure.  To force a 
	re-initialization nfb must be set to 0.
                                                                                                   
                                                                                                   
  History:
	02jul	ksl	Coding began
	0810	ksl	67 - Modified routine so that instead of exiting
			when there are more than NFB sets of data, it 
			creates a new set of data and assumes the oldest
			set can be discarded.  This was done primarily
			to accommodate some runs of balance.
                                                                                                   
 ************************************************************************/

int init_freebound_nfb;         /*Indicates the total number of freebound sets that
                                   could be used */

int
init_freebound (t1, t2, f1, f2)
     double t1, t2, f1, f2;
{
  double t;
  int i, j, nion;
  double ltmin, ltmax, dlt;
  double xinteg_fb ();
  int nput;


  if (nfb == 0)
  {
    if (t2 < t1)
    {
      Error ("init_freebound: t2(%g)<t1(%g)\n", t2, t1);
      exit (0);
    }

    ltmin = log10 (t1);
    ltmax = log10 (t2);
    dlt = (ltmax - ltmin) / (NTEMPS - 1);

    for (j = 0; j < NTEMPS; j++)
    {
      fb_t[j] = pow (10., ltmin + dlt * j);
    }

    Log ("init_freebound: Creating recombination coefficients\n");
    for (nion = 0; nion < nions; nion++)
    {

      for (j = 0; j < NTEMPS; j++)
      {
        t = fb_t[j];
        xnrecomb[nion][j] = xinteg_fb (t, 0.0, VERY_BIG, nion, FB_RATE);
        xninnerrecomb[nion][j] = xinteg_inner_fb (t, 0.0, VERY_BIG, nion, FB_RATE);
      }
    }
  }
  else if (fabs (fb_t[0] - t1) > 10. || fabs (fb_t[NTEMPS - 1] - t2) > 1000.)
  {
    Error ("init_freebound: Cannot initialize to new temps without resetting nfb");
    exit (0);

  }

/* Now check to see whether the freebound information has already
been calculated for these conditions, and if so simply return.
*/
  i = 0;
  while ((freebound[i].f1 != f1 || freebound[i].f2 != f2) && i < nfb)
    i++;

  if (i < nfb)
  {
    return (0);
  }

/* We have to calculate a new set of freebound data */
  if (i == NFB - 1)
  {
    /* We've filled all the available space in freebound so we start recycling elements, assuming that the latest
     * ones are still likelyt to be needed
     */
    nput = init_freebound_nfb % NFB;
    init_freebound_nfb++;

    Error ("init_freebound: Recycling freebound, storage for NFB (%d), need %d to avoid \n", NFB, init_freebound_nfb);

  }
  else
  {
    nput = init_freebound_nfb = nfb;
    nfb++;
  }


/* Having reached this point, a new set of fb emissivities
must be calculated.  Note that old information is not destroyed
unless nfb had been set to 0.  The new set is added to the old
on the assumption that the fb information will be reused.
*/


  Log ("init_freebound: Creating recombination emissivites between %e and %e\n", f1, f2);


  freebound[nput].f1 = f1;
  freebound[nput].f2 = f2;

  for (nion = 0; nion < nions; nion++)
  {
    for (j = 0; j < NTEMPS; j++)
    {                           //j covers the temps
      t = fb_t[j];
      freebound[nput].lum[nion][j] = xinteg_fb (t, f1, f2, nion, FB_FULL);
      freebound[nput].cool[nion][j] = xinteg_fb (t, f1, f2, nion, FB_REDUCED);
      freebound[nput].cool_inner[nion][j] = xinteg_inner_fb (t, f1, f2, nion, FB_REDUCED);

    }
  }

  return (0);
}




/**************************************************************************
                    Space Telescope Science Institute
                                                                                                   
                                                                                                   
  Synopsis: Return the recombination coefficient 
                                                                                                   
  Description:
                                                                                                   
  Arguments:  
                                                                                                   
                                                                                                   
  Returns:
                                                                                                   
  Notes:
                                                                                                   
                                                                                                   
  History:
		13sep	nsh	changed call to linterp to reflect new option
                                                                                                   
 ************************************************************************/


double
get_nrecomb (t, nion, mode)
     double t;
     int nion;
     int mode;
{
  int linterp ();
  double x;
  if (mode == 1)
    linterp (t, fb_t, xnrecomb[nion], NTEMPS, &x, 0);   //Interpolate in linear space
  else if (mode == 2)
    linterp (t, fb_t, xninnerrecomb[nion], NTEMPS, &x, 0);      //Interpolate in linear space
  else
  {
    Error ("Get_nrecomb - unkonwn mode %i", mode);
    exit (0);
  }
  return (x);
}


/* Return the specific emissivity due to recombination emission in an interval */

double
get_fb (t, nion, narray, fb_choice, mode)
     double t;
     int nion;
     int narray;
	 int fb_choice;
     int mode;
{
  int linterp ();
  double x;
  if (mode == OUTER_SHELL)
  {
	  if (fb_choice == FB_REDUCED)
    linterp (t, fb_t, &freebound[narray].cool[nion][0], NTEMPS, &x, 0);        //Interpolate in linear space
	  else if (fb_choice == FB_FULL)
    linterp (t, fb_t, &freebound[narray].lum[nion][0], NTEMPS, &x, 0);        //Interpolate in linear space
	  else
	  {
		  Error ("Get_fb - unexpected mode %i", mode);
		  exit(0);
		}
  }
  else if (mode == INNER_SHELL)
    linterp (t, fb_t, &freebound[narray].cool_inner[nion][0], NTEMPS, &x, 0);  //Interpolate in linear space

  else
  {
    Error ("Get_fb - unkonwn mode %i", mode);
    exit (0);
  }
  return (x);
}



/**************************************************************************
                    Space Telescope Science Institute
                                                                                                   
                                                                                                   
  Synopsis: xinteg_fb calculates the integrated emissivity of 
  an ion in the plasma.  
                                                                                                   
  Description:
                                                                                                   
  Arguments:  
                                                                                                   
                                                                                                   
  Returns:
                                                                                                   
  Notes:

  This routine is called by integ_fb.  It is not intended to be called 
  directly. 

  It's unusual nature is determined by the need to use a modififed 
  Numerical Recipes routine for integration of fb over a frequency range 
                                                                                                   
                                                                                                   
  History:
	02jul	ksl	Modified to elimate need to include information about the cell.
			At this point integ_fb yields answers per electron and per ion.
			Moved DENSITY_PHOT_MIN test out of integ_fb
	02jul	ksl	In attempt to store the fb coefficiencts this function
			renamed.  It actually calculates the emissivities, and integ_fb
			has become a kind of steering routine that normally, e.g. often
			reads the freebound arrays.
    04May   SS      Minor changes to exclude fb from macro atoms (which are treated elsewhere)
                                                                                                   
 ************************************************************************/


double
xinteg_fb (t, f1, f2, nion, fb_choice)
     double t;                  // The temperature at which to calculate the emissivity
     double f1, f2;             // The frequencies overwhich to integrate the emissivity
     int nion;                  // The ion for which the "specific emissivity is calculateed
     int fb_choice;             // 0=full, otherwise reduced
{
  int n;
  double fnu;
  double dnu;                   //NSH 140120 - a parameter to allow one to restrict the integration limits.
  double fthresh, fmax;
  double den_config ();
  int nmin, nmax;               // These are the limits over which number xsections we will use 
  double qromb ();


  dnu = 0.0;                    //Avoid compilation errors.

  if (-1 < nion && nion < nions)        //Get emissivity for this specific ion_number
  {
    if (ion[nion].phot_info > 0)        // topbase or hybrid
    {
      nmin = ion[nion].ntop_first;
      nmax = nmin + ion[nion].ntop;
    }
    else if (ion[nion].phot_info == 0)  // VFKY 
    {
      nmin = ion[nion].nxphot;
      nmax = nmin + 1;
    }
    else
      // the ion is a fullt ionized ion / doesn't have a cross-section, so return 0
      return (0.0);
  }
  else                          // Get the total emissivity
  {
    Error ("integ_fb: %d is unacceptable value of nion\n", nion);
    exit (0);
  }

  // Put information where it can be used by the integrating function
  fbt = t;
  fbfr = fb_choice;

  /* Limit the frequency range to one that is reasonable before integrating */

  if (f1 < 3e12)
    f1 = 3e12;                  // 10000 Angstroms
  if (f2 > 3e18)                // 110819 nsh increase upper limits to include  highly ionised ions that we are now seeing in x-ray illuminated nebulas.
    f2 = 3e18;                  // This is 1 Angstroms  - ksl
  if (f2 < f1)
    return (0.0);               /* Because there is nothing to integrate */

  fnu = 0.0;


  for (n = nmin; n < nmax; n++)
  {
    // loop over relevent Topbase or VFKY photoionzation x-sections
    fb_xtop = &phot_top[n];

    /* Adding an if statement here so that photoionization that's part of a macro atom is 
       not included here (these will be dealt with elsewhere). (SS, Apr04) */
    if (fb_xtop->macro_info == 0 || geo.macro_simple == 1 || geo.rt_mode == RT_MODE_2LEVEL)  //Macro atom check. (SS)
    {
      fthresh = fb_xtop->freq[0];
      fmax = fb_xtop->freq[fb_xtop->np - 1];    // Argues that this should be part of structure
      if (f1 > fthresh)
        fthresh = f1;
      if (f2 < fmax)
        fmax = f2;
      // Now calculate the emissivity as long as fmax exceeds xthreshold and there are ions to recombine
      if (fmax > fthresh)
      {
        //NSH 140120 - this is a test to ensure that the exponential will not go to zero in the integrations 
        dnu = 100.0 * (fbt / H_OVER_K);
        if (fthresh + dnu < fmax)
        {
          fmax = fthresh + dnu;
        }
        fnu += qromb (fb_topbase_partial, fthresh, fmax, 1.e-4);
      }
    }
  }


  return (fnu);
}


/**************************************************************************
                   Southampton University
                                                                                                   
                                                                                                   
  Synopsis: xinteg_inner_fb calculates the integrated fb emissivity of inner
  shell transitions in an ion at a given temperature
                                                                                                   
  Description: 
                                                                                                   
  Arguments:  
                                                                                                   
                                                                                                   
  Returns:
                                                                                                   
  Notes:

  This routine is  virtual copy of xinteg_fb but dones inner shell integrations.

  This routine is called by integ_fb.  It is not intended to be called 
  directly. 

  It's unusual nature is determined by the need to use a modififed 
  Numerical Recipes routine for integration of fb over a frequency range 
                                                                                                   
                                                                                                   
                                                                                                   
  History:
	16Jan NSH - coded
                                                                                                   
 ************************************************************************/


double
xinteg_inner_fb (t, f1, f2, nion, fb_choice)
     double t;                  // The temperature at which to calculate the emissivity
     double f1, f2;             // The frequencies overwhich to integrate the emissivity
     int nion;                  // The ion for which the "specific emissivity is calculateed
     int fb_choice;             // 0=full, otherwise reduced
{
  int n, nn;
  double fnu;
  double dnu;                   //NSH 140120 - a parameter to allow one to restrict the integration limits.
  double fthresh, fmax;
  double den_config ();
  double qromb ();


  dnu = 0.0;                    //Avoid compilation errors.
  fnu = 0.0;
  nn = -1;


  if (f1 < 3e12)
    f1 = 3e12;                  // 10000 Angstroms
  if (f2 > 3e18)                // 110819 nsh increase upper limits to include  highly ionised ions that we are now seeing in x-ray illuminated nebulas.
    f2 = 3e18;                  // This is 1 Angstroms  - ksl
  if (f2 < f1)
    return (0.0);               /* Because there is nothing to integrate */

  for (n = 0; n < n_inner_tot; n++)
  {
    if (inner_cross[n].nion == nion)
    {
      nn = n;
      fbt = t;
      fbfr = fb_choice;

      /* Limit the frequency range to one that is reasonable before integrating */



      // loop over relevent Topbase or VFKY photoionzation x-sections
      fb_xtop = &inner_cross[nn];

      /* Adding an if statement here so that photoionization that's part of a macro atom is 
         not included here (these will be dealt with elsewhere). (SS, Apr04) */
      if (fb_xtop->macro_info == 0 || geo.macro_simple == 1 || geo.rt_mode == RT_MODE_2LEVEL)        //Macro atom check. (SS)
      {
        fthresh = fb_xtop->freq[0];
        fmax = fb_xtop->freq[fb_xtop->np - 1];  // Argues that this should be part of structure
        if (f1 > fthresh)
          fthresh = f1;
        if (f2 < fmax)
          fmax = f2;

        // Now calculate the emissivity as long as fmax exceeds xthreshold and there are ions to recombine
        if (fmax > fthresh)
        {
          //NSH 140120 - this is a test to ensure that the exponential will not go to zero in the integrations 
          dnu = 100.0 * (fbt / H_OVER_K);
          if (fthresh + dnu < fmax)
          {
            fmax = fthresh + dnu;
          }
          fnu += qromb (fb_topbase_partial, fthresh, fmax, 1.e-4);
        }

      }
    }
  }



  return (fnu);
}



/***********************************************************
                                       Southampton University
                                                                                                                                      
 Synopsis:
        total_rr(nion, T)
                                                                                                                                      
Arguments:
        ion - ion for which we want a recombination rate - 
		this is the upper state, so the ion which is 
		doing the recombining, there is no rate for
		H1(ion0) but there is one for H(ion1)
	temperature - the temperature we want a rate for
Returns:
	rate - the total recombination rate for this ion 
		for this temperature
                                                                                                                                      
Description:
                                                                                                                                      
        This routine generates a total recombination rate for 
		a given ion at a given temperature using 
		badnell or shull type parameters. If these
		are not presnet, an error is produced but the code
		soldiers on with a value from the milne relation.
                                                                                                                                      
Notes:
	
                                                                                                                                      
History:
        12jul   nsh     73 -- Began coding
	24jul	nsh	73 -- Included the shull coefficients in the chianti database
	14aug	nsh	78b-- renamed - from bad_t_rr since we dont just use badnell data.
			Also rewritten to use the milne relation to get a value for the 
			recombination rate in the absence of data. This is all in preparation
			for the use of this routine to help populate a recombination rate matrix.
	
                                                                                                                                      
**************************************************************/

double
total_rrate (nion, T)
     int nion;
     double T;
{


  double rate;                  //The returned rate
  double rrA, rrB, rrT0, rrT1, rrC, rrT2;       //The parameters
  double term1, term2, term3;   //Some temporary parameters to make calculation simpler


  rate = 0.0;                   /* NSH 130605 to remove o3 compile error */


  if (ion[nion].total_rrflag == 1)      /*We have some kind of total radiative rate data */
  {
    if (total_rr[ion[nion].nxtotalrr].type == RRTYPE_BADNELL)
    {
      rrA = total_rr[ion[nion].nxtotalrr].params[0];
      rrB = total_rr[ion[nion].nxtotalrr].params[1];
      rrT0 = total_rr[ion[nion].nxtotalrr].params[2];
      rrT1 = total_rr[ion[nion].nxtotalrr].params[3];
      rrC = total_rr[ion[nion].nxtotalrr].params[4];
      rrT2 = total_rr[ion[nion].nxtotalrr].params[5];


      rrB = rrB + rrC * exp ((-1.0 * rrT2) / T);        //If C=0, this does nothing


      term1 = sqrt (T / rrT0);
      term2 = 1.0 + sqrt (T / rrT0);
      term2 = pow (term2, (1 - rrB));
      term3 = 1.0 + sqrt (T / rrT1);
      term3 = pow (term3, (1 + rrB));


      rate = pow ((term1 * term2 * term3), -1.0);
      rate *= rrA;
    }
    else if (total_rr[ion[nion].nxtotalrr].type == RRTYPE_SHULL)
    {
      rate = total_rr[ion[nion].nxtotalrr].params[0] * pow ((T / 1.0e4), -1.0 * total_rr[ion[nion].nxtotalrr].params[1]);
    }
    else
    {
      Error ("total_rrate: unknown parameter type for ion %i\n", nion);
      exit (0);                 /* NSH This is a serious problem! */
    }
  }
  else                          /*NSH 140812 - We dont have coefficients - in this case we can use xinteg_fb with mode 2 to use the milne relation to obtain a value for this - it is worth throwing an error though, since there rreally should be data for all ions. xinteg_fb
                                   is called with the lower ion in the pair, since it uses the photionization cross sectiuon of the lower ion */
  {
    Error ("total_rrate: No T_RR parameters for ion %i - using milne relation\n", nion);
    rate = xinteg_fb (T, 3e12, 3e18, nion - 1, 2);
  }


  return (rate);

}


/***********************************************************
                                       Southampton University
                                                                                                                                      
 Synopsis:
        gs_rr(nion,T)
                                                                                                                                      
Arguments:
        ion - ion for which we want a recombination rate - 
		this is the upper state, so the ion which is 
		doing the recombining, there is no rate for
		H1(ion0) but there is one for H(ion1)
	temperature - the temperature we want a rate for
Returns:
	rate - the resolved recombination rate for the ground
		state of this ion recombining into the GS of
		the lower ion.
		for this temperature
                                                                                                                                      
Description:
                                                                                                                                      
        This routine generates a recombination rate to the ground state for 
		a given ion at a given temperature using 
		badnell type parameters. If these parameters are not
		available for the given ion - the milne relation is
		used.   
                                                                                                                                      
Notes:
	
                                                                                                                                      
History:
        12jul   nsh     73 -- Began coding
  	14mar	nsh	77a-- Interpolaion now carried out in log space
	14aug	nsh	78b-- Renamed to gs_rr from bad_gs_rr and 
			rewritten to use the milne relation if badnell
			type paramerters are not available. This allows
			this code to be used to produce recombination
			rate coefficients for the matrix ionization scheme.

	
                                                                                                                                      
**************************************************************/

double
gs_rrate (nion, T)
     int nion;
     double T;
{
  double rate, drdt, dt;
  int i, imin, imax;
  double rates[BAD_GS_RR_PARAMS], temps[BAD_GS_RR_PARAMS];
  int ntmin;
  double fthresh, fmax, dnu;


  imin = imax = 0;              /* NSH 130605 to remove o3 compile error */


  //  if (ion[nion].bad_gs_rr_t_flag != 1 && ion[nion].bad_gs_rr_r_flag != 1)
  //    {
  //      Error ("bad_gs_rr: Insufficient GS_RR parameters for ion %i\n", nion);
  //      return (0);
  //    }

  if (ion[nion].bad_gs_rr_t_flag == 1 && ion[nion].bad_gs_rr_r_flag == 1)       //We have tabulated gs data

    //NSH force code to always use milne for a test REMOVE ME!!!
    //if (ion[nion].bad_gs_rr_t_flag == 100 && ion[nion].bad_gs_rr_r_flag == 100)       //We have tabulated gs data
  {
	  for (i = 0; i < BAD_GS_RR_PARAMS; i++)
    {
      rates[i] = bad_gs_rr[ion[nion].nxbadgsrr].rates[i];
      temps[i] = bad_gs_rr[ion[nion].nxbadgsrr].temps[i];
    }

    if (T < temps[0])           //we are below the range of GS data
    {
      Log_silent ("bad_gs_rr: Requested temp %e is below limit of data for ion %i(Tmin= %e)\n", T, nion, temps[0]);
      //      rate = rates[0];
      imax = 1;
      imin = 0;
    }

    else if (T >= temps[BAD_GS_RR_PARAMS - 1])  //we are above the range of GS data
    {
      Log_silent
        ("bad_gs_rr: Requested temp %e is above limit (%e) of data for ion %i\n",
         T, nion, bad_gs_rr[ion[nion].nxbadgsrr].temps[BAD_GS_RR_PARAMS - 1]);
      imax = BAD_GS_RR_PARAMS - 1;
      imin = BAD_GS_RR_PARAMS - 2;
      //We will try to extrapolate.

    }
    else                        //We must be within the range of tabulated data
    {
      for (i = 0; i < BAD_GS_RR_PARAMS - 1; i++)
      {
        if (temps[i] <= T && T < temps[i + 1])  //We have bracketed the correct temperature
        {
          imin = i;
          imax = i + 1;
        }
      }
      /* NSH 140313 - changed the following lines to interpolate in log space */
    }
    drdt = (log10 (rates[imax]) - log10 (rates[imin])) / (log10 (temps[imax]) - log10 (temps[imin]));
    dt = (log10 (T) - log10 (temps[imin]));
    rate = pow (10, (log10 (rates[imin]) + drdt * dt));
  }

  /* we will need to use the milne relation - 
     NB - this is different from using xinteg_fb because 
     that routine does recombination to all excited levels (at least for topbase ions).
   */
  else
  {
    rate = 0.0;                 /* NSH 130605 to remove o3 compile error */

    fbt = T;
    fbfr = FB_RATE;

    if (ion[nion - 1].phot_info > 0)    //topbase or hybrid
    {
      ntmin = ion[nion - 1].ntop_ground;
      fb_xtop = &phot_top[ntmin];
    }
    else if (ion[nion - 1].phot_info == 0)      //vfky 
    {
      fb_xtop = &phot_top[ion[nion - 1].nxphot];
    }

    fthresh = fb_xtop->freq[0];
    fmax = fb_xtop->freq[fb_xtop->np - 1];
    dnu = 100.0 * (fbt / H_OVER_K);

    if (fthresh + dnu < fmax)
    {
      fmax = fthresh + dnu;
    }

    rate = qromb (fb_topbase_partial, fthresh, fmax, 1e-5);
  }

  return (rate);
}



/* This routine takes an input array of doubles and sorts it into
 * numberical order.  It then eliminates duplicates.  
 * The input array is not destroyed
 */

int
sort_and_compress (array_in, array_out, npts)
     double *array_in, *array_out;
     int npts;
{
  double *values;
  int n, nfinal;
  int compare_doubles ();

  values = calloc (sizeof (double), npts);
  for (n = 0; n < npts; n++)
    {
      values[n] = array_in[n];
    }

  /* Sort the array in place */
  qsort (values, npts, sizeof (double), compare_doubles);
  

  array_out[0]=values[0]; //Copy the first jump into the output array  
  
  nfinal = 1;
  for (n = 1; n < npts; n++) //Loop over the remaining jumps in the array
    {
      if (values[n] > array_out[nfinal - 1])  //In the next point in the array is larger than the last one (i.e. not equal)
	{
	  array_out[nfinal] = values[n]; //Put the next point into the array
	  nfinal += 1;  //Increment the size of the array
	}
    }
	


  return (nfinal);
}

/* This routine just compares two double precision numbers and
 * returns 1 if a is greater than b, and 0 otherwise.  It is
 * used by qsort
 */



int
compare_doubles (const void *a, const void *b)
{
  if (*(double*)a > *(double*)b)
    return 1;
  else if (*(double*)a < *(double*)b)
    return -1;
  else
    return 0;
}
