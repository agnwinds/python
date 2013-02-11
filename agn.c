


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  These are routines needed to implement AGN into python

  Description:	

  Arguments:  


  Returns:

  Notes:


  History:
10oct	nsh	oded as part of initial effort to include a power law component
		to AGN

 ************************************************************************/

#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "atomic.h"
#include "python.h"


#include "log.h"

/***********************************************************
Space Telescope Science Institute

Synopsis: star_init (r, tstar, freqmin, freqmax, ioniz_or_final, f)
 Arguments:              
 Returns:
 Description:    
 
This routine calculates the luminosity of the star and the luminosity within the frequency boundaries. 
BB functions are assumed 
Notes:
History:
**************************************************************/



/* This is essentially a parallel routine that is set up for other types of sources.  It actually does not do very
 * much
 */

double
agn_init (r, lum, alpha, freqmin, freqmax, ioniz_or_final, f)
     double r, lum, alpha, freqmin, freqmax, *f;
     int ioniz_or_final;
{

  double t;
  double emit, emittance_bb (), emittance_continuum ();
  int spectype;


  if (ioniz_or_final == 1)
    spectype = geo.agn_spectype;	/* type for final spectrum */
  else
    spectype = geo.agn_ion_spectype;	/*type for ionization calculation */

  if (spectype >= 0)
    {
      /* Assume that if we simulate the continuum for an AGN that emit is the luminosity in a specific range */
      emit = emittance_continuum (spectype, freqmin, freqmax, lum, alpha);
      *f = emit;
    }
  else if (spectype == SPECTYPE_POW)
    {
      /* Emittance_pow actucally returns the specific luminosity directly */
      emit = emittance_pow (freqmin, freqmax, lum, alpha);
      *f = emit;
    }
  else if (spectype == SPECTYPE_BB)
    {
      t = alpha;		/* For this case the second variable is t */
      /* Note that emittane_bb is really the energy radiated per unit area  so
       * this must be multiplied by the area of the object */
      emit = emittance_bb (freqmin, freqmax, t);
      *f = emit;
      *f *= (4. * PI * r * r);
    }


  return (*f);			/* Return the luminosity    */
}


/* This routine returns the specific luminosity, that is the luminosity within the frequency interval */
double
emittance_pow (freqmin, freqmax,lum,alpha)
     double freqmin, freqmax,lum,alpha;
{
  double constant, emit;
  /* these are the frequencies over which the power law is defined - currently set to roughly what they 
     end up as when python is run normally with an accretion disk */

#define   XFREQMIN  1e14
#define   XFREQMAX  1e17

  /* first we need to calculate the constant for the power law function */
  constant =
    lum / (((pow (XFREQMAX, alpha + 1.)) - pow (XFREQMIN, alpha + 1.0)) /
	   (alpha + 1.0));

  /* now we need to work out the luminosity between our limited frequency range */
  /* we may need some checking routines to make sure that the requested frequency range is within the defined range,
     or it could default to zero outside the defined range */

  emit =
    constant * ((pow (freqmax, alpha + 1.0) - pow (freqmin, alpha + 1.0)) /
		(alpha + 1.0));

  return (emit);
}




int
photo_gen_agn (p, r, alpha, weight, f1, f2, spectype, istart, nphot)
     PhotPtr p;
     double r, alpha, weight;
     double f1, f2;		/* The freqency mininimum and maximum if a uniform distribution is selected */
     int spectype;		/*The spectrum type to generate: 0 is bb, 1 (or in fact anything but 0)
				   is uniform in frequency space */
     int istart, nphot;		/* Respecitively the starting point in p and the number of photons to generate */
{
  double freqmin, freqmax, dfreq, t;
  int i, iend;
  double dot ();
  double planck ();
  double plaw ();
  int randvec (), randvcos ();
  double zdisk ();
  if ((iend = istart + nphot) > NPHOT)
    {
      Error ("photo_gen_agn: iend %d > NPHOT %d\n", iend, NPHOT);
      exit (0);
    }
  if (f2 < f1)
    {
      Error
	("photo_gen_agn: Cannot generate photons if freqmax %g < freqmin %g\n",
	 f2, f1);
    }
  Log_silent ("photo_gen_agn creates nphot %5d photons from %5d to %5d \n",
	      nphot, istart, iend);
  freqmin = f1;
  freqmax = f2;
  dfreq = (freqmax - freqmin) / MAXRAND;
  r = (1. + EPSILON) * r;	/* Generate photons just outside the photosphere unnecessary for the AGN perhaps? */
  for (i = istart; i < iend; i++)
    {
      p[i].origin = PTYPE_AGN;	// For BL photons this is corrected in photon_gen 
      p[i].w = weight;
      p[i].istat = p[i].nscat = p[i].nrscat = 0;
      p[i].grid = 0;
      p[i].tau = 0.0;
      p[i].nres = -1;		// It's a continuum photon
      p[i].nnscat = 1;

      if (spectype == SPECTYPE_BB)
	{
	  t = alpha;
	  p[i].freq = planck (t, freqmin, freqmax);
	}
      else if (spectype == SPECTYPE_UNIFORM)
	{			/* Kurucz spectrum */
	  /*Produce a uniform distribution of frequencies */
	  p[i].freq = freqmin + rand () * dfreq;
	}
      else if (spectype == SPECTYPE_POW)	/* this is the call to the powerlaw routine we are most interested in */
	{
	  p[i].freq = get_rand_pow (freqmin, freqmax, alpha);
	}
      else
	{
	  p[i].freq =
	    one_continuum (spectype, t, geo.gstar, freqmin, freqmax);
	}

      if (p[i].freq < freqmin || freqmax < p[i].freq)
	{
	  Error_silent
	    ("photo_gen_agn: phot no. %d freq %g out of range %g %g\n",
	     i, p[i].freq, freqmin, freqmax);
	}

      randvec (p[i].x, r);

      /* Added by SS August 2004 for finite disk. */
      if (geo.disk_type == 2)
	{
	  while (fabs (p[i].x[2]) < zdisk (r))
	    {
	      randvec (p[i].x, r);
	    }


	  if (fabs (p[i].x[2]) < zdisk (r))
	    {
	      Error ("Photon_agn: agn photon %d in disk %g %g %g %g %g\n",
		     i, p[i].x[0], p[i].x[1], p[i].x[2], zdisk (r), r);
	      exit (0);
	    }
	}


      /* this last bit is the direction, might need a change */
      randvcos (p[i].lmn, p[i].x);
    }
  return (0);
}
