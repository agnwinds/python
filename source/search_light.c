


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

  	These are the routines that are needed to create photons 
	for the so-called search_light option that is a light 
	source for a situation where one wants to emulate a plane
	parallel option into ptyon

  Description:	

  Arguments: 		


  Returns:

  Notes:


  History:

  	1605	ksl	Work begun on a tranatlantic flight from Chicago
			to Tokyo

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

Synopsis: search_light_init (r, lum, alpha, freqmin, freqmax, ioniz_or_final, f)

 Arguments:              
 Returns:
 Description:    

 	Get the information needed to create a search_light source

Notes:
History:
**************************************************************/


int
search_light_init ()
{


  return (0);			/* Return the luminosity    */
}



/***********************************************************
Space Telescope Science Institute

Synopsis: 

Generate photons for a search light as defined bpy search_light_init()

 Arguments:              
 Returns:
 Description:    

 	Get the information needed to create a search_light source

Notes:
History:
**************************************************************/

int
photo_gen_search_light (p, r, alpha, weight, f1, f2, spectype, istart, nphot)
     PhotPtr p;
     double r, alpha, weight;
     double f1, f2;   /* The freqency mininimum and maximum if a uniform distribution is selected */
     int spectype;    /*The spectrum type to generate: 0 is bb, 1 (or in fact anything but 0)
           is uniform in frequency space */
     int istart, nphot;   /* Respecitively the starting point in p and the number of photons to generate */
{
  double freqmin, freqmax, dfreq, t;
  int i, iend;
  int n;
  double ftest;

  t = alpha;      /* NSH 130605, this is a slightly odd statmenent, put in to avoid 03 compilation 
           errors, but stemming from the odd behaviour that one can give the AGN a 
           themal spectral type, and the alpha parameter is then used to transmit 
           the temperature. */


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

  /* XXX - this line had been deleted from agn.c in domain, but it still exists in dev, so adding it back
   * as part of test of template_ionloop.pf.  It looks like agn.c in the two places have diverged */

  r = (1. + EPSILON) * r;       /* Generate photons just outside the photosphere unnecessary for the AGN perhaps? */

  /* Generate photons just outside the photosphere unnecessary for the AGN perhaps? */
  /* note this is only used in the spherical source */

  if (spectype == SPECTYPE_CL_TAB)  /* If we have a broken PL, we need to work out which alpha to use */
    {
      ftest = (freqmin + freqmax) / 2.0;
      for (n = 0; n < xband.nbands; n++)  /* work out what alpha and constant to use */
  {
    if (xband.f1[n] < ftest && xband.f2[n] > ftest)
      {
        alpha = xband.alpha[n];
      }
  }
    }
	
	


  for (i = istart; i < iend; i++)
    {
      p[i].origin = PTYPE_AGN;  // For BL photons this is corrected in photon_gen 
      p[i].w = weight;
      p[i].istat = p[i].nscat = p[i].nrscat = 0;
      p[i].grid = 0;
      p[i].tau = 0.0;
      p[i].nres = -1;   // It's a continuum photon
      p[i].nnscat = 1;

      if (spectype == SPECTYPE_BB)
  {
    p[i].freq = planck (t, freqmin, freqmax);
  }
      else if (spectype == SPECTYPE_UNIFORM)
  {     /* Kurucz spectrum */
    /*Produce a uniform distribution of frequencies */
    p[i].freq = freqmin + rand () * dfreq;
  }
      else if (spectype == SPECTYPE_POW)  /* this is the call to the powerlaw routine 
               we are most interested in */
  {
    p[i].freq = get_rand_pow (freqmin, freqmax, alpha);
  }
      else if (spectype == SPECTYPE_CL_TAB)
	{
	  p[i].freq = get_rand_pow (freqmin, freqmax, alpha);
	}
	else if (spectype == SPECTYPE_BREM)
	{
		p[i].freq = get_rand_brem(freqmin,freqmax);
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


    /* first option is for the original spherical X-ray source as used in e.g. Higginbottom+ 2013 */

      if (geo.pl_geometry == PL_GEOMETRY_SPHERE)
  {
    randvec (p[i].x, r);

    /* Added by SS August 2004 for finite disk. */
    if (geo.disk_type == DISK_VERTICALLY_EXTENDED)
      {
        /* JM XXX -- is this bit right? it seems to be that zdisk should use the x coordinate rather than
           magnitude of vector (r) */
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



      /* if we have a lamp post geometry then we should generate isotropic photons at a height
         above the disk plane */
      else if (geo.pl_geometry == PL_GEOMETRY_LAMP_POST)
  {
    /* x and y coordinates are 0 */
    p[i].x[0] = p[i].x[1] = 0.0;

    /* need to set the z coordinate to the lamp post height, but allow it to be above or below */
    if (rand () > MAXRAND / 2)
      {     /* Then the photon emerges in the upper hemisphere */
        p[i].x[2] = geo.lamp_post_height;
      }
    else
      {
        p[i].x[2] = -geo.lamp_post_height;
      }

    randvec (p[i].lmn, 1.0);  // lamp-post geometry is isotropic, so completely random vector
  }

    }

  return (0);
}
