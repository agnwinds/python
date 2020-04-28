
/***********************************************************/
/** @file  disk_init.c 
 * @author ksl
 * @date   April, 2020
 *
 * @brief  Primary routines for initializing the disk
 * as on a luminoisity weighted basis
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


#define STEPS 100000



/**********************************************************/
/**
 * @brief      calculates the total luminosity and the luminosity between freqqmin and freqmax
 * 	of the disk.  More importantly  divides the disk into annulae such that each
 * 	annulus contributes and equal amount to the lumionosity of the disk (within the frequency
 * 	limits).  Thus  initializes the structure "disk".
 *
 * @param [in] double  rmin   The minimum radius of the disk
 * @param [in] double  rmax   The maximum radius of the disk
 * @param [in] double  m   mass of central object
 * @param [in] double  mdot   mass accretion rate
 * @param [in] double  freqmin   The minimum frequency
 * @param [in] double  freqmax   The maximum frequency
 * @param [in] int  ioniz_or_final   A flag indicating whether this is an ionization or
 * a detailed spectral cycle (used to determine the spectral type to use)
 * @param [out] double *  ftot   The band limited luminosity in the freqency interval
 * @return     the total luminosity of the disk
 *
 * @details
 * This routine assumes the temperature distribution for the disk is
 * that of a simple Shakura-Sunyaev disk, and uses this to determine
 * the band limited luminosity of the disk.  Additionally, it divides
 * the disk in the rings of the same band-limited luminosity, so that
 * equal numbers of photons can be generated from each ring.  (The
 * reason the disk has to be initilaized mulitple times is because
 * the rings are different for different freqency intervals.)
 *
 * ### Notes ###
 * The information needed to generate photons from the disk is stored
 * in the disk structure.
 * The positional parameters x and v are at the edge of the ring,
 * but many of the other parameters (like temperature) are at the mid point.
 *
 *
 **********************************************************/

double
disk_init (rmin, rmax, m, mdot, freqmin, freqmax, ioniz_or_final, ftot)
     double rmin, rmax, m, mdot, freqmin, freqmax, *ftot;
     int ioniz_or_final;
{
  double t, tref, teff (), tdisk ();
  double log_g, gref, geff (), gdisk ();
  double dr, r;
  double logdr, logrmin, logrmax, logr;
  double f, ltot;
  double q1;
  int nrings, i, icheck;
  int spectype;
  double emit, emittance_bb (), emittance_continuum ();

  /* Calculate the reference temperature and luminosity of the disk */
  tref = tdisk (m, mdot, rmin);


  gref = gdisk (m, mdot, rmin);

  q_test_count = 0;
  /* Now compute the apparent luminosity of the disk.  This is not actually used
     to determine how annulae are set up.  It is just used to populate geo.ltot.
     It can change if photons hitting the disk are allowed to raise the temperature
   */

  logrmax = log (rmax);
  logrmin = log (rmin);
  logdr = (logrmax - logrmin) / STEPS;

  for (nrings = 0; nrings < NRINGS; nrings++)   //Initialise the structure
  {
    disk.nphot[nrings] = 0;
    disk.nphot[nrings] = 0;
    disk.r[nrings] = 0;
    disk.t[nrings] = 0;
    disk.nhit[nrings] = 0;
    disk.heat[nrings] = 0;
    disk.ave_freq[nrings] = 0;
    disk.w[nrings] = 0;
    disk.t_hit[nrings] = 0;
  }




  ltot = 0;

  for (logr = logrmin; logr < logrmax; logr += logdr)
  {
    r = exp (logr);
    dr = exp (logr + logdr) - r;
    t = teff (tref, (r + 0.5 * dr) / rmin);
    ltot += t * t * t * t * (2. * r + dr) * dr;
  }
  geo.lum_disk_init = ltot *= 2. * STEFAN_BOLTZMANN * PI;


  /* Now establish the type of spectrum to create */

  if (ioniz_or_final == 1)
    spectype = geo.disk_spectype;       /* type for final spectrum */
  else
    spectype = geo.disk_ion_spectype;   /*type for ionization calculation */

/* Next compute the band limited luminosity ftot */

/* The area of an annulus is  PI*((r+dr)**2-r**2) = PI * (2. * r +dr) * dr.
   The extra factor of two arises because the disk radiates from both of its sides.
   */

  q1 = 2. * PI;

  (*ftot) = 0;
  icheck = 0;


  for (logr = logrmin; logr < logrmax; logr += logdr)
  {
    r = exp (logr);
    dr = exp (logr + logdr) - r;
    t = teff (tref, (r + 0.5 * dr) / rmin);
    log_g = log10 (geff (gref, (r + 0.5 * dr) / rmin));

    if (spectype > -1)
    {                           // emittance from a continuum model
      emit = emittance_continuum (spectype, freqmin, freqmax, t, log_g);
    }
    else
    {
      emit = emittance_bb (freqmin, freqmax, t);

    }
    (*ftot) += emit * (2. * r + dr) * dr;
  }

  (*ftot) *= q1;



  /* If *ftot is 0 in this energy range then all the photons come elsewhere, e. g. the star or BL  */

  if ((*ftot) < EPSILON)
  {
    Log_silent ("disk_init: Warning! Disk does not radiate enough to matter in this wavelength range\n");
    return (ltot);
  }

  /* Now find the boundaries of the each annulus, which depends on the band limited flux.
     Note that disk.v is calculated at the boundaries, because vdisk() interporlates on
     the actual radius. */

  disk.r[0] = rmin;
  disk.v[0] = sqrt (GRAV * geo.mstar / rmin);
  nrings = 1;
  f = 0;

  i = 0;
  for (logr = logrmin; logr < logrmax; logr += logdr)
  {
    r = exp (logr);
    dr = exp (logr + logdr) - r;
    t = teff (tref, (r + 0.5 * dr) / rmin);
    log_g = log10 (geff (gref, (r + 0.5 * dr) / rmin));

    if (spectype > -1)
    {                           // continuum emittance
      emit = emittance_continuum (spectype, freqmin, freqmax, t, log_g);
    }
    else
    {
      emit = emittance_bb (freqmin, freqmax, t);
    }

    f += q1 * emit * (2. * r + dr) * dr;
    i++;
    /* EPSILON to assure that roundoffs don't affect result of if statement */
    if (f / (*ftot) * (NRINGS - 1) >= nrings)
    {
      if (r <= disk.r[nrings - 1])      //If the radius we have reached is smaller than or equal to the last assigned radius - we make a tiny annulus
      {
        r = disk.r[nrings - 1] * (1. + 1.e-10);
      }
      disk.r[nrings] = r;
      disk.v[nrings] = sqrt (GRAV * geo.mstar / r);
      nrings++;
      if (nrings >= NRINGS)
      {
//        Error_silent ("disk_init: Got to ftot %e at r %e < rmax %e. OK if freqs are high\n", f, r, rmax);             Not *really* an error, the error below deals with a *real* problem.
        break;
      }
    }
  }
  if (nrings < NRINGS - 1)
  {
    Error ("error: disk_init: Integration on setting r boundaries got %d nrings instead of %d\n", nrings, NRINGS - 1);
    Exit (0);
  }


  disk.r[NRINGS - 1] = exp (logrmax);
  disk.v[NRINGS - 1] = sqrt (GRAV * geo.mstar / disk.r[NRINGS - 1]);


  /* Now calculate the temperature and gravity of the annulae */

  for (nrings = 0; nrings < NRINGS - 1; nrings++)
  {
    r = 0.5 * (disk.r[nrings + 1] + disk.r[nrings]);
    disk.t[nrings] = teff (tref, r / rmin);
    disk.g[nrings] = geff (gref, r / rmin);
  }

  /* Wrap up by zerrowing other parameters */
  for (nrings = 0; nrings < NRINGS; nrings++)
  {
    disk.nphot[nrings] = 0;
    disk.nhit[nrings] = 0;
    disk.heat[nrings] = 0;
    disk.ave_freq[nrings] = 0;
    disk.w[nrings] = 0;
    disk.t_hit[nrings] = 0;
  }
  geo.lum_disk = ltot;
  return (ltot);
}


