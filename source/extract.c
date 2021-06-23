
/***********************************************************/
/** @file  extract.c
 * @author ksl
 * @date   April, 2018
 *
 * @brief  Routines for extracting photons during the 
 * spectral generation phase.
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


/**********************************************************/
/** 
 * @brief      A supervisory routine called to 
 * 	builds detailed spectra in the normal (extract) mode.
 *
 * @param [in] WindPtr  w   The entire wind
 * @param [in] PhotPtr  p   The photon to extract
 * @param [in] int  itype   An integer representing the type of photon 
 * for the purpose of being extracted.
 * @return     Always returns 0
 *
 * @details
 * extract is called when a photon begins its flight and every time that photon
 * scatters, unless the user has exercised the "live or die" option, in
 * which case it is not called.  
 *
 * extract carries out several preperatory steps for the extraction
 * and for each spectrum one wants to build, it calls extract_one,
 * where the actual incrementing of the spectrum is done.
 *
 * itype takes on the following values:
 * * PTYPE_STAR->the photon came for the star 
 * * PTYPE_BL->the photon came from the boundary layer
 * * PTYPE_DISK->the photon being redirected arose in the disk,
 * * PTYPE_WIND->the photon being redirected arose in the wind,
 * 
 * extract uses the types to prepare the photon for extraction, including
 * doppler shifting the photon if is of PTYPE_WIND or PTYPE_DISK.
 * 
 * Usually, Python constructs a spectrum of all photons, but there are
 * advanced options which allone to restrict the spectrum created to
 * those produced with a certain number of scatters or from photons 
 * that arise from the above or below the disk.  extract enforces
 * those choices before calling extract_one.
 * The parameters for this option all come in through
 * python.h, and are contained in the spectrum structure. 
 *
 * The basic choices, some of which can be used in tandom are as follows:
 *
 * *	If s[n].nscat>999; then one obtains the entire spectrum
 * *	if s[n].nscat is a positive number, then the spectrum is composed just 
 * 		of photons with that number of scatters
 * *	if s[n].nscat is a negative number, then the spectrum contains all photons
 * 		with >=  the absolute value of s[n].nscat
 * * 	s[n].top_bot=0 -> accept all photons
 * * 	s[n].top_bot    >0  -> accept photons last above the disk
 * * 	s[n].top_bot    <0  -> accept photons below the disk
 * *    One can also select photons that are to be extracted from a particular
 * spherical region.
 * ### Notes ###
 * 
 *
 **********************************************************/
double xlmn[3];                 // XTEST

int
extract (w, p, itype)
     WindPtr w;
     PhotPtr p;
     int itype;
{
  int n, mscat, mtopbot;
  struct photon pp, p_in, p_dummy;
  int extract_photon;
  double xdiff[3];
  double p_norm, tau_norm;
  double dvds_max;
  int ierr;


  ierr = check_frame (p, F_OBSERVER, "extract_start");
  if (ierr)
  {
    Error ("extract: check_frame failure at very start for itype %d\n", itype);
  }

  stuff_phot (p, &p_in);

  /* We increase weight to account for the number of scatters. This is done because in extract we multiply by the escape
     probability along a given direction, but we also need to divide the weight by the mean escape probability, which is
     equal to 1/nnscat.  See issue #710 for a more extended explanation of how the weight is renormalized stocahstically. */

  if (itype == PTYPE_WIND)
  {
    if (geo.scatter_mode == SCATTER_MODE_THERMAL && p_in.nres <= NLINES && p_in.nres > -1)
    {
      /* we normalised our rejection method by the escape probability along the vector of maximum velocity gradient.
         First find the sobolev optical depth along that vector. The -1 enforces calculation of the ion density 
       */

      dvds_max = get_dvds_max (&p_in);
      tau_norm = sobolev (&wmain[p_in.grid], p_in.x, -1.0, lin_ptr[p_in.nres], dvds_max);

      /* then turn into a probability */
      p_norm = p_escape_from_tau (tau_norm);

    }
    else
    {
      p_norm = 1.0;

      if (p_in.nnscat != 1)
        Error
          ("trans_phot: nnscat is %i for photon %i in scatter mode %i! nres %i NLINES %i\n",
           p_in.nnscat, p_in.np, geo.scatter_mode, p_in.nres, NLINES);
    }

    p_in.w *= p_in.nnscat / p_norm;
  }

  if (itype == PTYPE_WIND)
  {
    if ((ierr = observer_to_local_frame (&p_in, &p_in)))
      Error ("extract: wind photon not in observer frame %d\n", ierr);
  }
  if (itype == PTYPE_DISK)
  {
    if ((ierr = observer_to_local_frame_disk (&p_in, &p_in)))
      Error ("extract: disk photon not in observer frame %d\n", ierr);
  }

  /* At this point were are in a local frame for WIND and DISK photons, but the
     global frame for the central source and boundary layer */

  for (n = MSPEC; n < nspectra; n++)
  {
    /* If statement allows one to choose whether to construct the spectrum
       from all photons or just from photons that have scattered a specific number
       of times or in specific regions of the wind. A region is specified by a position
       and a radius. */

    extract_photon = TRUE;

    if ((mscat = xxspec[n].nscat) > 999 || p_in.nscat == mscat || (mscat < 0 && p_in.nscat >= (-mscat)))
      extract_photon = TRUE;
    else
      extract_photon = FALSE;

    if (extract_photon)
    {
      if ((mtopbot = xxspec[n].top_bot) == 0)
        extract_photon = TRUE;  // Then there are no positional parameters and we are done
      else if (mtopbot == -1 && p_in.x[2] < 0)
        extract_photon = TRUE;
      else if (mtopbot == 1 && p_in.x[2] > 0)
        extract_photon = TRUE;
      else if (mtopbot == 2)    // Then to count, the photom must originate within sn.r of sn.x
      {
        vsub (p_in.x, xxspec[n].x, xdiff);
        if (length (xdiff) > xxspec[n].r)
          extract_photon = FALSE;
      }
      else
        extract_photon = FALSE;
    }

    if (extract_photon)
    {
      /* Create a photon pp to use here and in extract_one, and send it in
       * the correct direction.  This assures we
       * have not modified p_in as part of extract.  Also, allow for aberration
       * of photons to assure that we are extracting at the correct angle
       * in the observer frame
       *
       * Note that a wind photon is in the local frame, and the only thing
       * we need to do is to figure out the extraction direction in the local
       * frame that will produce the directon we want
       *
       * A disk photon is in the observer frame
       */

      stuff_v (xxspec[n].lmn, xlmn);    //XTEST

      if ((rel_mode == REL_MODE_FULL || rel_mode == REL_MODE_SR_FREQ) && itype == PTYPE_WIND)
      {
        stuff_phot (&p_in, &p_dummy);
        p_dummy.frame = F_OBSERVER;
        stuff_v (xxspec[n].lmn, p_dummy.lmn);
        observer_to_local_frame (&p_dummy, &p_dummy);
        stuff_phot (&p_in, &pp);
        stuff_v (p_dummy.lmn, pp.lmn);
      }
      else if ((rel_mode == REL_MODE_FULL || rel_mode == REL_MODE_SR_FREQ) && itype == PTYPE_DISK)
      {
        stuff_phot (&p_in, &p_dummy);
        p_dummy.frame = F_OBSERVER;
        stuff_v (xxspec[n].lmn, p_dummy.lmn);
        observer_to_local_frame_disk (&p_dummy, &p_dummy);
        stuff_phot (&p_in, &pp);
        stuff_v (p_dummy.lmn, pp.lmn);
      }
      else
      {
        stuff_phot (&p_in, &pp);
        stuff_v (xxspec[n].lmn, pp.lmn);
      }

      /*At this point photons of type DISK or WIND are in the local frame, but others
         are in the global frame
       */

      if (modes.save_photons && 1180. < 2.997925e18 / pp.freq && 2.997925e18 / pp.freq < 1240.0)
      {
        save_photons (&pp, "Extract_start");
      }

      if (modes.save_extract_photons && 1545.0 < 2.997925e18 / pp.freq && 2.997925e18 / pp.freq < 1565.0)
      {
        save_extract_photons (n, p, &pp);
      }

      extract_one (w, &pp, itype, n);
    }
  }

  return (0);
}





/**********************************************************/
/** 
 * @brief      Extract a single photon along a single line of sight.
 *
 * @param [in] WindPtr  w   The entire wind
 * @param [in] PhotPtr  pp  The photon to be extracted
 * @param [in] int  itype   The type of photon (star, disk, wind, etc)
 * @param [in] int  nspec   the spectrum which will be incremented
 * @return     The photon status after translation
 *
 * @details
 * extract_one is analogous to the detailed portion of transphot except here the
 * basic point is to calculate the optical depth through the plasma in a certain
 * direction, and to increment the appropriate spectrum.  
 *
 * Unlike trans_phot, this routine also checks to see if 
 * whether the photon hits the secondary star, if one exists.
 *
 * ### Notes ###
 * The logic behind the weighting of the photons is described in Christian Knigge's thesis in
 * section 2.3.4.  According to equation 2.19
 * 	Pc/Pw=12 cos(theta)*(1+b cos(theta)/(3+2b) where b=1.5 corresponds to the
 * Eddington approximation.
 * 
 * In Python, and in extract and transphot in particular, tau generally refers to the tau associated
 * with scattering processes, and the weight contains the effect of dimunition of the energy of
 * the photon bundle due to pure absorption processes.  So, in extract, we add pp->w * exp(-tau)
 * to the spectrum.
 *
 * Note that both linearly and logarithmically spaced spectra are produced.
 *
 **********************************************************/

int
extract_one (w, pp, itype, nspec)
     WindPtr w;
     PhotPtr pp;
     int itype, nspec;

{
  int istat, nres;
  struct photon pstart;
  struct photon pdummy, pdummy_orig;
  double weight_min;
  int icell;
  int k, k1;
  double x[3];
  double tau;
  double zz;
  double dvds;
  double lfreqmin, lfreqmax, ldfreq;
  int ishell;
  double normal[3];
  int ierr;


  weight_min = EPSILON * pp->w;
  istat = P_INWIND;
  tau = 0;

  // We want pstart to be in the observer frame

  if (itype == PTYPE_WIND)
  {
    ierr = local_to_observer_frame (pp, &pstart);
    if (ierr)
      Error ("extract_one: pp of type WIND not in local frame %d\n", ierr);
  }
  else if (itype == PTYPE_DISK)
  {
    ierr = local_to_observer_frame_disk (pp, &pstart);
    if (ierr)
      Error ("extract_one: pp of type DISK not in local frame %d\n", ierr);
  }
  else
  {
    ierr = check_frame (pp, F_OBSERVER, "extract_one: start");
    if (ierr)
    {
      Error ("extract_one: check_frame_failure for itype %d\n", itype);
    }
    stuff_phot (pp, &pstart);
  }

  /*
   * At this stage, we need to transform the photon back into the observer
   * frame
   */

  if (itype == PTYPE_DISK)
  {
    if ((ierr = local_to_observer_frame_disk (pp, pp)))
      Error ("extract_one: disk photon not in local frame");
  }
  if (itype == PTYPE_WIND)
  {
    stuff_phot (pp, &pdummy);

    if ((ierr = local_to_observer_frame (pp, pp)))
      Error ("extract_one: wind photon not in local frame\n");

    if (pp->x[2] * pstart.x[2] < 0)
    {
      Error ("Extract_one: Went through xz plane on local2observer frame\n");
      Error ("Extract_one: start %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
             pstart.x[0], pstart.x[1], pstart.x[2], pstart.lmn[0], pstart.lmn[1], pstart.lmn[2]);
      Error ("Extract_one:   was %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
             pdummy.x[0], pdummy.x[1], pdummy.x[2], pdummy.lmn[0], pdummy.lmn[1], pdummy.lmn[2]);
      Error ("Extract_one:    is %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
             pp->x[0], pp->x[1], pp->x[2], pp->lmn[0], pp->lmn[1], pp->lmn[2]);
    }

//OLD    if (run_ztest == FALSE)
//OLD    {
//OLD      if ((ierr = reposition (pp)))
//OLD      {
//OLD        Error ("extract_one: reposition returned error %d\n", ierr);    //BEFORE
//OLD      }                         //AFTER
//OLD    }
    pstart.lmn[0] = pp->lmn[0];
    pstart.lmn[1] = pp->lmn[1];
    pstart.lmn[2] = pp->lmn[2];
  }

  /* Re-weight the photons. Note that photons have already been frequency
   * shifted prior to entering extract. For disk and central object photons,
   * see Eqn 2.19 Knigge's thesis
   */

  if (itype == PTYPE_STAR || itype == PTYPE_BL || itype == PTYPE_AGN)
  {
    stuff_v (pp->x, x);
    renorm (x, 1.);
    zz = fabs (dot (x, xxspec[nspec].lmn));
    pp->w *= zz * (2.0 + 3.0 * zz);
  }
  else if (itype == PTYPE_DISK)
  {
    zz = fabs (xxspec[nspec].lmn[2]);
    pp->w *= zz * (2.0 + 3.0 * zz);
  }
  else if (pp->nres > -1 && pp->nres < NLINES)
  {

    /* It was a wind photon.  In this case, what we do depends
       on whether it is a photon which arose via line radiation or
       some other process.

       If geo.scatter_mode==SCATTER_MODE_ISOTROPIC then there is no need
       to reweight.  This is the isotropic assumption.  Otherwise, one
       needs to reweight

       //OLD       Once the photon is reweighted, reposition it so it does not interact ith
       //OLD       the same resonance a second time.
     */

    if (geo.scatter_mode == SCATTER_MODE_THERMAL)
    {
      dvds = dvwind_ds_cmf (pp);
      ishell = pp->grid;
      tau = sobolev (&w[ishell], pp->x, -1.0, lin_ptr[pp->nres], dvds);
      if (tau > 0.0)
        pp->w *= p_escape_from_tau (tau);
      tau = 0.0;
    }

    /* XXX - It is unclear why reposition needs to be here, but at present this
     * produces  better agreement with live or die than below */
//OLD    if (run_ztest == TRUE)
//OLD    {
//OLD      /* Fudge to be able to reposition photon while it is in the local frame */
//OLD      pp->frame = F_OBSERVER;
//OLD      if ((ierr = reposition (pp)))
//OLD        Error ("extract_one: reposition returned error %d\n", ierr);    //BEFORE
//OLD      pp->frame = F_LOCAL;
//OLD    }
//HOLD    if (pp->x[2] * pstart.x[2] < 0)
//HOLD    {
//HOLD      Error ("Extract_one: Went through xz plane on reposition\n");
//HOLD      Error ("Extract_one: start %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
//HOLD             pstart.x[0], pstart.x[1], pstart.x[2], pstart.lmn[0], pstart.lmn[1], pstart.lmn[2]);
//HOLD      Error ("Extract_one:    is %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
//HOLD             pp->x[0], pp->x[1], pp->x[2], pp->lmn[0], pp->lmn[1], pp->lmn[2]);
//HOLD    }
  }

  stuff_phot (pp, &pdummy);     // Actually not clear to me why we do/need this dummy photon at this stage

  if (tau > TAU_MAX)
  {
    istat = P_ABSORB;           /* Check to see if tau already too large */
    Error ("extract: tau should not be large\n");
  }
  else if (geo.binary == TRUE)
  {
    istat = hit_secondary (pp); /* Check to see if it hit secondary */
//HOLD    if (istat)
//HOLD    {
//HOLD      Log ("extract: Hit secondary %3d %5d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
//HOLD           geo.pcycle, pp->np, pp->x[0], pp->x[1], pp->x[2], pp->lmn[0], pp->lmn[1], pp->lmn[2]);
//HOLD    }
  }

//HOLD  if (pp->x[2] * pstart.x[2] < 0)
//HOLD  {
//HOLD    Error ("Extract_one: Went through xz plane after hit secondary \n");
//HOLD    Error ("Extract_one: start %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
//HOLD           pstart.x[0], pstart.x[1], pstart.x[2], pstart.lmn[0], pstart.lmn[1], pstart.lmn[2]);
//HOLD    Error ("Extract_one:   was %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
//HOLD           pdummy.x[0], pdummy.x[1], pdummy.x[2], pdummy.lmn[0], pdummy.lmn[1], pdummy.lmn[2]);
//HOLD    Error ("Extract_one:    is %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
//HOLD           pp->x[0], pp->x[1], pp->x[2], pp->lmn[0], pp->lmn[1], pp->lmn[2]);
//HOLD  }


  if (modes.save_photons && 1180. < 2.997925e18 / pp->freq && 2.997925e18 / pp->freq < 1240.0)
  {
    save_photons (pp, "AfterRepoObs");
  }

  double xdot, xxang;
  xdot = dot (pp->lmn, xlmn);
  xxang = acos (xdot) * 57.29578;
  if (xxang > 0.1)
  {

    Log ("gotcha type %d nphot %8d lmn  %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e not equal theta %.6f\n",
         itype, pp->np, xlmn[0], xlmn[1], xlmn[2], pp->lmn[0], pp->lmn[1], pp->lmn[2], xxang);

  }

//HOLD  pp->lmn[0] = xlmn[0];
//HOLD  pp->lmn[1] = xlmn[1];
//HOLD  pp->lmn[2] = xlmn[2];





/* Preserve the starting position of the photon so one can use this to determine whether the
 * photon encountered the disk or star as it tried to exist the wind.
 */



  if (modes.save_photons && 1180. < 2.997925e18 / pp->freq && 2.997925e18 / pp->freq < 1240.0)
  {
    save_photons (pp, "BeforeExtract");
  }

  tau = 0;
  pp->ds = 0;
  icell = 0;
/* Now we can actually extract the reweighted photon */

//HOLD  istat = walls (pp, &pstart, normal);
  stuff_phot (pp, &pdummy_orig);

  ierr = check_frame (pp, F_OBSERVER, "extract_one: photon not in observer frame at start");
  ierr = walls (&pdummy_orig, &pstart, normal);
  if (pdummy.istat != pp->istat)
  {
    Error ("extract_one: Surprising state change made by walls %d _> %d\n", pp->istat, pdummy.istat);
  }

  while (istat == P_INWIND)
  {
    istat = translate (w, pp, 20., &tau, &nres);
    icell++;

    stuff_phot (pp, &pdummy);
    istat = walls (pp, &pstart, normal);
    if (istat == -1)
    {

      Error ("Extract_one: Abnormal return from translate (icell %d) of phot no %5d\n", icell, pp->np);
      Error ("Extract_one: start %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
             pstart.x[0], pstart.x[1], pstart.x[2], pstart.lmn[0], pstart.lmn[1], pstart.lmn[2]);
      Error ("Extract_one:  orig %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
             pdummy_orig.x[0], pdummy_orig.x[1], pdummy_orig.x[2], pdummy_orig.lmn[0], pdummy_orig.lmn[1], pdummy_orig.lmn[2]);
      Error ("Extract_one:   was %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
             pdummy.x[0], pdummy.x[1], pdummy.x[2], pdummy.lmn[0], pdummy.lmn[1], pdummy.lmn[2]);
      Error ("Extract_one:    is %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
             pp->x[0], pp->x[1], pp->x[2], pp->lmn[0], pp->lmn[1], pp->lmn[2]);
      break;
    }

    if (pp->w < weight_min)
    {
      istat = P_ABSORB;         /*This photon was absorbed within the wind */
      break;
    }

    if (istat == P_HIT_STAR)
    {                           /* It was absorbed in the photosphere */
      break;
    }
    if (istat == P_HIT_DISK)
    {                           /* It was absorbed in the disk */
      break;
    }
    if (istat == P_SCAT)
    {                           /* Cause the photon to scatter and reinitilize */
      break;
    }
  }

  if (modes.save_photons && 1180. < 2.997925e18 / pp->freq && 2.997925e18 / pp->freq < 1240.0)
  {
    save_photons (pp, "AfterExtract");
  }

  if (istat == P_ESCAPE)
  {

    if (!(0 <= tau && tau < 1.e4))
      Error_silent ("Warning: extract_one: ignoring very high tau  %8.2e at %g\n", tau, pp->freq);
    else
    {
      k = (int) ((pp->freq - xxspec[nspec].freqmin) / xxspec[nspec].dfreq);

      /* Force the frequency to be in range of that recorded in the spectrum */

      if (k < 0)
        k = 0;
      else if (k > NWAVE_EXTRACT - 1)
        k = NWAVE_EXTRACT - 1;


      lfreqmin = log10 (xxspec[nspec].freqmin);
      lfreqmax = log10 (xxspec[nspec].freqmax);
      ldfreq = (lfreqmax - lfreqmin) / NWAVE_EXTRACT;

      /* find out where we are in log space */
      k1 = (int) ((log10 (pp->freq) - log10 (xxspec[nspec].freqmin)) / ldfreq);
      if (k1 < 0)
      {
        k1 = 0;
      }
      if (k1 > NWAVE_EXTRACT - 1)
      {
        k1 = NWAVE_EXTRACT - 1;
      }

      /* Increment the spectrum.  Note that the photon weight has not been diminished
       * by its passage through th wind, even though it may have encounterd a number
       * of resonance, and so the weight must be reduced by tau
       */

      xxspec[nspec].f[k] += pp->w * exp (-(tau));
      xxspec[nspec].lf[k1] += pp->w * exp (-(tau));



      /* If this photon was a wind photon, then also increment the "reflected" spectrum */
      if (pp->origin == PTYPE_WIND || pp->origin == PTYPE_WIND_MATOM || pp->nscat > 0)
      {

        xxspec[nspec].f_wind[k] += pp->w * exp (-(tau));
        xxspec[nspec].lf_wind[k1] += pp->w * exp (-(tau));

      }


      /* Records the total distance travelled by extracted photon if in reverberation mode */
      if (geo.reverb != REV_NONE)
      {
        if (pstart.nscat > 0 || pstart.origin > 9 || (pstart.nres > -1 && pstart.nres < nlines))
        {                       //If this photon has scattered, been reprocessed, or originated in the wind it's important
          pstart.w = pp->w * exp (-(tau));
          stuff_v (xxspec[nspec].lmn, pstart.lmn);
          delay_dump_single (&pstart, nspec);
        }
      }


    }

  }


  if (istat > -1 && istat < 9)
    xxspec[nspec].nphot[istat]++;
  else
    Error
      ("Extract: Abnormal photon %d %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e\n",
       istat, pp->x[0], pp->x[1], pp->x[2], pp->lmn[0], pp->lmn[1], pp->lmn[2]);

  return (istat);
}
