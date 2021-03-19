/***********************************************************/
/** @file  matom.c
 * @author SS,JM,SWM
 * @date   February, 2004
 * @brief  Macro-atom functions
 *
 * File containing Macro-atom functions.
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "atomic.h"
#include "python.h"

/* External variables for use in matom */

double jprbs_known[NLEVELS_MACRO][2 * (NBBJUMPS + NBFJUMPS)], eprbs_known[NLEVELS_MACRO][2 * (NBBJUMPS + NBFJUMPS)];
double pjnorm_known[NLEVELS_MACRO], penorm_known[NLEVELS_MACRO];
int prbs_known[NLEVELS_MACRO];
int matom_cell = -1;
int matom_z = -1;
int matom_cycle = -1;

/**********************************************************/
/** 
 * @brief The core of the implementation of Macro Atoms in python
 *
 * @param [in,out]  PhotPtr p   the packet at the point of activation and deactivation
 * @param [in,out]  int nres    the process which activates and deactivates the Macro Atom
 * @param [out]  int escape  flag to tell us whether the matom de-activation
 *                             is via an r-packet (1) or a k-packet (0)
 * @return returns the number of level jumps before exiting the routine, or -1
 * if an error was encountered.
 *
 * Matom is the core of the implementation of Leon Lucy's Macro Atom
 * approach to dealing with radiation-matter interactions. It is called whenever
 * a photon packet activates a macro atom. As input it takes the process by
 * which activation occurred (identified by the label "nres") which allow it 
 * to deduce the level that has been excited. It then calculates all the 
 * Macro Atom jumping/deactivation probabilities following Lucy and so deduces
 * the process by which deactivation occurs. 
 * 
 * The two possibilities are that the macro-atom deactivates by an r-packet
 * (radiative deactivations) or a k-packet (putting the energy  into the
 * thermal pool) as indicated by the variable escape.
 * 
 * At the point where the routine exits, nres idenifies the process of
 * deactivation and the packet information is also updated 
 * 
 *
 * ###Notes###
 * 
 * in macro atom mode, nres = -1 indicates electron scattering, 
 * nres = -2 indicates ff, and nres > NLINES indicates bound-free. 
 * [nres == NLINES is never used. Note also that NLINES is the *max* number of lines, whereas nlines
 * is the *actual* number of lines. So, actually, it's not just nres = NLINES that's never used, but 
 * the entire range of nlines <= nres <= NLINES]
 * 
 * 
 * 
***********************************************************/

int
matom (p, nres, escape)
     PhotPtr p;
     int *nres;
     int *escape;
{
  struct lines *line_ptr;
  struct topbase_phot *cont_ptr;
  int uplvl, uplvl_old;
  int icheck;
  double jprbs[2 * (NBBJUMPS + NBFJUMPS)];
  double eprbs[NBBJUMPS + NBFJUMPS];
  double pjnorm, penorm;
  double threshold, run_tot;
  double sp_rec_rate;
  int n, m;
  int njumps;
  int nbbd, nbbu, nbfd, nbfu;
  double t_e, ne;
  double bb_cont, choice, bf_cont;
  WindPtr one;
  double rad_rate, coll_rate, lower_density, density_ratio;
  PlasmaPtr xplasma;
  MacroPtr mplasma;
  int z;
//OLD  double jprbs_known[NLEVELS_MACRO][2 * (NBBJUMPS + NBFJUMPS)], eprbs_known[NLEVELS_MACRO][2 * (NBBJUMPS + NBFJUMPS)];
//OLD  double pjnorm_known[NLEVELS_MACRO], penorm_known[NLEVELS_MACRO];
//OLD  int prbs_known[NLEVELS_MACRO];


//OLD  for (n = 0; n < NLEVELS_MACRO; n++)
//OLD  {
//OLD    prbs_known[n] = -1;         //flag all as unknown
//OLD  }


  one = &wmain[p->grid];
  xplasma = &plasmamain[one->nplasma];
  check_plasma (xplasma, "matom");

  mplasma = &macromain[xplasma->nplasma];


  t_e = xplasma->t_e;
  ne = xplasma->ne;

  /* these are used later for stimulated recomb */
  lower_density = density_ratio = 0.0;


  /* The first step is to identify the configuration that has been excited.
   * If *nres < NLINES the level will have been excited by a bb transioion
   * but if it is greater than this, it will bave been produced by photoionization.
   */

  uplvl = 0;
  if (*nres < NLINES)
  {
    uplvl = lin_ptr[*nres]->nconfigu;
    z = lin_ptr[*nres]->z;

  }
  else if (*nres > NLINES)
  {
    uplvl = phot_top[*nres - NLINES - 1].uplev;
    z = phot_top[*nres - NLINES - 1].z;
  }
  else
  {
    Error ("matom: upper level not identified. nres = %d in photon %d of cycle %d/%d in thread %d\n",
           *nres, p->np, geo.wcycle, geo.pcycle, rank_global);
    *escape = TRUE;
    p->istat = P_ERROR_MATOM;
    return (-1);
  }

  if (z != matom_z || p->grid != matom_cell || geo.wcycle != matom_cycle)
  {
    for (n = 0; n < NLEVELS_MACRO; n++)
    {
      prbs_known[n] = FALSE;    //flag all as unknown
    }
    matom_z = z;
    matom_cell = p->grid;
    matom_cycle = geo.wcycle;
  }
//OLD  else
//OLD  {
//OLD    Log ("Using old for photon in grid %d and element %e\n", p->grid, z);
//OLD  }

  /* Now follows the main loop to govern the macro atom jumps. Keeps jumping until
     an emission occurs (at which point it breaks). */


  for (njumps = 0; njumps < MAXJUMPS; njumps++)
  {
    /*  The excited configuration is now known. Now compute all the probabilities of deactivation
       /jumping from this configuration. Then choose one. */

    nbbd = config[uplvl].n_bbd_jump;    //store these for easy access -- number of bb downward jumps
    nbbu = config[uplvl].n_bbu_jump;    // number of bb upward jump from this configuration
    nbfd = config[uplvl].n_bfd_jump;    // number of bf downward jumps from this transition
    nbfu = config[uplvl].n_bfu_jump;    // number of bf upward jumps from this transiion

    if (prbs_known[uplvl] == FALSE)
    {

      m = 0;                    //m counts the total number of possible ways to leave the level
      pjnorm = 0.0;             //stores the total jump probability
      penorm = 0.0;             //stores the total emission probability

      for (n = 0; n < nbbd + nbfd; n++)
      {
        eprbs_known[uplvl][n] = eprbs[n] = 0;   //stores the individual emission probabilities SS
        jprbs_known[uplvl][n] = jprbs[n] = 0;   //stores the individual jump probabilities SS
      }
      for (n = nbbd + nbfd; n < nbbd + nbfd + nbbu + nbfu; n++)
      {
        jprbs_known[uplvl][n] = jprbs[n] = 0;   /*slots for upward jumps */
      }


      /* bb */

      /* First downward jumps. (I.e. those that have emission probabilities. */

      /* For bound-bound decays the jump probability is A-coeff * escape-probability * energy */
      /* At present the escape probability is only approximated (p_escape). This should be improved. */
      /* The collisional contribution to both the jumping and deactivation probabilities are now added (SS, Apr04) */

      for (n = 0; n < nbbd; n++)
      {
        line_ptr = &line[config[uplvl].bbd_jump[n]];

        rad_rate = (a21 (line_ptr) * p_escape (line_ptr, xplasma));
        coll_rate = ne * q21 (line_ptr, t_e);


        bb_cont = rad_rate + coll_rate;
        jprbs_known[uplvl][m] = jprbs[m] = bb_cont * config[line_ptr->nconfigl].ex;     //energy of lower state

        eprbs_known[uplvl][m] = eprbs[m] = bb_cont * (config[uplvl].ex - config[line[config[uplvl].bbd_jump[n]].nconfigl].ex);  //energy difference

        if (jprbs[m] < 0.)      //test (can be deleted eventually SS)
        {
          Error ("Negative probability (matom, 1). Abort.");
          *escape = TRUE;
          p->istat = P_ERROR_MATOM;
          return (-1);
        }
        if (eprbs[m] < 0.)      //test (can be deleted eventually SS)
        {
          Error ("Negative probability (matom, 2). Abort.");
          *escape = TRUE;
          p->istat = P_ERROR_MATOM;
          return (-1);
        }

        pjnorm += jprbs[m];
        penorm += eprbs[m];
        m++;
      }

      /* bf */
      for (n = 0; n < nbfd; n++)
      {

        cont_ptr = &phot_top[config[uplvl].bfd_jump[n]];        //pointer to continuum
        if (n < 25)
        {
          sp_rec_rate = mplasma->recomb_sp[config[uplvl].bfd_indx_first + n];   //need this twice so store it
          bf_cont = (sp_rec_rate + q_recomb (cont_ptr, t_e) * ne) * ne;
        }
        else
        {
          bf_cont = 0.0;
        }

        jprbs_known[uplvl][m] = jprbs[m] = bf_cont * config[phot_top[config[uplvl].bfd_jump[n]].nlev].ex;       //energy of lower state
        eprbs_known[uplvl][m] = eprbs[m] = bf_cont * (config[uplvl].ex - config[phot_top[config[uplvl].bfd_jump[n]].nlev].ex);  //energy difference
//OLD        if (jprbs[m] < 0.)      //test (can be deleted eventually SS)
//OLD        {
//OLD          Error ("Negative probability (matom, 3). Abort.");
//OLD          *escape = 1;
//OLD          p->istat = P_ERROR_MATOM;
//OLD          return (0);
//OLD        }
//OLD        if (eprbs[m] < 0.)      //test (can be deleted eventually SS)
//OLD        {
//OLD          Error ("Negative probability (matom, 4). Abort.");
//OLD        }
        pjnorm += jprbs[m];
        penorm += eprbs[m];
        m++;
      }

      /* Now upwards jumps. */

      /* bb */
      /* For bound-bound excitation the jump probability is B-coeff times Jbar with a correction 
         for stimulated emission. To avoid the need for recalculation all the time, the code will
         be designed to include the stimulated correction in Jbar - i.e. the stimulated correction
         factor will NOT be included here. (SS) */
      /* There is no emission probability for upwards transitions. */
      /* Collisional contribution to jumping probability added. (SS,Apr04) */

      for (n = 0; n < nbbu; n++)
      {
        line_ptr = &line[config[uplvl].bbu_jump[n]];
        rad_rate = (b12 (line_ptr) * mplasma->jbar_old[config[uplvl].bbu_indx_first + n]);

        coll_rate = ne * q12 (line_ptr, t_e);   // this is multiplied by ne below


        jprbs_known[uplvl][m] = jprbs[m] = (rad_rate + coll_rate) * config[uplvl].ex;   //energy of lower state



//OLD        if (jprbs[m] < 0.)      //test (can be deleted eventually SS)
//OLD        {
//OLD          Error ("Negative probability (matom, 5). Abort.");
//OLD          *escape = 1;
//OLD          p->istat = P_ERROR_MATOM;
//OLD          return (0);
//OLD        }
        pjnorm += jprbs[m];
        m++;
      }

      /* bf */
      for (n = 0; n < nbfu; n++)
      {
        /* For bf ionization the jump probability is just gamma * energy
           gamma is the photoionisation rate. Stimulated recombination also included. */
        cont_ptr = &phot_top[config[uplvl].bfu_jump[n]];        //pointer to continuum

        /* first let us take care of the situation where the lower level is zero or close to zero */
        lower_density = den_config (xplasma, cont_ptr->nlev);
        if (lower_density >= DENSITY_PHOT_MIN)
        {
          density_ratio = den_config (xplasma, cont_ptr->uplev) / lower_density;
        }
        else
          density_ratio = 0.0;

        jprbs_known[uplvl][m] = jprbs[m] = (mplasma->gamma_old[config[uplvl].bfu_indx_first + n] - (mplasma->alpha_st_old[config[uplvl].bfu_indx_first + n] * xplasma->ne * density_ratio) + (q_ioniz (cont_ptr, t_e) * ne)) * config[uplvl].ex;        //energy of lower state

        /* this error condition can happen in unconverged hot cells where T_R >> T_E.
           for the moment we set to 0 and hope spontaneous recombiantion takes care of things */
        /* note that we check and report this in check_stimulated_recomb() in estimators.c once a cycle */
        if (jprbs[m] < 0.)      //test (can be deleted eventually SS)
        {
          //Error ("Negative probability (matom, 6). Abort?\n");
          jprbs_known[uplvl][m] = jprbs[m] = 0.0;

        }
        pjnorm += jprbs[m];
        m++;
      }

      pjnorm_known[uplvl] = pjnorm;
      penorm_known[uplvl] = penorm;
      prbs_known[uplvl] = TRUE;

    }

    /* Probabilities of jumping (j) and emission (e) are now known. The normalisation factor pnorm has
       also been recorded. the integer m now gives the total number of possibilities too. 
       now select what happens next. Start by choosing the random threshold value at which the
       event will occur. */



    if ((pjnorm_known[uplvl] + penorm_known[uplvl]) <= 0.0)
    {
      Error ("matom: macro atom level has no way out: uplvl %d pj %g pe %g t_e %.3g  ne %.3g\n", uplvl, pjnorm_known[uplvl],
             penorm_known[uplvl], t_e, ne);
      Error ("matom: macro atom level has no way out: z %d istate %d nion %d ilv %d nbfu %d nbfd %d nbbu %d nbbd %d\n", config[uplvl].z,
             config[uplvl].istate, config[uplvl].nion, config[uplvl].ilv, nbfu, nbfd, nbbu, nbbd);
      *escape = TRUE;
      p->istat = P_ERROR_MATOM;
      return (-1);
    }

    threshold = random_number (0.0, 1.0);
    if (((pjnorm_known[uplvl] / (pjnorm_known[uplvl] + penorm_known[uplvl])) < threshold) || (pjnorm_known[uplvl] == 0))
      break;                    // An emission occurs and so we leave the for loop.

    uplvl_old = uplvl;

// Continue on if a jump has occured 

    /* Use a running total to decide where event occurs. */


    threshold = random_number (0.0, 1.0);
    threshold = threshold * pjnorm_known[uplvl_old];

    n = 0;
    run_tot = 0;
    while (run_tot < threshold)
    {
      run_tot += jprbs_known[uplvl_old][n];
      n++;
    }
    // This added to prevent case where theshold is essentially 0. It is already checked that the jumping probability is not zero
    if (n > 0)
    {
      n = n - 1;
    }

    /* n now identifies the jump that occurs - now set the new level. */
    icheck = 0;
    if (n < nbbd)
    {                           /* bb downwards jump */
      uplvl = line[config[uplvl].bbd_jump[n]].nconfigl;
      icheck = 1;
    }
    else if (n < (nbbd + nbfd))
    {                           /* bf downwards jump */
      uplvl = phot_top[config[uplvl].bfd_jump[n - nbbd]].nlev;
      icheck = 2;
    }
    else if (n < (nbbd + nbfd + nbbu))
    {                           /* bb upwards jump */
      uplvl = line[config[uplvl].bbu_jump[n - nbbd - nbfd]].nconfigu;
      icheck = 3;
    }
    else if (n < (nbbd + nbfd + nbbu + nbfu))
    {                           /* bf upwards jump */
      uplvl = phot_top[config[uplvl].bfu_jump[n - nbbd - nbfd - nbbu]].uplev;
      icheck = 4;
    }
    else
    {
      Error ("Trying to jump but nowhere to go! Matom. Abort");
      *escape = TRUE;
      p->istat = P_ERROR_MATOM;
      return (-1);
    }

  }

  /* Check that the the number of level jumps has not exceeded the maxium value */


  if (njumps == MAXJUMPS)
  {
    Error ("Matom: jumped %d times with no emission for photon %d from upper level %d  pjnorm %e penorm %e Abort.\n", MAXJUMPS, p->np,
           uplvl, pjnorm_known[uplvl], penorm_known[uplvl]);
    *escape = TRUE;
    p->istat = P_ERROR_MATOM;
    return (-1);
  }


  /* At this point we know the macro actom has deactivated and we know the level
   * at which the deactivation occurs.  We still need to decide the exact procees
   * by which the macro actom deactivates.  
   */

  run_tot = 0;
  n = 0;

  threshold = random_number (0.0, 1.0);

  threshold = threshold * penorm_known[uplvl];  //normalise to total emission prob.

  while (run_tot < threshold)
  {
    run_tot += eprbs_known[uplvl][n];
    n++;
  }
  n = n - 1;
  /* n now identifies the jump that occurs - now set nres for the return value. */
  if (n < nbbd)
  {                             /* bb downwards jump */
    /* With collisions included we need to decide whether the deactivation is radiative (r-packet)
       or collisional (k-packet). Get a random number and then use the ratio of the collisional
       emission probability to the (already known) collisional+radiative probability to decide whether
       collisional or radiative deactivation occurs. */
    choice = random_number (0.0, 1.0);


    line_ptr = &line[config[uplvl].bbd_jump[n]];        //pointer for the bb transition

    rad_rate = a21 (line_ptr) * p_escape (line_ptr, xplasma);

    coll_rate = ne * q21 (line_ptr, t_e);

    if (choice > (coll_rate / (rad_rate + coll_rate)))
    {
      /* It's a r-packet (radiative) deactivation */
      *escape = TRUE;
      *nres = line[config[uplvl].bbd_jump[n]].where_in_list;
      p->freq = line[config[uplvl].bbd_jump[n]].freq;
    }
    else
    {
      /* It's a k-packet (collisional deactivation. */
      *escape = FALSE;
    }
  }
  else if (n < (nbbd + nbfd))
  {                             /* bf downwards jump */
    /* With collisional recombination included we need to decide whether to deactivate
       radiatively or make a k-packet. */

    cont_ptr = &phot_top[config[uplvl].bfd_jump[n - nbbd]];

    rad_rate = mplasma->recomb_sp[config[uplvl].bfd_indx_first + n - nbbd];     //again using recomb_sp rather than alpha_sp (SS July 04)
    coll_rate = ne * q_recomb (cont_ptr, t_e);

    choice = random_number (0.0, 1.0);

    if (choice > (coll_rate / (rad_rate + coll_rate)))
    {
      /*It's a r-packet (radiative) deactivation */
      *escape = TRUE;
      *nres = config[uplvl].bfd_jump[n - nbbd] + NLINES + 1;
      p->freq = matom_select_bf_freq (one, config[uplvl].bfd_jump[n - nbbd]);

    }
    else
    {
      /* It's a k-packet (collisional) deactivation */
      *escape = FALSE;
    }

  }
  else
  {
    Error ("Trying to emitt from Macro Atom but no available route (matom). Abort.");
    *escape = TRUE;
    p->istat = P_ERROR_MATOM;
    return (-1);
  }

//  Log ("Leaving phot %d in cell %d after % d jumps via %d\n", p->np, p->grid, njumps, *escape);

  return (njumps);
}





/************************************
**  
* @brief the b12 Einstein coefficient.
*
* @param struct lines line_ptr line pointer to calculate
* 
* ###Notes### 
* ksl OK, Probably should be moved to lines.c for consistency, but that can be done
* later.
***********************************/

#define B12_CONSTANT 5.01983e25

struct lines *b12_line_ptr;
double b12_a;

double
b12 (line_ptr)
     struct lines *line_ptr;
{
  double freq;

  if (b12_line_ptr != line_ptr)
  {
    freq = line_ptr->freq;
    b12_a = B12_CONSTANT * line_ptr->f / freq;
    b12_line_ptr = line_ptr;
  }

  return (b12_a);
}

/************************************************************/
/* As for similar routines in recomb.c, in order to use the integrator the 
   following external structures are used (SS)*/
/* This relates to the alpha_sp routines at the end of this file */

struct topbase_phot *cont_ext_ptr;      //continuum pointer passed externally
double temp_ext;                //temperature passed externally
int temp_choice;                //choice of type of calcualation for alpha_sp

/*****************************************************************************/




/**********************************************************/
/**  
 * @brief the matom estimator for the spontaneous recombination rate.
 * 
 * The rate is given by 
 * 
 *    (4 pi /c2) (gu/gl) (h2/2 pi m k T)^(3/2) 
 * times the integral of   a(nu) nu2 exp [(chi- h nu)/kT].
 * 
 * ###Notes###
 * 04jul30	ksl	Modified so that one does not need to have multiple versions
 * 		of the code depending on how whether the integrand is 
 * 		multiplied by 1, f/fthresh, or f/fthresh-1.  This was
 * 		done eliminate alpha_sp_e as a separate set of routines
 * 		and to assure that bf rates are positive 
 * 			ichoice = 0   --> spontanous recombination
 * 			ichoice = 1   --> energy weighted recombination 
 * 			ichoice = 2   --> the difference between energy_weighted
 * 					and spontaneous
 * 
 * 	06may	ksl	57+ -- Modified to use plasma structure
 *  Energy weighted means that the integrand has an extra factor nu/nu_threshold
 *  The difference case is (nu-nu_threshold)/nu_threhold 
***********************************************************/
#define ALPHA_SP_CONSTANT 5.79618e-36

double
alpha_sp (cont_ptr, xplasma, ichoice)
     struct topbase_phot *cont_ptr;
     PlasmaPtr xplasma;
     int ichoice;
{
  double alpha_sp_value;
  double fthresh, flast;

  temp_choice = ichoice;
  temp_ext = xplasma->t_e;      //external for use in alph_sp_integrand
  cont_ext_ptr = cont_ptr;      //"
  fthresh = cont_ptr->freq[0];  //first frequency in list
  flast = cont_ptr->freq[cont_ptr->np - 1];     //last frequency in list
  if ((H_OVER_K * (flast - fthresh) / temp_ext) > ALPHA_MATOM_NUMAX_LIMIT)
  {
    //flast is currently very far into the exponential tail: so reduce flast to limit value of h nu / k T.
    flast = fthresh + temp_ext * ALPHA_MATOM_NUMAX_LIMIT / H_OVER_K;
  }
  // alpha_sp_value = qromb (alpha_sp_integrand, fthresh, flast, 1e-4);
  alpha_sp_value = num_int (alpha_sp_integrand, fthresh, flast, 1e-4);

  /* The lines above evaluate the integral in alpha_sp. Now we just want to multiply 
     through by the appropriate constant. */
  if (cont_ptr->macro_info == TRUE && geo.macro_simple == FALSE)
  {
    alpha_sp_value = alpha_sp_value * config[cont_ptr->nlev].g / config[cont_ptr->uplev].g * pow (xplasma->t_e, -1.5);
  }
  else                          //case for simple element
  {
    alpha_sp_value = alpha_sp_value * config[cont_ptr->nlev].g / ion[cont_ptr->nion + 1].g * pow (xplasma->t_e, -1.5);  //g for next ion up used
  }

  alpha_sp_value = alpha_sp_value * ALPHA_SP_CONSTANT;

  return (alpha_sp_value);
}

/**********************************************************/
/**  
 *  @brief the matom estimator for the spontaneous recombination rate.
 * 
 * The rate is given by 
 * 
 *    (4 pi /c2) (gu/gl) (h2/2 pi m k T)^(3/2) 
 * times the integral of   a(nu) nu2 exp [(chi- h nu)/kT].
 * 
 * ###Notes###
 * 04jul30	ksl	Modified so that one does not need to have multiple versions
 * 		of the code depending on how whether the integrand is 
 * 		multiplied by 1, f/fthresh, or f/fthresh-1.  This was
 * 		done eliminate alpha_sp_e as a separate set of routines
 * 		and to assure that bf rates are positive 
 * 			ichoice = 0   --> spontanous recombination
 * 			ichoice = 1   --> energy weighted recombination 
 * 			ichoice = 2   --> the difference between energy_weighted
 * 					and spontaneous
 * 
 * 	06may	ksl	57+ -- Modified to use plasma structure
***********************************************************/
#define ALPHA_SP_CONSTANT 5.79618e-36

double
scaled_alpha_sp_integral_band_limited (cont_ptr, xplasma, ichoice, fmin, fmax)
     struct topbase_phot *cont_ptr;
     PlasmaPtr xplasma;
     int ichoice;
     double fmin, fmax;
{
  double alpha_sp_value;
  double fthresh, flast;

  temp_choice = ichoice;
  temp_ext = xplasma->t_e;      //external for use in alph_sp_integrand
  cont_ext_ptr = cont_ptr;      //"
  fthresh = cont_ptr->freq[0];  //first frequency in list
  flast = cont_ptr->freq[cont_ptr->np - 1];     //last frequency in list
  if ((H_OVER_K * (flast - fthresh) / temp_ext) > ALPHA_MATOM_NUMAX_LIMIT)
  {
    //flast is currently very far into the exponential tail: so reduce flast to limit value of h nu / k T.
    flast = fthresh + temp_ext * ALPHA_MATOM_NUMAX_LIMIT / H_OVER_K;
  }
  if (flast < fmax)
  {
    fmax = flast;
  }
  // alpha_sp_value = qromb (alpha_sp_integrand, fmin, fmax, 1e-4);
  alpha_sp_value = num_int (alpha_sp_integrand, fmin, fmax, 1e-4);

  return (alpha_sp_value);
}



/**********************************************************/
/** 
 *  @brief This returns the integrand for alpha_sp at a chosen frequency - 
 *
 * Note:
 *  The peculiar variable *params is needed by the gsl routine
 *  that does the integration.
***********************************************************/

double
alpha_sp_integrand (double freq, void *params)
{
  double fthresh;
  double x;
  double integrand;
  double tt;

  fthresh = cont_ext_ptr->freq[0];
  tt = temp_ext;

  if (freq < fthresh)
    return (0.0);               // No recombination at frequencies lower than the threshold freq occur

  x = sigma_phot (cont_ext_ptr, freq);  //this is the cross-section
  integrand = x * freq * freq * exp (H_OVER_K * (fthresh - freq) / tt);


  if (temp_choice == 1)
    return (integrand * freq / fthresh);        //energy weighted case
  if (temp_choice == 2)
    return (integrand * (freq - fthresh) / fthresh);    // difference case
  return (integrand);           //spontaneous case
}




/**********************************************************/
/** 
 * @brief deals with the elimination of k-packets.
 *
 * Whenever a k-packet is created it is
 * immediately destroyed insitu by this routine. At output "nres" identified the process
 * that destroys the k-packet and the packet information has been updated in the same
 * way as at the end of matom
 *
 * @param [in,out]  PhotPtr p   the packet at the point of activation and deactivation
 * @param [in,out]  int nres    the process which activates and deactivates the Macro Atom
 * @param [out]  int escape  to tell us whether the matom de-activation
 *                             is via an r-packet (1 TRUE) or a k-packet (0 FALSE)
 * @param [in] int mode         switch which allows photon to be deactivated by a non-radiative
 * term.  (non_zero is true)
 * @return 0
 *
 * ###Notes###
************************************************************/

int
kpkt (p, nres, escape, mode)
     PhotPtr p;
     int *nres;
     int *escape;
     int mode;
{

  int i;
  double cooling_adiabatic;
  double cooling_normalisation;
  double destruction_choice;
  double electron_temperature;
  double upweight_factor;
  WindPtr one;
  PlasmaPtr xplasma;
  MacroPtr mplasma;
  double freqmin, freqmax;


  /* Idea is to calculated the cooling
     terms for all the processes and then choose one at random to destroy the k-packet and
     turn it back into a photon bundle. 
     The routine considers bound-free, collision excitation and ff
     emission. */

  /* The loop here runs over all the bf processes, regardless of whether they are macro or simple
     ions. The reason that we can do this is that in the macro atom method stimulated recombination
     has already been considered before creating the k-packet (by the use of gamma-twiddle) in "scatter"
     and so we only have spontaneous recombination to worry about here for ALL cases. */


  one = &wmain[p->grid];
  xplasma = &plasmamain[one->nplasma];
  check_plasma (xplasma, "kpkt");
  mplasma = &macromain[xplasma->nplasma];

  electron_temperature = xplasma->t_e;

  /* Set maximum and minimum frequency limits.  See #187. We need band limits for free free packet
     generation (see call to one_ff below). A bandaid is applied if there is not enough
     separation between freqmin and max */

  freqmin = xband.f1[0];
  freqmax = ALPHA_FF * xplasma->t_e / H_OVER_K;
  if (freqmax < 1.1 * freqmin)
  {
    freqmax = 1.1 * freqmin;
  }

  /* If the kpkt destruction rates for this cell are not known they are calculated here.  This happens
   * every time the wind is updated. Note that most of the cooling rates do not include a factor
   * of ne, since one power of ne applies to all of the cooling rates and we are only concerned
   * with relative cooling rates */

  if (mplasma->kpkt_rates_known != TRUE)
  {
    fill_kpkt_rates (xplasma, escape, p);
  }

/* This is the end of the cooling rate calculation, which is done only once for each cell
   and once for each cycle.  

   The next little section deals whith handling adiabatic cooling and shock heating.
   */

  cooling_normalisation = mplasma->cooling_normalisation - mplasma->cooling_adiabatic;
  cooling_adiabatic = 0.0;

  if (mode == KPKT_MODE_ALL)
  {
    if (KPKT_NET_HEAT_MODE && geo.nonthermal)
    {
      if (xplasma->cool_adiabatic > xplasma->heat_shock)
      {
        cooling_adiabatic = (xplasma->cool_adiabatic - xplasma->heat_shock) / xplasma->vol / xplasma->ne;
      }
    }
    else
    {
      cooling_adiabatic = mplasma->cooling_adiabatic;
    }
  }
  cooling_normalisation += cooling_adiabatic;


  /* The cooling rates for the recombination and collisional processes are now known. 

     Choose which process destroys the k-packet with a random number. */

  destruction_choice = random_number (0.0, 1.0) * cooling_normalisation;

  /* This logic of what follows may not be obvious.  For choosing the basic
   * process, we just look to see if the destruction choice is less than bf, bf+bb, bf+bb+ff
   * etc, but inside the bhe "basic_choices", we iteratively reduce! the destruction
   * choice until we get to one that is less than the cooling associated with a specific
   * tranistion. If we do not find such a transition, within for example the bf if statement
   * we drop all the way down to the Error at the end.
   */


  if (destruction_choice < mplasma->cooling_bftot)
  {                             //destruction by bf

    /* JM 1503 -- we used to loop over ntop_phot here, 
       but we should really loop over the tabulated Verner Xsections too
       see #86, #141 */
    for (i = 0; i < nphot_total; i++)
    {
      if (destruction_choice < mplasma->cooling_bf[i])
      {

        if (i > nphot_total - 1)
        {
          Error ("kpkt (matom.c): trying to destroy k-packet in unknown process. Abort.\n");
          *escape = TRUE;
          p->istat = P_ERROR_MATOM;
          return (0);
        }


        *nres = i + NLINES + 1;
        *escape = TRUE;

        p->freq = matom_select_bf_freq (one, i);

        /* if the cross-section corresponds to a simple ion (macro_info == FALSE)
           or if we are treating all ions as simple, then adopt the total emissivity
           approach to choosing photon weights - this means we 
           multipy down the photon weight by a factor nu/(nu-nu_0)
           and we force a kpkt to be created */
#if BF_SIMPLE_EMISSIVITY_APPROACH
        if (phot_top[i].macro_info == FALSE || geo.macro_simple == TRUE)
        {
          upweight_factor = xplasma->recomb_simple_upweight[i];
          p->w *= upweight_factor;

          /* record the amount of energy being extracted from the simple ion ionization pool */
          xplasma->bf_simple_ionpool_out += p->w - (p->w / upweight_factor);
        }
#endif

        return (0);
      }
      else
      {
        destruction_choice = destruction_choice - mplasma->cooling_bf[i];
      }
    }
  }
  else if (destruction_choice < (mplasma->cooling_bftot + mplasma->cooling_bbtot))
  {
    /* a collisional destruction has occurred and so, if the line is  associated with
       a macro atom, it  must be excited.
     */

    destruction_choice = destruction_choice - mplasma->cooling_bftot;
    for (i = 0; i < nlines; i++)
    {
      if (destruction_choice < mplasma->cooling_bb[i])
      {
        *nres = line[i].where_in_list;
        if (line[i].macro_info == TRUE && geo.macro_simple == FALSE)
        {
          *escape = FALSE;
        }
        else
        {
          *escape = TRUE;
          p->freq = line[i].freq;
        }
        return (0);
      }
      else
      {
        destruction_choice = destruction_choice - mplasma->cooling_bb[i];
      }
    }
  }

  else if (destruction_choice < (mplasma->cooling_bftot + mplasma->cooling_bbtot + mplasma->cooling_ff))
  {
    /* If reached this point, it is a FF destruction event */
    /* consult issues #187, #492 regarding free-free */
    *escape = TRUE;
    *nres = -2;
    p->freq = one_ff (one, freqmin, freqmax);
    return (0);
  }
  else if (destruction_choice < (mplasma->cooling_bftot + mplasma->cooling_bbtot + mplasma->cooling_ff + mplasma->cooling_ff_lofreq))
  {
    /*this is ff at a frequency that is so low frequency that it is not worth tracking further */
    *escape = TRUE;
    *nres = -2;
    p->istat = P_LOFREQ_FF;
    return (0);
  }


  else if (destruction_choice <
           (mplasma->cooling_bftot + mplasma->cooling_bbtot + mplasma->cooling_ff + mplasma->cooling_ff_lofreq + cooling_adiabatic))
  {
    /* It is a k-packat that is destroyed by adiabatic cooling */

    if (geo.adiabatic == 0 || mode != KPKT_MODE_ALL)
    {
      Error ("kpkt: Destroying kpkt by adiabatic cooling even though it is turned off.\n");
    }
    *escape = TRUE;
    *nres = -2;
    p->istat = P_ADIABATIC;

    return (0);
  }


  else
  {
    /* It is a k-packed destroyed by collisional ionization in a macro atom. */
    destruction_choice =
      destruction_choice - mplasma->cooling_bftot - mplasma->cooling_bbtot - mplasma->cooling_ff - mplasma->cooling_ff_lofreq -
      cooling_adiabatic;

    for (i = 0; i < nphot_total; i++)
    {
      if (destruction_choice < mplasma->cooling_bf_col[i])
      {
        if (i > nphot_total - 1)
        {
          Error ("kpkt: trying to destroy k-packet in unknown process. Abort.\n");
          *escape = TRUE;
          p->istat = P_ERROR_MATOM;
          return (0);
        }

        *nres = i + NLINES + 1;
        *escape = FALSE;

        return (0);
      }
      else
      {
        destruction_choice = destruction_choice - mplasma->cooling_bf_col[i];
      }
    }
  }



  Error ("kpkt: Failed to select a destruction process in kpkt. Abort.\n");
  Error
    ("kpkt: choice %8.4e norm %8.4e cooling_bftot %g, cooling_bbtot %g, cooling_ff %g, cooling_ff_lofreq %g, cooling_bf_coltot %g cooling_adiabatic %g cooling_adiabatic %g\n",
     destruction_choice, cooling_normalisation, mplasma->cooling_bftot, mplasma->cooling_bbtot, mplasma->cooling_ff,
     mplasma->cooling_ff_lofreq, mplasma->cooling_bf_coltot, mplasma->cooling_adiabatic, cooling_adiabatic);

  *escape = TRUE;
  p->istat = P_ERROR_MATOM;

  return (0);

}



/************************************************************
 ** 
 * @brief routine for dealing with bound-bound "simple ions" within the hybrid macro-atom framework
 *
 * @param [in,out] PhotPtr p   the packet at the point of activation
 * @param [in,out] int nres    the process which activates the Macro Atom
 * @param [in,out] int escape  in principle this tells us whether de-activation is
 *                             via an r-packet or a k-packet. For this routine at the 
 *                             moment only r-packets are possible so it always returns
 *                             escape = TRUE
 *
 * @detail fake_matom_bb is the macro atom routine that deals with line events involving
 * simple ions (i.e. ions for which a full macro atom treatment is not employed.
 * When this routine is called a simple line has absorbed a packet. This routine 
 * creates a fake two-level macro atom and determines whether the packet energy
 * is simply re-emitted in the line or is thermalised. If it is thermalised it
 * turns into a k-packet and the appropriate routine is called. 
 *
************************************************************/

int
fake_matom_bb (p, nres, escape)
     PhotPtr p;
     int *nres;
     int *escape;
{
  double kprb, rprb;
  WindPtr one;
  PlasmaPtr xplasma;
  struct lines *line_ptr;
  double electron_temperature;
  double normalisation;
  double choice;

  one = &wmain[p->grid];
  xplasma = &plasmamain[one->nplasma];

  line_ptr = lin_ptr[*nres];    //record the line pointer
  electron_temperature = xplasma->t_e;

  /* Upon calling we know that the upper level of our fake two level macro
     atom is excited. Since it's only two-levels there are no jumping probabilities
     and the only decision that needs to be made is whether it de-excites by an r-packet
     in the line or a collisionally produced k-packet. To decide this we compute the
     two emission probabilities. rprb is the probability of the r-packet and kprb
     for the k-packet. */

  /* Since both k- and r-packet emission have the same delta energy associated with them then
     we can ignore the energy factor that is otherwise needed in the de-activation probabilities.
     The upper level population is also factored out. (SS) */

  rprb = a21 (line_ptr) * p_escape (line_ptr, xplasma);
  // Einstein-A * escape probability

  kprb = q21 (line_ptr, electron_temperature) * xplasma->ne * (1. - exp (-H_OVER_K * line_ptr->freq / electron_temperature));
  // Collision de-excitation rate coeff. * electron density 
  /* The extra factor of (1 - exp (-h nu / k T)) is explained in KSL's notes on Python. It appears because we are
     using the "scattering fraction" formalism for simple ions. (SS May 04). */

  normalisation = kprb + rprb;

  kprb = kprb / normalisation;

  rprb = rprb / normalisation;

  /* Now just use a random number to decide what happens. */

  choice = random_number (0.0, 1.0);


  /* If "choice" is less than rprb then we have chosen a radiative decay - for this fake macro atom there 
     is only one line so there's nothing to do - the energy is re-radiated in the line and that's it. We
     don't need to change nres. Otherwise we've chosen a collisional destruction and so we need to get
     kpkt to deal with the resulting k-packet and give us a new value of nres following thermalisation and
     re-emission of the energy. */

  if (choice < rprb)
  {
    *escape = TRUE;             //it's an r-packet de-activation
    p->freq = line_ptr->freq;
  }
  else
  {
    /* Return the instruction to macro_gov to
       make a k-packet. */

    *escape = FALSE;

  }


  return (0);

}


/************************************************************
 ** 
 * @brief routine for dealing with bound-free "simple ions" within the hybrid macro-atom framework
 *
 * @param [in,out] PhotPtr p   the packet at the point of activation
 * @param [in,out] int nres    the process which activates the Macro Atom
 * @param [in,out] int escape  in principle this tells us whether de-activation is
 *                             via an r-packet or a k-packet. For this routine at the 
 *                             moment only r-packets are possible so it always returns
 *                             escape = TRUE
 *
 * @details fake_matom_bf is the macro atom routine that deals with photoionisation 
 * events involving simple ions (i.e. ions for which a full macro atom treatment 
 * is not employed). When this routine is called a photoionisation has absorbed a packet. 
 * The idea of this routine is to deal with the subsequenct
 * by creating a fake two-level atom. However, in the absense of collisional recombination 
 * (or something similar) there's only one thing that can happen - radiative 
 * recombination. Therefore there's no need to do anything here unless
 * collisional recombination is to be introduced at some point in the
 * future. All this routine does for now is choose a new frequency for the 
 * emitted r-packet.
 * 
************************************************************/

int
fake_matom_bf (p, nres, escape)
     PhotPtr p;
     int *nres;
     int *escape;
{
  WindPtr one;
  PlasmaPtr xplasma;

  one = &wmain[p->grid];
  xplasma = &plasmamain[one->nplasma];

  *escape = TRUE;

  p->freq = matom_select_bf_freq (one, *nres - NLINES - 1);


  return (0);

}




/**********************************************************/
/** 
 * @brief a scaled down version of matom which deals with deactivation only for spectral cycles.
 *
 * @param [in]     WindPtr w   the ptr to the structure defining the wind
 * @param [in]     int upper   the upper level that we deactivate from
 * @param [in,out]  PhotPtr p   the packet at the point of activation and deactivation
 * @param [in,out]  int nres    the process by which deactivation occurs
 * @return 0
 *
 * emit_matom is a scaled down version of matom which deals with the emission due
 * to deactivating macro atoms in the detailed spectrum part of the calculation.
 * 
***********************************************************/

int
emit_matom (w, p, nres, upper, fmin, fmax)
     WindPtr w;
     PhotPtr p;
     int *nres;
     int upper;
     double fmin, fmax;
{
  struct lines *line_ptr;
  struct topbase_phot *cont_ptr;
  int uplvl;
  double eprbs[NBBJUMPS + NBFJUMPS];
  double penorm;
  double threshold, run_tot;
  double sp_rec_rate;
  int n, m;
  int nbbd, nbfd;
  double t_e, ne;
  double bb_cont;
  WindPtr one;
  PlasmaPtr xplasma;


  one = &w[p->grid];
  xplasma = &plasmamain[one->nplasma];

  t_e = xplasma->t_e;
  ne = xplasma->ne;


  /* The first step is to identify the configuration that has been excited. */

  uplvl = upper;

  /* Unlike the proper matom routine, we don't want any jumps, only deactivations here. */


  /*  The excited configuration is now known. Now compute all the probabilities of deactivation
     /jumping from this configuration. Then choose one. */


  nbbd = config[uplvl].n_bbd_jump;      //store these for easy access -- number of bb downward jumps
  nbfd = config[uplvl].n_bfd_jump;      // number of bf downared jumps from this transition


  // Start by setting everything to 0

  m = 0;                        //m counts the total number of possible ways to leave the level
  penorm = 0.0;                 //stores the total emission probability

  for (n = 0; n < nbbd + nbfd; n++)
  {
    eprbs[n] = 0;               //stores the individual emission probabilities SS
  }


  /* bb */

  /* First downward jumps. */

  for (n = 0; n < nbbd; n++)
  {
    line_ptr = &line[config[uplvl].bbd_jump[n]];
    /* Since we are only interested in making an r-packet here we can (a) ignore collisional
       deactivation and (b) ignore lines outside the frequency range of interest. */
    if ((line_ptr->freq > fmin) && (line_ptr->freq < fmax))     // correct range
    {
      bb_cont = (a21 (line_ptr) * p_escape (line_ptr, xplasma));

      eprbs[m] = bb_cont * (config[uplvl].ex - config[line[config[uplvl].bbd_jump[n]].nconfigl].ex);    //energy difference

      if (eprbs[m] < 0.)        //test (can be deleted eventually SS)
      {
        Error ("Negative probability (matom, 2). Abort.");
        Exit (0);
      }

      penorm += eprbs[m];
    }
    m++;
  }

  /* bf */
  for (n = 0; n < nbfd; n++)
  {
    cont_ptr = &phot_top[config[uplvl].bfd_jump[n]];    //pointer to continuum

    /* If the edge is above the frequency range we are interested in then we need not consider this
       bf process. */
    if (cont_ptr->freq[0] < fmax)       //means that it may contribute
    {
      sp_rec_rate = alpha_sp (cont_ptr, xplasma, 0);
      eprbs[m] = sp_rec_rate * ne * (config[uplvl].ex - config[phot_top[config[uplvl].bfd_jump[n]].nlev].ex);   //energy difference
      if (eprbs[m] < 0.)        //test (can be deleted eventually SS)
      {
        Error ("Negative probability (matom, 4). Abort.");
        Exit (0);
      }
      penorm += eprbs[m];
    }
    m++;
  }



  /* No need for upwards jumps - they don't provide emission mechanisms. */

  /* Probabilities of emission (e) are now known. The normalisation factor pnorm has
     also been recorded. the integer m now gives the total number of possibilities too. 
     now select what happens next. Start by choosing the random threshold value at which the
     event will occur. */

  threshold = random_number (0.0, 1.0);


  run_tot = 0;
  n = 0;
  threshold = threshold * penorm;       //normalise to total emission prob.
  while (run_tot < threshold)
  {
    run_tot += eprbs[n];
    n++;
  }
  n = n - 1;
  /* n now identifies the jump that occurs - now set nres for the return value. */
  if (n < nbbd)
  {                             /* bb downwards */
    line_ptr = &line[config[uplvl].bbd_jump[n]];
    *nres = line[config[uplvl].bbd_jump[n]].where_in_list;
    p->freq = line[config[uplvl].bbd_jump[n]].freq;
  }
  else if (n < (nbbd + nbfd))
  {                             /* bf downwards jump */
    *nres = config[uplvl].bfd_jump[n - nbbd] + NLINES + 1;
    /* continuua are indicated by nres > NLINES */

    p->freq = matom_select_bf_freq (one, config[uplvl].bfd_jump[n - nbbd]);


    /*Currently this assumed hydrogenic shape cross-section - Improve */
  }
  else
  {
    Error ("Trying to emit from Macro Atom but no available route (emit_matom). Abort.");
    Exit (0);
  }
  /* The frequency and the value of nres have been set correctly. All done. */
  return (0);
}

/**********************************************************/
/** 
 * @brief Prob. of cell emitting in a given line
 *
 * @param [in] one            Pointer to cell of interest
 * @param [in] line_ptr_emit  Pointer to line of interest
 * @return     Probability of line emission
 *
 * Given a cell and a line, calculates the probabiltiy that
 * that cell will emit in that line.
 *
 * ###Notes###
 * 6/15 - Written by SWM
***********************************************************/
double
matom_emit_in_line_prob (WindPtr one, struct lines *line_ptr_emit)
{
  double eprbs, eprbs_line;
  struct lines *line_ptr;
  struct topbase_phot *cont_ptr;
  int uplvl;
  double alpha_sp ();
  double penorm, freqmax, freqmin;
  double sp_rec_rate;
  int n;
  int nbbd, nbfd;
  double ne;
  double bb_cont;
  PlasmaPtr xplasma;

  line_ptr = NULL;              // Avoid compiler warning
  xplasma = &plasmamain[one->nplasma];
  ne = xplasma->ne;


  //We know the line we're de-exciting into, so find the upper level of this line
  uplvl = line_ptr_emit->nconfigu;


  //Now find the number of jumps available from this level
  nbbd = config[uplvl].n_bbd_jump;      //store these for easy access -- number of bb downward jumps
  nbfd = config[uplvl].n_bfd_jump;      // number of bf downared jumps from this transition

  penorm = 0.0;                 //stores the total emission probability
  eprbs_line = 0.0;

  // Set frequency range to search to be the spectral range
  freqmin = VLIGHT / (geo.swavemax * 1e-8);
  freqmax = VLIGHT / (geo.swavemin * 1e-8);

  /* bb */
  /* First downward jumps. */
  for (n = 0; n < nbbd; n++)
  {
    line_ptr = &line[config[uplvl].bbd_jump[n]];
    /* Since we are only interested in making an r-packet here we can (a) ignore collisional
       deactivation and (b) ignore lines outside the frequency range of interest. */
    if ((line_ptr->freq > freqmin) && (line_ptr->freq < freqmax))
    {
      bb_cont = (a21 (line_ptr) * p_escape (line_ptr, xplasma));
      eprbs = bb_cont * (config[uplvl].ex - config[line[config[uplvl].bbd_jump[n]].nconfigl].ex);       //energy difference
      penorm += eprbs;
      if (line_ptr == line_ptr_emit)
      {
        eprbs_line = eprbs;
      }
    }
  }

  if (eprbs_line == 0.0)
  {
    Error ("matom_emit_in_line_prob: Line frequency %g lies outside spectral range %g-%g!\n", line_ptr->freq, geo.sfmin, geo.sfmax);
    return (-1.0);
  }

  /* bf */
  /* There should be no bf jumps for Ha/Hb but included for potential use for other lines */
  for (n = 0; n < nbfd; n++)
  {
    cont_ptr = &phot_top[config[uplvl].bfd_jump[n]];

    /* If the edge is above the frequency range we are interested in then we need not consider this
       bf process. */
    if (cont_ptr->freq[0] < freqmax)    //means that it may contribute
    {
      sp_rec_rate = alpha_sp (cont_ptr, xplasma, 0);
      eprbs = sp_rec_rate * ne * (config[uplvl].ex - config[phot_top[config[uplvl].bfd_jump[n]].nlev].ex);      //energy difference
      penorm += eprbs;
    }
  }
  //We now have the total probability of all avenues the matom could de-excite into

  //Return probability that the matom in this cell de-excites into this line
  return (eprbs_line / penorm);
}
