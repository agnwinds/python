#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_blas.h>
#include "my_linalg.h"


/**********************************************************/
/** 
 * @brief calculate the matrix of probabilities for the accelerated macro-atom scheme 
 *
 * @param [in] PlasmaPtr  xplasma  
 * @param [in,out] double **matom_matrix 
 *        the 2D matrix array we will populate with normalised probabilities 
 *
 * @details
 * given an activation state, this routine calculates the probability that a packet
 * will deactivate from a given state. Let's suppose that the \f$(i,j)\f$-th element of matrix 
 * \f$Q\f$ contains the **jumping** probability from state i to state j, and the diagonals 
 * \f$(i,i)\f$ of matrix \f$R\f$ contain the emission probabilities from each state, then it 
 * can be shown that (see short notes from Vogl, or 
 * <a href="Ergon et al. 2018">https://ui.adsabs.harvard.edu/abs/2018A%26A...620A.156E</a>)
 * the quantity we want is then \f$B = N R\f$, where \f$N = (I - Q)^{-1}\f$, where \f$I\f$
 * is the identity matrix. This routine does this calculation.
 *
 **********************************************************/

void
calc_matom_matrix (xplasma, matom_matrix)
     PlasmaPtr xplasma;
     double **matom_matrix;
{
  MacroPtr mplasma;
  double t_e, ne;
  int nbbd, nbbu, nbfd, nbfu;
  int uplvl, target_level;
  double Qcont;
  struct lines *line_ptr;
  struct topbase_phot *cont_ptr;
  double rad_rate, coll_rate;
  int n, i, nn, mm;
  double Qcont_kpkt, bb_cont, sp_rec_rate, bf_cont, lower_density, density_ratio;
  double kpacket_to_rpacket_rate, norm, Rcont;
  double *a_data;
  mplasma = &macromain[xplasma->nplasma];       //telling us where in the matom structure we are

  t_e = xplasma->t_e;           //electron temperature 
  ne = xplasma->ne;             //electron number density

  /* allocate arrays for matrices and normalsiations
     we want to include all macro-atom levels + 1 kpacket levels */
  int nrows = nlevels_macro + 1;
  double **R_matrix = (double **) calloc (sizeof (double *), nrows);
  double **Q_matrix = (double **) calloc (sizeof (double *), nrows);
  double *Q_norm = (double *) calloc (sizeof (double), nrows);

  for (i = 0; i < nrows; i++)
  {
    R_matrix[i] = (double *) calloc (sizeof (double), nrows);
    Q_matrix[i] = (double *) calloc (sizeof (double), nrows);
  }

  /* initialise everything to zero */
  norm = 0.0;
  for (uplvl = 0; uplvl < nrows; uplvl++)
  {
    Q_norm[uplvl] = 0.0;
    for (target_level = 0; target_level < nrows; target_level++)
    {
      Q_matrix[uplvl][target_level] = 0.0;
      R_matrix[uplvl][target_level] = 0.0;
      matom_matrix[uplvl][target_level] = 0.0;
    }
  }

  /* loop over all macro-atom levels and populate the rate matrix */
  for (uplvl = 0; uplvl < nlevels_macro; uplvl++)
  {
    nbbd = config[uplvl].n_bbd_jump;    //store these for easy access -- number of bb downward jumps
    nbbu = config[uplvl].n_bbu_jump;    // number of bb upward jump from this configuration
    nbfd = config[uplvl].n_bfd_jump;    // number of bf downward jumps from this transition
    nbfu = config[uplvl].n_bfu_jump;    // number of bf upward jumps from this transiion


    /* bound-bound */
    for (n = 0; n < nbbd; n++)
    {

      line_ptr = &line[config[uplvl].bbd_jump[n]];

      rad_rate = (a21 (line_ptr) * p_escape (line_ptr, xplasma));
      coll_rate = q21 (line_ptr, t_e);  // this is multiplied by ne below

      if (coll_rate < 0)
      {
        coll_rate = 0;
      }

      bb_cont = rad_rate + (coll_rate * ne);

      target_level = line_ptr->nconfigl;

      //internal jump to another macro atom level
      Q_matrix[uplvl][target_level] += Qcont = bb_cont * config[target_level].ex;       //energy of lower state

      //jump to the k-packet pool (we used to call this "deactivation")
      Q_matrix[uplvl][nlevels_macro] += Qcont_kpkt = (coll_rate * ne) * (config[uplvl].ex - config[target_level].ex);   //energy of lower state

      //deactivation back to r-packet
      R_matrix[uplvl][uplvl] += Rcont = rad_rate * (config[uplvl].ex - config[target_level].ex);        //energy difference

      Q_norm[uplvl] += Qcont + Qcont_kpkt + Rcont;
    }

    /* bound-free */
    for (n = 0; n < nbfd; n++)
    {

      cont_ptr = &phot_top[config[uplvl].bfd_jump[n]];  //pointer to continuum
      if (n < 25)               //?
      {
        sp_rec_rate = mplasma->recomb_sp[config[uplvl].bfd_indx_first + n];     //need this twice so store it
        bf_cont = (sp_rec_rate + q_recomb (cont_ptr, t_e) * ne) * ne;
      }
      else
      {
        bf_cont = sp_rec_rate = 0.0;
      }


      target_level = phot_top[config[uplvl].bfd_jump[n]].nlev;

      if (bf_cont > 0.0)
      {

        //internal jump to another macro atom level
        Q_matrix[uplvl][target_level] += Qcont = bf_cont * config[target_level].ex;     //energy of lower state

        //jump to the k-packet pool (we used to call this "deactivation")
        Q_matrix[uplvl][nlevels_macro] += Qcont_kpkt = q_recomb (cont_ptr, t_e) * ne * ne * (config[uplvl].ex - config[target_level].ex);       //energy difference

        //deactivation back to r-packet
        R_matrix[uplvl][uplvl] += Rcont = ne * sp_rec_rate * (config[uplvl].ex - config[target_level].ex);   //energy difference

        Q_norm[uplvl] += Qcont + Qcont_kpkt + Rcont;
      }
    }

    /* Now upwards jumps. */

    /* bound-bound */
    for (n = 0; n < nbbu; n++)
    {
      line_ptr = &line[config[uplvl].bbu_jump[n]];
      rad_rate = (b12 (line_ptr) * mplasma->jbar_old[config[uplvl].bbu_indx_first + n]);

      coll_rate = q12 (line_ptr, t_e);  // this is multiplied by ne below

      if (coll_rate < 0)
      {
        coll_rate = 0;
      }
      target_level = line[config[uplvl].bbu_jump[n]].nconfigu;
      Q_matrix[uplvl][target_level] += Qcont = ((rad_rate) + (coll_rate * ne)) * config[uplvl].ex;      //energy of lower state


      Q_norm[uplvl] += Qcont;
    }

    /* bound-free */
    for (n = 0; n < nbfu; n++)
    {
      /* For bf ionization the jump probability is just gamma * energy
         gamma is the photoionisation rate. Stimulated recombination also included. */
      cont_ptr = &phot_top[config[uplvl].bfu_jump[n]];  //pointer to continuum

      /* first let us take care of the situation where the lower level is zero or close to zero */
      lower_density = den_config (xplasma, cont_ptr->nlev);
      if (lower_density >= DENSITY_PHOT_MIN)
      {
        density_ratio = den_config (xplasma, cont_ptr->uplev) / lower_density;
      }
      else
        density_ratio = 0.0;

      target_level = phot_top[config[uplvl].bfu_jump[n]].uplev;
      Qcont = (mplasma->gamma_old[config[uplvl].bfu_indx_first + n] - (mplasma->alpha_st_old[config[uplvl].bfu_indx_first + n] * xplasma->ne * density_ratio) + (q_ioniz (cont_ptr, t_e) * ne)) * config[uplvl].ex;     //energy of lower state

      /* this error condition can happen in unconverged hot cells where T_R >> T_E.
         for the moment we set to 0 and hope spontaneous recombiantion takes care of things */
      /* note that we check and report this in check_stimulated_recomb() in estimators.c once a cycle */
      if (Qcont < 0.)           //test (can be deleted eventually SS)
      {
        //Error ("Negative probability (matom, 6). Abort?\n");
        Qcont = 0.0;

      }
      Q_matrix[uplvl][target_level] += Qcont;
      Q_norm[uplvl] += Qcont;
    }


  }

  /* Now need to do k-packet processes */

  int escape_dummy = 0;
  int istat_dummy = 0;
  fill_kpkt_rates (xplasma, escape_dummy, istat_dummy);
  /* Cooling due to collisional transitions in lines and collision ionization [for macro atoms] constitute internal transitions from the k-packet pool to macro atom states. */
  kpacket_to_rpacket_rate = 0.0;        // keep track of rate for kpacket_to_rpacket channel

  for (i = 0; i < nlines; i++)
  {
    if (line[i].macro_info == 1 && geo.macro_simple == 0)       //line is for a macro atom
    {
      target_level = line[i].nconfigu;
      Q_matrix[nlevels_macro][target_level] += Qcont = mplasma->cooling_bb[i];
      Q_norm[nlevels_macro] += Qcont;
    }
    else
    { 
      /* XXX - ask stuart about this! */
      /* the idea here is that if it is a simple line then it *must* create an r-packet originally?? */
      kpacket_to_rpacket_rate += mplasma->cooling_bb[i];
    }
  }

  for (i = 0; i < nphot_total; i++)
  {
    if (phot_top[i].macro_info == 1 && geo.macro_simple == 0)   //part of macro atom
    {
      target_level = phot_top[i].uplev;
      Q_matrix[nlevels_macro][target_level] += Qcont = mplasma->cooling_bf[i];
      Q_norm[nlevels_macro] += Qcont;

    }
    else
    {
      /* XXX - ask stuart about this! */
      kpacket_to_rpacket_rate += mplasma->cooling_bf_col[i];
    }
  }

  /* Cooling due to other processes will always contribute to k-packet -> r-packet channel */

  kpacket_to_rpacket_rate += mplasma->cooling_bftot;
  kpacket_to_rpacket_rate += mplasma->cooling_adiabatic;
  kpacket_to_rpacket_rate += mplasma->cooling_ff + mplasma->cooling_ff_lofreq;
  R_matrix[nlevels_macro][nlevels_macro] += Rcont = kpacket_to_rpacket_rate;
  Q_norm[nlevels_macro] += Rcont;

  /* end kpacket */

  /* now in one step, we multiply by the identity matrix and normalise the probabilities
     this means that what is now stored in Q_matrix is no longer Q, but N=(I - Q) using Vogl 
     notation. We check that Q_norm is 0, because some states (ground states) can have 0 
     jumping probabilities and so zero normalisation to */
  for (uplvl = 0; uplvl < nrows; uplvl++)
  {
    for (target_level = 0; target_level < nrows; target_level++)
    {
      if (Q_norm[uplvl] > 0.0)
      {
        Q_matrix[uplvl][target_level] = -1. * Q_matrix[uplvl][target_level] / Q_norm[uplvl];
      }
    }
    Q_matrix[uplvl][uplvl] = 1.0;

    if (Q_norm[uplvl] > 0.0)
    {
      R_matrix[uplvl][uplvl] = R_matrix[uplvl][uplvl] / Q_norm[uplvl];
    }
  }


  /* Check normalisation of the matrix. the diagonals of R, minus N which is now stored in 
     Q_matrix, should be zero */
  for (uplvl = 0; uplvl < nrows; uplvl++)
  {
    norm = R_matrix[uplvl][uplvl];
    for (target_level = 0; target_level < nlevels_macro + 1; target_level++)
    {
      norm -= Q_matrix[uplvl][target_level];
    }

    /* throw an error if this normalisation is not zero */
    /* note that the ground state is a special case here (improve error check) */
    if ( (abs(norm) > 1e-15 && uplvl > 0) || sane_check(norm))
      Error ("calc_matom_matrix: matom accelerator matrix has bad normalisation for level %d: %8.4e\n", 
              norm, uplvl);
  }

  /* This next line produces an array of the correct size to hold the rate matrix */
  a_data = (double *) calloc (sizeof (double), nrows * nrows);

  /* We now copy our rate matrix into the prepared matrix */
  for (mm = 0; mm < nrows; mm++)
  {
    for (nn = 0; nn < nrows; nn++)
    {
      a_data[mm * nrows + nn] = Q_matrix[mm][nn];
    }
  }

  /* now get ready for the matrix operations. first let's assign variables for use with GSL */
  gsl_matrix_view N;
  gsl_matrix *inverse_matrix;
  gsl_matrix *output;
  gsl_permutation *p, *pp;
  int ierr, s;
  
  /* create a view into the array we just created */
  N = gsl_matrix_view_array (a_data, nrows, nrows);

  /* permuations are special structures that contain integers 0 to nrows-1, which can
   * be manipulated */
  p = gsl_permutation_alloc (nrows);
  inverse_matrix = gsl_matrix_alloc (nrows, nrows);

  /* first we do the LU decomposition */
  ierr = gsl_linalg_LU_decomp (&N.matrix, p, &s);

  /* from the decomposition, get the inverse of the Q matrix */
  ierr = gsl_linalg_LU_invert (&N.matrix, p, inverse_matrix);

  /* We now copy our rate matrix into the prepared matrix */
  for (mm = 0; mm < nrows; mm++)
  {
    for (nn = 0; nn < nrows; nn++)
    {
      /* copy the matrix from the gsl object to the array we return to user */
      /* in Christian Vogl's notation this is doing his equation 3: B = (N * R)
         where N is the inverse matrix we have just calculated. */
      /* the reason this matrix multiplication is so simple here is because R is a diagonal matrix */
      matom_matrix[mm][nn] = gsl_matrix_get (inverse_matrix, mm, nn) * R_matrix[nn][nn];
    }
  }
 
  /* free memory */
  gsl_permutation_free (p);
  gsl_matrix_free (inverse_matrix);
  free (a_data);
  free (R_matrix);
  free (Q_matrix);
}


/**********************************************************/
/** 
 * @brief calculate the cooling rates for the conversion of k-packets.
 *
 * @param [in] PlasmaPtr  xplasma  
 * @param [in,out] double **matom_matrix 
 *        the 2D matrix array we will populate with normalised probabilities 
 *
 *
 * @details
 * 
 *
 **********************************************************/

int
fill_kpkt_rates (xplasma, escape, istat)
     PlasmaPtr xplasma;
     int *escape;
     int *istat;
{

  int i;
  int ulvl;
  double cooling_bf[nphot_total];
  double cooling_bf_col[nphot_total];   //collisional cooling in bf transitions
  double cooling_bb[NLINES];
  double cooling_adiabatic;
  struct topbase_phot *cont_ptr;
  struct lines *line_ptr;
  double cooling_normalisation;
  double destruction_choice;
  double electron_temperature;
  double cooling_bbtot, cooling_bftot, cooling_bf_coltot;
  double lower_density, upper_density;
  double cooling_ff, upweight_factor;
  WindPtr one;

  MacroPtr mplasma;

  double coll_rate, rad_rate;
  double freqmin, freqmax;

  one = &wmain[xplasma->nwind];


  /* Idea is to calculated the cooling
     terms for all the processes and then choose one at random to destroy the k-packet and
     turn it back into a photon bundle. 
     The routine considers bound-free, collision excitation and ff
     emission. */

  /* The loop here runs over all the bf processes, regardless of whether they are macro or simple
     ions. The reason that we can do this is that in the macro atom method stimulated recombination
     has already been considered before creating the k-packet (by the use of gamma-twiddle) in "scatter"
     and so we only have spontaneous recombination to worry about here for ALL cases. */

  mplasma = &macromain[xplasma->nplasma];

  electron_temperature = xplasma->t_e;

  /* JM 1511 -- Fix for issue 187. We need band limits for free free packet
     generation (see call to one_ff below) */
  freqmin = xband.f1[0];
  freqmax = ALPHA_FF * xplasma->t_e / H_OVER_K;

  /* ksl This is a Bandaid for when the temperatures are very low */
  /* in this case cooling_ff should be low compared to cooling_ff_lofreq anyway */
  if (freqmax < 1.1 * freqmin)
  {
    freqmax = 1.1 * freqmin;
  }

  /* ksl 091108 - If the kpkt destruction rates for this cell are not known they are calculated here.  This happens
   * every time the wind is updated */

  if (mplasma->kpkt_rates_known != 1)
  {
    cooling_normalisation = 0.0;
    cooling_bftot = 0.0;
    cooling_bbtot = 0.0;
    cooling_ff = 0.0;
    cooling_bf_coltot = 0.0;

    /* JM 1503 -- we used to loop over ntop_phot here, 
       but we should really loop over the tabulated Verner Xsections too
       see #86, #141 */
    for (i = 0; i < nphot_total; i++)
    {
      cont_ptr = &phot_top[i];
      ulvl = cont_ptr->uplev;

      if (cont_ptr->macro_info == 1 && geo.macro_simple == 0)
      {
        upper_density = den_config (xplasma, ulvl);
        /* SS July 04 - for macro atoms the recombination coefficients are stored so use the
           stored values rather than recompue them. */
        cooling_bf[i] = mplasma->cooling_bf[i] =
          upper_density * H * cont_ptr->freq[0] * (mplasma->recomb_sp_e[config[ulvl].bfd_indx_first + cont_ptr->down_index]);
        // _sp_e is defined as the difference 
      }
      else
      {
        upper_density = xplasma->density[cont_ptr->nion + 1];

        cooling_bf[i] = mplasma->cooling_bf[i] = upper_density * H * cont_ptr->freq[0] * (xplasma->recomb_simple[i]);
      }


      /* Note that the electron density is not included here -- all cooling rates scale
         with the electron density. */
      if (cooling_bf[i] < 0)
      {
        Error ("kpkt: bf cooling rate negative. Density was %g\n", upper_density);
        Error ("alpha_sp(cont_ptr, xplasma,2) %g \n", alpha_sp (cont_ptr, xplasma, 2));
        Error ("i, ulvl, nphot_total, nion %d %d %d %d\n", i, ulvl, nphot_total, cont_ptr->nion);
        Error ("nlev, z, istate %d %d %d \n", cont_ptr->nlev, cont_ptr->z, cont_ptr->istate);
        Error ("freq[0] %g\n", cont_ptr->freq[0]);
        cooling_bf[i] = mplasma->cooling_bf[i] = 0.0;
      }
      else
      {
        cooling_bftot += cooling_bf[i];
      }

      cooling_normalisation += cooling_bf[i];

      if (cont_ptr->macro_info == 1 && geo.macro_simple == 0)
      {
        /* Include collisional ionization as a cooling term in macro atoms. Don't include
           for simple ions for now.  SS */

        lower_density = den_config (xplasma, cont_ptr->nlev);
        cooling_bf_col[i] = mplasma->cooling_bf_col[i] = lower_density * H * cont_ptr->freq[0] * q_ioniz (cont_ptr, electron_temperature);

        cooling_bf_coltot += cooling_bf_col[i];

        cooling_normalisation += cooling_bf_col[i];

      }



    }

    /* end of loop over nphot_total */

    for (i = 0; i < nlines; i++)
    {
      line_ptr = &line[i];
      if (line_ptr->macro_info == 1 && geo.macro_simple == 0)
      {                         //It's a macro atom line and so the density of the upper level is stored
        cooling_bb[i] = mplasma->cooling_bb[i] =
          den_config (xplasma, line_ptr->nconfigl) * q12 (line_ptr, electron_temperature) * line_ptr->freq * H;

        /* Note that the electron density is not included here -- all cooling rates scale
           with the electron density so I've factored it out. */
      }
      else
      {                         //It's a simple line. Get the upper level density using two_level_atom

        two_level_atom (line_ptr, xplasma, &lower_density, &upper_density);

        /* the collisional rate is multiplied by ne later */
        coll_rate = q21 (line_ptr, electron_temperature) * (1. - exp (-H_OVER_K * line_ptr->freq / electron_temperature));

        cooling_bb[i] =
          (lower_density * line_ptr->gu / line_ptr->gl -
           upper_density) * coll_rate / (exp (H_OVER_K * line_ptr->freq / electron_temperature) - 1.) * line_ptr->freq * H;

        rad_rate = a21 (line_ptr) * p_escape (line_ptr, xplasma);

        /* Now multiply by the scattering probability - i.e. we are only going to consider bb cooling when
           the photon actually escapes - we don't to waste time by exciting a two-level macro atom only so that
           it makes another k-packet for us! (SS May 04) */

        cooling_bb[i] *= rad_rate / (rad_rate + (coll_rate * xplasma->ne));
        mplasma->cooling_bb[i] = cooling_bb[i];
      }

      if (cooling_bb[i] < 0)
      {
        cooling_bb[i] = mplasma->cooling_bb[i] = 0.0;
      }
      else
      {
        cooling_bbtot += cooling_bb[i];
      }
      cooling_normalisation += cooling_bb[i];
    }

    /* end of loop over nlines  */


    /* 57+ -- This might be modified later since we "know" that xplasma cannot be for a grid with zero
       volume.  Recall however that vol is part of the windPtr */
    if (one->vol > 0)
    {
      cooling_ff = mplasma->cooling_ff = total_free (one, xplasma->t_e, freqmin, freqmax) / xplasma->vol / xplasma->ne; // JM 1411 - changed to use filled volume
      cooling_ff += mplasma->cooling_ff_lofreq = total_free (one, xplasma->t_e, 0.0, freqmin) / xplasma->vol / xplasma->ne;
    }
    else
    {
      /* SS June 04 - This should never happen, but sometimes it does. I think it is because of 
         photons leaking from one cell to another due to the push-through-distance. It is sufficiently
         rare (~1 photon in a complete run of the code) that I'm not worrying about it for now but it does
         indicate a real problem somewhere. */

      /* SS Nov 09: actually I've not seen this problem for a long
         time. Don't recall that we ever actually fixed it,
         however. Perhaps the improved volume calculations
         removed it? We delete this whole "else" if we're sure
         volumes are never zero. */

      cooling_ff = mplasma->cooling_ff = mplasma->cooling_ff_lofreq = 0.0;
      Error ("kpkt: A scattering event in cell %d with vol = 0???\n", one->nwind);
      //Diagnostic      return(-1);  //57g -- Cannot diagnose with an exit
      *escape = 1;
      *istat = P_ERROR_MATOM;
      return (0);
    }


    if (cooling_ff < 0)
    {
      Error ("kpkt: ff cooling rate negative. Abort.");
      *escape = 1;
      *istat = P_ERROR_MATOM;
      return (0);
    }
    else
    {
      cooling_normalisation += cooling_ff;
    }


    /* JM -- 1310 -- we now want to add adiabatic cooling as another way of destroying kpkts
       this should have already been calculated and stored in the plasma structure. Note that 
       adiabatic cooling does not depend on type of macro atom excited */

    /* note the units here- we divide the total luminosity of the cell by volume and ne to give cooling rate */

    cooling_adiabatic = xplasma->cool_adiabatic / xplasma->vol / xplasma->ne;   // JM 1411 - changed to use filled volume


    if (geo.adiabatic == 0 && cooling_adiabatic > 0.0)
    {
      Error ("Adiabatic cooling turned off, but non zero in cell %d", xplasma->nplasma);
    }


    /* JM 1302 -- Negative adiabatic coooling- this used to happen due to issue #70, where we incorrectly calculated dvdy, 
       but this is now resolved. Now it should only happen for cellspartly in wind, because we don't treat these very well.
       Now, if cooling_adiabatic < 0 then set it to zero to avoid runs exiting for part in wind cells. */
    if (cooling_adiabatic < 0)
    {
      Error ("kpkt: Adiabatic cooling negative! Major problem if inwind (%d) == 0\n", one->inwind);
      Log ("kpkt: Setting adiabatic kpkt destruction probability to zero for this matom.\n");
      cooling_adiabatic = 0.0;
    }

    /* When we generate photons in the wind, from photon_gen we need to prevent deactivation by non-radiative cooling
     * terms.  If mode is True we include adiabatic cooling
     * this is now dealt with by setting cooling_adiabatic to 0 
     */

    cooling_normalisation += cooling_adiabatic;

    mplasma->cooling_bbtot = cooling_bbtot;
    mplasma->cooling_bftot = cooling_bftot;
    mplasma->cooling_bf_coltot = cooling_bf_coltot;
    mplasma->cooling_adiabatic = cooling_adiabatic;
    mplasma->cooling_normalisation = cooling_normalisation;
    mplasma->kpkt_rates_known = 1;

  }

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
 *
 * ###Notes###
***********************************************************/

double
f_matom_emit_accelerate (w, p, nres, upper, fmin, fmax)
     WindPtr w;
     PhotPtr p;
     int *nres;
     int upper;
     double fmin, fmax;
{
  struct lines *line_ptr;
  struct topbase_phot *cont_ptr;
  int uplvl;
  double eprbs[NBBJUMPS + NBFJUMPS], eprbs_band[NBBJUMPS + NBFJUMPS];
  double penorm, penorm_band;
  double threshold, run_tot;
  double sp_rec_rate;
  int n, m;
  int nbbd, nbfd;
  double t_e, ne;
  double bb_cont;
  WindPtr one;
  PlasmaPtr xplasma;
  double flast, fthresh, bf_int_full, bf_int_inrange;


  one = &w[p->grid];            //This is to identify the grip cell in which we are
  xplasma = &plasmamain[one->nplasma];

  t_e = xplasma->t_e;           //electron temperature 
  ne = xplasma->ne;             //electron number density


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
  penorm_band = 0.0;

  for (n = 0; n < nbbd + nbfd; n++)
  {
    eprbs[n] = 0;               //stores the individual emission probabilities SS
    eprbs_band[n] = 0;
  }
  /* Finished zeroing. */


  /* bb */

  /* First downward jumps. */

  for (n = 0; n < nbbd; n++)
  {
    line_ptr = &line[config[uplvl].bbd_jump[n]];
    /* Since we are only interested in making an r-packet here we can (a) ignore collisional
       deactivation and (b) ignore lines outside the frequency range of interest. */
    bb_cont = (a21 (line_ptr) * p_escape (line_ptr, xplasma));
    
    eprbs[n] = bb_cont * (config[uplvl].ex - config[line[config[uplvl].bbd_jump[n]].nconfigl].ex);    //energy difference
    
    if (eprbs[n] < 0.)        //test (can be deleted eventually SS)
      {
        Error ("Negative probability (matom, 2). Abort.");
        Exit (0);
      }

    penorm += eprbs[n];
    if ((line_ptr->freq > fmin) && (line_ptr->freq < fmax))   // correct range
      {
	eprbs_band[m] = eprbs[n];
	penorm_band += eprbs_band[m];
	m++;
      }
  }

  /* bf */
  for (n = 0; n < nbfd; n++)
  {
    cont_ptr = &phot_top[config[uplvl].bfd_jump[n]];    //pointer to continuum

    /* If the edge is above the frequency range we are interested in then we need not consider this
       bf process. */
    sp_rec_rate = alpha_sp (cont_ptr, xplasma, 0);
    eprbs[n] = sp_rec_rate * ne * (config[uplvl].ex - config[phot_top[config[uplvl].bfd_jump[n]].nlev].ex);   //energy difference
    if (eprbs[n] < 0.)        //test (can be deleted eventually SS)
      {
        Error ("Negative probability (matom, 4). Abort.");
        Exit (0);
      }
    penorm += eprbs[n];
    if (cont_ptr->freq[0] < fmax)  //means that it may contribute
      {
	eprbs_band[m] = eprbs[n];
	fthresh = cont_ptr->freq[0];  //first frequency in list
	flast = cont_ptr->freq[cont_ptr->np - 1];     //last frequency in list
	bf_int_full = qromb (alpha_sp_integrand, fthresh, flast, 1e-4);
	if (fthresh < fmin && flast > fmax)
	  {
	  bf_int_inrange = qromb (alpha_sp_integrand, fmin, fmax, 1e-4);
	  }
	else if (fthresh < fmin && flast < fmax)
	  {
	  bf_int_inrange = qromb (alpha_sp_integrand, fmin, flast, 1e-4);
	  }
	else if (fthresh > fmin && flast > fmax)
	  {
	    bf_int_inrange = qromb (alpha_sp_integrand, fthresh, flast, 1e-4);
	  }
	else if (fthresh > fmin && flast < fmax)
	  {
	    bf_int_inrange = qromb (alpha_sp_integrand, fthresh, fmax, 1e-4);
	  }
	else
	  {
	    Error("Something wrong here: f_matom_emit_accelerate broke the law!");
	  }
	penorm_band += eprbs_band[m] * bf_int_inrange  / bf_int_full;	
	m++;
      }
  }


  
  return (penorm_band / penorm);
}

/* The frequency and the value of nres have been set correctly. All done. */
