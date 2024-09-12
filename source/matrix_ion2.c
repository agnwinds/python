
/***********************************************************/
/** @file  matrix_ion.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief  Contains routines which compute the ionization state using a matrix inversion technique
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <float.h>
#include "atomic.h"
#include "python.h"


/**********************************************************/
/**
 * @brief      A matrix solver for the ionization state in a cell
 *
 * @param [in,out] PlasmaPtr  xplasma   The plasma cell we are working on
 * @param [in] int  mode   What type of model to use for J_nu - 1=fitted model 2=blackbody
 * @return     0 if successful
 *
 * The abundances contained in xplasma are updated witht he results of
 * the calculation.
 *
 * @details
 * modes:
 *   1 - use a power law and or exponential model for J
 *   2 - use a dilute BB model for J defined by t_r and w
 * The general structure of matrix_ion_populations is as follows:
 * First we compute all of the rates that we have data for linking all ionization stages
 * We make an initial guess at the electron density
 * We attempt to solve the ionization rates matrix, and produce a new guess at the electron density
 * We then calculate new rates, re-solve and proceed until the electron density converges.
 * We solve the matrix equation A x = b for the vector x, where A is a square rate matrix linking
 * each state with each other state (often containing zeros if the ions arent linked) b is a vector
 * containing the total density of each element arranged in a certain way (see later) and
 * x is a vector containing the relative ion populations. We invert A to get x=bA^-1
 *
 * ### Notes ###
 * Uses a relative abundance scheme, in order to reduce large number issues
 *
 * Various parameters for the calculation, and in particular the t_e are passed
 * via the PlasmaPtr
 *
 **********************************************************/

int
matrix_ion_populations2 (xplasma, mode)
     PlasmaPtr xplasma;
     int mode;

{
  double elem_dens[NELEMENTS];  //The density of each element
  int nn, mm, nrows;
  double rate_matrix[nions][nions];     //The rate matrix that we are going to try and solve
  double newden[NIONS];         //A temporary array to hold our intermediate solutions
  double nh, nh1, nh2, t_e;
  double xne, xxne, xxxne;      //Various stores for intermediate guesses at electron density
  double b_temp[nions];         //The b matrix
  double *b_data, *a_data;      //These arrays are allocated later and sent to the matrix solver
  double *populations;          //This array is allocated later and is retrieved from the matrix solver
  int matrix_err, niterate;     //counters for errors and the number of iterations we have tried to get a converged electron density
  double xnew;
  int xion[nions];              // This array keeps track of what ion is in each line
//OLD  int xelem[nions];             // This array keeps track of the element for each ion
  double pi_rates[nions];       //photoionization rate coefficients
  double rr_rates[nions];       //radiative recombination rate coefficients
  double inner_rates[n_inner_tot];      //This array contains the rates for each of the inner shells. Where they go to requires the electron yield array

  nh1 = nh2 = 0;

  /* Copy some quantities from the cell into local variables */

  nh = xplasma->rho * rho2nh;   // The number density of hydrogen ions - computed from density
  t_e = xplasma->t_e;           // The electron temperature in the cell - used for collisional processes


  /* We now calculate the total abundances for each element to allow us to use fractional abundances */

  /* Belt and braces, we zero our array */
  for (mm = 0; mm < ion[nions - 1].z + 1; mm++)
  {
    elem_dens[mm] = 0.0;
  }

  /* Now we populate the elemental abundance array */
  for (mm = 0; mm < nions; mm++)
  {
    elem_dens[ion[mm].z] = elem_dens[ion[mm].z] + xplasma->density[mm];
  }

  /* Dielectronic recombination, collisional ionization coefficients, three body recombination and
     charge_exchange rate coefficients depend only on electron temperature, calculate them now -
     they will not change they are all stored in global arrays */

  compute_dr_coeffs (t_e);
  compute_di_coeffs (t_e);
  compute_qrecomb_coeffs (t_e);
  compute_ch_ex_coeffs (t_e);

  /* In the following loop, over all ions in the simulation, we compute the radiative recombination rates, and photionization
     rates OUT OF each ionization stage. The PI rates are calculated either using the modelled mean intensity in a cell, or
     using the dilute blackbody approximation, depending on which mode we are in. At the same time, we copy the ion densities
     from the plasma structure into a local array. We will only overwrite the numbers in the structure if we believe the
     results are an improvement on what is there. */

  for (mm = 0; mm < nions; mm++)
  {
    newden[mm] = xplasma->density[mm] / elem_dens[ion[mm].z];   // newden is our local fractional density array
    xion[mm] = mm;              // xion is an array we use to track which ion is in which row of the matrix
    if (mm != ele[ion[mm].nelem].firstion)      // We can recombine since we are not in the first ionization stage
    {
      rr_rates[mm] = total_rrate (mm, xplasma->t_e);    // radiative recombination rates
    }
    if (ion[mm].istate != ele[ion[mm].nelem].istate_max)        // we can photoionize, since we are not in the highest ionization state
    {
      if (mode == NEBULARMODE_MATRIX_BB)
      {
        pi_rates[mm] = calc_pi_rate (mm, xplasma, 2, 1);        // PI rate for the BB model
      }
      else if (mode == NEBULARMODE_MATRIX_SPECTRALMODEL)
      {
        pi_rates[mm] = calc_pi_rate (mm, xplasma, 1, 1);        // PI rate for an explicit spectral model
      }
      else if (mode == NEBULARMODE_MATRIX_ESTIMATORS)
      {
        pi_rates[mm] = xplasma->ioniz[mm] / xplasma->density[mm];       // PI rate logged during the photon passage
      }
      else
      {
        Error ("matrix_ion_populations: Unknown mode %d\n", mode);
        Exit (0);
      }
    }

/* ksl: This whole lope is not actually used, and produces warnings in gcc11 */
//OLD    for (nn = 0; nn < nelements; nn++)
//OLD    {
//OLD      if (ion[mm].z == ele[nn].z)
//OLD      {
//OLD        xelem[mm] = nn;         /* xelem logs which element each row in the arrays refers to. This is important because we need
//OLD                                   to know what the total density will be for a group of rows all representing the same
//OLD                                   element. */
//OLD      }
//OLD    }
  }


  /* The next loop generates the inner shell ionization rates, if they are present in the atomic data and
     we wish to compute auger ionizaion rates. This only computes the rates out of each ion, we also need to
     consult the electron yield array if we are to compute the change in electron number */

  if (geo.auger_ionization == 1)
  {
    for (mm = 0; mm < n_inner_tot; mm++)
    {
      if (mode == NEBULARMODE_MATRIX_BB)
      {
        inner_rates[mm] = calc_pi_rate (mm, xplasma, 2, 2);
      }
      else if (mode == NEBULARMODE_MATRIX_SPECTRALMODEL)
      {
        inner_rates[mm] = calc_pi_rate (mm, xplasma, 1, 2);
      }
      else if (mode == NEBULARMODE_MATRIX_ESTIMATORS)   //We are using estimators - so we will need to reorder the rates from freq order to cross section order to match the electron yields. This takes time, so we only do it if we need to.
      {
        for (nn = 0; nn < n_inner_tot; nn++)
        {
          if (inner_cross[mm].nion == inner_cross_ptr[nn]->nion && inner_cross[mm].freq[0] == inner_cross_ptr[nn]->freq[0])     //Check for a match
          {
            inner_rates[mm] = xplasma->inner_ioniz[nn] / xplasma->density[inner_cross_ptr[nn]->nion];
          }
        }
      }
    }
  }

  /* This next line sets the partition function for each ion. This has always been the place here python calculates the
     partition functions and sets the level densities for each ion. It needs to be done, or other parts of the code which rely
     on sensible level populations don't work properly. In the case of the dilute blackbody, the code works well, however we do
     not currently (v78) have a procedure to calucate the levels for a spectral model case. We therefore call partition
     functions with mode 4 - this is a special mode which forces the ions to be in the ground state. This is reasonable for a
     radiation dominated plasma, since any excitations into higher states will quickly be followed by radative de-excitation.
     We should really do better though... */

  if (mode == NEBULARMODE_MATRIX_BB)
  {
    partition_functions (xplasma, NEBULARMODE_ML93);    // We use t_r and the radiative weight
  }
  else if (mode == NEBULARMODE_MATRIX_SPECTRALMODEL || mode == NEBULARMODE_MATRIX_ESTIMATORS)
  {
    partition_functions (xplasma, NEBULARMODE_LTE_GROUND);      // Set to ground state
  }

  /* Next we need to obtain an initial guess for the electron density. In the past this has been done by calculating the
     hydrogen density directly from the Saha equation at the current electron temperature. In initial testing of this mode -
     this seemed to be a little unstable. At the moment, we use the last ionization state to compute n_e. For the first
     iteraion, this will just be calculated in LTE from the initial temperature structure of the wind, which will give almost
     the same result as the original procedure, or for successive calculations, it should be a better guess. I've leftin the
     original code, commented out...  */

  xne = xxne = xxxne = get_ne (xplasma->density);       //Even though the abundances are fractional, we need the real electron density


  /* xne is the current working number xxne */


  /* We are now going to iterate on the electron density - MAXITERATIONS is set in python.h and is currently (78) set to 200.
     We would normally expect to converge much fater than this */

  niterate = 0;
  while (niterate < MAXITERATIONS)
  {

    /* In order to compute charge exchange reactions, we need the number density of neutral and ionized hydrogen
       the current ethod of matrix solving does not allow is to link ions multiplicatively so we need to compute
       these outside the loop. If hydrogen were not dominant - this could cause iddues since the ionization(recombination)
       of a metal or helium via this process would cause an identical recombination(ionization) of hydrogen. The loop below
       is a bit belt and braces, since one would almost always expect hydrogen to be ion=0 and 1  */

    for (nn = 0; nn < nions; nn++)
    {
      if (ion[nn].z == 1 && ion[nn].istate == 1)
      {
        nh1 = newden[nn] * elem_dens[ion[nn].z];
      }
      if (ion[nn].z == 1 && ion[nn].istate == 2)
      {
        nh2 = newden[nn] * elem_dens[ion[nn].z];
      }
    }


    populate_ion_rate_matrix (rate_matrix, pi_rates, inner_rates, rr_rates, b_temp, xne, nh1, nh2);


    /* The array is now fully populated, and we can begin the process of solving it */

    nrows = nions;              /* This is a placeholder, we may end up removing rows and columns that have no rates (or very
                                   low rates) connecting them to other rows. This may improve stability but will need to be
                                   done carefully */

    /* Here we solve the matrix equation M x = b, where x is our vector containing level populations as a fraction w.r.t the
       whole element. The actual LU decomposition- the process of obtaining a solution - is done by the routine solve_matrix() */

    a_data = (double *) calloc (sizeof (double), nrows * nrows);
    populations = (double *) calloc (nrows, sizeof (double));

    /* This b_data column matrix is the total number density for each element, placed into the row which relates to the neutral
       ion. This matches the row in the rate matrix which is just 1 1 1 1 for all stages. NB, we could have chosen any line for
       this. */

    b_data = (double *) calloc (nrows, sizeof (double));

    /* We now copy our rate matrix into the prepared matrix */
    for (mm = 0; mm < nrows; mm++)
    {
      for (nn = 0; nn < nrows; nn++)
      {
        a_data[mm * nrows + nn] = rate_matrix[mm][nn];  /* row-major */
      }
    }

    for (nn = 0; nn < nrows; nn++)
    {
      b_data[nn] = b_temp[nn];
    }

    matrix_err = solve_matrix (a_data, b_data, nrows, populations, xplasma->nplasma);

    if (matrix_err)
    {
      Error ("matrix_ion_populations: %s\n", get_matrix_error_string (matrix_err));
    }

    /* free memory */
    free (a_data);
    free (b_data);

    if (matrix_err == 4)
    {
      free (populations);
      return (-1);
    }

    /* Calculate level populations for macro-atoms */
    if (geo.macro_ioniz_mode == MACRO_IONIZ_MODE_ESTIMATORS)
    {
      int mp_err = macro_pops (xplasma, xne);
      if (mp_err != EXIT_SUCCESS)
      {
        return -1;
      }
    }

    /* We now have the populations of all the ions stored in the matrix populations. We copy this data into the newden array
       which will temperarily store all the populations. We wont copy this to the plasma structure until we are sure thatwe
       have made things better. We just loop over all ions, and copy. The complexity in the loop is to future proof us against
       the possibility that there are some ions that are not included in the matrix scheme becuse there is no route into or
       out of it. */

    for (nn = 0; nn < nions; nn++)
    {
      newden[nn] = 0.0;

      /* if the ion is being treated by macro_pops then use the populations just computed */
      if ((ion[nn].macro_info == TRUE) && (geo.macro_simple == FALSE) && (geo.macro_ioniz_mode == MACRO_IONIZ_MODE_ESTIMATORS)
          && (modes.no_macro_pops_for_ions == FALSE))
      {
        newden[nn] = xplasma->density[nn] / elem_dens[ion[nn].z];
      }

      /* if the ion is "simple" then find it's calculated ionization state in populations array */
      else
      {

        for (mm = 0; mm < nrows; mm++)  // inner loop over the elements of the population array
        {
          if (xion[mm] == nn)   // if this element contains the population of the ion is question
          {
            newden[nn] = populations[mm];       // get the population
          }
        }
      }

      if (newden[nn] < DENSITY_MIN)     // this wil also capture the case where population doesnt have a value for this ion
        newden[nn] = DENSITY_MIN;
    }
    free (populations);


/* We need to get the 'true' new electron density so we need to do a little loop here to compute it */


    xnew = 0.0;
    for (nn = 0; nn < nions; nn++)
    {
      xnew += newden[nn] * (ion[nn].istate - 1) * elem_dens[ion[nn].z];
    }



    if (xnew < DENSITY_MIN)
      xnew = DENSITY_MIN;       /* fudge to keep a floor on ne */


    if (fabs ((xne - xnew) / (xnew)) < FRACTIONAL_ERROR || xnew < 1.e-6)        /* We have converged, or have the situation where we
                                                                                   have a neutral plasma */
    {
      break;                    /* Break out of the while loop - we have finished our iterations */
    }
    xne = xxxne = (xnew + xne) / 2.;    /* New value of ne */

    niterate++;


    if (niterate == MAXITERATIONS)
    {
      Error ("matrix_ion_populations: failed to converge for cell %i t %e nh %e xnew %e\n", xplasma->nplasma, t_e, nh, xnew);


      for (nn = 0; nn < geo.nxfreq; nn++)
      {
        Log
          ("numin= %e (%e) numax= %e (%e) Model= %2d PL_log_w= %e PL_alpha= %e Exp_w= %e EXP_temp= %e\n",
           xplasma->fmin_mod[nn], geo.xfreq[nn], xplasma->fmax_mod[nn],
           geo.xfreq[nn + 1], xplasma->spec_mod_type[nn],
           xplasma->pl_log_w[nn], xplasma->pl_alpha[nn], xplasma->exp_w[nn], xplasma->exp_temp[nn]);
      }

      Error ("matrix_ion_populations: xxne %e theta %e\n", xxne);

      return (-1);              /* If we get to MAXITERATIONS, we return without copying the new populations into plasma */
    }
  }                             /* This is the end of the iteration loop */


  xplasma->ne = xnew;
  for (nn = 0; nn < nions; nn++)
  {
    /* If statement added here to suppress interference with macro populations */
    if (ion[nn].macro_info == FALSE || geo.macro_ioniz_mode == MACRO_IONIZ_MODE_NO_ESTIMATORS || geo.macro_simple == TRUE)
    {
      xplasma->density[nn] = newden[nn] * elem_dens[ion[nn].z]; //We return to absolute densities here
    }
    if ((sane_check (xplasma->density[nn])) || (xplasma->density[nn] < 0.0))
      Error ("matrix_ion_populations: ion %i has population %8.4e in cell %i\n", nn, xplasma->density[nn], xplasma->nplasma);
  }

  xplasma->ne = get_ne (xplasma->density);

  if (n_charge_exchange > 0)
  {
    xplasma->heat_ch_ex = ch_ex_heat (&wmain[xplasma->nwind], xplasma->t_e);    //Compute the charge exchange heating
    xplasma->heat_tot += xplasma->heat_ch_ex;
  }

  /*We now need to populate level densities in order to later calculate line emission (for example).
     We call partition functions to do this. At present, we do not have a method for accurately computing
     level populations for modelled mean intensities. We get around this by setting all level populations
     other than ground state to zero - this is a pretty good situation for very dilute radiation fields,
     which is what we are normally looking at, however one should really do better in the long term.
   */

  partition_functions (xplasma, NEBULARMODE_LTE_GROUND);

  return (0);
}

