/***********************************************************/
/** @file  photo_gen_matom.c
 * @author ksl, ss, jm
 * @date   January, 2018
 *
 * @brief functions for calculating emissivities and generating photons from macro-atoms and k-packets.
 *   during the spectral cycles. The actual functions which do the jumps inside an activated 
 *  macro-atom are in matom.c. This is partly done to prevent overly long files (JM1504)
 *
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


/**********************************************************/
/** 
 * @brief      returns the specific band-limited luminosity in macro-atoms
 *
 * @param [in] int  mode   vvariable which controls whether or not we need to compute the
 *            emissivities (CALCULATE_MATOM_EMISSIVITIES) or use stored ones
 *            because we are restarting a spectral cycle (USE_STORED_MATOM_EMISSIVITIES)
 *            see #define statements in python.h and code in xdefine_phot().
 * @return double lum  The energy radiated by the deactivation of macro atoms in the wind in the 
 *            wavelength range required for the specrum calculation.
 *
 * @details
 * this routine calculates the luminosity in the band needed for the computation of the
 * spectrum. It gets the total energy radiated by the deactivation of macro atoms in the
 * required wavelength range. This can be a slow process, as there is no priori way 
 * (at least with our method) to work out where a photon is going to come out when a 
 * macro-atom is generated
 *
 * ### Notes ###
 * Consult Matthews thesis. 
 *
 **********************************************************/

double
get_matom_f (mode)
     int mode;
{
  int n, m, mm;
  double lum;
  double level_emit_doub[NLEVELS_MACRO], kpkt_emit_doub;
  int n_tries, n_tries_local;
  double norm;
  int which_out;
  int i, j;
  int my_nmin, my_nmax;         //These variables are used even if not in parallel mode


  if (mode == USE_STORED_MATOM_EMISSIVITIES)
  {
    Log ("geo.pcycle (%i) > 0, so using saved level emissivities for macro atoms\n", geo.pcycle);
  }

  else                          // we need to compute the emissivities
  {
#ifdef MPI_ON
    int num_mpi_cells, num_mpi_extra, position, ndo, n_mpi, num_comm, n_mpi2;
    int size_of_commbuffer;
    char *commbuffer;

    /* the commbuffer needs to communicate 2 variables and the number of macor levels, 
       plus the variable for how many cells each thread is doing */
    size_of_commbuffer = 8 * (3 + nlevels_macro) * (floor (NPLASMA / np_mpi_global) + 1);

    commbuffer = (char *) malloc (size_of_commbuffer * sizeof (char));
#endif

    /* if we are using the accelerated macro-atom scheme then we want to allocate an array 
       for the macro-atom probabilities and various other quantities */
#if (ACCELERATED_MACRO == TRUE)
    PlasmaPtr xplasma;
    int nrows = nlevels_macro + 1;
    double **matom_matrix = (double **) calloc (sizeof (double *), nrows);

    for (i = 0; i < nrows; i++)
    {
      matom_matrix[i] = (double *) calloc (sizeof (double), nrows);
    }
#else
    /* these variables are only used in the non-accelerated scheme */
    struct photon ppp;
    int ss, level_emit[NLEVELS_MACRO], kpkt_emit;
    double contribution;
    int nres, esc_ptr;
#endif


    /* add the non-radiative k-packet heating to the kpkt_abs quantity */
    get_kpkt_heating_f ();

    which_out = 0;
    n_tries = 5000000;
    geo.matom_radiation = 0;
    n_tries_local = 0;



    /* zero all the emissivity counters and check absorbed quantities */
    norm = 0;
    for (n = 0; n < NPLASMA; n++)
    {
      for (m = 0; m < nlevels_macro; m++)
      {
        norm += macromain[n].matom_abs[m];
        macromain[n].matom_emiss[m] = 0.0;
        if (sane_check (macromain[n].matom_abs[m]))
          Error ("matom_abs is %8.4e in matom %i level %i\n", macromain[n].matom_abs[m], n, m);
      }
      norm += plasmamain[n].kpkt_abs;
      plasmamain[n].kpkt_emiss = 0.0;
      if (sane_check (plasmamain[n].kpkt_abs))
        Error ("kpkt_abs is %8.4e in matom %i\n", plasmamain[n].kpkt_abs, n);
    }



    Log ("Calculating macro-atom and k-packet emissivities- this might take a while...\n");
    Log ("Number of macro-atom levels: %d\n", nlevels_macro);


    /* For MPI parallelisation, the following loop will be distributed over multiple tasks. 
       Note that the mynmim and mynmax variables are still used even without MPI on */
    my_nmin = 0;
    my_nmax = NPLASMA;

#ifdef MPI_ON

    num_mpi_cells = floor (NPLASMA / np_mpi_global);    // divide the cells between the threads
    num_mpi_extra = NPLASMA - (np_mpi_global * num_mpi_cells);  // the remainder from the above division

    /* this next loop splits the cells up between the threads. All threads with 
       rank_global<num_mpi_extra deal with one extra cell to account for the remainder */
    if (rank_global < num_mpi_extra)
    {
      my_nmin = rank_global * (num_mpi_cells + 1);
      my_nmax = (rank_global + 1) * (num_mpi_cells + 1);
    }
    else
    {
      my_nmin = num_mpi_extra * (num_mpi_cells + 1) + (rank_global - num_mpi_extra) * (num_mpi_cells);
      my_nmax = num_mpi_extra * (num_mpi_cells + 1) + (rank_global - num_mpi_extra + 1) * (num_mpi_cells);
    }
    ndo = my_nmax - my_nmin;

    Log_parallel ("Thread %d is calculating macro atom emissivities for macro atoms %d to %d\n", rank_global, my_nmin, my_nmax);

#endif


    for (n = my_nmin; n < my_nmax; n++)
    {

      /* JM 1309 -- these lines are just log statements which track progress, as this section
         can take a long time */
#ifdef MPI_ON
      if (n % 50 == 0)
        Log
          ("Thread %d is calculating  macro atom emissivity for macro atom %7d of %7d or %6.3f per cent\n",
           rank_global, n, my_nmax, n * 100. / my_nmax);
#else
      if (n % 50 == 0)
        Log ("Calculating macro atom emissivity for macro atom %7d of %7d or %6.3f per cent\n", n, my_nmax, n * 100. / my_nmax);
#endif

#if (ACCELERATED_MACRO == TRUE)

      /* use the accelerated macro-atom scheme */
      xplasma = &plasmamain[n];
      calc_matom_matrix (xplasma, matom_matrix);
      /* before we calculate the emissivities we need to know what fraction of the energy 
         from each level comes out in the frequency band we care about */

      for (i = 0; i < nlevels_macro; i++)
      {
        level_emit_doub[i] = f_matom_emit_accelerate (xplasma, i, geo.sfmin, geo.sfmax);
      }

      /* do the same for k-packets */
      kpkt_emit_doub = f_kpkt_emit_accelerate (xplasma, geo.sfmin, geo.sfmax);

      /* Now use the matrix to calculate the fraction of the absorbed energy that comes out in a given level */
      for (i = 0; i < nlevels_macro; i++)
      {
        for (j = 0; j < nlevels_macro; j++)
        {
          macromain[n].matom_emiss[j] += macromain[n].matom_abs[i] * matom_matrix[i][j];
        }
      }

      /* do the same for the thermal pool. we also normalise by banded_emiss_frac here */
      for (j = 0; j < nlevels_macro; j++)
      {
        macromain[n].matom_emiss[j] += plasmamain[n].kpkt_abs * matom_matrix[nlevels_macro][j];
        macromain[n].matom_emiss[j] *= (1.0 * level_emit_doub[j]);
      }
      plasmamain[n].kpkt_emiss += plasmamain[n].kpkt_abs * matom_matrix[nlevels_macro][nlevels_macro];
      plasmamain[n].kpkt_emiss *= (1.0 * kpkt_emit_doub);
    }
#else
      for (m = 0; m < nlevels_macro + 1; m++)
      {
        if ((m == nlevels_macro && plasmamain[n].kpkt_abs > 0) || (m < nlevels_macro && macromain[n].matom_abs[m] > 0))
        {
          if (m < nlevels_macro)
          {
            if (macromain[n].matom_abs[m] > 0)
            {
              n_tries_local = (n_tries * macromain[n].matom_abs[m] / norm) + 10;
            }
            else
            {
              n_tries_local = 0;
            }
          }
          else if (m == nlevels_macro)
          {
            if (plasmamain[n].kpkt_abs > 0)
            {
              n_tries_local = (n_tries * plasmamain[n].kpkt_abs / norm) + 10;
            }
            else
            {
              n_tries_local = 0;
            }
          }
          /* We know "matom_abs" is the amount of energy absorbed by
             each macro atom level in each cell. We now want to determine what fraction of
             that energy re-appears in the frequency range we want and from which macro atom
             level it is re-emitted (or if it appears via a k-packet). */

          for (mm = 0; mm < NLEVELS_MACRO; mm++)
          {
            level_emit[mm] = 0;
          }
          kpkt_emit = 0;
          if (n_tries_local > 0)
          {
            for (ss = 0; ss < n_tries_local; ss++)
            {
              if (m < nlevels_macro)
              {
                /* Dealing with excitation of a macro atom level. */
                /* First find a suitable transition in which to excite the macro atom (to get it
                   in the correct starting level. */

                /* As Knox pointed out in an earlier comment
                   here, the next loop isn't the best way to
                   do this - however I'm not sure we can
                   just replace nlines->nlines_macro since
                   we don't know for sure that the
                   macro_lines are the first lines in the
                   array. If this is really a bottle neck we
                   should set up a quicker way to identify
                   the level - the purpose is just to set
                   the macro atom state to the appropriate
                   upper level so all we need is some way to
                   tell it which state it should be. */

                nres = 0;
                while ((nres < nlines))
                {
                  if (lin_ptr[nres]->nconfigu != m)
                  {
                    nres += 1;
                  }
                  else
                  {
                    break;
                  }
                }


                if (nres > (nlines - 1))
                {
                  nres = NLINES + 1;
                  while (nres < (NLINES + 1 + nphot_total))
                  {
                    if (phot_top[nres - NLINES - 1].uplev != m)
                    {
                      nres += 1;
                    }
                    else
                    {
                      break;
                    }
                  }
                }

                if (nres > NLINES + nphot_total)
                {
                  Error ("Problem in get_matom_f (1). Abort. nres is %d, NLINES %d, nphot_total %d m %d %8.4e\n",
                         nres, NLINES, nphot_total, m, macromain[n].matom_abs[m]);
                  Exit (0);
                }

                ppp.nres = nres;
                ppp.grid = plasmamain[n].nwind;
                ppp.w = 0;
                /* this needs to be initialised because we set to istat to P_ADIABATIC
                   for adiabatic destruction */
                ppp.istat = P_INWIND;

                macro_gov (&ppp, &nres, 1, &which_out);


                /* Now a macro atom has been excited and followed until an r-packet is made. Now, if that 
                   r-packet is in the correct frequency range we record it. If not, we throw it away. */
              }
              else if (m == nlevels_macro)
              {
                /* kpkt case. */

                nres = -2;      //will do - just need
                //something that will tigger kpkt

                ppp.nres = nres;
                ppp.grid = plasmamain[n].nwind;
                ppp.w = 0;
                /* this needs to be initialised because we set to istat to P_ADIABATIC
                   for adiabatic destruction */
                ppp.istat = P_INWIND;

                macro_gov (&ppp, &nres, 2, &which_out);

                /* We have an r-packet back again. */
              }


              if (ppp.freq > geo.sfmin && ppp.freq < geo.sfmax && ppp.istat != P_ADIABATIC)
              {
                if (which_out == 1)
                {
                  if (nres < 0)
                  {
                    Error ("Negative out from matom?? Abort.\n");
                    Exit (0);
                  }

                  /* It was a macro atom de-activation. */
                  if (nres < nlines)
                  {             /* It was a line. */
                    level_emit[lin_ptr[nres]->nconfigu] += 1;
                  }
                  else
                  {
                    level_emit[phot_top[nres - NLINES - 1].uplev] += 1;
                  }
                }
                else if (which_out == 2)
                {
                  /* It was a k-packet de-activation. */
                  kpkt_emit += 1;
                }
                else
                {
                  Error ("Packet didn't emerge from matom or kpkt??? Abort. \n");
                  Exit (0);
                }
              }
            }


            /* Now we've done all the runs for this level de-activating so we can add the contributions
               to the level emissivities, the k-packet emissivity and the total luminosity. */

            for (mm = 0; mm < nlevels_macro; mm++)
            {
              contribution = 0;
              if (m < nlevels_macro)
              {
                macromain[n].matom_emiss[mm] += contribution = level_emit[mm] * macromain[n].matom_abs[m] / n_tries_local;
              }
              else if (m == nlevels_macro)
              {
                macromain[n].matom_emiss[mm] += contribution = level_emit[mm] * plasmamain[n].kpkt_abs / n_tries_local;
              }
            }

            if (m < nlevels_macro)
            {
              plasmamain[n].kpkt_emiss += kpkt_emit * macromain[n].matom_abs[m] / n_tries_local;
            }
            else if (m == nlevels_macro)
            {
              plasmamain[n].kpkt_emiss += kpkt_emit * plasmamain[n].kpkt_abs / n_tries_local;
            }
          }
        }
      }
    }

#endif

    /*This is the end of the update loop that is parallelised. We now need to exchange data between the tasks.
       This is done much the same way as in wind_update */
#ifdef MPI_ON

    /* JM Add an MPI Barrier here */
    MPI_Barrier (MPI_COMM_WORLD);

    for (n_mpi = 0; n_mpi < np_mpi_global; n_mpi++)
    {
      /* here we loop over the number of threads. If the thread is this thread then we pack the macromain information
         into memory and then broadcast it to ther other threads */
      position = 0;

      if (rank_global == n_mpi)
      {

        Log ("MPI task %d is working on matoms %d to max %d (total size %d).\n", rank_global, my_nmin, my_nmax, NPLASMA);

        MPI_Pack (&ndo, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        for (n = my_nmin; n < my_nmax; n++)
        {

          /* pack the number of the cell, and the kpkt and macro atom emissivites for that cell // */
          MPI_Pack (&n, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
          MPI_Pack (&plasmamain[n].kpkt_emiss, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
          MPI_Pack (macromain[n].matom_emiss, nlevels_macro, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);

        }

        Log_silent ("MPI task %d broadcasting matom emissivity information.\n", rank_global);
      }


      /* Set MPI_Barriers and broadcast information to other threads */
      MPI_Barrier (MPI_COMM_WORLD);
      MPI_Bcast (commbuffer, size_of_commbuffer, MPI_PACKED, n_mpi, MPI_COMM_WORLD);
      MPI_Barrier (MPI_COMM_WORLD);
      Log_parallel ("MPI task %d survived broadcasting matom emissivity information.\n", rank_global);



      position = 0;

      /* If not this thread then we unpack the macromain information from the other threads */

      if (rank_global != n_mpi)
      {
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &num_comm, 1, MPI_INT, MPI_COMM_WORLD);
        for (n_mpi2 = 0; n_mpi2 < num_comm; n_mpi2++)
        {

          /* unpack the number of the cell, and the kpkt and macro atom emissivites for that cell */
          MPI_Unpack (commbuffer, size_of_commbuffer, &position, &n, 1, MPI_INT, MPI_COMM_WORLD);
          MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].kpkt_emiss, 1, MPI_DOUBLE, MPI_COMM_WORLD);
          MPI_Unpack (commbuffer, size_of_commbuffer, &position, macromain[n].matom_emiss, nlevels_macro, MPI_DOUBLE, MPI_COMM_WORLD);
        }
      }
    }                           // end of parallelised section

    /* add an MPI Barrier after unpacking stage */
    MPI_Barrier (MPI_COMM_WORLD);
#endif

#ifdef MPI_ON
    free (commbuffer);
#endif
#if (ACCELERATED_MACRO == TRUE)
    for (i = 0; i < nrows; i++)
    {
      free (matom_matrix[i]);
    }

    free (matom_matrix);
#endif

  }                             // end of if loop which controls whether to compute the emissivities or not 


  /* this next loop just calculates lum to be the correct summed value in parallel mode */
  /* if mode == USE_STORED_MATOM_EMISSIVITIES this is all this routine does */
  lum = 0.0;                    // need to rezero, fixes segfault bug #59 

  for (n = 0; n < NPLASMA; n++)
  {

    for (mm = 0; mm < nlevels_macro; mm++)
    {
      lum += macromain[n].matom_emiss[mm];
    }
  }

  geo.matom_radiation = 1;

  return (lum);
}

/* All done. */


