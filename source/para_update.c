
/***********************************************************/
/** @file  para_update.c
 * @author ksl, jm
 * @date   January, 2018
 *
 * @brief  routines for communicating MC estimators and spectra between MPI threads.
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

/**********************************************************/
/**
 * @brief helper routine for splitting up tasks in MPI
 * @param   [in]      int   rank       processor rank (typically set from rank_global)
 * @param   [in]      int   ntotal     total number of tasks, e.g. NPLASMA
 * @param   [in]       int   nproc      total number of MPI processors
 * @param   [in,out]  int   *my_nmax   pointer to integer value of first task
 * @param   [in,out]  int   *my_nmax   pointer to integer value of final task
 * @return            int   ndo        number of tasks this thread is working on

 * @details  For a process with ntotal tasks,
 * this routine calculates which thread will be given each task.
 * typically ntotal is NPLASMA and the thread is splitting up wind cells.
 * The routine deals with remainders by distributing the remainder over the
 * threads if the cells do not divide evenly by thread
 **********************************************************/

int
get_parallel_nrange (rank, ntotal, nproc, my_nmin, my_nmax)
     int rank;
     int ntotal;
     int nproc;
     int *my_nmin;
     int *my_nmax;
{
  /* divide the cells between the threads */
  int ndo;
  int num_mpi_cells = floor (ntotal / nproc);

  /* the remainder from the above division */
  int num_mpi_extra = ntotal - (nproc * num_mpi_cells);

  /* this section distributes the remainder over the threads if the cells
     do not divide evenly by thread */
  if (rank < num_mpi_extra)
  {
    *my_nmin = rank_global * (num_mpi_cells + 1);
    *my_nmax = (rank_global + 1) * (num_mpi_cells + 1);
  }
  else
  {
    *my_nmin = num_mpi_extra * (num_mpi_cells + 1) + (rank - num_mpi_extra) * (num_mpi_cells);
    *my_nmax = num_mpi_extra * (num_mpi_cells + 1) + (rank - num_mpi_extra + 1) * (num_mpi_cells);
  }
  ndo = *my_nmax - *my_nmin;

  return ndo;
}


/**********************************************************/
/** 
 * @brief      communicates the MC estimators between tasks
 *
 * @details
 * communicates the MC estimators between tasks relating to 
 * spectral models, heating and cooling and cell diagnostics like IP. 
 * In the case of some variables, the quantities are maxima and minima so the 
 * flag MPI_MAX or MPI_MIN is used in MPI_Reduce. For summed
 * quantities like heating we use MPI_SUM.
 *  
 * This routine should only do anything if the MPI_ON flag was present 
 * in compilation. It communicates all the information
 * required for the spectral model ionization scheme, and 
 * also heating and cooling quantities in cells.
 **********************************************************/

int
communicate_estimators_para ()
{
#ifdef MPI_ON                   // these routines should only be called anyway in parallel but we need these to compile

  int mpi_i, mpi_j;
  double *maxfreqhelper, *maxfreqhelper2;
  /*NSH 131213 the next line introduces new helper arrays for the max and min frequencies in bands */
  double *maxbandfreqhelper, *maxbandfreqhelper2, *minbandfreqhelper, *minbandfreqhelper2;
  double *redhelper, *redhelper2, *qdisk_helper, *qdisk_helper2;
  double *ion_helper, *ion_helper2;
  double *inner_ion_helper, *inner_ion_helper2;
  double *flux_helper, *flux_helper2;

//OLD  int nspec, size_of_commbuffer;
  int size_of_commbuffer;

  int *iredhelper, *iredhelper2, *iqdisk_helper, *iqdisk_helper2;
  // int size_of_helpers;
  int plasma_double_helpers, plasma_int_helpers;

  /* The size of the helper array for doubles. We transmit 10 numbers
     for each cell, plus three arrays, each of length NXBANDS */

  plasma_double_helpers = (39 + 3 * NXBANDS) * NPLASMA;

  /* The size of the helper array for integers. We transmit 7 numbers
     for each cell, plus one array of length NXBANDS */
  plasma_int_helpers = (7 + NXBANDS) * NPLASMA;


  maxfreqhelper = calloc (sizeof (double), NPLASMA);
  maxfreqhelper2 = calloc (sizeof (double), NPLASMA);
  /* NSH 131213 - allocate memory for the band limited max and min frequencies */
  maxbandfreqhelper = calloc (sizeof (double), NPLASMA * NXBANDS);
  maxbandfreqhelper2 = calloc (sizeof (double), NPLASMA * NXBANDS);
  minbandfreqhelper = calloc (sizeof (double), NPLASMA * NXBANDS);
  minbandfreqhelper2 = calloc (sizeof (double), NPLASMA * NXBANDS);
  redhelper = calloc (sizeof (double), plasma_double_helpers);
  redhelper2 = calloc (sizeof (double), plasma_double_helpers);

  ion_helper = calloc (sizeof (double), NPLASMA * nions);
  ion_helper2 = calloc (sizeof (double), NPLASMA * nions);
  inner_ion_helper = calloc (sizeof (double), NPLASMA * n_inner_tot);
  inner_ion_helper2 = calloc (sizeof (double), NPLASMA * n_inner_tot);
  /* JM -- added routine to average the qdisk quantities. The 2 is because
     we only have two doubles to worry about (heat and ave_freq) and
     two integers (nhit and nphot) */
  qdisk_helper = calloc (sizeof (double), NRINGS * 2);
  qdisk_helper2 = calloc (sizeof (double), NRINGS * 2);

  flux_helper = calloc (sizeof (double), NPLASMA * NFLUX_ANGLES * 3);
  flux_helper2 = calloc (sizeof (double), NPLASMA * NFLUX_ANGLES * 3);

  // the following blocks gather all the estimators to the zeroth (Master) thread


  for (mpi_i = 0; mpi_i < NPLASMA; mpi_i++)
  {
    maxfreqhelper[mpi_i] = plasmamain[mpi_i].max_freq;
    redhelper[mpi_i] = plasmamain[mpi_i].j / np_mpi_global;
    redhelper[mpi_i + NPLASMA] = plasmamain[mpi_i].ave_freq / np_mpi_global;
    redhelper[mpi_i + 2 * NPLASMA] = plasmamain[mpi_i].cool_tot / np_mpi_global;
    redhelper[mpi_i + 3 * NPLASMA] = plasmamain[mpi_i].heat_tot / np_mpi_global;
    redhelper[mpi_i + 4 * NPLASMA] = plasmamain[mpi_i].heat_lines / np_mpi_global;
    redhelper[mpi_i + 5 * NPLASMA] = plasmamain[mpi_i].heat_ff / np_mpi_global;
    redhelper[mpi_i + 6 * NPLASMA] = plasmamain[mpi_i].heat_comp / np_mpi_global;
    redhelper[mpi_i + 7 * NPLASMA] = plasmamain[mpi_i].heat_ind_comp / np_mpi_global;
    redhelper[mpi_i + 8 * NPLASMA] = plasmamain[mpi_i].heat_photo / np_mpi_global;
    redhelper[mpi_i + 9 * NPLASMA] = plasmamain[mpi_i].ip / np_mpi_global;
    redhelper[mpi_i + 10 * NPLASMA] = plasmamain[mpi_i].j_direct / np_mpi_global;
    redhelper[mpi_i + 11 * NPLASMA] = plasmamain[mpi_i].j_scatt / np_mpi_global;
    redhelper[mpi_i + 12 * NPLASMA] = plasmamain[mpi_i].ip_direct / np_mpi_global;
    redhelper[mpi_i + 13 * NPLASMA] = plasmamain[mpi_i].ip_scatt / np_mpi_global;
    redhelper[mpi_i + 14 * NPLASMA] = plasmamain[mpi_i].heat_auger / np_mpi_global;
    redhelper[mpi_i + 15 * NPLASMA] = plasmamain[mpi_i].rad_force_es[0] / np_mpi_global;
    redhelper[mpi_i + 16 * NPLASMA] = plasmamain[mpi_i].rad_force_es[1] / np_mpi_global;
    redhelper[mpi_i + 17 * NPLASMA] = plasmamain[mpi_i].rad_force_es[2] / np_mpi_global;
    redhelper[mpi_i + 18 * NPLASMA] = plasmamain[mpi_i].rad_force_es[3] / np_mpi_global;
    redhelper[mpi_i + 19 * NPLASMA] = plasmamain[mpi_i].F_vis[0] / np_mpi_global;
    redhelper[mpi_i + 20 * NPLASMA] = plasmamain[mpi_i].F_vis[1] / np_mpi_global;
    redhelper[mpi_i + 21 * NPLASMA] = plasmamain[mpi_i].F_vis[2] / np_mpi_global;
    redhelper[mpi_i + 22 * NPLASMA] = plasmamain[mpi_i].F_vis[3] / np_mpi_global;
    redhelper[mpi_i + 23 * NPLASMA] = plasmamain[mpi_i].F_UV[0] / np_mpi_global;
    redhelper[mpi_i + 24 * NPLASMA] = plasmamain[mpi_i].F_UV[1] / np_mpi_global;
    redhelper[mpi_i + 25 * NPLASMA] = plasmamain[mpi_i].F_UV[2] / np_mpi_global;
    redhelper[mpi_i + 26 * NPLASMA] = plasmamain[mpi_i].F_UV[3] / np_mpi_global;
    redhelper[mpi_i + 27 * NPLASMA] = plasmamain[mpi_i].F_Xray[0] / np_mpi_global;
    redhelper[mpi_i + 28 * NPLASMA] = plasmamain[mpi_i].F_Xray[1] / np_mpi_global;
    redhelper[mpi_i + 29 * NPLASMA] = plasmamain[mpi_i].F_Xray[2] / np_mpi_global;
    redhelper[mpi_i + 30 * NPLASMA] = plasmamain[mpi_i].F_Xray[3] / np_mpi_global;
    redhelper[mpi_i + 31 * NPLASMA] = plasmamain[mpi_i].rad_force_bf[0] / np_mpi_global;
    redhelper[mpi_i + 32 * NPLASMA] = plasmamain[mpi_i].rad_force_bf[1] / np_mpi_global;
    redhelper[mpi_i + 33 * NPLASMA] = plasmamain[mpi_i].rad_force_bf[2] / np_mpi_global;
    redhelper[mpi_i + 34 * NPLASMA] = plasmamain[mpi_i].rad_force_bf[3] / np_mpi_global;
    redhelper[mpi_i + 35 * NPLASMA] = plasmamain[mpi_i].rad_force_ff[0] / np_mpi_global;
    redhelper[mpi_i + 36 * NPLASMA] = plasmamain[mpi_i].rad_force_ff[1] / np_mpi_global;
    redhelper[mpi_i + 37 * NPLASMA] = plasmamain[mpi_i].rad_force_ff[2] / np_mpi_global;
    redhelper[mpi_i + 38 * NPLASMA] = plasmamain[mpi_i].rad_force_ff[3] / np_mpi_global;

    for (mpi_j = 0; mpi_j < NXBANDS; mpi_j++)
    {
      redhelper[mpi_i + (39 + mpi_j) * NPLASMA] = plasmamain[mpi_i].xj[mpi_j] / np_mpi_global;
      redhelper[mpi_i + (39 + NXBANDS + mpi_j) * NPLASMA] = plasmamain[mpi_i].xave_freq[mpi_j] / np_mpi_global;
      redhelper[mpi_i + (39 + 2 * NXBANDS + mpi_j) * NPLASMA] = plasmamain[mpi_i].xsd_freq[mpi_j] / np_mpi_global;
      /* 131213 NSH populate the band limited min and max frequency arrays */
      maxbandfreqhelper[mpi_i * NXBANDS + mpi_j] = plasmamain[mpi_i].fmax[mpi_j];
      minbandfreqhelper[mpi_i * NXBANDS + mpi_j] = plasmamain[mpi_i].fmin[mpi_j];
    }
    for (mpi_j = 0; mpi_j < nions; mpi_j++)
    {
      ion_helper[mpi_i * nions + mpi_j] = plasmamain[mpi_i].ioniz[mpi_j] / np_mpi_global;
    }
    for (mpi_j = 0; mpi_j < n_inner_tot; mpi_j++)
    {
      inner_ion_helper[mpi_i * n_inner_tot + mpi_j] = plasmamain[mpi_i].inner_ioniz[mpi_j] / np_mpi_global;
    }
    for (mpi_j = 0; mpi_j < NFLUX_ANGLES; mpi_j++)
    {
      flux_helper[mpi_i * (3 * NFLUX_ANGLES) + mpi_j] = plasmamain[mpi_i].F_UV_ang_x[mpi_j] / np_mpi_global;
      flux_helper[mpi_i * (3 * NFLUX_ANGLES) + NFLUX_ANGLES + mpi_j] = plasmamain[mpi_i].F_UV_ang_y[mpi_j] / np_mpi_global;
      flux_helper[mpi_i * (3 * NFLUX_ANGLES) + 2 * NFLUX_ANGLES + mpi_j] = plasmamain[mpi_i].F_UV_ang_z[mpi_j] / np_mpi_global;
    }
  }

  for (mpi_i = 0; mpi_i < NRINGS; mpi_i++)
  {
    qdisk_helper[mpi_i] = qdisk.heat[mpi_i] / np_mpi_global;
    qdisk_helper[mpi_i + NRINGS] = qdisk.ave_freq[mpi_i] / np_mpi_global;
  }


  /* 131213 NSH communicate the min and max band frequencies these use MPI_MIN or MPI_MAX */
  MPI_Reduce (minbandfreqhelper, minbandfreqhelper2, NPLASMA * NXBANDS, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce (maxbandfreqhelper, maxbandfreqhelper2, NPLASMA * NXBANDS, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce (maxfreqhelper, maxfreqhelper2, NPLASMA, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce (redhelper, redhelper2, plasma_double_helpers, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (flux_helper, flux_helper2, NPLASMA * 3 * NFLUX_ANGLES, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Reduce (ion_helper, ion_helper2, NPLASMA * nions, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (inner_ion_helper, inner_ion_helper2, NPLASMA * n_inner_tot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  /* JM 1607 -- sum up the qdisk values */
  MPI_Reduce (qdisk_helper, qdisk_helper2, 2 * NRINGS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Bcast (redhelper2, plasma_double_helpers, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (maxfreqhelper2, NPLASMA, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  /* 131213 NSH Send out the global min and max band limited frequencies to all threads */
  MPI_Bcast (minbandfreqhelper2, NPLASMA * NXBANDS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (maxbandfreqhelper2, NPLASMA * NXBANDS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (flux_helper2, NPLASMA * 3 * NFLUX_ANGLES, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (ion_helper2, NPLASMA * nions, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (inner_ion_helper2, NPLASMA * n_inner_tot, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  /* JM 1607 -- send out the qdisk values to all threads */
  MPI_Bcast (qdisk_helper2, NRINGS, MPI_DOUBLE, 0, MPI_COMM_WORLD);


  for (mpi_i = 0; mpi_i < NPLASMA; mpi_i++)
  {
    plasmamain[mpi_i].max_freq = maxfreqhelper2[mpi_i];
    plasmamain[mpi_i].j = redhelper2[mpi_i];
    plasmamain[mpi_i].ave_freq = redhelper2[mpi_i + NPLASMA];
    plasmamain[mpi_i].cool_tot = redhelper2[mpi_i + 2 * NPLASMA];
    plasmamain[mpi_i].heat_tot = redhelper2[mpi_i + 3 * NPLASMA];
    plasmamain[mpi_i].heat_lines = redhelper2[mpi_i + 4 * NPLASMA];
    plasmamain[mpi_i].heat_ff = redhelper2[mpi_i + 5 * NPLASMA];
    plasmamain[mpi_i].heat_comp = redhelper2[mpi_i + 6 * NPLASMA];
    plasmamain[mpi_i].heat_ind_comp = redhelper2[mpi_i + 7 * NPLASMA];
    plasmamain[mpi_i].heat_photo = redhelper2[mpi_i + 8 * NPLASMA];
    plasmamain[mpi_i].ip = redhelper2[mpi_i + 9 * NPLASMA];
    plasmamain[mpi_i].j_direct = redhelper2[mpi_i + 10 * NPLASMA];
    plasmamain[mpi_i].j_scatt = redhelper2[mpi_i + 11 * NPLASMA];
    plasmamain[mpi_i].ip_direct = redhelper2[mpi_i + 12 * NPLASMA];
    plasmamain[mpi_i].ip_scatt = redhelper2[mpi_i + 13 * NPLASMA];
    plasmamain[mpi_i].heat_auger = redhelper2[mpi_i + 14 * NPLASMA];
    plasmamain[mpi_i].rad_force_es[0] = redhelper2[mpi_i + 15 * NPLASMA];
    plasmamain[mpi_i].rad_force_es[1] = redhelper2[mpi_i + 16 * NPLASMA];
    plasmamain[mpi_i].rad_force_es[2] = redhelper2[mpi_i + 17 * NPLASMA];
    plasmamain[mpi_i].rad_force_es[3] = redhelper2[mpi_i + 18 * NPLASMA];
    plasmamain[mpi_i].F_vis[0] = redhelper2[mpi_i + 19 * NPLASMA];
    plasmamain[mpi_i].F_vis[1] = redhelper2[mpi_i + 20 * NPLASMA];
    plasmamain[mpi_i].F_vis[2] = redhelper2[mpi_i + 21 * NPLASMA];
    plasmamain[mpi_i].F_vis[3] = redhelper2[mpi_i + 22 * NPLASMA];
    plasmamain[mpi_i].F_UV[0] = redhelper2[mpi_i + 23 * NPLASMA];
    plasmamain[mpi_i].F_UV[1] = redhelper2[mpi_i + 24 * NPLASMA];
    plasmamain[mpi_i].F_UV[2] = redhelper2[mpi_i + 25 * NPLASMA];
    plasmamain[mpi_i].F_UV[3] = redhelper2[mpi_i + 26 * NPLASMA];
    plasmamain[mpi_i].F_Xray[0] = redhelper2[mpi_i + 27 * NPLASMA];
    plasmamain[mpi_i].F_Xray[1] = redhelper2[mpi_i + 28 * NPLASMA];
    plasmamain[mpi_i].F_Xray[2] = redhelper2[mpi_i + 29 * NPLASMA];
    plasmamain[mpi_i].F_Xray[3] = redhelper2[mpi_i + 30 * NPLASMA];
    plasmamain[mpi_i].rad_force_bf[0] = redhelper2[mpi_i + 31 * NPLASMA];
    plasmamain[mpi_i].rad_force_bf[1] = redhelper2[mpi_i + 32 * NPLASMA];
    plasmamain[mpi_i].rad_force_bf[2] = redhelper2[mpi_i + 33 * NPLASMA];
    plasmamain[mpi_i].rad_force_bf[3] = redhelper2[mpi_i + 34 * NPLASMA];
    plasmamain[mpi_i].rad_force_ff[0] = redhelper2[mpi_i + 35 * NPLASMA];
    plasmamain[mpi_i].rad_force_ff[1] = redhelper2[mpi_i + 36 * NPLASMA];
    plasmamain[mpi_i].rad_force_ff[2] = redhelper2[mpi_i + 37 * NPLASMA];
    plasmamain[mpi_i].rad_force_ff[3] = redhelper2[mpi_i + 38 * NPLASMA];

    for (mpi_j = 0; mpi_j < NXBANDS; mpi_j++)
    {
      plasmamain[mpi_i].xj[mpi_j] = redhelper2[mpi_i + (39 + mpi_j) * NPLASMA];
      plasmamain[mpi_i].xave_freq[mpi_j] = redhelper2[mpi_i + (39 + NXBANDS + mpi_j) * NPLASMA];
      plasmamain[mpi_i].xsd_freq[mpi_j] = redhelper2[mpi_i + (39 + NXBANDS * 2 + mpi_j) * NPLASMA];

      /* 131213 NSH And unpack the min and max banded frequencies to the plasma array */
      plasmamain[mpi_i].fmax[mpi_j] = maxbandfreqhelper2[mpi_i * NXBANDS + mpi_j];
      plasmamain[mpi_i].fmin[mpi_j] = minbandfreqhelper2[mpi_i * NXBANDS + mpi_j];
    }
    for (mpi_j = 0; mpi_j < nions; mpi_j++)
    {
      plasmamain[mpi_i].ioniz[mpi_j] = ion_helper2[mpi_i * nions + mpi_j];
    }
    for (mpi_j = 0; mpi_j < n_inner_tot; mpi_j++)
    {
      plasmamain[mpi_i].inner_ioniz[mpi_j] = inner_ion_helper2[mpi_i * n_inner_tot + mpi_j];
    }
    for (mpi_j = 0; mpi_j < NFLUX_ANGLES; mpi_j++)
    {
      plasmamain[mpi_i].F_UV_ang_x[mpi_j] = flux_helper2[mpi_i * (3 * NFLUX_ANGLES) + mpi_j];
      plasmamain[mpi_i].F_UV_ang_y[mpi_j] = flux_helper2[mpi_i * (3 * NFLUX_ANGLES) + NFLUX_ANGLES + mpi_j];
      plasmamain[mpi_i].F_UV_ang_z[mpi_j] = flux_helper2[mpi_i * (3 * NFLUX_ANGLES) + 2 * NFLUX_ANGLES + mpi_j];
    }
  }


  for (mpi_i = 0; mpi_i < NRINGS; mpi_i++)
  {
    qdisk.heat[mpi_i] = qdisk_helper2[mpi_i];
    qdisk.ave_freq[mpi_i] = qdisk_helper2[mpi_i + NRINGS];
  }

  /* now we've done all the doubles so we can free their helper arrays */
  free (qdisk_helper);
  free (qdisk_helper2);
  free (redhelper);
  free (redhelper2);
  free (maxfreqhelper);
  free (maxfreqhelper2);
  free (maxbandfreqhelper);
  free (maxbandfreqhelper2);
  free (minbandfreqhelper);
  free (minbandfreqhelper2);

  free (ion_helper);
  free (ion_helper2);

  free (inner_ion_helper);
  free (inner_ion_helper2);
  /* allocate the integer helper arrays, set a barrier, then do all the integers. */
  iqdisk_helper = calloc (sizeof (int), NRINGS * 2);
  iqdisk_helper2 = calloc (sizeof (int), NRINGS * 2);
  iredhelper = calloc (sizeof (int), plasma_int_helpers);
  iredhelper2 = calloc (sizeof (int), plasma_int_helpers);

  for (mpi_i = 0; mpi_i < NPLASMA; mpi_i++)
  {
    iredhelper[mpi_i] = plasmamain[mpi_i].ntot;
    iredhelper[mpi_i + NPLASMA] = plasmamain[mpi_i].ntot_star;
    iredhelper[mpi_i + 2 * NPLASMA] = plasmamain[mpi_i].ntot_bl;
    iredhelper[mpi_i + 3 * NPLASMA] = plasmamain[mpi_i].ntot_disk;
    iredhelper[mpi_i + 4 * NPLASMA] = plasmamain[mpi_i].ntot_wind;
    iredhelper[mpi_i + 5 * NPLASMA] = plasmamain[mpi_i].ntot_agn;
    iredhelper[mpi_i + 6 * NPLASMA] = plasmamain[mpi_i].nioniz;

    for (mpi_j = 0; mpi_j < NXBANDS; mpi_j++)
    {
      iredhelper[mpi_i + (7 + mpi_j) * NPLASMA] = plasmamain[mpi_i].nxtot[mpi_j];
    }
  }

  for (mpi_i = 0; mpi_i < NRINGS; mpi_i++)
  {
    iqdisk_helper[mpi_i] = qdisk.nphot[mpi_i];
    iqdisk_helper[mpi_i + NRINGS] = qdisk.nhit[mpi_i];
  }

  MPI_Reduce (iredhelper, iredhelper2, plasma_int_helpers, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (iqdisk_helper, iqdisk_helper2, 2 * NRINGS, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Bcast (iredhelper2, plasma_int_helpers, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (iqdisk_helper2, NRINGS, MPI_INT, 0, MPI_COMM_WORLD);

  for (mpi_i = 0; mpi_i < NPLASMA; mpi_i++)
  {
    plasmamain[mpi_i].ntot = iredhelper2[mpi_i];
    plasmamain[mpi_i].ntot_star = iredhelper2[mpi_i + NPLASMA];
    plasmamain[mpi_i].ntot_bl = iredhelper2[mpi_i + 2 * NPLASMA];
    plasmamain[mpi_i].ntot_disk = iredhelper2[mpi_i + 3 * NPLASMA];
    plasmamain[mpi_i].ntot_wind = iredhelper2[mpi_i + 4 * NPLASMA];
    plasmamain[mpi_i].ntot_agn = iredhelper2[mpi_i + 5 * NPLASMA];
    plasmamain[mpi_i].nioniz = iredhelper2[mpi_i + 6 * NPLASMA];

    for (mpi_j = 0; mpi_j < NXBANDS; mpi_j++)
    {
      plasmamain[mpi_i].nxtot[mpi_j] = iredhelper2[mpi_i + (7 + mpi_j) * NPLASMA];
    }
  }

  for (mpi_i = 0; mpi_i < NRINGS; mpi_i++)
  {
    qdisk.nphot[mpi_i] = iqdisk_helper2[mpi_i];
    qdisk.nhit[mpi_i] = iqdisk_helper2[mpi_i + NRINGS];
  }

  free (iredhelper);
  free (iredhelper2);
  free (iqdisk_helper);
  free (iqdisk_helper2);


/* Now during ionization cycles, process the cell spectra */

  /* The size of the commbuffers need to be the number of spectra x the length of each */

  if (geo.ioniz_or_extract == CYCLE_IONIZ)
  {
    size_of_commbuffer = NPLASMA * NBINS_IN_CELL_SPEC;
//OLD    nspec = NPLASMA;

    redhelper = calloc (sizeof (double), size_of_commbuffer);
    redhelper2 = calloc (sizeof (double), size_of_commbuffer);

    for (mpi_i = 0; mpi_i < NBINS_IN_CELL_SPEC; mpi_i++)
    {
      for (mpi_j = 0; mpi_j < NPLASMA; mpi_j++)
      {
        redhelper[mpi_i * NPLASMA + mpi_j] = plasmamain[mpi_j].cell_spec_flux[mpi_i] / np_mpi_global;

      }
    }

    MPI_Reduce (redhelper, redhelper2, size_of_commbuffer, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast (redhelper2, size_of_commbuffer, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (mpi_i = 0; mpi_i < NBINS_IN_CELL_SPEC; mpi_i++)
    {
      for (mpi_j = 0; mpi_j < NPLASMA; mpi_j++)
      {
        plasmamain[mpi_j].cell_spec_flux[mpi_i] = redhelper2[mpi_i * NPLASMA + mpi_j];

      }
    }

    free (redhelper);
    free (redhelper2);
  }


#endif
  return (0);
}


/**********************************************************/
/**
 * @brief sum up the synthetic and cell spectra between threads.
 *
 * @details
 * sum up the synthetic spectra between threads. Does an
 * MPI_Reduce then an MPI_Bcast for each element of the
 * linear and log spectra arrays (xxspec)
 *
 **********************************************************/


int
gather_spectra_para ()
{
#ifdef MPI_ON
  double *redhelper, *redhelper2;
  int mpi_i, mpi_j;
  int size_of_commbuffer, nspec;

  if (geo.ioniz_or_extract == CYCLE_EXTRACT)
  {
    nspec = MSPEC + geo.nangles;
  }
  else
  {
    nspec = MSPEC;
  }

  size_of_commbuffer = 4 * nspec * NWAVE_MAX;   //we need space for all 4 separate spectra we are normalizing

  redhelper = calloc (sizeof (double), size_of_commbuffer);
  redhelper2 = calloc (sizeof (double), size_of_commbuffer);

  for (mpi_i = 0; mpi_i < NWAVE_MAX; mpi_i++)
  {
    for (mpi_j = 0; mpi_j < nspec; mpi_j++)
    {
      redhelper[mpi_i * nspec + mpi_j] = xxspec[mpi_j].f[mpi_i] / np_mpi_global;
      redhelper[mpi_i * nspec + mpi_j + (NWAVE_MAX * nspec)] = xxspec[mpi_j].lf[mpi_i] / np_mpi_global;
      redhelper[mpi_i * nspec + mpi_j + (2 * NWAVE_MAX * nspec)] = xxspec[mpi_j].f_wind[mpi_i] / np_mpi_global;
      redhelper[mpi_i * nspec + mpi_j + (3 * NWAVE_MAX * nspec)] = xxspec[mpi_j].lf_wind[mpi_i] / np_mpi_global;
    }
  }

  MPI_Reduce (redhelper, redhelper2, size_of_commbuffer, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Bcast (redhelper2, size_of_commbuffer, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  for (mpi_i = 0; mpi_i < NWAVE_MAX; mpi_i++)
  {
    for (mpi_j = 0; mpi_j < nspec; mpi_j++)
    {
      xxspec[mpi_j].f[mpi_i] = redhelper2[mpi_i * nspec + mpi_j];
      xxspec[mpi_j].lf[mpi_i] = redhelper2[mpi_i * nspec + mpi_j + (NWAVE_MAX * nspec)];
      xxspec[mpi_j].f_wind[mpi_i] = redhelper2[mpi_i * nspec + mpi_j + (2 * NWAVE_MAX * nspec)];
      xxspec[mpi_j].lf_wind[mpi_i] = redhelper2[mpi_i * nspec + mpi_j + (3 * NWAVE_MAX * nspec)];
    }
  }

  free (redhelper);
  free (redhelper2);
#endif

  return (0);
}





/**********************************************************/
/**
 * @brief
 *
 * @details averages the macro-atom estimators between tasks using MPI_Reduce.
 *   It should only be called if the MPI_ON flag was present
 *   in compilation, and returns 0 immediately if no macro atom levels.
 *   This should probably be improved by working out exactly
 *   what is needed in simple-ion only mode.
 *
 **********************************************************/

int
communicate_matom_estimators_para ()
{
#ifdef MPI_ON                   // these routines should only be called anyway in parallel but we need these to compile

  int n, mpi_i;
  double *gamma_helper, *alpha_helper;
  double *level_helper, *cell_helper, *jbar_helper;
  double *gamma_helper2, *alpha_helper2;
  double *level_helper2, *cell_helper2, *jbar_helper2;
  double *cooling_bf_helper, *cooling_bb_helper;
  double *cooling_bf_helper2, *cooling_bb_helper2;

  if (nlevels_macro == 0 && geo.nmacro == 0)
  {
    /* in this case no space would have been allocated for macro-atom estimators */
    Log ("No need to communicate matom estimators as no macro-atoms!\n");
    return (0);
  }



  /* allocate helper arrays for the estimators we want to communicate */
  /* the sizes of these arrays should match the allocation in calloc_estimators in gridwind.c */
  /* we need two arrays for each set of variables. Note that we stick all estimators of
     the same size in the same helper array */
  /* could put some error conditions here to check memory allocation worked */
  jbar_helper = calloc (sizeof (double), NPLASMA * size_Jbar_est);
  gamma_helper = calloc (sizeof (double), NPLASMA * 4 * size_gamma_est);
  alpha_helper = calloc (sizeof (double), NPLASMA * 2 * size_alpha_est);
  level_helper = calloc (sizeof (double), NPLASMA * nlevels_macro);
  cell_helper = calloc (sizeof (double), 8 * NPLASMA);
  cooling_bf_helper = calloc (sizeof (double), NPLASMA * 2 * nphot_total);
  cooling_bb_helper = calloc (sizeof (double), NPLASMA * nlines);

  jbar_helper2 = calloc (sizeof (double), NPLASMA * size_Jbar_est);
  gamma_helper2 = calloc (sizeof (double), NPLASMA * 4 * size_gamma_est);
  alpha_helper2 = calloc (sizeof (double), NPLASMA * 2 * size_alpha_est);
  level_helper2 = calloc (sizeof (double), NPLASMA * nlevels_macro);
  cell_helper2 = calloc (sizeof (double), 8 * NPLASMA);
  cooling_bf_helper2 = calloc (sizeof (double), NPLASMA * 2 * nphot_total);
  cooling_bb_helper2 = calloc (sizeof (double), NPLASMA * nlines);


  /* now we loop through each cell and copy the values of our variables
     into our helper arrays */
  for (mpi_i = 0; mpi_i < NPLASMA; mpi_i++)
  {
    /* one kpkt_abs quantity per cell */
    cell_helper[mpi_i] = plasmamain[mpi_i].kpkt_abs / np_mpi_global;

    /* each of the cooling sums and normalisations also have one quantity per cell */
    cell_helper[mpi_i + NPLASMA] = macromain[mpi_i].cooling_normalisation / np_mpi_global;
    cell_helper[mpi_i + 2 * NPLASMA] = macromain[mpi_i].cooling_bftot / np_mpi_global;
    cell_helper[mpi_i + 3 * NPLASMA] = macromain[mpi_i].cooling_bf_coltot / np_mpi_global;
    cell_helper[mpi_i + 4 * NPLASMA] = macromain[mpi_i].cooling_bbtot / np_mpi_global;
    cell_helper[mpi_i + 5 * NPLASMA] = macromain[mpi_i].cooling_ff / np_mpi_global;
    cell_helper[mpi_i + 6 * NPLASMA] = macromain[mpi_i].cooling_ff_lofreq / np_mpi_global;
    cell_helper[mpi_i + 7 * NPLASMA] = macromain[mpi_i].cooling_adiabatic / np_mpi_global;



    for (n = 0; n < nlevels_macro; n++)
    {
      level_helper[mpi_i + (n * NPLASMA)] = macromain[mpi_i].matom_abs[n] / np_mpi_global;
    }

    for (n = 0; n < size_Jbar_est; n++)
    {
      jbar_helper[mpi_i + (n * NPLASMA)] = macromain[mpi_i].jbar[n] / np_mpi_global;
    }

    for (n = 0; n < size_gamma_est; n++)
    {
      gamma_helper[mpi_i + (n * NPLASMA)] = macromain[mpi_i].alpha_st[n] / np_mpi_global;
      gamma_helper[mpi_i + ((n + size_gamma_est) * NPLASMA)] = macromain[mpi_i].alpha_st_e[n] / np_mpi_global;
      gamma_helper[mpi_i + ((n + 2 * size_gamma_est) * NPLASMA)] = macromain[mpi_i].gamma[n] / np_mpi_global;
      gamma_helper[mpi_i + ((n + 3 * size_gamma_est) * NPLASMA)] = macromain[mpi_i].gamma_e[n] / np_mpi_global;
    }

    for (n = 0; n < size_alpha_est; n++)
    {
      alpha_helper[mpi_i + (n * NPLASMA)] = macromain[mpi_i].recomb_sp[n] / np_mpi_global;
      alpha_helper[mpi_i + ((n + size_alpha_est) * NPLASMA)] = macromain[mpi_i].recomb_sp_e[n] / np_mpi_global;
    }

    for (n = 0; n < nphot_total; n++)
    {
      cooling_bf_helper[mpi_i + (n * NPLASMA)] = macromain[mpi_i].cooling_bf[n] / np_mpi_global;
      cooling_bf_helper[mpi_i + ((n + nphot_total) * NPLASMA)] = macromain[mpi_i].cooling_bf_col[n] / np_mpi_global;
    }

    for (n = 0; n < nlines; n++)
    {
      cooling_bb_helper[mpi_i + (n * NPLASMA)] = macromain[mpi_i].cooling_bb[n] / np_mpi_global;
    }
  }

  /* because in the above loop we have already divided by number of processes, we can now do a sum
     with MPI_Reduce, passing it MPI_SUM as an argument. This will give us the mean across threads */
  MPI_Reduce (cell_helper, cell_helper2, NPLASMA * 8, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (level_helper, level_helper2, NPLASMA * nlevels_macro, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (jbar_helper, jbar_helper2, NPLASMA * size_Jbar_est, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (gamma_helper, gamma_helper2, NPLASMA * 4 * size_gamma_est, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (alpha_helper, alpha_helper2, NPLASMA * 2 * size_alpha_est, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (cooling_bf_helper, cooling_bf_helper2, NPLASMA * 2 * nphot_total, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (cooling_bb_helper, cooling_bb_helper2, NPLASMA * nlines, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Bcast (cell_helper2, NPLASMA * 8, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (level_helper2, NPLASMA * nlevels_macro, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (jbar_helper2, NPLASMA * size_Jbar_est, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (gamma_helper2, NPLASMA * 4 * size_gamma_est, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (alpha_helper2, NPLASMA * 2 * size_alpha_est, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (cooling_bf_helper2, NPLASMA * 2 * nphot_total, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (cooling_bb_helper2, NPLASMA * nlines, MPI_DOUBLE, 0, MPI_COMM_WORLD);




  /* We now need to copy these reduced variables to the plasma structure in each thread */


  for (mpi_i = 0; mpi_i < NPLASMA; mpi_i++)
  {
    /* one kpkt_abs quantity per cell */
    plasmamain[mpi_i].kpkt_abs = cell_helper2[mpi_i];

    /* each of the cooling sums and normalisations also have one quantity per cell */
    macromain[mpi_i].cooling_normalisation = cell_helper2[mpi_i + NPLASMA];
    macromain[mpi_i].cooling_bftot = cell_helper2[mpi_i + 2 * NPLASMA];
    macromain[mpi_i].cooling_bf_coltot = cell_helper2[mpi_i + 3 * NPLASMA];
    macromain[mpi_i].cooling_bbtot = cell_helper2[mpi_i + 4 * NPLASMA];
    macromain[mpi_i].cooling_ff = cell_helper2[mpi_i + 5 * NPLASMA];
    macromain[mpi_i].cooling_ff_lofreq = cell_helper2[mpi_i + 6 * NPLASMA];
    macromain[mpi_i].cooling_adiabatic = cell_helper2[mpi_i + 7 * NPLASMA];


    for (n = 0; n < nlevels_macro; n++)
    {
      macromain[mpi_i].matom_abs[n] = level_helper2[mpi_i + (n * NPLASMA)];
    }

    for (n = 0; n < size_Jbar_est; n++)
    {
      macromain[mpi_i].jbar[n] = jbar_helper2[mpi_i + (n * NPLASMA)];
    }

    for (n = 0; n < size_gamma_est; n++)
    {
      macromain[mpi_i].alpha_st[n] = gamma_helper2[mpi_i + (n * NPLASMA)];
      macromain[mpi_i].alpha_st_e[n] = gamma_helper2[mpi_i + ((n + size_gamma_est) * NPLASMA)] / np_mpi_global;
      macromain[mpi_i].gamma[n] = gamma_helper2[mpi_i + ((n + 2 * size_gamma_est) * NPLASMA)];
      macromain[mpi_i].gamma_e[n] = gamma_helper2[mpi_i + ((n + 3 * size_gamma_est) * NPLASMA)];
    }

    for (n = 0; n < size_alpha_est; n++)
    {
      macromain[mpi_i].recomb_sp[n] = alpha_helper2[mpi_i + (n * NPLASMA)];
      macromain[mpi_i].recomb_sp_e[n] = alpha_helper2[mpi_i + ((n + size_alpha_est) * NPLASMA)];
    }

    for (n = 0; n < nphot_total; n++)
    {
      macromain[mpi_i].cooling_bf[n] = cooling_bf_helper2[mpi_i + (n * NPLASMA)];
      macromain[mpi_i].cooling_bf_col[n] = cooling_bf_helper2[mpi_i + ((n + nphot_total) * NPLASMA)];
    }

    for (n = 0; n < nlines; n++)
    {
      macromain[mpi_i].cooling_bb[n] = cooling_bb_helper2[mpi_i + (n * NPLASMA)];
    }
  }

  free (cell_helper);
  free (level_helper);
  free (jbar_helper);
  free (gamma_helper);
  free (alpha_helper);
  free (cooling_bf_helper);
  free (cooling_bb_helper);

  free (cell_helper2);
  free (level_helper2);
  free (jbar_helper2);
  free (gamma_helper2);
  free (alpha_helper2);
  free (cooling_bf_helper2);
  free (cooling_bb_helper2);
#endif


  return (0);
}

/**********************************************************/
/**
 * @brief communicates the macro-atom B matrices between threads
 *
 *
 * @details communicates the macro-atom B matrices between threads 
 * using MPI_Reduce. Only does anything if MPI_ON flag is on, 
 * and should only be called if geo.rt_mode == RT_MODE_MACRO 
 * and nlevels_macro > 0
 *
 **********************************************************/

int
communicate_matom_matrices ()
{
#ifdef MPI_ON
  int size_of_commbuffer, nrows, n_mpi, n_mpi2, num_comm;
  int my_nmax, my_nmin, ndo, n, position, i;
  char *commbuffer;
  ndo = get_parallel_nrange (rank_global, NPLASMA, np_mpi_global, &my_nmin, &my_nmax);

  nrows = nlevels_macro + 1;
  size_of_commbuffer = 8.0 * ((nrows * nrows) + 2) * (floor ((double) NPLASMA / np_mpi_global) + 1);
  commbuffer = (char *) malloc (size_of_commbuffer * sizeof (char));

  for (n_mpi = 0; n_mpi < np_mpi_global; n_mpi++)
  {
    position = 0;

    if (rank_global == n_mpi)
    {
      MPI_Pack (&ndo, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
      for (n = my_nmin; n < my_nmax; n++)
      {
        MPI_Pack (&n, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);

        /* we only communicate the matrix if it is being stored in this cell */
        if (macromain[n].store_matom_matrix == TRUE)
        {
          for (i = 0; i < nrows; i++)
          {
            MPI_Pack (macromain[n].matom_matrix[i], nrows, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
          }
        }
      }
    }
    MPI_Bcast (commbuffer, size_of_commbuffer, MPI_PACKED, n_mpi, MPI_COMM_WORLD);

    position = 0;

    if (rank_global != n_mpi)
    {
      MPI_Unpack (commbuffer, size_of_commbuffer, &position, &num_comm, 1, MPI_INT, MPI_COMM_WORLD);
      for (n_mpi2 = 0; n_mpi2 < num_comm; n_mpi2++)
      {
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &n, 1, MPI_INT, MPI_COMM_WORLD);

        /* we only communicate the matrix if it is being stored in this cell */
        if (macromain[n].store_matom_matrix == TRUE)
        {
          for (i = 0; i < nrows; i++)
          {
            MPI_Unpack (commbuffer, size_of_commbuffer, &position, macromain[n].matom_matrix[i], nrows, MPI_DOUBLE, MPI_COMM_WORLD);
          }
        }
      }
    }
  }

  free (commbuffer);
#endif
  return (0);
}

/**********************************************************/
/**
 * @brief
 *
 * @details
 *
 * ### Notes ###
 *
 **********************************************************/

int
communicate_plasma_cells (int n_start_rank, int n_stop_rank)
{
#ifdef MPI_ON
  int n;
  int position, n_mpi, num_comm, n_mpi2;
  int size_of_commbuffer, size_of_specbuffer;
  char *commbuffer;

  int ndo = n_stop_rank - n_start_rank;

  /* JM 1409 -- Added for issue #110 to ensure correct reporting in parallel */
//  int nmax_r_temp, nmax_e_temp;
//  double dt_e_temp, dt_r_temp;

  /* The commbuffer needs to be larger enough to pack all variables in MPI_Pack and MPI_Unpack routines
   * The cmombuffer is currently sized to be the minimum requred.  Therefor when variables are added, the
   * size must be increased.
   *
   * The cell spectra are handled separately
   */

  int n_cells_max = ceil ((double) NPLASMA / np_mpi_global);

  size_of_commbuffer = 8 * (n_inner_tot + 10 * nions + nlte_levels + 3 * nphot_total + 15 * NXBANDS + 133 * NFLUX_ANGLES) * n_cells_max;
  size_of_commbuffer += 2 * 4 * n_cells_max * N_PHOT_PROC;      /* Add space for storing number of bf transitions, up and down */
  size_of_specbuffer = 8 * n_cells_max * NBINS_IN_CELL_SPEC;
  size_of_commbuffer += size_of_specbuffer;
  commbuffer = malloc (size_of_commbuffer);

  Log ("Beginning communication for plasma cells\n");

  for (n_mpi = 0; n_mpi < np_mpi_global; n_mpi++)
  {
    position = 0;

    if (rank_global == n_mpi)
    {
      Log ("MPI task %d is working on cells %d to max %d (total size %d).\n", rank_global, n_start_rank, n_stop_rank, NPLASMA);
      MPI_Pack (&ndo, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
      for (n = n_start_rank; n < n_stop_rank; n++)
      {
        MPI_Pack (&n, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].nwind, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].nplasma, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ne, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].rho, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].vol, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].xgamma, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].density, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].partition, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].levden, nlte_levels, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].kappa_ff_factor, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].nscat_es, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].recomb_simple, nphot_total, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].recomb_simple_upweight, nphot_total, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].kpkt_emiss, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].kpkt_abs, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].kbf_use, nphot_total, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].kbf_nuse, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].t_r, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].t_r_old, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].t_e, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].t_e_old, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].dt_e, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].dt_e_old, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_tot, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].abs_tot, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_tot_old, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_lines, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_ff, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_comp, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_ind_comp, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_lines_macro, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_photo_macro, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_photo, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_auger, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].abs_photo, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].abs_auger, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_z, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].w, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ntot, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ntot_star, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ntot_bl, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ntot_disk, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ntot_wind, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ntot_agn, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].mean_ds, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].n_ds, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].nrad, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].nioniz, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].ioniz, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].inner_ioniz, n_inner_tot, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].recomb, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].scatters, nions, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].xscatters, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].heat_ion, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].heat_inner_ion, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].cool_rr_ion, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].lum_rr_ion, nions, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].j, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].j_direct, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].j_scatt, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ave_freq, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_tot, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].xj, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].xave_freq, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].xsd_freq, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].nxtot, NXBANDS, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].F_vis, 4, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].F_UV, 4, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].F_Xray, 4, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].max_freq, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_lines, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_ff, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_adiabatic, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].comp_nujnu, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_comp, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_dr, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_ch_ex, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_di, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_rr, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_rr, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_rr_metals, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_rr_metals, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_tot, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_tot_old, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_tot_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_lines_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_ff_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_adiabatic_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_comp_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_dr_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_di_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_rr_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_rr_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].cool_rr_metals_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].lum_tot_ioniz, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].heat_shock, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].dmo_dt, 3, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].rad_force_es, 4, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].rad_force_ff, 4, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].rad_force_bf, 4, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].gain, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].converge_t_r, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].converge_t_e, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].converge_hc, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].trcheck, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].techeck, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].hccheck, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].converge_whole, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].converging, 1, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].spec_mod_type, NXBANDS, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].pl_alpha, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].pl_log_w, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].exp_temp, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].exp_w, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].fmin_mod, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].fmax_mod, NXBANDS, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ip, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ip_direct, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].ip_scatt, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].xi, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].bf_simple_ionpool_in, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].bf_simple_ionpool_out, 1, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);

        MPI_Pack (&plasmamain[n].n_bf_in, N_PHOT_PROC, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n].n_bf_out, N_PHOT_PROC, MPI_INT, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);


        MPI_Pack (plasmamain[n].cell_spec_flux, NBINS_IN_CELL_SPEC, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].F_UV_ang_x, NFLUX_ANGLES, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].F_UV_ang_y, NFLUX_ANGLES, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
        MPI_Pack (plasmamain[n].F_UV_ang_z, NFLUX_ANGLES, MPI_DOUBLE, commbuffer, size_of_commbuffer, &position, MPI_COMM_WORLD);
      }
    }

    MPI_Bcast (commbuffer, size_of_commbuffer, MPI_PACKED, n_mpi, MPI_COMM_WORLD);

    position = 0;

    if (rank_global != n_mpi)
    {
      MPI_Unpack (commbuffer, size_of_commbuffer, &position, &num_comm, 1, MPI_INT, MPI_COMM_WORLD);
      for (n_mpi2 = 0; n_mpi2 < num_comm; n_mpi2++)
      {
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &n, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].nwind, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].nplasma, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ne, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].rho, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].vol, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].xgamma, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].density, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].partition, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].levden, nlte_levels, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].kappa_ff_factor, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].nscat_es, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].recomb_simple, nphot_total, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].recomb_simple_upweight, nphot_total, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].kpkt_emiss, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].kpkt_abs, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].kbf_use, nphot_total, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].kbf_nuse, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].t_r, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].t_r_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].t_e, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].t_e_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].dt_e, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].dt_e_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].abs_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_tot_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_lines, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_ff, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_comp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_ind_comp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_lines_macro, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_photo_macro, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_photo, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_auger, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].abs_photo, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].abs_auger, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_z, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].w, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ntot, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ntot_star, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ntot_bl, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ntot_disk, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ntot_wind, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ntot_agn, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].mean_ds, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].n_ds, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].nrad, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].nioniz, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].ioniz, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].inner_ioniz, n_inner_tot, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].recomb, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].scatters, nions, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].xscatters, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].heat_ion, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].heat_inner_ion, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].cool_rr_ion, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].lum_rr_ion, nions, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].j, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].j_direct, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].j_scatt, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ave_freq, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].xj, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].xave_freq, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].xsd_freq, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].nxtot, NXBANDS, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].F_vis, 4, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].F_UV, 4, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].F_Xray, 4, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].max_freq, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_lines, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_ff, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_adiabatic, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].comp_nujnu, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_comp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_dr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_ch_ex, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_di, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_rr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_rr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_rr_metals, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_rr_metals, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_tot_old, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_tot_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_lines_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_ff_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_adiabatic_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_comp_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_dr_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_di_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_rr_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_rr_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].cool_rr_metals_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].lum_tot_ioniz, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].heat_shock, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].dmo_dt, 3, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].rad_force_es, 4, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].rad_force_ff, 4, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].rad_force_bf, 4, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].gain, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].converge_t_r, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].converge_t_e, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].converge_hc, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].trcheck, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].techeck, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].hccheck, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].converge_whole, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].converging, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].spec_mod_type, NXBANDS, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].pl_alpha, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].pl_log_w, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].exp_temp, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].exp_w, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].fmin_mod, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].fmax_mod, NXBANDS, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ip, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ip_direct, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].ip_scatt, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].xi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].bf_simple_ionpool_in, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].bf_simple_ionpool_out, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].n_bf_in, N_PHOT_PROC, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, &plasmamain[n].n_bf_out, N_PHOT_PROC, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].cell_spec_flux, NBINS_IN_CELL_SPEC, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].F_UV_ang_x, NFLUX_ANGLES, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].F_UV_ang_y, NFLUX_ANGLES, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (commbuffer, size_of_commbuffer, &position, plasmamain[n].F_UV_ang_z, NFLUX_ANGLES, MPI_DOUBLE, MPI_COMM_WORLD);
      }
    }
  }

  free (commbuffer);
  Log ("Ending communication for plasma cells\n");
#endif
  return EXIT_SUCCESS;
}

/**********************************************************/
/**
 * @brief
 *
 * @details
 *
 * ### Notes ###
 *
 **********************************************************/

void
communicate_wind_luminosity (const int n_start, const int n_stop, const int n_cells_rank)
{
#ifdef MPI_ON
  int n_plasma;
  int size_doubles;
  int size_ints;
  int current_rank;

  const int n_cells_max = ceil ((double) NPLASMA / np_mpi_global);
  const int num_ints = 1 + 1 * n_cells_max;
  const int num_doubles = 4 * n_cells_max;

  MPI_Pack_size (num_doubles, MPI_DOUBLE, MPI_COMM_WORLD, &size_doubles);
  MPI_Pack_size (num_ints, MPI_INT, MPI_COMM_WORLD, &size_ints);
  const int comm_buffer_size = size_doubles + size_ints;
  char *comm_buffer = malloc (comm_buffer_size);        // comm_buffer_size is already in bytes

  for (current_rank = 0; current_rank < np_mpi_global; ++current_rank)
  {
    if (rank_global == current_rank)
    {
      int pack_position = 0;
      MPI_Pack (&n_cells_rank, 1, MPI_INT, comm_buffer, comm_buffer_size, &pack_position, MPI_COMM_WORLD);
      for (n_plasma = n_start; n_plasma < n_stop; ++n_plasma)
      {
        MPI_Pack (&n_plasma, 1, MPI_INT, comm_buffer, comm_buffer_size, &pack_position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].lum_tot, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &pack_position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].lum_lines, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &pack_position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].lum_rr, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &pack_position, MPI_COMM_WORLD);
        MPI_Pack (&plasmamain[n_plasma].lum_ff, 1, MPI_DOUBLE, comm_buffer, comm_buffer_size, &pack_position, MPI_COMM_WORLD);
      }
    }

    MPI_Bcast (comm_buffer, comm_buffer_size, MPI_PACKED, current_rank, MPI_COMM_WORLD);

    if (rank_global != current_rank)
    {
      int num_cells;
      int unpack_position = 0;
      MPI_Unpack (comm_buffer, comm_buffer_size, &unpack_position, &num_cells, 1, MPI_INT, MPI_COMM_WORLD);
      for (n_plasma = 0; n_plasma < num_cells; ++n_plasma)
      {
        int cell;
        MPI_Unpack (comm_buffer, comm_buffer_size, &unpack_position, &cell, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &unpack_position, &plasmamain[cell].lum_tot, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &unpack_position, &plasmamain[cell].lum_lines, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &unpack_position, &plasmamain[cell].lum_rr, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &unpack_position, &plasmamain[cell].lum_ff, 1, MPI_DOUBLE, MPI_COMM_WORLD);

      }
    }
  }
  free (comm_buffer);
#endif
}


/**********************************************************/
/**
 * @brief
 *
 * @details
 *
 * ### Notes ###
 *
 **********************************************************/

void
communicate_alpha_sp (const int n_start, const int n_stop, const int n_cells_rank)
{
#ifdef MPI_ON
  int i;
  int int_size;
  int double_size;
  int current_rank;

  const int n_cells_max = ceil ((double) NPLASMA / np_mpi_global);
  MPI_Pack_size (1 + n_cells_max, MPI_INT, MPI_COMM_WORLD, &int_size);
  MPI_Pack_size (n_cells_max * (2 * size_alpha_est + 2 * nphot_total), MPI_DOUBLE, MPI_COMM_WORLD, &double_size);
  int comm_buffer_size = double_size + int_size;
  char *comm_buffer = malloc (comm_buffer_size);        // comm_buffer_size is already in bytes


  for (current_rank = 0; current_rank < np_mpi_global; ++current_rank)
  {
    if (rank_global == current_rank)
    {
      int pack_position = 0;
      MPI_Pack (&n_cells_rank, 1, MPI_INT, comm_buffer, comm_buffer_size, &pack_position, MPI_COMM_WORLD);      // how many cells to unpack
      for (i = n_start; i < n_stop; ++i)
      {
        MPI_Pack (&i, 1, MPI_INT, comm_buffer, comm_buffer_size, &pack_position, MPI_COMM_WORLD);       // which cell we're working on
        if (nlevels_macro > 0)
        {
          MPI_Pack (macromain[i].recomb_sp, size_alpha_est, MPI_DOUBLE, comm_buffer, comm_buffer_size, &pack_position, MPI_COMM_WORLD);
          MPI_Pack (macromain[i].recomb_sp_e, size_alpha_est, MPI_DOUBLE, comm_buffer, comm_buffer_size, &pack_position, MPI_COMM_WORLD);
        }
        if (nphot_total > 0)
        {
          MPI_Pack (plasmamain[i].recomb_simple, nphot_total, MPI_DOUBLE, comm_buffer, comm_buffer_size, &pack_position, MPI_COMM_WORLD);
          MPI_Pack (plasmamain[i].recomb_simple_upweight, nphot_total, MPI_DOUBLE, comm_buffer, comm_buffer_size,
                    &pack_position, MPI_COMM_WORLD);
        }
      }
    }

    // TODO: this should be non-blocking eventually
    MPI_Bcast (comm_buffer, comm_buffer_size, MPI_PACKED, current_rank, MPI_COMM_WORLD);

    if (rank_global != current_rank)
    {
      int unpack_position = 0;
      int n_cells_to_do;
      MPI_Unpack (comm_buffer, comm_buffer_size, &unpack_position, &n_cells_to_do, 1, MPI_INT, MPI_COMM_WORLD);
      for (i = 0; i < n_cells_to_do; ++i)
      {
        int current_cell;
        MPI_Unpack (comm_buffer, comm_buffer_size, &unpack_position, &current_cell, 1, MPI_INT, MPI_COMM_WORLD);

        if (nlevels_macro > 0)
        {
          MPI_Unpack (comm_buffer, comm_buffer_size, &unpack_position, macromain[current_cell].recomb_sp, size_alpha_est,
                      MPI_DOUBLE, MPI_COMM_WORLD);
          MPI_Unpack (comm_buffer, comm_buffer_size, &unpack_position, macromain[current_cell].recomb_sp_e,
                      size_alpha_est, MPI_DOUBLE, MPI_COMM_WORLD);
        }
        if (nphot_total > 0)
        {
          MPI_Unpack (comm_buffer, comm_buffer_size, &unpack_position, plasmamain[current_cell].recomb_simple,
                      nphot_total, MPI_DOUBLE, MPI_COMM_WORLD);
          MPI_Unpack (comm_buffer, comm_buffer_size, &unpack_position, plasmamain[current_cell].recomb_simple_upweight,
                      nphot_total, MPI_DOUBLE, MPI_COMM_WORLD);
        }
      }
    }
  }

  free (comm_buffer);
#endif
}


/**********************************************************/
/**
 * @brief
 *
 * @details
 *
 * ### Notes ###
 *
 **********************************************************/

int
communicate_macro_cells (int n_start, int n_stop)
{
#ifdef MPI_ON
  int n_plasma;
  int current_rank;
  int i;
  int size_ints;
  int size_doubles;

  const int n_do = n_stop - n_start;
  const int n_cells_max = ceil ((double) NPLASMA / np_mpi_global);
  const int num_ints = 1 + n_cells_max + 2 * n_cells_max;
  const int num_doubles = ((6 * size_gamma_est) + (2 * size_Jbar_est)) * n_cells_max;
  MPI_Pack_size (num_ints, MPI_INT, MPI_COMM_WORLD, &size_ints);
  MPI_Pack_size (num_doubles, MPI_DOUBLE, MPI_COMM_WORLD, &size_doubles);
  const int comm_buffer_size = size_ints + size_doubles;
  char *comm_buffer = malloc (comm_buffer_size);

  Log ("Beginning communication for macro cells\n");

  for (current_rank = 0; current_rank < np_mpi_global; ++current_rank)
  {
    if (rank_global == current_rank)
    {
      int pack_position = 0;
      MPI_Pack (&n_do, 1, MPI_INT, comm_buffer, comm_buffer_size, &pack_position, MPI_COMM_WORLD);
      for (n_plasma = n_start; n_plasma < n_stop; ++n_plasma)
      {
        MPI_Pack (&n_plasma, 1, MPI_INT, comm_buffer, comm_buffer_size, &pack_position, MPI_COMM_WORLD);
        MPI_Pack (&macromain[n_plasma].kpkt_rates_known, 1, MPI_INT, comm_buffer, comm_buffer_size, &pack_position, MPI_COMM_WORLD);
        MPI_Pack (&macromain[n_plasma].matrix_rates_known, 1, MPI_INT, comm_buffer, comm_buffer_size, &pack_position, MPI_COMM_WORLD);

        // gamma
        MPI_Pack (macromain[n_plasma].gamma, size_gamma_est, MPI_DOUBLE, comm_buffer, comm_buffer_size, &pack_position, MPI_COMM_WORLD);
        MPI_Pack (macromain[n_plasma].gamma_old, size_gamma_est, MPI_DOUBLE, comm_buffer, comm_buffer_size, &pack_position, MPI_COMM_WORLD);

        // gamma_e
        MPI_Pack (macromain[n_plasma].gamma_e, size_gamma_est, MPI_DOUBLE, comm_buffer, comm_buffer_size, &pack_position, MPI_COMM_WORLD);
        MPI_Pack (macromain[n_plasma].gamma_e_old, size_gamma_est, MPI_DOUBLE, comm_buffer, comm_buffer_size,
                  &pack_position, MPI_COMM_WORLD);

        // alpha_st
        MPI_Pack (macromain[n_plasma].alpha_st, size_gamma_est, MPI_DOUBLE, comm_buffer, comm_buffer_size, &pack_position, MPI_COMM_WORLD);
        MPI_Pack (macromain[n_plasma].alpha_st_old, size_gamma_est, MPI_DOUBLE, comm_buffer, comm_buffer_size,
                  &pack_position, MPI_COMM_WORLD);

        // jbar
        MPI_Pack (macromain[n_plasma].jbar, size_Jbar_est, MPI_DOUBLE, comm_buffer, comm_buffer_size, &pack_position, MPI_COMM_WORLD);
        MPI_Pack (macromain[n_plasma].jbar_old, size_Jbar_est, MPI_DOUBLE, comm_buffer, comm_buffer_size, &pack_position, MPI_COMM_WORLD);
      }
    }

    MPI_Bcast (comm_buffer, comm_buffer_size, MPI_PACKED, current_rank, MPI_COMM_WORLD);

    if (rank_global != current_rank)
    {
      int unpack_position = 0;
      int n_cells_to_do;
      MPI_Unpack (comm_buffer, comm_buffer_size, &unpack_position, &n_cells_to_do, 1, MPI_INT, MPI_COMM_WORLD);
      for (i = 0; i < n_cells_to_do; ++i)
      {
        int current_cell;
        MPI_Unpack (comm_buffer, comm_buffer_size, &unpack_position, &current_cell, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &unpack_position, &macromain[current_cell].kpkt_rates_known, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &unpack_position, &macromain[current_cell].matrix_rates_known, 1, MPI_INT,
                    MPI_COMM_WORLD);

        // gamma
        MPI_Unpack (comm_buffer, comm_buffer_size, &unpack_position, macromain[current_cell].gamma, size_gamma_est, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &unpack_position, macromain[current_cell].gamma_old, size_gamma_est, MPI_DOUBLE,
                    MPI_COMM_WORLD);

        // gamma_e
        MPI_Unpack (comm_buffer, comm_buffer_size, &unpack_position, macromain[current_cell].gamma_e, size_gamma_est, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &unpack_position, macromain[current_cell].gamma_e_old, size_gamma_est, MPI_DOUBLE,
                    MPI_COMM_WORLD);

        // alpha_st
        MPI_Unpack (comm_buffer, comm_buffer_size, &unpack_position, macromain[current_cell].alpha_st, size_gamma_est, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &unpack_position, macromain[current_cell].alpha_st_old, size_gamma_est, MPI_DOUBLE,
                    MPI_COMM_WORLD);

        // jbar
        MPI_Unpack (comm_buffer, comm_buffer_size, &unpack_position, macromain[current_cell].jbar, size_Jbar_est, MPI_DOUBLE,
                    MPI_COMM_WORLD);
        MPI_Unpack (comm_buffer, comm_buffer_size, &unpack_position, macromain[current_cell].jbar_old, size_Jbar_est, MPI_DOUBLE,
                    MPI_COMM_WORLD);
      }
    }
  }

  free (comm_buffer);
  Log ("Ending communication for macro cells\n");
#endif
  return EXIT_SUCCESS;
}
