/***********************************************************/
/** @file  communicate_wind.c
 * @author EJP
 * @date   December 2023
 *
 * @brief Functions for communicating wind properties
 *
 * @TODO: as much as this as possible should use non-blocking communication
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
 * @brief sum up the synthetic and cell spectra between threads.
 *
 * @details
 * sum up the synthetic spectra between threads. Does an
 * MPI_Reduce then an MPI_Bcast for each element of the
 * linear and log spectra arrays (xxspec)
 *
 **********************************************************/

int
gather_extracted_spectrum (void)
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
