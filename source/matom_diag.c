
/***********************************************************/
/** @file  matom_diag.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief  diagnostic macro-atom printouts.
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <errno.h>

#include "atomic.h"
#include "python.h"


/**********************************************************/
/** 
 * @brief      a routine which reports the matom level and kpkt emissivities summed
 * 	over all cells. It is called in define_phot() in the spectral cycles, after get_matom_f()
 *
 * @return     simply logs the level emissivites and kpkt emissivites in macro atom mode.
 *
 * @details
 * prints out the level emissivities in macro atoms, summed over each cell.
 *     	Should only be used in macro atom mode in the spectral cycles
 *     	(when geo.matom_radiation = 1).
 *
 * ### Notes ###
 *
 **********************************************************/

int
matom_emiss_report ()
{

  int n, m;
  double emiss_sum, abs_sum;

  /* Cycle over macro atom levels and log emissivities */

  for (m = 0; m < nlevels_macro; m++)
  {
    emiss_sum = 0.0;
    abs_sum = 0.0;
    for (n = 0; n < NPLASMA; n++)
    {
      emiss_sum += macromain[n].matom_emiss[m];
      abs_sum += macromain[n].matom_abs[m];
    }

    Log ("Macro Atom level emissivities (summed): z %2d i %2d macro %2d n %2d matom_abs %8.4e matom_emiss %8.4e\n",
         config[m].z, config[m].istate, config[m].macro_info, m, abs_sum, emiss_sum);
  }



  /* Log kpkt emissivities as well */
  emiss_sum = 0.0;
  abs_sum = 0.0;

  for (n = 0; n < NPLASMA; n++)
  {
    emiss_sum += plasmamain[n].kpkt_emiss;
    abs_sum += plasmamain[n].kpkt_abs;

  }

  Log ("Kpkt emissivities (summed over cells): kpkt_abs %8.4e kpkt_emiss %8.4e\n", abs_sum, emiss_sum);

  /* Log totals */
  Log ("Totals: f_matom %le f_kpkt %le\n", geo.f_matom, geo.f_kpkt);

  return (0);
}

/**********************************************************/
/**
 * @brief
 *
 * @param
 *
 * @return
 *
 * @details
 *
 **********************************************************/

int directory_init = FALSE;

void
create_matrix_diag_folder (void)
{
  int err;

  err = mkdir ("matrix_output", 0777);
  if (err)
  {
    if (errno != EEXIST)
    {
      perror ("create_matrix_diag_folder");
      Exit (1);
    }
  }
}

/**********************************************************/
/**
 * @brief
 *
 * @param
 *
 * @return
 *
 * @details
 *
 **********************************************************/

void
write_1d_matrix_to_file (char *filename, double *matrix, int len_i)
{
  int i;
  int err;
  FILE *fp;

  if (directory_init == FALSE)
  {
    create_matrix_diag_folder ();
    directory_init = TRUE;
  }

  fp = fopen (filename, "w");
  if (fp == NULL)
  {
    Error ("Unable to open %s to write 1D matrix to file\n", filename);
    perror ("write_1d_matrix_to_file");
    Exit (1);
  }

  for (i = 0; i < len_i; ++i)
  {
    fprintf (fp, "%-+15.6g", matrix[i]);
  }

  err = fclose (fp);
  if (err)
  {
    perror ("write_1d_matrix_to_file");
    Exit (1);
  }
}

/**********************************************************/
/**
 * @brief
 *
 * @param
 *
 * @return
 *
 * @details
 *
 **********************************************************/

void
write_flat_2d_matrix_to_file (char *filename, double *matrix, int len_i, int len_j)
{
  int i, j;
  int err;
  FILE *fp;

  if (directory_init == FALSE)
  {
    create_matrix_diag_folder ();
    directory_init = TRUE;
  }

  fp = fopen (filename, "w");
  if (fp == NULL)
  {
    Error ("Unable to open %s to write 1D matrix to file\n", filename);
    perror ("write_flat_2d_matrix_to_file");
    Exit (1);
  }

  for (i = 0; i < len_i; ++i)
  {
    for (j = 0; j < len_j; j++)
    {
      fprintf (fp, "%-+15.6g", matrix[i * len_j + j]);
    }
    fprintf (fp, "\n");
  }

  err = fclose (fp);
  if (err)
  {
    perror ("write_flat_2d_matrix_to_file");
    Exit (1);
  }
}


/**********************************************************/
/**
 * @brief
 *
 * @param
 *
 * @return
 *
 * @details
 *
 **********************************************************/

void
write_2d_matrix_to_file (char *filename, double matrix[NLEVELS_MACRO][NLEVELS_MACRO], int len_i, int len_j)
{
  int i, j;
  int err;
  FILE *fp;

  if (directory_init == FALSE)
  {
    create_matrix_diag_folder ();
    directory_init = TRUE;
  }

  fp = fopen (filename, "w");
  if (fp == NULL)
  {
    Error ("Unable to open %s to write 2D matrix to file\n", filename);
    perror ("write_2d_matrix_to_file");
    Exit (1);
  }

  for (i = 0; i < len_i; ++i)
  {
    for (j = 0; j < len_j; j++)
    {
      fprintf (fp, "%-+15.6g", matrix[i][j]);
    }
    fprintf (fp, "\n");
  }

  err = fclose (fp);
  if (err)
  {
    perror ("write_2d_matrix_to_file");
    Exit (1);
  }
}
