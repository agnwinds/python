/** ************************************************************************* */
/**
 * @file     py_optical_depth_output.c
 * @author   Edward Parkinson
 * @date     May 2021
 *
 * @brief
 *
 * ************************************************************************** */

#include <stdio.h>
#include <stdlib.h>

#include "atomic.h"
#include "python.h"
#include "py_optical_depth.h"

/* ************************************************************************** */
/**
 * @brief  Print the various optical depths calculated using this routine
 *
 * @details
 *
 * Prints the different optical depths for each angle and optical depth.
 * Historically, this used to also print to file. But I found that it was mostly
 * useless to do this.
 *
 * ************************************************************************** */

void
print_optical_depths (SightLines_t * inclinations, int n_inclinations, Edges_t edges[], int n_edges, double *optical_depth_values,
                      double *column_density_values)
{
  int i, j;
  int c_linelen;
  char str[LINELENGTH];
  const int MAX_COL = 120;

  printf ("\nOptical depths along the defined line of sights for domain %i:\n\n", N_DOMAIN);

  for (i = 0; i < n_inclinations; i++)
  {
    if (COLUMN_MODE == COLUMN_MODE_RHO)
    {
      printf ("%-8s: Mass column density     : %3.2e cm^-2\n", inclinations[i].name, column_density_values[i]);
      printf ("%-8s: Hydrogen column density : %3.2e cm^-2\n", "", column_density_values[i] * rho2nh);
    }
    else
    {
      printf ("%-8s: %s %i column density    : %3.2e cm^-2\n", inclinations[i].name, ele[ion[COLUMN_MODE_ION_NUMBER].nelem].name,
              ion[COLUMN_MODE_ION_NUMBER].istate, column_density_values[i]);
    }

    c_linelen = 0;
    for (j = 0; j < n_edges; j++)
    {
      c_linelen += snprintf (str, LINELENGTH, "tau_%-9s: %3.2e  ", edges[j].name, optical_depth_values[i * n_edges + j]);
      if (c_linelen > MAX_COL)
      {
        c_linelen = 0;
        printf ("\n");
      }
      printf ("%s", str);
    }
    printf ("\n\n");
  }
}

/* ************************************************************************* */
/**
 * @brief           Write the various optical depth spectra to file
 *
 * @details
 *
 * Simply write the optical depth spectra to the file named root.tau_spec.diag.
 * This file will be located in the diag folder.
 *
 * ************************************************************************** */

void
write_optical_depth_spectrum (SightLines_t * inclinations, int n_inclinations, double *tau_spectrum, double freq_min, double d_freq)
{
  int i, j;
  double c_wavelength, c_frequency;
  char filename[LINELENGTH];
  FILE *fp;

  int len = snprintf (filename, LINELENGTH, "%s.spec_tau", files.root);
  if (len < 0)
  {
    exit (EXIT_FAILURE);
  }

  fp = fopen (filename, "w");
  if (fp == NULL)
  {
    printf ("write_optical_depth_spectrum: uh oh, could not open optical depth spectrum output file\n");
    exit (EXIT_FAILURE);
  }

  /*
   * Write out the spec_tau header
   */

  fprintf (fp, "%-12s %-12s ", "Freq.", "Lambda");
  for (j = 0; j < n_inclinations; j++)
    fprintf (fp, "%-12s ", inclinations[j].name);
  fprintf (fp, "\n");

  /*
   * Write out the spectrum for each inclination angle, omne line at a time
   */

  c_frequency = freq_min;
  for (i = 0; i < N_FREQ_BINS; i++)
  {
    c_wavelength = VLIGHT / c_frequency / ANGSTROM;
    fprintf (fp, "%-12e %-12e ", c_frequency, c_wavelength);
    for (j = 0; j < n_inclinations; j++)
      fprintf (fp, "%-12e ", tau_spectrum[j * N_FREQ_BINS + i]);
    fprintf (fp, "\n");
    c_frequency += d_freq;
  }

  if (fclose (fp))
    printf ("write_optical_depth_spectrum: uh oh, could not close optical depth spectrum output file\n");
}

/* ************************************************************************* */
/**
 * @brief
 *
 * @details
 *
 * ************************************************************************** */

void
write_photosphere_location_to_file (Positions_t * positions, int n_inclinations)
{
  int i;

  char filename[LINELENGTH];
  FILE *fp;

  int len = snprintf (filename, LINELENGTH, "%s.photosphere.txt", files.root);
  if (len < 0)
  {
    exit (EXIT_FAILURE);
  }

  fp = fopen (filename, "w");

  fprintf (fp, "# TAU DEPTH = %f\n", TAU_DEPTH);
  fprintf (fp, "# x y z\n");

  for (i = 0; i < n_inclinations; i++)
  {
    fprintf (fp, "%e %e %e\n", positions[i].x, positions[i].y, positions[i].z);
  }

  fclose (fp);
}
