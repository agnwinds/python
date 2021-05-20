/** ************************************************************************* */
/**
 * @file     py_optical_depth_util.c
 * @author   Edward Parkinson
 * @date     May 2021
 *
 * @brief
 *
 * ************************************************************************** */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


#include "atomic.h"
#include "python.h"
#include "py_optical_depth.h"

/* ************************************************************************* */
/**
 * @brief
 *
 * @details
 *
 * ************************************************************************** */

SightLines_t *
outward_initialize_2d_model_angles (int *n_angles)
{
  int i;
  long mem_req;
  const double default_phase = 0.5;
  const double default_angles[] = { 0.0, 10.0, 30.0, 45.0, 60.0, 75.0, 85.0, 90.0 };
  const int n_default_angles = sizeof default_angles / sizeof default_angles[0];
  SightLines_t *inclinations = NULL;

  /*
   * Use the angles specified for by the user for spectrum generation, this
   * requires for xxspec to be initialised. Otherwise use the pre-defined angles
   * above.
   */

  if (xxspec != NULL && geo.nangles > 0)
  {
    *n_angles = geo.nangles;
    inclinations = calloc (geo.nangles, sizeof *inclinations);

    if (inclinations == NULL)
    {
      mem_req = geo.nangles * (int) sizeof *inclinations;
      printf ("outward_initialize_2d_model_angles: cannot allocate %ld bytes for observers array\n", mem_req);
      exit (EXIT_FAILURE);
    }

    for (i = MSPEC; i < MSPEC + geo.nangles; i++)
    {
      strcpy (inclinations[i - MSPEC].name, xxspec[i].name);
      stuff_v (xxspec[i].lmn, inclinations[i - MSPEC].lmn);
    }
  }
  else
  {
    printf ("outward_initialize_2d_model_angles: No spec.save file has been found, using a default set of inclination angles\n");
    *n_angles = n_default_angles;
    inclinations = calloc (n_default_angles, sizeof *inclinations);

    if (inclinations == NULL)
    {
      mem_req = n_default_angles * (int) sizeof *inclinations;
      printf ("initialize_inclination_angles: cannot allocate %ld bytes for observers array\n", mem_req);
      exit (EXIT_FAILURE);
    }

    for (i = 0; i < n_default_angles; i++)
    {
      sprintf (inclinations[i].name, "A%02.0fP%04.2f", default_angles[i], default_phase);
      inclinations[i].lmn[0] = sin (default_angles[i] / RADIAN) * cos (-default_phase * 360.0 / RADIAN);
      inclinations[i].lmn[1] = sin (default_angles[i] / RADIAN) * sin (-default_phase * 360.0 / RADIAN);
      inclinations[i].lmn[2] = cos (default_angles[i] / RADIAN);
    }
  }

  return inclinations;
}

/* ************************************************************************* */
/**
 * @brief
 *
 * @details
 *
 * ************************************************************************** */

SightLines_t *
outward_initialize_1d_model_angles (int *n_angles)
{
  const int n_default_angles = 1;
  const double default_angle = 45.0;
  const double default_phase = 0.5;
  SightLines_t *inclinations = NULL;

  *n_angles = n_default_angles;
  inclinations = calloc (n_default_angles, sizeof (SightLines_t));

  if (inclinations == NULL)
  {
    printf ("outward_initialize_1d_model_angles: unable to allocate %ld bytes for observers array\n", sizeof *inclinations);
    exit (EXIT_FAILURE);
  }

  sprintf (inclinations[0].name, "A%02.0fP%04.2f", default_angle, default_phase);
  inclinations[0].lmn[1] = sin (default_angle / RADIAN) * sin (-default_phase * 360.0 / RADIAN);
  inclinations[0].lmn[0] = sin (default_angle / RADIAN) * cos (-default_phase * 360.0 / RADIAN);
  inclinations[0].lmn[2] = cos (default_angle / RADIAN);


  return inclinations;
}

/* ************************************************************************* */
/**
 * @brief
 *
 * @details
 *
 *
 * ************************************************************************** */

SightLines_t *
photosphere_initialize_angles (int *n_angles)
{
  int i;
  const double default_phase = 1;
  const int n_default_angles = 180;
  const double d_theta = 90.0 / (double) n_default_angles;
  SightLines_t *inclinations = NULL;

  *n_angles = n_default_angles;
  inclinations = calloc (n_default_angles, sizeof (SightLines_t));

  if (inclinations == NULL)
  {
    printf ("photosphere_initialize_angles: unable to allocate memory for sightlines array\n");
    exit (EXIT_FAILURE);
  }

  for (i = 0; i < n_default_angles; i++)
  {
    double default_angle = i * d_theta;
    sprintf (inclinations[i].name, "A%02.0fP%04.2f", default_angle, default_phase);
    inclinations[i].lmn[0] = sin (default_angle / RADIAN) * cos (-default_phase * 360.0 / RADIAN);
    inclinations[i].lmn[1] = sin (default_angle / RADIAN) * sin (-default_phase * 360.0 / RADIAN);
    inclinations[i].lmn[2] = cos (default_angle / RADIAN);
  }

  return inclinations;
}

/* ************************************************************************* */
/**
 * @brief
 *
 * @details
 *
 *
 * ************************************************************************** */

SightLines_t *
initialize_inclination_angles (int *n_angles)
{
  SightLines_t *inclinations;

  if (MODE == RUN_MODE_OUTWARD)
  {
    if (zdom[N_DOMAIN].coord_type == SPHERICAL)
    {
      inclinations = outward_initialize_1d_model_angles (n_angles);
    }
    else
    {
      inclinations = outward_initialize_2d_model_angles (n_angles);
    }
  }
  else
  {
    inclinations = photosphere_initialize_angles (n_angles);
  }

  return inclinations;
}

/* ************************************************************************* */
/**
 * @brief  Generate a photon packet with a given frequency nu for use with the
 *         optical depth diagnostic routines.
 *
 * @param[in,out]  p_out  The photon packet to initialise
 * @param[in]  freq  The frequency of the photon packet
 *
 * @return  EXIT_SUCCESS or EXIT_FAILURE
 *
 * @details
 *
 * This routine assumes that the user wishes for the photon to be generated from
 * the surface of the central source, taking into account the current direction
 * of the sight-line being extracted. It's currently not possible, without a
 * tedious re-write, to place the photon at the very origin of the grid when
 * there is a central source because of how we check boundaries.
 *
 * Note that photons are initialised with a weight of f_tot as photons are
 * required to have weight, but since functions do not care about the weight of
 * the photon, it is set to something large to make sure it does not get
 * destroyed by accident somehow.
 *
 * ************************************************************************** */

int
create_photon (PhotPtr p_out, double freq, double *lmn)
{
  int n;

  if (freq < 0)
  {
    printf ("create_photon: photon can't be created with negative frequency\n");
    return EXIT_FAILURE;
  }

  p_out->freq = p_out->freq_orig = freq;
  p_out->origin = p_out->origin_orig = PTYPE_DISK;
  p_out->istat = P_INWIND;
  p_out->w = p_out->w_orig = geo.f_tot;
  p_out->tau = 0.0;
  p_out->frame = F_OBSERVER;
  p_out->x[0] = p_out->x[1] = p_out->x[2] = 0.0;
  stuff_v (lmn, p_out->lmn);

  if (MODE == RUN_MODE_OUTWARD)
  {
    move_phot (p_out, geo.rstar + DFUDGE);
  }
  else
  {
    move_phot (p_out, zdom[N_DOMAIN].rmax - DFUDGE);
    printf ("p_out->x [%e, %e, %e] length(x) = %e\n", p_out->x[0], p_out->x[1], p_out->x[2], length (p_out->x));
    for (n = 0; n < 3; ++n)
    {
      p_out->lmn[n] *= -1.0;
    }
  }

  return EXIT_SUCCESS;
}
