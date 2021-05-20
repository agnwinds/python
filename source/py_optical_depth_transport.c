/** ************************************************************************* */
/**
 * @file     py_optical_depth_transport.c
 * @author   Edward Parkinson
 * @date     April 2019
 *
 * @brief
 *
 * ************************************************************************** */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "atomic.h"
#include "python.h"
#include "py_optical_depth.h"

const double MAXDIFF = VCHECK / VLIGHT; // For linear velocity requirement for photon transport

/* ************************************************************************* */
/**
 * @brief  Calculate the total optical depth a photon experiences across the
 *         cell of distance SMAX_FRAC * smax.

 * @param[in]  photon  The photon packet
 * @param[in,out]  *c_column_density  The column density the photon has moved
 *                                    through
 * @param[in,out]  *c_optical_depth  The optical depth experienced by the photon
 *
 * @return p_istat  The current photon status or EXIT_FAILURE on failure.
 *
 * @details
 *
 * This function is concerned with finding the opacity of the photon's current
 * cell, the distance the photon can move in the cell and hence it increments
 * the optical depth tau a photon has experienced as it moves through the wind.
 *
 * ************************************************************************** */

int
integrate_tau_across_cell (PhotPtr photon, double *c_column_density, double *c_optical_depth)
{
  int p_istat;
  int n_domain, n_plasma;
  double kappa_total, density;
  double smax, diff;
  double freq_inner, freq_outer, mean_freq;
  WindPtr c_wind_cell;
  PlasmaPtr c_plasma_cell;
  struct photon p_start, p_stop, p_now;

  photon->grid = where_in_grid (wmain[photon->grid].ndom, photon->x);
  if (photon->grid < 0)
  {
    printf ("integrate_tau_across_cell: photon is not in a grid cell\n");
    return EXIT_FAILURE;
  }

  /*
   * Determine where in the plasma grid the cell is. This is required so a
   * column density can be calculated. Also calculate the maximum distance the
   * photon can traverse across the cell
   */

  c_wind_cell = &wmain[photon->grid];
  n_domain = c_wind_cell->ndom;
  n_plasma = wmain[photon->grid].nplasma;
  c_plasma_cell = &plasmamain[n_plasma];

  if (COLUMN_MODE == COLUMN_MODE_RHO)
  {
    density = c_plasma_cell->rho;
  }
  else
  {
    density = c_plasma_cell->density[COLUMN_MODE_ION_NUMBER];
  }

  smax = smax_in_cell (photon);
  if (smax < 0)
  {
    printf ("integrate_tau_across_cell: smax %e < 0 in cell %d\n", smax, photon->grid);
    return EXIT_FAILURE;
  }

  /*
   * Transform the photon into the CMF and create a photon in the CMF at the
   * end of the path of length smax
   */

  observer_to_local_frame (photon, &p_start);
  stuff_phot (photon, &p_stop);
  move_phot (&p_stop, smax);
  observer_to_local_frame (&p_stop, &p_stop);

  /* At this point p_start and p_stop are in the local frame
   * at the and p_stop is at the maximum distance it can
   * travel. We want to check that the frequency shift is
   * not too great along the path that a linear approximation
   * to the change in frequency is not reasonable
   */

  while (smax > DFUDGE)
  {
    stuff_phot (photon, &p_now);
    move_phot (&p_now, smax * 0.5);
    observer_to_local_frame (&p_now, &p_now);
    diff = fabs (p_now.freq - 0.5 * (p_start.freq + p_stop.freq)) / p_start.freq;
    if (diff < MAXDIFF)
      break;
    stuff_phot (&p_now, &p_stop);
    smax *= 0.5;
  }

  freq_inner = p_start.freq;
  freq_outer = p_stop.freq;
  mean_freq = 0.5 * (freq_inner + freq_outer);

  /*
   * Now we can finally calculate the opacity due to all the continuum
   * processes. In macro-atom mode, we need to calculate the continuum opacity
   * using kappa_bf and kappa_ff using the macro treatment. For simple mode, we
   * can simply use radiation which **SHOULD** return the continuum opacity
   * as well, plus something from induced Compton heating. In either cases,
   * we still then need to add the optical depth from electron scattering at
   * the end.
   */

  kappa_total = 0;

  if (MODE != RUN_MODE_PHOTOSPHERE)
  {
    if (geo.rt_mode == RT_MODE_2LEVEL)
    {
      kappa_total += radiation (photon, smax);
    }
    else                        // macro atom case
    {
      if (c_wind_cell->vol > 0)
      {
        kappa_total += kappa_bf (c_plasma_cell, freq_inner, 0);
        kappa_total += kappa_ff (c_plasma_cell, freq_inner);
      }
    }
  }

  kappa_total += klein_nishina (mean_freq) * c_plasma_cell->ne * zdom[n_domain].fill;

  /*
   * Increment the optical depth and column density variables and move the
   * photon to the edge of the cell
   */

  *c_column_density += smax * density;
  *c_optical_depth += smax * kappa_total;
  move_phot (photon, smax);
  p_istat = photon->istat;

  return p_istat;
}

/* ************************************************************************* */
/**
 * @brief           Extract the optical depth the photon packet porig must
 *                  travel through to reach the observer.
 *
 * @param[in]  photon  The photon packet to extract
 * @param[out]  *c_column_density  The column depth of the extracted photon angle
 * @param[out]  *c_optical_depth  The optical depth from photon origin to
 *                                the observer
 *
 * @return  EXIT_SUCCESS or EXIT_FAILURE
 *
 * @details
 *
 * The direction of the observer is set in porig.lmn and should be set prior
 * to passing a photon to this function. If any values of lmn are negative, the
 * photon will translate in the negative direction. However, have no fear as this
 * is normal and is fine due to the assumed symmetry of models in Python.
 *
 * ************************************************************************** */

int
integrate_tau_across_wind (PhotPtr photon, double *c_column_density, double *c_optical_depth)
{
  int err;
  int n_dom, where;
  enum istat_enum p_istat;
  const int max_translate_in_space = 10;
  int n_in_space;
  double norm[3];
  struct photon p_extract;

  p_istat = P_INWIND;           // assume photon is in wind for initialisation reasons
  stuff_phot (photon, &p_extract);

  n_in_space = 0;
  while (p_istat == P_INWIND)
  {
    where = where_in_wind (p_extract.x, &n_dom);

    if (where < 0)
    {
      translate_in_space (&p_extract);
      if (++n_in_space > max_translate_in_space)
      {
        printf ("integrate_tau_across_wind: something has gone wrong as this photon has translated in space %d times\n", n_in_space);
        return EXIT_FAILURE;
      }
    }
    else if ((p_extract.grid = where_in_grid (n_dom, p_extract.x)) >= 0)
    {
      err = integrate_tau_across_cell (&p_extract, c_column_density, c_optical_depth);
      if (err)
        return EXIT_FAILURE;
    }
    else
    {
      printf ("integrate_tau_across_wind: photon in unknown location grid stat %i\n", p_extract.grid);
      return EXIT_FAILURE;
    }

    p_istat = walls (&p_extract, photon, norm);

    if (MODE == RUN_MODE_PHOTOSPHERE)
    {
      if (*c_optical_depth >= TAU_DEPTH)
      {
        // printf ("Reached optical depth limit of %e %e\n", *c_optical_depth, TAU_DEPTH);
        // printf ("p_extract->x [%e, %e, %e]\n", p_extract.x[0], p_extract.x[1], p_extract.x[2]);
        // printf ("\n");
        break;
      }
    }
  }

  /*
   * If we are in RUN_MODE_PHOTOSPHERE, then we shouldn't care about hitting the
   * star or disc, since we are aiming for the origin of the system
   */

  if (MODE == RUN_MODE_OUTWARD)
  {
    if (p_istat == P_HIT_STAR || p_istat == P_HIT_DISK)
    {
      printf ("integrate_tau_across_wind: photon hit central source or disk incorrectly istat = %i\n", p_istat);
      return EXIT_FAILURE;
    }
  }
  else
  {
    stuff_phot (&p_extract, photon);
    if (p_istat == P_HIT_DISK)
    {
      printf ("integrate_tau_across_wind: the photon hit the disk whilst in RUN_MODE_PHOTOSPHERE when it should hit the central source\n");
      return EXIT_FAILURE;
    }
    if (p_istat == P_HIT_STAR)
    {
      printf ("Hit the star, as we should\n");
    }
  }

  return EXIT_SUCCESS;
}

/* ************************************************************************* */
/**
 * @brief  Create spectra of tau vs lambda for each observer angle
 *
 * @details
 *
 * This is the main function which will generate the optical depth spectra for
 * each observer angle in xxspec. The algorithm is similar to extract and the
 * tau_diag algorithm which this function is called in.
 *
 * A photon is generated at the central source of the model and is extracted
 * from this location towards the observer where it escapes, where integrate_tau_across_wind
 * returns the integrated optical depth along its path to escape. This is done
 * for a range of photon frequencies to find how optical depth changes with
 * frequency.
 *
 * This processes can take some time compared to tau_evalulate_photo_edges. But,
 * since N_FREQ_BINS photons are being generated for each spectrum and the fact
 * that these photons do not interact, the spectra does not actually take that
 * long to complete.
 *
 * ************************************************************************** */

void
create_optical_depth_spectrum (void)
{
  int i, j;
  int err;
  double *tau_spectrum;
  double c_optical_depth, c_column_density;
  double c_frequency, freq_min, freq_max, d_freq;
  struct photon photon;

  int n_inclinations;
  SightLines_t *inclinations = initialize_inclination_angles (&n_inclinations);

  printf ("Creating optical depth spectra:\n");

  tau_spectrum = calloc (n_inclinations * N_FREQ_BINS, sizeof *tau_spectrum);
  if (tau_spectrum == NULL)
  {
    printf ("create_optical_depth_spectrum: cannot allocate %lu bytes for tau_spectrum\n",
            n_inclinations * N_FREQ_BINS * sizeof *tau_spectrum);
    exit (EXIT_FAILURE);
  }

  /*
   * Define the limits of the spectra in frequency space. If xxpsec is NULL,
   * then the frequency range will be over a default 100 - 10,000 Angstrom
   * band.
   */

  if ((geo.nangles == 0 && xxspec == NULL) || (geo.swavemax == 0 && geo.swavemin == 0))
  {
    printf ("create_optical_depth_spectrum: xxspec is uninitialized, defaulting spectral wavelength range to 100 - 10,000 Angstrom\n");
    freq_min = VLIGHT / (10000 * ANGSTROM);
    freq_max = VLIGHT / (100 * ANGSTROM);
  }
  else
  {
    freq_min = VLIGHT / (geo.swavemax * ANGSTROM);
    freq_max = VLIGHT / (geo.swavemin * ANGSTROM);
    if (sane_check (freq_min))
    {
      freq_min = VLIGHT / (10000 * ANGSTROM);
      printf ("create_optical_depth_spectrum: freq_min has an invalid value setting to %e\n", freq_min);
    }
    if (sane_check (freq_max))
    {
      freq_max = VLIGHT / (100 * ANGSTROM);
      printf ("create_optical_depth_spectrum: freq_min has an invalid value setting to %e\n", freq_max);
    }
  }

  d_freq = (freq_max - freq_min) / N_FREQ_BINS;
  kbf_need (freq_min, freq_max);

  /*
   * Now create the optical depth spectra for each inclination
   */

  for (i = 0; i < n_inclinations; i++)
  {
    printf ("  - Creating spectrum: %s\n", inclinations[i].name);
    c_frequency = freq_min;
    for (j = 0; j < N_FREQ_BINS; j++)
    {
      c_optical_depth = 0.0;
      c_column_density = 0.0;

      err = create_photon (&photon, c_frequency, inclinations[i].lmn);
      if (err == EXIT_FAILURE)
      {
        printf ("create_optical_depth_spectrum: skipping photon of frequency %e when creating spectrum\n", c_frequency);
        continue;
      }
      photon.np = j;
      err = integrate_tau_across_wind (&photon, &c_column_density, &c_optical_depth);
      if (err == EXIT_FAILURE)
        continue;
      tau_spectrum[i * N_FREQ_BINS + j] = c_optical_depth;
      c_frequency += d_freq;
    }
  }

  write_optical_depth_spectrum (inclinations, n_inclinations, tau_spectrum, freq_min, d_freq);
  free (tau_spectrum);
  free (inclinations);
}

/* ************************************************************************* */
/**
 * @brief Calculate the optical depth for various optical depth edges and
 *        extract the column density.
 *
 * @details
 *
 * This is the main function which will control the procedure for calculating
 * various diagnostic numbers for the optical depth's experienced in the current
 * model. Namely, this function aims to show the total integrated optical depth
 * to each observer angle using (originally) the following optical depths:
 *
 *  - Lymann edge
 *  - Balmer edge
 *  - Helium II edge
 *
 * Once these integrated optical depths have been calculated for each angle, a
 * spectrum of optical depth vs wavelength is created for each angle.
 *
 * The aim of these diagnostic numbers it to provide some sort of quick metric
 * on the optical thickness of the current model.
 *
 * ************************************************************************** */

void
evaluate_photoionization_edges (void)
{
  int i, j;
  int err;
  double c_frequency, c_optical_depth, c_column_density;
  double *optical_depth_values = NULL, *column_density_values = NULL;
  struct photon photon;

  Edges_t edges[] = {
    {"HLymanEdge", 3.387485e+15},
    {"HBalmerEdge", 8.293014e+14},
    {"HeI24eVEdge", 5.9483e+15},
    {"HeII54eVEdge", 1.394384e+16}
  };

  const int n_edges = sizeof edges / sizeof edges[0];

  int n_inclinations;
  SightLines_t *inclinations = initialize_inclination_angles (&n_inclinations);

  optical_depth_values = calloc (n_inclinations * n_edges, sizeof *optical_depth_values);
  if (optical_depth_values == NULL)
  {
    printf ("evaluate_photoionization_edges: cannot allocate %lu bytes for optical_depths\n",
            n_inclinations * n_edges * sizeof *optical_depth_values);
    exit (EXIT_FAILURE);
  }

  column_density_values = calloc (n_inclinations, sizeof *column_density_values);
  if (column_density_values == NULL)
  {
    printf ("evaluate_photoionization_edges: cannot allocate %lu bytes for column_densities\n",
            n_inclinations * sizeof *column_density_values);
    exit (EXIT_FAILURE);
  }

  /*
   * Now extract the optical depths and mass column densities. We loop over
   * each PI edge for each inclination angle.
   */

  for (i = 0; i < n_inclinations; i++)
  {
    for (j = 0; j < n_edges; j++)
    {
      c_optical_depth = 0.0;
      c_column_density = 0.0;
      c_frequency = edges[j].freq;

      err = create_photon (&photon, c_frequency, inclinations[i].lmn);
      if (err == EXIT_FAILURE)
      {
        printf ("evaluate_photoionization_edges: skipping photon of frequency %e\n", c_frequency);
        continue;
      }

      err = integrate_tau_across_wind (&photon, &c_column_density, &c_optical_depth);
      if (err == EXIT_FAILURE)
        continue;               // do not throw extra warning when one is already thrown in integrate_tau_across_wind

      optical_depth_values[i * n_edges + j] = c_optical_depth;
      column_density_values[i] = c_column_density;
    }
  }

  print_optical_depths (inclinations, n_inclinations, edges, n_edges, optical_depth_values, column_density_values);
  free (inclinations);
  free (optical_depth_values);
  free (column_density_values);
}

/* ************************************************************************* */
/**
 * @brief
 *
 * @details
 *
 * ************************************************************************** */

void
find_photosphere (void)
{
  int i, err;
  struct photon photon;
  int n_inclinations;
  SightLines_t *inclinations;

  inclinations = initialize_inclination_angles (&n_inclinations);

  typedef struct Positions_s
  {
    double x, y, z;
  } Positions_t;

  Positions_t *positions = calloc (n_inclinations, sizeof (Positions_t));
  if (positions == NULL)
  {
    printf ("Unable to allocate memory for the positions array\n");
    exit (EXIT_FAILURE);
  }

  double tau = 0;
  double nh = 0;
  // const double test_freq = VLIGHT / (100 * ANGSTROM);
  const double test_freq = 8e14;

  Positions_t *original = calloc (n_inclinations, sizeof (Positions_t));


  FILE *fp3 = fopen ("lmn_directions.txt", "w");

  for (i = 0; i < n_inclinations; i++)
  {
    err = create_photon (&photon, test_freq, inclinations[i].lmn);
    if (err)
      printf ("find_photosphere: error creating photon\n");

    original[i].x = photon.x[0];
    original[i].y = photon.x[1];
    original[i].z = photon.x[2];

    err = integrate_tau_across_wind (&photon, &tau, &nh);
    if (err)
      printf ("find_photosphere: something went wrong with tau integration");

    positions[i].x = photon.x[0];
    positions[i].y = photon.x[1];
    positions[i].z = photon.x[2];

    fprintf (fp3, "%e %e %e\n", photon.lmn[0], photon.lmn[1], photon.lmn[2]);

  }

  fclose (fp3);

  FILE *fp = fopen ("locations.txt", "w");
  FILE *fp2 = fopen ("original_locations.txt", "w");

  for (i = 0; i < n_inclinations; i++)
  {
    printf ("angle %s", inclinations[i].name);
    printf (" %e %e %e\n", positions[i].x, positions[i].y, positions[i].z);

    printf ("angle %s", inclinations[i].name);
    printf (" %e %e %e\n\n", original[i].x, original[i].y, original[i].z);

    fprintf (fp, "%e %e %e\n", positions[i].x, positions[i].y, positions[i].z);
    fprintf (fp2, "%e %e %e\n", original[i].x, original[i].y, original[i].z);
  }

  fclose (fp);
  fclose (fp2);
}

/* ************************************************************************* */
/**
 * @brief  Main control function for create optical depth diagnostics.
 *
 * @details
 *
 * This function is the main steering function for generating the optical depth
 * diagnostics.
 *
 * ************************************************************************** */

void
control_program (void)
{
  if (MODE == RUN_MODE_OUTWARD)
  {
    evaluate_photoionization_edges ();
    create_optical_depth_spectrum ();
  }
  else if (MODE == RUN_MODE_PHOTOSPHERE)
  {
    find_photosphere ();
  }
  else
  {
    printf ("Mode %d is an unknown run mode, not sure how you got here so exiting the program\n", MODE);
    exit (EXIT_FAILURE);
  }
}
