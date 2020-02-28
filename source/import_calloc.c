/* ************************************************************************** */
/**
 * @file  import_calloc.c
 * @author EJP
 * @date   Feb 2019
 *
 * @brief    Functions for memory management for the import structures.
 *
 * ************************************************************************** */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"
#include "import.h"


/* ************************************************************************** */
/**
 * @brief   Allocate memory for the import_model structures.
 *
 * @param[in]   int coord_type     The coordinate system of the imported grid
 *
 * @details
 *
 * This function is simply a wrapper around a lot of calloc () calls to allocate
 * memory for the different arrays in the structure. It was opted to do this way
 * rather than static memory to avoid storing the raw imported grid details
 * when it isn't required.
 *
 * NDIM_MAX (or NDIM_MAX * NDIM_MAX) elements are allocated for each array, as
 * at this point we do not know how many elements comprise the grid. Note that
 *
 * ************************************************************************** */

void
calloc_import (int coord_type)
{
  if (coord_type == SPHERICAL)
  {
    import_model_1d.element = calloc (sizeof *import_model_1d.element, NDIM_MAX);
    import_model_1d.r = calloc (sizeof *import_model_1d.r, NDIM_MAX);
    import_model_1d.v_r = calloc (sizeof *import_model_1d.v_r, NDIM_MAX);
    import_model_1d.mass_rho = calloc (sizeof *import_model_1d.mass_rho, NDIM_MAX);
    import_model_1d.t_r = calloc (sizeof *import_model_1d.t_r, NDIM_MAX);
    import_model_1d.t_e = calloc (sizeof *import_model_1d.t_e, NDIM_MAX);
  }
  else if (coord_type == CYLIND || coord_type == RTHETA)
  {
    import_model_2d.i = calloc (sizeof *import_model_2d.i, NDIM_MAX * NDIM_MAX);
    import_model_2d.j = calloc (sizeof *import_model_2d.j, NDIM_MAX * NDIM_MAX);
    import_model_2d.inwind = calloc (sizeof *import_model_2d.inwind, NDIM_MAX * NDIM_MAX);
    import_model_2d.v_x = calloc (sizeof *import_model_2d.v_x, NDIM_MAX * NDIM_MAX);
    import_model_2d.v_y = calloc (sizeof *import_model_2d.v_y, NDIM_MAX * NDIM_MAX);
    import_model_2d.v_z = calloc (sizeof *import_model_2d.v_z, NDIM_MAX * NDIM_MAX);
    import_model_2d.mass_rho = calloc (sizeof *import_model_2d.mass_rho, NDIM_MAX * NDIM_MAX);
    import_model_2d.t_r = calloc (sizeof *import_model_2d.t_r, NDIM_MAX * NDIM_MAX);
    import_model_2d.wind_x = calloc (sizeof *import_model_2d.wind_x, NDIM_MAX * NDIM_MAX);
    import_model_2d.wind_z = calloc (sizeof *import_model_2d.wind_z, NDIM_MAX * NDIM_MAX);
    import_model_2d.wind_midx = calloc (sizeof *import_model_2d.wind_midx, NDIM_MAX * NDIM_MAX);
    import_model_2d.wind_midz = calloc (sizeof *import_model_2d.wind_midz, NDIM_MAX * NDIM_MAX);

    if (coord_type == CYLIND)
    {
      import_model_2d.x = calloc (sizeof *import_model_2d.x, NDIM_MAX * NDIM_MAX);
      import_model_2d.z = calloc (sizeof *import_model_2d.z, NDIM_MAX * NDIM_MAX);
    }
    else
    {
      import_model_2d.r = calloc (sizeof *import_model_2d.r, NDIM_MAX * NDIM_MAX);
      import_model_2d.theta = calloc (sizeof *import_model_2d.theta, NDIM_MAX * NDIM_MAX);
    }
  }
  else
  {
    Error ("%s: %i: Unknown coord_type %i\n", __FILE__, __LINE__, coord_type);
    Exit (1);
  }
}



/* ************************************************************************** */
/**
 * @brief    Free the memory allocated for the import structures.
 *
 * @param[in]   int coord_type     The coordinate system of the imported grid
 *
 * @details
 *
 * When the raw imported grid is no longer needed as the relevant data has been
 * assigned to the wind or domain structure, then we should free the memory
 * so as not to have extra data hanging around.
 *
 * ************************************************************************** */

void
free_import (int coord_type)
{
  if (coord_type == SPHERICAL)
  {
    free (import_model_1d.element);
    free (import_model_1d.r);
    free (import_model_1d.v_r);
    free (import_model_1d.mass_rho);
    free (import_model_1d.t_r);
    free (import_model_1d.t_e);
  }
  else if (coord_type == CYLIND || coord_type == RTHETA)
  {
    free (import_model_2d.i);
    free (import_model_2d.j);
    free (import_model_2d.inwind);
    free (import_model_2d.v_x);
    free (import_model_2d.v_y);
    free (import_model_2d.v_z);
    free (import_model_2d.mass_rho);
    free (import_model_2d.t_r);
    free (import_model_2d.wind_x);
    free (import_model_2d.wind_z);
    free (import_model_2d.wind_midx);
    free (import_model_2d.wind_midz);
    free (import_model_2d.x);
    free (import_model_2d.z);
    free (import_model_2d.r);
    free (import_model_2d.theta);
  }
  else
  {
    Error ("%s: %i: Unknown coord_type %i\n", __FILE__, __LINE__, coord_type);
    Exit (1);
  }
}
