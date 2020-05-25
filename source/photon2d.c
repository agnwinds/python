
/***********************************************************/
/** @file  photon2d.c
 * @author ksl
 * @date   April, 2018
 *
 * @brief  This file contains the routines which translate
 * photons, both inside and outside of the domains.
 *
 * Except for translating the photon the routines do as
 * litle as possible w.r.t modifying any of the other
 * information in the photon structure.
 *
 * There is a steering routine translate which calle
 * either translate_in_wind or translate_in_space, depending
 * on whether a photon bundle is in a wind domain or in space.
 * The routines translate the photon until hits something,
 * either the edge of the wind, or the edge of a cell, or
 * a resonance and then return with a status noting what
 * has happened.
 *
 * ### Notes ###
 *
 * Some of these routines pass the WindPtr for the entire
 * wind which should not be necessary and had been removed
 * in favor of using wmain.  It would probably be a good idea
 * to that here.
 *
 * These routines are called from trans_phot.  One could argue
 * that the filename is misnamed for that reason.
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


/**********************************************************/
/**
 * @brief      a steering routine that either calls translate_in_space or translate_in_wind  depending upon the
 * 	current location of the photon.
 *
 * @param [in] WindPtr  w   A pointer to the wind
 * @param [in, out] PhotPtr  pp   A photon
 * @param [in] double  tau_scat   the depth at which the photon should scatter
 * @param [ out] double *  tau   The optical depth associated with a resonanant scatter when
 * this is the reason the photon has stopped
 * @param [in, out] int *  nres   The number of the resonance if the photon stopped.
 * due to reaching the scattering optical depth
 * @return     A status that states what caused the photon to stp as it did
 *
 * @details
 *
 * The routine takes the current position of the photon and determines
 * where it is.  Depending on this it calls either translate_in_space
 * or translate_in_wind
 *
 * On exit, the position of the photon will have been updated.
 *
 * ### Notes ###
 * Translate controls the flow of the photon through one grid cell.  In Python,
 * there are additional possibilities since there are wind free regions.
 *
 * The routine first calls where_in_wind which returns the correct domain for a given position, or domain 0 if the
 * position is not in the wind of a domain.  The photon structure does not have the domain
 * number directly encoded, but it can be obtained from the grid number, which where in_grid updates
 *
 *
 * This construction is
 * used so that the calling routines for translate (trans_phot and extract) can be the same
 * in the 1 and 2d versions.
 *
 **********************************************************/

int
translate (w, pp, tau_scat, tau, nres)
     WindPtr w;                 //w here refers to entire wind, not a single element
     PhotPtr pp;
     double tau_scat;
     double *tau;
     int *nres;
{
  int istat;
  int ndomain;



  if (where_in_wind (pp->x, &ndomain) < 0)      //If the return is negative, this means we are outside the wind
  {
    istat = translate_in_space (pp);    //And so we should translate in space

  }
  else if ((pp->grid = where_in_grid (ndomain, pp->x)) >= 0)
  {
    istat = translate_in_wind (w, pp, tau_scat, tau, nres);
  }
  else
  {
    istat = pp->istat = P_ERROR;        /* It's not in the wind and it's not in the grid.  Bummer! */
    Error ("translate: Found photon that was not in wind or grid, istat %i\n", where_in_wind (pp->x, &ndomain));
  }

  return (istat);
}




/**********************************************************/
/**
 * @brief      translates the photon from its current position to the
 * 	edge of the wind.
 *
 * @param [in, out] PhotPtr  pp   A photon
 * @return     A status flag indication why the photon stopped
 *
 * @details
 * The routine translate the photon (which is assumed to be
 * in the wind, but not in an active wind cell) to the closest
 * boundary.
 *
 * The photon position is updated in the process.
 *
 *
 * ### Notes ###
 * @bug  There are questions in the comments about why an additonal
 * check is needed as to whther the photon as hit the star. It seems
 * superfluous so someone should check whether this addtional check
 * can be removed.
 *
 **********************************************************/

int
translate_in_space (pp)
     PhotPtr pp;
{
  double ds, delta, s, smax;
  int ndom, ndom_next;
  struct photon ptest;
  int ifail;

  ds = ds_to_wind (pp, &ndom);

  /* For IMPORT, although we have reached the edge of the wind, we may be in a cell that is
   * not really in the wind, so we have to address this situtation here.  The first problem
   * we have though is that we do not know what domain we are in.*/

  if (ndom >= 0 && zdom[ndom].wind_type == IMPORT)
  {
    stuff_phot (pp, &ptest);
    move_phot (&ptest, ds + DFUDGE);    /* So now ptest is at the edge of the wind as defined by the boundary
                                           From here on we should be in the grid  */
    ds += DFUDGE;               //Fix for Bug #592 - we need to keep track of the little DFUDGE we moved the test photon
    /* XXX this is a test.  We check at the start whether we are in the grid */

    if ((ifail = where_in_grid (ndom, ptest.x)) < 0)
    {
    }

    /* XXX this ends the test */



    /* Note there is a possiblity that we reach the other side 
     * of the grid without actually encoutering a
     * wind cell
     */


    if (where_in_wind (ptest.x, &ndom_next) < 0)
    {

      smax = ds_to_wind (&ptest, &ndom_next);   // This is the maximum distance can go in this domain
      s = 0;
      while (s < smax && where_in_wind (ptest.x, &ndom_next) < 0)
      {
        if ((delta = ds_in_cell (ndom, &ptest)) > 0)
        {
          move_phot (&ptest, delta + DFUDGE);
          s += delta + DFUDGE;  // The distance the photon has moved
        }
        else
        {
          break;
        }
      }
      /* So at this point we have either gotten out of the domain or we have found a cell that
       * is actually in the wind or we encoutered the error above
       */
      if (s > 0)
      {
        ds += s - DFUDGE;       /* We are going to add DFUDGE back later */
      }
    }

  }

  move_phot (pp, ds + DFUDGE);


  return (pp->istat);
}




/**********************************************************/
/**
 * @brief      calculates the photon pathlength to the edge of the wind.
 *
 * @param [in, out] PhotPtr  pp   A photon
 * @param [in, out] int *  ndom_current   The current domain
 * @return     The distance to the nearest boundary of the wind and the domain for which
 * 	the boudary applies.
 *
 * @details
 * Python defines the boundaries of the wind in term of the intersection
 * of a biconical flow (defined by an inner and outer windcone) and
 * an inner and outer radius centered on the star.	 If the user is interested
 * in a biconical flow, he/she will most likely set the radii so that they
 * do not truncate the wind.  Similarly, by choosing parameters for the windcones
 * properly, one can assure that one has a completely spherical flow.  However,
 * the user may also perversely set the parameters to that both the conical
 * flow boundaries and the min and maximum radius come into play.
 *
 * Usually, one would
 * efine the two types of The two types of boundaries are
 * sually defined so that one includes the other and so the cones apply when one is
 * ealing with a CV type wind, while the inner and outer radius applies for a spherical
 * odel.  However this routine does not require this to be the case, since it just
 * alculates where the edges are.
 *
 * n any event, if you are inside the wind already ds_to_wind calculates the distance to the edge of the wind.
 * If you are outside, It will also be to the nearest edge.
 *
 * The routine distinguishes between  two basic cases.  Either the photon is already in the wind
 * in which case one looks for the nearest boundary, or the photon is not currently
 * in the wind and you have to translate it until it reaches the the boundary (or
 * VERY_BIG)
 *
 * ### Notes ###
 * There is no guarantee that you will still be in the region defined by the grid.
 *
 * @bug 1802 -ksl - At present this routine for imported models this routine only deals with
 * cylindrical models.  Additionally for imported models we skip all of the
 * of the wind_cones.  This is inefficient, and needs to be corrected for
 * rtheta and spherical models which can easily be handled using wind cones.
 *
 **********************************************************/

double
ds_to_wind (pp, ndom_current)
     PhotPtr pp;
     int *ndom_current;
{
  struct photon ptest, qtest;
  double ds, x, rho, z;
  int ndom;

  stuff_phot (pp, &ptest);

  /* First calculated the distance to the edge of the of
     all of the "computatational domain */

  ds = ds_to_sphere (geo.rmax, &ptest);
  *ndom_current = (-1);
  xxxbound = BOUND_NONE;

  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {
    if ((zdom[ndom].wind_type != IMPORT || (zdom[ndom].wind_type == IMPORT && zdom[ndom].coord_type != CYLIND)) &&
        zdom[ndom].wind_type != CORONA)
    {
      /* Check if the photon hits the inner or outer radius of the wind */
      if ((x = ds_to_sphere (zdom[ndom].rmax, &ptest)) < ds)
      {
        ds = x;
        *ndom_current = ndom;
        xxxbound = BOUND_RMIN;
      }

      if ((x = ds_to_sphere (zdom[ndom].rmin, &ptest)) < ds)
      {
        ds = x;
        *ndom_current = ndom;
        xxxbound = BOUND_RMAX;
      }

      /* Check if the photon hits the inner or outer windcone */

      if ((x = ds_to_cone (&zdom[ndom].windcone[0], &ptest)) < ds)
      {
        ds = x;
        *ndom_current = ndom;
        xxxbound = BOUND_INNER_CONE;
      }
      if ((x = ds_to_cone (&zdom[ndom].windcone[1], &ptest)) < ds)
      {
        ds = x;
        *ndom_current = ndom;
        xxxbound = BOUND_OUTER_CONE;
      }
    }

    /* For this rectangular region we check whether we are in side the grid,
     * which should effectively.  For * an imported region file we may not be
     * inside the wind, since some cells may be empty
     */

    else if (zdom[ndom].wind_type == CORONA || (zdom[ndom].wind_type == IMPORT && zdom[ndom].coord_type == CYLIND))
    {


      x = ds_to_plane (&zdom[ndom].windplane[0], &ptest);
      if (x > 0 && x < ds)
      {
        stuff_phot (pp, &qtest);
        move_phot (&qtest, x);
        rho = sqrt (qtest.x[0] * qtest.x[0] + qtest.x[1] * qtest.x[1]);
        if (zdom[ndom].wind_rho_min <= rho && rho <= zdom[ndom].wind_rho_min)
        {

          ds = x;
          *ndom_current = ndom;
          xxxbound = BOUND_ZMIN;
        }
      }
      x = ds_to_plane (&zdom[ndom].windplane[1], &ptest);
      if (x > 0 && x < ds)
      {
        stuff_phot (pp, &qtest);
        move_phot (&qtest, x);
        rho = sqrt (qtest.x[0] * qtest.x[0] + qtest.x[1] * qtest.x[1]);
        if (zdom[ndom].wind_rho_min <= rho && rho <= zdom[ndom].wind_rho_min)
        {

          ds = x;
          *ndom_current = ndom;
          xxxbound = BOUND_ZMAX;
        }
      }

      x = ds_to_cylinder (zdom[ndom].wind_rho_min, &ptest);
      if (x > 0 && x < ds)
      {
        stuff_phot (pp, &qtest);
        move_phot (&qtest, x);
        z = fabs (qtest.x[2]);
        if (zdom[ndom].zmin <= z && z <= zdom[ndom].zmax)

        {

          ds = x;
          *ndom_current = ndom;
          xxxbound = BOUND_INNER_RHO;
        }
      }

      x = ds_to_cylinder (zdom[ndom].wind_rho_max, &ptest);
      if (x > 0 && x < ds)
      {
        stuff_phot (pp, &qtest);
        move_phot (&qtest, x);
        z = fabs (qtest.x[2]);
        if (zdom[ndom].zmin <= z && z <= zdom[ndom].zmax)
        {

          ds = x;
          *ndom_current = ndom;
          xxxbound = BOUND_OUTER_RHO;
        }
      }

    }
    else
    {
      Error ("ds_to_wind:Do not know how to deal with this combination of coordinate type %d and wind_type %d\n", zdom[ndom].coord_type,
             zdom[ndom].wind_type);
      Exit (0);

    }

  }


  return (ds);
}


/** Added because there were cases where the number
 * of photons passing through a cell with a neglible volume was becoming
 * too large and stopping the program.  This is a bandaide since if this
 * is occurring a lot we should be doing a better job at calculating the
 * volume
 */
int neglible_vol_count = 0;
int translate_in_wind_failure = 0;

/**********************************************************/
/**
 * @brief      translates the photon within a single cell in the wind.
 *
 * @param [in] WindPtr  w   The entire wind
 * @param [in, out] PhotPtr  p   A photon
 * @param [in] double  tau_scat   The depth at which the photon will scatter
 * @param [out] double *  tau   The tau of a resonance
 * @param [out] int *  nres   The resonaance which caused the photon to stop
 * @return     A status indicated whether the photon has stopped for a scattering
 * even of for some other teason
 *
 * @details
 * It calculates and updates the final position of the photon p, the optical depth tau after
 * 	having translated the photon, and, if the photon is scattered, the number of the resonance
 * 	responsible for scattering (or absorbing) the photon bundle.
 *
 * 	These last two quantities are calculated in ds_calculate and simply passed back.
 *
 * 	translate_in_wind returns that status directly to the calling routine.
 *
 * ### Notes ###
 *
 * In addition to moving the photon, the routine also updates some values
 * having to do with the radiation filed both in and out of macro atom
 * mode.
 *
 **********************************************************/
int
translate_in_wind (w, p, tau_scat, tau, nres)
     WindPtr w;                 //w here refers to entire wind, not a single element
     PhotPtr p;
     double tau_scat, *tau;
     int *nres;


{

  int n;
  double smax, s, ds_current;
  int istat;
  int nplasma;
  int ndom, ndom_current;
  int inwind;

  WindPtr one;
  PlasmaPtr xplasma;


/* First verify that the photon is in the grid, and if not
return and record an error */

  if ((p->grid = n = where_in_grid (wmain[p->grid].ndom, p->x)) < 0)
  {
//OLD    if (translate_in_wind_failure < 1000)
//OLD    {
//OLD     if (modes.save_photons)
//OLD       {
//OLD         save_photons (p, "NotInGrid_translate_in_wind");
//OLD       }
//OLD    }
    return (n);                 /* Photon was not in grid */
  }
/* Assign the pointers for the cell containing the photon */

  one = &wmain[n];              /* one is the grid cell where the photon is */
  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];
  ndom = one->ndom;
  inwind = one->inwind;



/* Calculate the maximum distance the photon can travel in the cell */

  if ((smax = ds_in_cell (ndom, p)) < 0)
  {
    return ((int) smax);
  }
  if (one->inwind == W_PART_INWIND)
  {                             /* The cell is partially in the wind */
    s = ds_to_wind (p, &ndom_current);  /* smax is set to be the distance to edge of the wind */
    if (s < smax)
      smax = s;
    s = ds_to_disk (p, 0);      /* the 0 implies ds_to_disk can not return a negative distance */
    if (s > 0 && s < smax)
      smax = s;
  }
  else if (one->inwind == W_IGNORE)
  {
    smax += one->dfudge;
    move_phot (p, smax);
    return (p->istat);

  }
  else if (one->inwind == W_NOT_INWIND)
  {                             /* The cell is not in the wind at all */

    Error ("translate_in_wind: Grid cell %d of photon is not in wind, moving photon %.2e\n", n, smax);
    Error ("translate_in_wind: photon %d position: x %g y %g z %g\n", p->np, p->x[0], p->x[1], p->x[2]);
    move_phot (p, smax);
    return (p->istat);

  }

  if (modes.save_photons)
  {
    Diag ("smax  %10.3e tau_scat %10.3e tau %10.3e\n", smax, tau_scat, *tau);
  }


/* At this point we now know how far the photon can travel in it's current grid cell */

  smax += one->dfudge;          /* dfudge is to force the photon through the cell boundaries. */

/* Set limits the distance a photon can travel.  There are
a good many photons which travel more than this distance without this
limitation, at least in the standard 30 x 30 instantiation.  It does
make small differences in the structure of lines in some cases.
The choice of SMAX_FRAC can affect execution time.*/

  if (smax > SMAX_FRAC * length (p->x))
  {
    smax = SMAX_FRAC * length (p->x);
  }

  /* We now determine whether scattering prevents the photon from reaching the far edge of
     the cell.  calculate_ds calculates whether there are scatterings and makes use of the
     current position of the photon and the position of the photon at the far edge of the
     shell.  It needs a "trial photon at the maximum distance however */


/* Note that ds_current does not alter p in any way */

  ds_current = calculate_ds (w, p, tau_scat, tau, nres, smax, &istat);

  if (p->nres < 0)
    xplasma->nscat_es++;
  if (p->nres > 0)
    xplasma->nscat_res++;



/* OK now we increment the radiation field in the cell, translate the photon and wrap
   things up If the photon is going to scatter in this cell, radiation also reduces
   the weight of the photon due to continuum absorption, e.g. free free */

  if (geo.rt_mode == RT_MODE_MACRO)
  {                             // Macro-method
    /* In the macro-method, b-f and other continuum processes do not reduce the photon
       weight, but are treated as as scattering processes.  Therefore most of what was in
       subroutine radiation can be avoided.
     */

    one = &w[p->grid];
    nplasma = one->nplasma;
    xplasma = &plasmamain[nplasma];

    if (geo.ioniz_or_extract == 1)
    {
      xplasma->ntot++;          // EP 11-19: Moved so only increments during ionisation cycles

      /* For an ionization cycle */
      bf_estimators_increment (one, p, ds_current);

      /*photon weight times distance in the shell is proportional to the mean intensity */
      xplasma->j += p->w * ds_current;

      /* frequency weighted by the weights and distance in the shell.  See eqn 2 ML93 */
      xplasma->ave_freq += p->freq * p->w * ds_current;

    }
  }
  else
  {
    radiation (p, ds_current);
  }

  move_phot (p, ds_current);

  p->nres = (*nres);

  return (p->istat = istat);


}




/**********************************************************/
/**
 * @brief      calculates the distance photon can travel within the cell
 * 	that it is currently in.
 *
 * @param [in, out] int  ndom   The current domain
 * @param [in, out] PhotPtr  p   A photon
 * @return     A distance indicating how far the photon can travel
 *
 * @details
 * The routine is basically just a steering routine
 * and calls ds_in_whatever for various coordinate systems
 *
 * ### Notes ###
 *
 * The other routines are contained in routines like cylindrical.c,
 * rtheta.c etc.
 *
 **********************************************************/

double
ds_in_cell (ndom, p)
     int ndom;
     PhotPtr p;

{

  int n;
  double smax;

  WindPtr one;


/* First verify that the photon is in the grid, and if not
return and record an error */

  if ((p->grid = n = where_in_grid (ndom, p->x)) < 0)
  {
//OLD      if (modes.save_photons)
//OLD   {
//OLD     save_photons (p, "NotInGrid_ds_in_cell");
//OLD   }
    return (n);
  }

/* Assign the pointers for the cell containing the photon */

  one = &wmain[n];              /* one is the grid cell where the photon is */

/* Calculate the maximum distance the photon can travel in the cell */

  if (zdom[ndom].coord_type == CYLIND)
  {
    smax = cylind_ds_in_cell (ndom, p); // maximum distance the photon can travel in a cell
  }
  else if (zdom[ndom].coord_type == RTHETA)
  {
    smax = rtheta_ds_in_cell (ndom, p);
  }
  else if (zdom[ndom].coord_type == SPHERICAL)
  {
    smax = spherical_ds_in_cell (ndom, p);
  }
  else if (zdom[ndom].coord_type == CYLVAR)
  {
    smax = cylvar_ds_in_cell (ndom, p);
  }
  else
  {
    Error ("ds_in_cell: Don't know how to find ds_in_cell in this coord system %d\n", zdom[ndom].coord_type);
    Exit (0);
    return (0);
  }

  return (smax);
}
