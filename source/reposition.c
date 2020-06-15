
/***********************************************************/
/** @file  reposition.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief  Reposition a photon so it goes across a cell boundary
 *
 * ### Programming Notes ###
 *
 * This is a very short function and may not justify being in its own file. The
 * (only) function reposition could be moved into extract.c or trans_phot.c.
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


/**********************************************************/
/**
 * @brief      p) attempts to assure that a photon is not scattered
 * 	a second time inappropriately by the same transition
 *
 * @param [in,out] PhotPtr  p   A photons
 * @return    Normally returns 0, but returns a negative number
 * if p is not in the wind in the domain it is supposed to be in
 *
 * @details
 * For resonant scatters, the routine moves the photon by
 * a distance dfudge. For non-resonant scattering the routine
 * simply returns.
 *
 * ### Notes ###
 *
 **********************************************************/

int
reposition (PhotPtr p)
{
  int n;
  double s, s_disk, s_star;
  int hit_disk;

  if (p->nres < 0)
    return (0);                 /* Do nothing for non-resonant scatters */

  if ((p->grid = n = where_in_grid (wmain[p->grid].ndom, p->x)) < 0)
  {
    Error ("reposition: Photon not in grid when routine entered %d \n", n);
    return (n);                 /* Photon was not in wind */
  }

  s = wmain[p->grid].dfudge;

  if (geo.disk_type != DISK_NONE)
  {
    s_disk = ds_to_disk (p, 1, &hit_disk);      // Allow negative values
    if (s_disk > 0 && s_disk < s)
    {
      s = 0.1 * s_disk;
    }
  }
  s_star = ds_to_sphere (geo.rstar, p);
  if (s_star > 0 && s_star < s)
  {
    s = 0.1 * s_star;
  }

  move_phot (p, s);

  return (0);
}

/* The next routine has been removed because it should not longer be necessary
   with the modifcations to repostion above.  It should be noted that this routine
   as written would fail with vertically extended disks, so if it is needd then
   it needs to be modified to reflect this.  A better approach prbably is to work
   on reposition above.  ksl - 200514
*/

//OLD /* ************************************************************************* */
//OLD /**
//OLD  * @brief           Reposition a photon which was lost due to dfudge pushing
//OLD  *                  the photon into the disk
//OLD  *
//OLD  * @param[in,out]   PhotPtr   p     The photon to be repositioned
//OLD  *
//OLD  * @return          void
//OLD  *
//OLD  * @details
//OLD  *
//OLD  * For resonant scatters, this function will push the photon a distance away
//OLD  * from the location where it previously interacted to assure that the photon
//OLD  * does not interact with the same resonance twice.
//OLD  *
//OLD  * Created in response to issue #584 on the GitHub repository. The purpose of
//OLD  * this function is to calculate the distance to the surface of the accretion
//OLD  * disc and to make it some distance towards the disc instead of dfudge, which
//OLD  * previously pushed the photon through the disc plane accidentally.
//OLD  *
//OLD  * ************************************************************************** */

//OLD void
//OLD reposition_lost_disk_photon (PhotPtr p)
//OLD {
//OLD   double smax;
//OLD 
//OLD   if (p->nres < 0)
//OLD     return;                     /* Do nothing for non-resonant scatters */
//OLD 
//OLD   if ((p->grid = where_in_grid (wmain[p->grid].ndom, p->x)) < 0)
//OLD   {
//OLD     Error ("%s:%s(%i): Photon not in grid\n", __FILE__, __func__, __LINE__);
//OLD     return;                     /* Photon was not in wind */
//OLD   }
//OLD 
//OLD   smax = -p->x[2] / p->lmn[2] * 0.999;  // Move some distance toward the disc
//OLD   move_phot (p, smax);
//OLD }
