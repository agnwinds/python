
/***********************************************************/
/** @file  reposition.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief  Reposition a photon so it goes across a cell boundary
 *
 * @bug This is too short a routine to be in a separate file.  It
 * should be incorporated into trans_phot or extract which conatin the
 * calling routines.
 ***********************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


/**********************************************************/
/** @name      reposition
 * @brief      p) attempts to assure that a photon is not scattered 
 * 	a second time inappropriately by the same transition
 *
 * @param [in,out] PhotPtr  p   A photons
 * @return    Normally returns 0, but returns a negative number
 * if p is not in the wind in the domain it is supposed to be in
 *
 * @details
 * For resonant scatters, the routine moves the photon by 
 * a distance DFUDGE.  For non-resonant scattering the routine 
 * simple returns
 *
 * ### Notes ###
 * @bug This seems like an error.  We have calculated cell specific
 * values of dfudge, but this routine used the global DFUDGE to
 * move the photon though a cell wall.
 *
 **********************************************************/

int
reposition (p)
     PhotPtr p;

{

  int n;


  if (p->nres < 0)
  {                             // Do nothing for non-resonant scatters
    return (0);
  }

  if ((p->grid = n = where_in_grid (wmain[p->grid].ndom, p->x)) < 0)
  {
    Error ("reposition: Photon not in grid when routine entered %d \n", n);
    return (n);                 /* Photon was not in wind */
  }


  move_phot (p, DFUDGE);
  return (0);


}
