/***********************************************************/
/** @file  frame.c
 * @author ksl
 * @date   May, 2020
 *
 * @brief
 * These routines have to do with switching from the observer 
 * to the local frame and back
 *
 * ###Notes###
 *
 * These routines are written so that if p and p_obs point 
 * to the same variable then everything occurs in place
 *
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


/**********************************************************/
/**
 * @brief  Check that a photon is in the desired frame.
 *
 * @param [in] PhotPtr  p      The photon whose frame is to be checked
 * @param [in] frame    frame  The desired frame of the photon
 * @param [in] char    *msg    A comment to print as an error message if the
 *                             photon is in the incorrect frame
 * @return 0 (False) if the photons is in the desired frame, 1 (True) if the
 *  photon is NOT in the desired frame.
 *
 * @details
 *
 * The purpose of this function is to avoid situations where photons are being
 * transformed incorrectly between the local and observer frame.
 *
 * If the photon is in the incorrect frame and save_photons mode is enabled,
 * the photon will be dumped to an external file.
 *
 **********************************************************/

int ncheck_frame = 0;

int
check_frame (p, desired_frame, msg)
     PhotPtr p;
     enum frame desired_frame;
     char *msg;
{
  if (p->frame == desired_frame)
  {
    return (0);
  }
  else if (ncheck_frame < 100)
  {
    Error ("check_frame: %s\n", msg);
    ncheck_frame++;

    if (modes.save_photons)
      save_photons (p, "PhotonInIncorrectFrame");

    return (1);
  }
  else
  {
    return (0);
  }
}



/**********************************************************/
/**
 * @brief      carries out the transformation of a all the quantities
 *      in a photon structure
 *      from the observer (or global) frame to the local (or co-moving)
 *      frame
 *
 * @param [in] PhotPtr  p_in   The photon in the observer frame
 * @param [out] PhotPtr  p_out   The photon in the local frame
 *
 * @return    The routine returns 0, unless there was an error
 *
 *
 * @details
 * The routine copies all of the quantities contain in the p_in to
 * p_out, obtains the velocity at that point in the wind, and then
 * performs an in place transformation from the * global to the local
 * frame of the photon allowing for special
 * relativity.
 *
 *
 * ### Notes ###
 *
 * It p_in and p_out are the same then the transformation will
 * be performed in place.
 *
 **********************************************************/

int
observer_to_local_frame (p_in, p_out)
     PhotPtr p_in, p_out;
{
  WindPtr one;
  int ndom;
  double f;
  double x;
  double v[3], vel;
  double gamma;
  int i, ierr;
  char msg[LINELENGTH];


  sprintf (msg, "observer_to_local_frame: Photon (%d) of type (%d) not in observer frame", p_in->np, p_in->istat);

  ierr = check_frame (p_in, F_OBSERVER, msg);

  /* Initialize the output photon */
  stuff_phot (p_in, p_out);

  /* Calculate the local velocity of the wind at this position */
  one = &wmain[p_in->grid];
  ndom = one->ndom;
  vwind_xyz (ndom, p_in, v);

  vel = dot (p_in->lmn, v);

  if (geo.rel_mode == REL_MODE_LINEAR)
  {
    f = p_out->freq = p_in->freq * (1. - vel / VLIGHT);
    return (ierr);
  }



  // beta=length(v)/VLIGHT;

  gamma = sqrt (1 - (dot (v, v) / (VLIGHT * VLIGHT)));

  f = p_out->freq = p_in->freq * gamma * (1. - vel / VLIGHT);

  x = gamma / VLIGHT * (1.0 - (gamma * vel / ((gamma + 1) * VLIGHT)));

  for (i = 0; i < 3; i++)
  {
    p_out->lmn[i] = f / p_in->freq * (p_in->lmn[i] - x * v[i]);
  }

  p_out->w *= (f / p_in->freq);
  p_out->frame = F_LOCAL;

  return (ierr);
}



/**********************************************************/
/**
 * @brief      carries out the transformation of all the quantities
 *      in a photon structure 
 *      from the local (or co-moving) frame to the observer (or global)
 *      frame
 *
 * @param [in] PhotPtr  p_in   The photon in the local frame                   
 * @param [out] PhotPtr  p_out   The photon in the global frame                  
 *
 * @return    The routine returns 0 unless there was an error       
 *
 *
 * @details
 * The routine copies all of the quantities contain in the p_in to 
 * p_out, obtains the velocity at that point in the wind, and then 
 * performs an in place transformation from the local to the observer
 * frame of the photon allowing for special 
 * relativity.
 *
 *
 * ### Notes ###
 *
 * It p_in and p_out are the same then the transformation will
 * be performed in place.
 *
 **********************************************************/


int
local_to_observer_frame (p_in, p_out)
     PhotPtr p_in, p_out;
{
  WindPtr one;
  int ndom;
  double f;
  double x;
  double v[3], vel;
  double gamma;
  int i;
  char msg[LINELENGTH];
  int ierr;


  sprintf (msg, "local_to_observer_frame: Photon (%d) of type (%d) not_in_local_frame", p_in->np, p_in->istat);


  ierr = check_frame (p_in, F_LOCAL, msg);

  /* Initialize the output photon */
  stuff_phot (p_in, p_out);

  /* Calculate the local velocity of the wind at this position */
  one = &wmain[p_in->grid];
  ndom = one->ndom;
  vwind_xyz (ndom, p_in, v);

  vel = dot (p_in->lmn, v);
  if (geo.rel_mode == REL_MODE_LINEAR)
  {
    f = p_out->freq = p_in->freq / (1. - vel / VLIGHT);
    return (ierr);
  }
  gamma = sqrt (1 - (dot (v, v) / (VLIGHT * VLIGHT)));
  f = p_out->freq = p_in->freq * gamma * (1. + vel / VLIGHT);


/* Need to worry about sign changes, etc. here */
  x = gamma / VLIGHT * (1.0 + (gamma * vel / ((gamma + 1) * VLIGHT)));

  for (i = 0; i < 3; i++)
  {
    p_out->lmn[i] = f / p_in->freq * (p_in->lmn[i] + x * v[i]);
  }

  p_out->w *= (f / p_in->freq);
  p_out->frame = F_OBSERVER;

  return (ierr);
}



/**********************************************************/
/**
 * @brief      carries out the transformation of a all the quantities
 *      in a photon emitted from the surface of the disk
 *      from the local (or co-moving) frame to the observer (or global)
 *      frame
 *
 * @param [in] PhotPtr  p_in   The photon in the local frame                   
 * @param [out] PhotPtr  p_out   The photon in the global frame                  
 *
 * @return    The routine returns 0 unless the was an error         
 *
 *
 * @details
 * The routine copies all of the quantities contain in the p_in to 
 * p_out, obtains the velocity at that point in the wind, and then 
 * performs an in place transformation from the local to the observer
 * frame of the photon allowing for special 
 * relativity.
 *
 *
 * ### Notes ###
 *
 * It p_in and p_out are the same then the transformation will
 * be performed in place.  This is normally the way the 
 * routine is used.   One could have omitted p_out in this case
 * but to make the interface coisistent, both photon structures
 * are included.  
 *
 **********************************************************/

int
local_to_observer_frame_disk (p_in, p_out)
     PhotPtr p_in, p_out;
{
  double f;
  double x;
  double v[3], vel;
  double gamma;
  int i;
  char msg[LINELENGTH];
  int ierr;

  sprintf (msg, "local_to_observer_frame_disk: Photon (%d) of type (%d) not in local frame", p_in->np, p_in->istat);
  ierr = check_frame (p_in, F_LOCAL, msg);

  /* Initialize the output photon */
  stuff_phot (p_in, p_out);


  /* Calculate the local velocity of the disk at this position */
  vdisk (p_in->x, v);
  vel = dot (p_in->lmn, v);


  if (geo.rel_mode == REL_MODE_LINEAR)
  {
    f = p_out->freq = p_in->freq / (1. - vel / VLIGHT);
    return (ierr);
  }

  gamma = sqrt (1 - (dot (v, v) / (VLIGHT * VLIGHT)));
  f = p_out->freq = p_in->freq * gamma * (1. + vel / VLIGHT);

/* Need to worry about sign changes, etc. here */
  x = gamma / VLIGHT * (1.0 + (gamma * vel / ((gamma + 1) * VLIGHT)));

  for (i = 0; i < 3; i++)
  {
    p_out->lmn[i] = f / p_in->freq * (p_in->lmn[i] + x * v[i]);
  }

  p_out->w *= (f / p_in->freq);
  p_out->frame = F_OBSERVER;

  return (ierr);
}


/**********************************************************/
/**
 * @brief      calculate the doppler  shift given the direction of the incoming
 * 	and outgoing photons and the local  velocity of the wind
 *
 * @param [in] PhotPtr  pin   The pre-scattered  photon (used for its direction)
 * @param [in,out] PhotPtr  pout   the scattered  photon (used for its direction)
 * @param [in] double  v[]   the velocity of the wind where the scatter occurred
 * @param [in] int  nres   either the number of the scatter
 * @return    The routine returns 0 unless there was an error
 *
 * pout->freq is updated
 *
 * @details
 * Given the incoming and the outgoing photon direction,
 * doppler calculates the outgoing frequency (to first order in beta).
 * There are two basic cases, depending on whether it was a resonant
 * or a nonresonant scatter.
 *
 * ### Notes ###
 * Called from scatter in resonate.c for a resonant transition  
 * and in extract for disk or wind photons, because you need to 
 * direct photons along a specific line of sight.  
 *
 **********************************************************/

/* Notes
   doppler is currently called in these 3 places

extract.c:159:        doppler (p, &pp, -1);   This is just to Doppler shift a disk photon into the observer frame
extract.c:169:        doppler (p, &pp, pp.nres);      This similar for a wind photon, the implication being that when
                extract is called the photon has not been put in the local frame
resonate.c:1271:    doppler (&pold, p, *nres)

*/

int
doppler (p_in, p_out, nres)
     PhotPtr p_in, p_out;
     int nres;

{
  WindPtr one;
  int ndom;
  double v[3];

  /* Calculate the local velocity of the wind at this position */
  one = &wmain[p_in->grid];
  ndom = one->ndom;
  vwind_xyz (ndom, p_in, v);




  if (nres == -1)               //Electron scattering (SS)
  {                             /*It was a non-resonant scatter */
    p_out->freq = p_in->freq * (1 - dot (v, p_in->lmn) / VLIGHT) / (1 - dot (v, p_out->lmn) / VLIGHT);

  }
  else if (nres > -1 && nres < nlines)
  {                             /* It was a resonant scatter. */
    p_out->freq = lin_ptr[nres]->freq / (1. - dot (v, p_out->lmn) / VLIGHT);
  }
  else if ((nres > NLINES && nres < NLINES + nphot_total + 1) || nres == -2)
    /* It was continuum emission - new comoving frequency has been chosen by
       the matom/kpkt routine, but now need to convert in the same way
       as for lines (SS) */
  {
    /*
       If a 2-level atom run, one should never arrive here.
       Just do a check that all is well - this can be removed eventually (SS)
     */
    if (geo.rt_mode == RT_MODE_2LEVEL)
    {
      Error ("doppler: Not using macro atoms but trying to deexcite one? Abort.\n");
      Exit (0);
    }
    p_out->freq = p_out->freq / (1. - dot (v, p_out->lmn) / VLIGHT);
  }
/* Now do one final check that nothing is awry.  This is another
 * check added by SS that should probably be deleted or done before this point.
 * I have made it fatal so that we will pay attention to it if it occurs. ksl */

  else
  {
    Error ("doppler: nres %d > NLINES + nphot_total %d\n", nres, NLINES + nphot_total);
    Exit (0);
  }

  p_out->frame = F_OBSERVER;

  return (0);

}
