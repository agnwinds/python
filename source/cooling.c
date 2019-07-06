
/***********************************************************/
/** @file  cooling.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief  File containg all (most) cooling mechanisms
 *
 ***********************************************************/




#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"



/**********************************************************/
/**
 * @brief      calculate the total cooling of a cell at a temperature t
 *
 * @param [in,out] PlasmaPtr  xxxplasma   A cell in the wind
 * @param [in] double  t   The temperature where cooling is calculated
 * @return     The total amount of cooling in a plasma cell,
 *
 * @details
 * The routine calculates cooling of a wind cell, storing the results
 * of various different processes in various parameters of the plasma
 * cell.
 *
 * ### Notes ###
 * @bug There appear to be redundant calls here to xtotal_emission (which
 * include calls to total_fb, and to total_fb in this routine.  This makes
 * it difficult to understand whether something is being counted twice.
 *
 **********************************************************/
double
cooling (xxxplasma, t)
     PlasmaPtr xxxplasma;
     double t;
{

  xxxplasma->t_e = t;


  if (geo.adiabatic)
  {
    if (wmain[xxxplasma->nwind].div_v >= 0.0)
    {
      /* This is the case where we have adiabatic cooling - we want to retain the old behaviour,
         so we use the 'test' temperature to compute it. If div_v is less than zero, we don't do
         anything here, and so the existing value of adiabatic cooling is used - this was computed
         in wind_updates2d before the call to ion_abundances. */
      xxxplasma->cool_adiabatic = adiabatic_cooling (&wmain[xxxplasma->nwind], t);
    }
  }

  else
  {
    xxxplasma->cool_adiabatic = 0.0;
  }


  /*81c - nsh - we now treat DR cooling as a recombinational process - still unsure as to how to treat emission, so at the moment
     it remains here */

  xxxplasma->cool_dr = total_fb (&wmain[xxxplasma->nwind], t, 0, VERY_BIG, FB_REDUCED, INNER_SHELL);

  /* 78b - nsh adding this line in next to calculate direct ionization cooling without generating photons */

  xxxplasma->cool_di = total_di (&wmain[xxxplasma->nwind], t);

  /* 70g compton cooling calculated here to avoid generating photons */

  xxxplasma->cool_comp = total_comp (&wmain[xxxplasma->nwind], t);

  /* we now call xtotal emission which computes the cooling rates for processes which can, in principle, make photons. */

  xxxplasma->cool_tot =
    xxxplasma->cool_adiabatic + xxxplasma->cool_dr + xxxplasma->cool_di +
    xxxplasma->cool_comp + xtotal_emission (&wmain[xxxplasma->nwind], 0., VERY_BIG);

  return (xxxplasma->cool_tot);
}






/**********************************************************/
/**
 * @brief      calculate the cooling of a single cell
 * 			due to free free, bound - boound and recombination
 * 			processes.
 *
 * @param [in out] WindPtr  one   A singel wind cell
 * @param [in] double  f1   The minimum frequecny
 * @param [in] double  f2   The maximum frequecny
 * @return     The cooling rate
 *
 * @details
 * xtotal_emission gives the total enery loss due to photons.  It does
 * 	not include other coooling sources, e. g. adiabatic expansion.
 *
 * ### Notes ###
 *
 * It returns the total cooling, but also stores the luminosity due
 * to various types of emssion, e.g ff and  lines into the
 * Plasms cells. The fb cooling calculated here is *not* equal to
 * the fb lumniosity and so this value is stored in cool_rr.
 *
 * @bug The call to this routine was changed when PlasmaPtrs
 * were introduced, but it appears that the various routines
 * that were called were not changed. This needs to be fixed for
 * consistency
 *
 **********************************************************/

double
xtotal_emission (one, f1, f2)
     WindPtr one;               /* WindPtr to a specific cell in the wind */
     double f1, f2;             /* The minimum and maximum frequency over which the emission is
                                   integrated */
{
  double t_e;
  int nplasma;
  double cooling;
  PlasmaPtr xplasma;

  cooling = 0.0;
  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];

  t_e = xplasma->t_e;           // Change so calls to total emission are simpler

  if (f2 < f1)
  {
    xplasma->cool_tot = xplasma->lum_lines = xplasma->lum_ff = xplasma->cool_rr = 0;    //NSH 1108 Zero the new cool_comp variable NSH 1101 - removed
  }
  else
  {
    if (geo.rt_mode == RT_MODE_MACRO)   //Switch for macro atoms (SS)
    {
      xplasma->cool_rr = total_fb_matoms (xplasma, t_e, f1, f2) + total_fb (one, t_e, f1, f2, FB_REDUCED, OUTER_SHELL); //outer shellrecombinations
      //The first term here is the fb cooling due to macro ions and the second gives
      //the fb cooling due to simple ions.
      //total_fb has been modified to exclude recombinations treated using macro atoms.
      //Note: This the fb_matom call makes no use of f1 or f2. They are passed for
      //now in case they should be used in the future. But they could
      //also be removed.
      // (SS)
      cooling = xplasma->cool_rr;
      xplasma->lum_lines = total_bb_cooling (xplasma, t_e);
      cooling += xplasma->lum_lines;
      /* total_bb_cooling gives the total cooling rate due to bb transisions whether they
         are macro atoms or simple ions. */
      xplasma->lum_ff = total_free (one, t_e, f1, f2);
      cooling += xplasma->lum_ff;


    }
    else                        //default (non-macro atoms) (SS)
    {
      /*The line cooling is equal to the line emission */
      cooling = xplasma->lum_lines = total_line_emission (one, f1, f2);
      /* The free free cooling is equal to the free free emission */
      cooling += xplasma->lum_ff = total_free (one, t_e, f1, f2);
      /*The free bound cooling is equal to the recomb rate x the electron energy - the boinding energy - this is computed
         with the FB_REDUCED switch */
      cooling += xplasma->cool_rr = total_fb (one, t_e, f1, f2, FB_REDUCED, OUTER_SHELL);       //outer shell recombinations


    }
  }


  return (cooling);


}






/**********************************************************/
/**
 * @brief      determines the amount of
 * 	adiabatic cooling in a cell, in units of luminosity.
 *
 * @param [in] WindPtr  one   pointer to wind cell
 * @param [in] double  t   electron temperature
 * @return The adiabiatic cooling (or heating) of a cell
 *
 * @details
 * Adiabatic cooling is the amount of PdV work done by
 * a fluid element, per unit time dt. Thus it is equal to P dV/dt.
 * dV/dt is given by * the volume * div v.  As is evident from
 * this what we call adiabatic cooling can heat the wind if div
 * v is negative.
 *
 * ### Notes ###
 *
 * In calculating the adiabatic cooling we include the pressure
 * from all particles.  Adiabatic coolling due to the radiation
 * pressure is not considered.
 *
 * The routine does not populate xplasma->cool_adiabatic.
 *
 * Note also that this function should only be called
 * if geo.adiabatic == 1, in which case it populates
 * xplasma->cool_adiabatic. This is used in heating and cooling
 * balance. We also use it as a potential destruction choice for
 * kpkts in which case the kpkt is thrown away by setting its istat
 * to P_ADIABATIC.
 *
 * @bug The statements here about not calling this routine are
 * odd.  There is no reason that the routine cannot be called
 * any time since it does not populate cool_adiabatic.  And
 * what really shold happen if we do not want to use cool_adiabatic
 * is that we should check when we want to calculate the total cooling
 * that it is not done.
 *
 **********************************************************/

double
adiabatic_cooling (one, t)
     WindPtr one;
     double t;
{
  double cooling;
  int nplasma, nion;
  double nparticles;
  PlasmaPtr xplasma;

  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];

  nparticles = xplasma->ne;

  /* loop over all ions as they all contribute to the pressure */
  for (nion = 0; nion < nions; nion++)
  {
    nparticles += xplasma->density[nion];
  }

  cooling = nparticles * BOLTZMANN * t * xplasma->vol * one->div_v;

  return (cooling);
}


/**********************************************************/
/**
 * @brief      determines the amount of
 * 	shock related heating n a cell, in units of luminosity.
 *
 * @param [in] WindPtr  one   pointer to wind cell
 * @return The amount os shock-heating in a cell according
 * to a specific formula.  
 *
 * @details
 *
 * Calculates the shock-related heating according to
 * a simple formula provided by Lee Hartmann.  The
 * formula is supposed to represent
 *
 *
 * H = C * (r/rstar)**-4 
 *
 * where c is given by
 *
 *
 * H_tot/(4PI rstar**3)
 *
 * and H_tot is
 *
 * eta*Mdot_wind*v_infty**3
 *
 * that is a fraction of the kinetic energy of the wind
 *
 *
 *
 * ### Notes ###
 *
 * This is implemented for the FU Ori project and may not be
 * appropriate in other situations.
 *
 *
 * The heating is multiplied by the volume of the plasma in
 * the cell, and so the units are ergs/s.  
 *
 * This is implemented analagously to adiabatic cooling
 *
 **********************************************************/


double
shock_heating (one)
     WindPtr one;
{
  int nplasma;
  double x, r;
  PlasmaPtr xplasma;

  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];

  r = length (one->xcen) / geo.rstar;

  x = geo.shock_factor / (r * r * r * r);

  x *= xplasma->vol;

  return (x);
}

/**********************************************************/
/**
 * @brief      calculate the cooling rate of the entire
 * wind between freqencies f1 and f2 (including non-radiative processes
 *
 * @return     The total cooling rate of the wind
 *
 * Various parameters in geo having to do with the total cooling
 * for the for specific processes are also populated.
 *
 * @details
 * This routine cycles through all of the wind cells and
 * calculates the cooling rate in each cell.
 *
 * ### Notes ###
 * This is not the same as the luminosty of the wind, because it
 * includes non-radiative processes, such as adiabatic and Compton
 * cooling.
 *
 * Non-radiative process include adiabatic cooling, Compton cooling,
 * dielectronic recombination, direct ionization, and in some cases
 * non-thermal or shock heating)
 *
 * Adiabatic cooling is summed separately for cells where the
 * adiabatic cooling is positive (tends to cool the wind) and
 * negative (tends to heat the wind.)
 *
 **********************************************************/

double
wind_cooling ()
{
  double cool, lum_lines, cool_rr, lum_ff, cool_comp, cool_dr, cool_di, cool_adiab, heat_adiab;
  double nonthermal;
  int n;
  double x;
  int nplasma;


  cool = lum_lines = cool_rr = lum_ff = cool_comp = 0;
  cool_dr = cool_di = cool_adiab = heat_adiab = 0;      //1108 NSH Zero the new counter 1109 including DR counter 1408 and the DI counter
  nonthermal = 0;

  for (n = 0; n < NDIM2; n++)
  {

    if (wmain[n].vol > 0.0)
    {
      nplasma = wmain[n].nplasma;
      cool += x = cooling (&plasmamain[nplasma], plasmamain[nplasma].t_e);
      lum_lines += plasmamain[nplasma].lum_lines;
      cool_rr += plasmamain[nplasma].cool_rr;
      lum_ff += plasmamain[nplasma].lum_ff;
      cool_comp += plasmamain[nplasma].cool_comp;       //1108 NSH Increment the compton luminosity for that cell.
      cool_dr += plasmamain[nplasma].cool_dr;   //1109 NSH Increment the DR luminosity for the cell.
      cool_di += plasmamain[nplasma].cool_di;   //1408 NSH Increment the DI luminosity for the cell.

      if (geo.adiabatic)        //Caculate the total adiatbaic heating/cooling separating these into two variables
      {

        if (plasmamain[nplasma].cool_adiabatic >= 0.0)
        {
          cool_adiab += plasmamain[nplasma].cool_adiabatic;
        }
        else
        {
          heat_adiab += plasmamain[nplasma].cool_adiabatic;
        }
      }

      else
      {
        cool_adiab = 0.0;
      }

      /* Caculate the non-thermal heating (for FU Ori models with extra wind heating */
      if (geo.nonthermal)
      {
        nonthermal += plasmamain[nplasma].heat_shock;
      }



      if (x < 0)
        Error ("wind_cooling: xtotal emission %8.4e is < 0!\n", x);

      if (recipes_error != 0)
        Error ("wind_cooling: Received recipes error on cell %d\n", n);
    }
  }


  /* Store the results into the geo structure */


  geo.lum_lines = lum_lines;
  geo.cool_rr = cool_rr;
  geo.lum_ff = lum_ff;
  geo.cool_comp = cool_comp;    //The compton luminosity
  geo.cool_dr = cool_dr;        //the DR luminosity
  geo.cool_di = cool_di;        //the DI luminosity 
  geo.cool_adiabatic = cool_adiab;
  geo.heat_adiabatic = heat_adiab;
  geo.heat_shock = nonthermal;


  return (cool);
}
