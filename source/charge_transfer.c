
/***********************************************************/
/** @file  charge_transfer.c
 * @author nsh
 * @date   May, 2020
 *
 * @brief  charge_transfer contains routines to compute charge transfer rates
 *
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "atomic.h"
#include "python.h"
#include "recipes.h"


/**********************************************************/
/**
 * @brief      computes the rates due to charge exchange
 *
 * @param [in] double  T  temperature in K of the cell in question
 * @return     0 on success
 *
 * @details
 * 
 *
 * ### Notes ###
 * ??? NOTES ???
 *
 **********************************************************/

int
compute_ch_trans_coeffs (double T)
{
  int n, m;
  double temp;
  double T_used;
/*  for (n = 0; n < n_charge_exchange; n++)
  {
    if (T >= charge_exchange[n].tmin && T <= charge_exchange[n].tmax)
    {
      temp = charge_exchange[n].c * exp (charge_exchange[n].d * T / 1e4);
      temp = temp + 1.0;
      temp = temp * pow (T / 1e4, charge_exchange[n].b);
      charge_exchange_rates[n] = temp * charge_exchange[n].a * 1e-9;
      printf ("BOOM %i %i %e\n", ion[charge_exchange[n].nion2].z, ion[charge_exchange[n].nion2].istate, charge_exchange_rates[n]);
    }
    else
    {
      charge_exchange_rates[n] = 0.0;   //Set the rate to zero
    }
  }*/

  for (n = 0; n < nions; n++)
  {
    if (ion[n].z != 1 && ion[n].istate > 1)     //No rates for hydrogen or neutral ions
    {
      if (ion[n].n_ch_ex < 0)   //We dont have a proper rate, so use the approximation
      {
        charge_exchange_recomb_rates[n] = 1.92e-9 * (ion[n].istate - 1);        //Approxmation - multiplied by charge 
      }
      else
      {
        m = ion[n].n_ch_ex;
        if (T < charge_exchange[m].tmin)
          T_used = charge_exchange[m].tmin;
        else if (T > charge_exchange[m].tmax)
          T_used = charge_exchange[m].tmax;
        else
          T_used = T;
//        if (T >= charge_exchange[m].tmin && T <= charge_exchange[m].tmax)

//        {
        temp = charge_exchange[m].c * exp (charge_exchange[m].d * T_used / 1e4);
        temp = temp + 1.0;
        temp = temp * pow (T_used / 1e4, charge_exchange[m].b);
        charge_exchange_recomb_rates[n] = temp * charge_exchange[m].a * 1e-9;
//          printf ("IN  TEMP %e %e %e ", T, charge_exchange[m].tmin, charge_exchange[m].tmax);
//        }
        //       else
//        {
//          charge_exchange_recomb_rates[n] = 0.0;        //Set the rate to zero
//          printf ("OUT TEMP %e %e %e ", T, charge_exchange[m].tmin, charge_exchange[m].tmax);

//        }
      }
    }
    else
    {
      charge_exchange_recomb_rates[n] = 0.0;    //Set the rate to zero
    }
    printf ("BLAH %i %i %i %e\n", ion[n].z, ion[n].istate, ion[n].n_ch_ex, charge_exchange_recomb_rates[n]);
  }

  //Now we make ionization rates - do it in a slightly different way
  for (n = 0; n < n_charge_exchange; n++)
  {
    if (ion[charge_exchange[n].nion2].z == 1)   //The recombining ion is hydrogen so we have an ionization rate to make
    {
//      if (T >= charge_exchange[n].tmin && T <= charge_exchange[n].tmax)
//      {
      if (T < charge_exchange[n].tmin)
        T_used = charge_exchange[n].tmin;
      else if (T > charge_exchange[n].tmax)
        T_used = charge_exchange[n].tmax;
      else
        T_used = T;

      temp = charge_exchange[n].c * exp (charge_exchange[n].d * T_used / 1e4);
      temp = temp + 1.0;
      temp = temp * pow (T_used / 1e4, charge_exchange[n].b);
      charge_exchange_ioniz_rates[n] = temp * charge_exchange[n].a * 1e-9;
      charge_exchange_ioniz_rates[n] = charge_exchange_ioniz_rates[n] * exp (-1. * charge_exchange[n].delta_e_ovr_k * 1e4 / T);
//      }
//      else
//      {
//        charge_exchange_ioniz_rates[n] = 0.0;   //Set the rate to zero
//      }
    }
    else
    {
      charge_exchange_ioniz_rates[n] = 0.0;     //Set the rate to zero          
    }
    printf ("BLAH ioniz %i %i %e %e %e %e %e\n", ion[charge_exchange[n].nion1].z, ion[charge_exchange[n].nion1].istate,
            charge_exchange_ioniz_rates[n], charge_exchange[n].delta_e_ovr_k, T, exp (-1. * charge_exchange[n].delta_e_ovr_k * 1e4 / T),
            temp * charge_exchange[n].a * 1e-9);
  }



//  exit (0);
  return (0);
}



double
total_ch_ex (one, t_e)
     WindPtr one;               // Pointer to the current wind cell - we need the cell volume, this is not in the plasma structure
     double t_e;                //Current electron temperature of the cell

{
  double x;                     //The returned variable
  int nplasma;                  //The cell number in the plasma array
  PlasmaPtr xplasma;            //pointer to the relevant cell in the plasma structure
  int n;                        //loop pointers


  nplasma = one->nplasma;       //Get the correct plasma cell related to this wind cell
  xplasma = &plasmamain[nplasma];       //copy the plasma structure for that cell to local variable
  x = 0;                        //zero the luminosity


  compute_ch_trans_coeffs (t_e);        //Calculate the DI rate coefficients for this cell


  for (n = 0; n < n_charge_exchange; n++)
  {
    if (charge_exchange[n].energy_defect > 0)   //This is a heating process  
    {
      x += 0.0;
    }
    else                        //Multiply the rate coeff by the density of electrons and ions, by the ionization potential (in eV), the volume and convert to ergs
    {
      x +=
        xplasma->vol * charge_exchange_ioniz_rates[n] * xplasma->density[charge_exchange[n].nion1] *
        xplasma->density[charge_exchange[n].nion2] * -1.0 * charge_exchange[n].energy_defect;
    }
  }
  return (0.0);
}

double
ch_ex_heat (one, t_e)
     WindPtr one;               // Pointer to the current wind cell - we need the cell volume, this is not in the plasma structure
     double t_e;                //Current electron temperature of the cell

{
  double x;                     //The returned variable
  int nplasma;                  //The cell number in the plasma array
  PlasmaPtr xplasma;            //pointer to the relevant cell in the plasma structure
  int n;                        //loop pointers
  double nh1, nh2;


  nplasma = one->nplasma;       //Get the correct plasma cell related to this wind cell
  xplasma = &plasmamain[nplasma];       //copy the plasma structure for that cell to local variable
  x = 0;                        //zero the luminosity


  compute_ch_trans_coeffs (t_e);        //Calculate the DI rate coefficients for this cell

  for (n = 0; n < nions; n++)
  {
    if (ion[n].z == 1 && ion[n].istate == 1)
    {
      nh1 = xplasma->density[n];
    }
    if (ion[n].z == 1 && ion[n].istate == 2)
    {
      nh2 = xplasma->density[n];
    }
  }


  for (n = 0; n < nions; n++)
  {
    if (ion[n].z != 1 && ion[n].istate > 1)     //No rates for hydrogen or neutral ions
    {
      if (ion[n].n_ch_ex < 0)   //We dont have a proper rate, so use the approximation
      {
        x += xplasma->vol * charge_exchange_recomb_rates[n] * nh1 * xplasma->density[n] * 2.86 * (ion[n].istate - 1) * EV2ERGS;
        printf ("BOOM1 z %i istate %i rate %e defect %e nh1 %e density %e x %e\n", ion[n].z, ion[n].istate, charge_exchange_recomb_rates[n],
                2.86 * (ion[n].istate - 1) * EV2ERGS, nh1, xplasma->density[n], x);
      }
      else
      {
        x += xplasma->vol * charge_exchange_recomb_rates[n] * nh1 * xplasma->density[n] * charge_exchange[ion[n].n_ch_ex].energy_defect;
        printf ("BOOM2 z %i istate %i rate %e defect %e nh1 %e density %e x %e\n", ion[n].z, ion[n].istate, charge_exchange_recomb_rates[n],
                charge_exchange[ion[n].n_ch_ex].energy_defect / EV2ERGS, nh1, xplasma->density[n], x);
      }
    }
  }
  for (n = 0; n < n_charge_exchange; n++)
    if (ion[charge_exchange[n].nion2].z == 1)   //A hydrogen recomb - metal ionization rate
    {
      x +=
        xplasma->vol * charge_exchange_ioniz_rates[n] * nh2 * xplasma->density[charge_exchange[n].nion1] * charge_exchange[n].energy_defect;
      printf ("BOOM3 z %i istate %i rate %e defect %e nh2 %e density %e x %e\n", ion[charge_exchange[n].nion1].z,
              ion[charge_exchange[n].nion1].istate, charge_exchange_ioniz_rates[n], charge_exchange[n].energy_defect / EV2ERGS, nh2,
              xplasma->density[charge_exchange[n].nion1], x);
    }
  return (x);
}
