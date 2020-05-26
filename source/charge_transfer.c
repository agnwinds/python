
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
    int n;
    double temp;
    for (n=0;n<n_charge_exchange;n++)
    {
        if (T>=charge_exchange[n].tmin && T<=charge_exchange[n].tmax)
        {
            temp=charge_exchange[n].c*exp(charge_exchange[n].d*T/1e4);
            temp=temp+1.0;
            temp=temp*pow(T/1e4,charge_exchange[n].b);
            charge_exchange_rates[n]=temp*charge_exchange[n].a*1e-9;
        }
        else
        {
            charge_exchange_rates[n]=0.0; //Set the rate to zero
        }
    }
    return(0);
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


  compute_ch_trans_coeffs (t_e);      //Calculate the DI rate coefficients for this cell


  for (n = 0; n < n_charge_exchange; n++)
  {
    if (charge_exchange[n].energy_defect>0) //This is a heating process  
    {
      x += 0.0;                 
    }
    else                        //Multiply the rate coeff by the density of electrons and ions, by the ionization potential (in eV), the volume and convert to ergs
    {
      x += xplasma->vol * charge_exchange_rates[n] * xplasma->density[charge_exchange[n].nion1] * xplasma->density[charge_exchange[n].nion2] * -1.0 *charge_exchange[n].energy_defect;
    }
  }
  return (x);
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


  nplasma = one->nplasma;       //Get the correct plasma cell related to this wind cell
  xplasma = &plasmamain[nplasma];       //copy the plasma structure for that cell to local variable
  x = 0;                        //zero the luminosity


  compute_ch_trans_coeffs (t_e);      //Calculate the DI rate coefficients for this cell


  for (n = 0; n < n_charge_exchange; n++)
  {
    if (charge_exchange[n].energy_defect<0) //This is a heating process  
    {
      x += 0.0;                 
    }
    else                        //Multiply the rate coeff by the density of electrons and ions, by the ionization potential (in eV), the volume and convert to ergs
    {
      x += xplasma->vol * charge_exchange_rates[n] * xplasma->density[charge_exchange[n].nion1] * xplasma->density[charge_exchange[n].nion2] *charge_exchange[n].energy_defect;
    }
  }
  return (x);
}


