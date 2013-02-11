/**************************************************************************
                    Southampton University
                                                                                                   
                                                                                                   
  Synopsis:

The routines in this file all have to do with dielectronic recombination.


  Description:

The first routine is compute_dr_coeffs
It takes as its input a temperature, and it populates the fd_coeff array.
This is an initial test, and it will be called every time a cell needs a set
of DR coefficients. It could well be quicker in the future to generate one,
temperature split array, and then sample it. Lets see...



                                                                                                   
  Arguments:  
                                                                                                   
                                                                                                   
  Returns:
                                                                                                   
  Notes:
                                                                                                   
                                                                                                   

                                                                                                   
  History:
	11aug	nsh	Began work

                                                                                                   
 ************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "atomic.h"
#include "python.h"
#include "recipes.h"




/**************************************************************************
                    Space Telescope Science Institute
                                                                                                   
                                                                                                   
  Synopsis: compute_dr_coeffs returns the volumetric dielectronic rate
 	coefficients for a given temperature. 
                                                                                                   
  Description:
                                                                                                   
  Arguments:  
	temperature
                                                                                                   
                                                                                                   
  Returns:
	nothing, but populates the array dr_coeffs
                                                                                                   
  Notes:
 
                                                                                                   
                                                                                                   
                                                                                                   
  History:
	11sep	nsh	Written as part of python70 effort to incorporate DR.
                                                                                                   
 ************************************************************************/

int
compute_dr_coeffs (temp)
     double temp;
{
int n,n1,n2;
for (n=0;n<nions;n++)
	{	
	if (ion[n].drflag==0)
		{
		dr_coeffs[n]=0.0;
//		printf ("DDDDDD ion %i has no DR coefficients\n",n);
		}
	else
		{	
		n1=ion[n].nxdrecomb;
		dr_coeffs[n]=0.0;
		for (n2=0;n2<drecomb[n1].nparam;n2++)
			{
			dr_coeffs[n]+=(drecomb[n1].c[n2]*exp(-1*(drecomb[n1].e[n2]/temp)));
			}
		dr_coeffs[n]*=pow(temp,-1.5);
//		printf ("DDDDDD ion %i has DR coefficient %e\n",n,dr_coeffs[n]);
		}
	}
return (0);
}


/**************************************************************************
                    Space Telescope Science Institute
                                                                                                   
                                                                                                   
  Synopsis: total_dr calculates the total luminosity from DR. 
                                                                                                   
  Description:
                                                                                                   
  Arguments:  
	pointer to grid cell we are interested in	
	temperature
	
                                                                                                   
                                                                                                   
  Returns:
  	the total luminosity of this cell due to dielectronic recombinaions
                                                                                                   
  Notes:
 
                                                                                                   
                                                                                                   
                                                                                                   
  History:
	11sep	nsh	Written as part of python70 effort to incorporate DR. Initally we are just doing a ROM calculation by multiplying the volumetric rate be the ion density, the electron density and mean eelectron energy.
                                                                                                   
 ************************************************************************/

   double total_dr(one,t_e)
	WindPtr one; 	// Pointer to the current wind cell - we need the cell volume, this is not in the plasma structure
	double t_e;   	//Current electron temperature of the cell
{
	double x;   	//The returned variable
	double meanv, meanke;  //The mean velocity and kinetic energy of electrons in the cell
	int nplasma; 	//The cell number in the plasma array
	PlasmaPtr xplasma;   //pointer to the relevant cell in the plasma structure
	int n;  //loop pointers


	nplasma=one->nplasma;  //Get the correct plasma cell related to this wind cell
	xplasma=&plasmamain[nplasma];   //copy the plasma structure for that cell to local variable
	x=0; //zero the luminosity

	meanv=pow((2*BOLTZMANN*t_e/MELEC),0.5);
	meanke=0.5*MELEC*meanv*meanv;

for (n=0;n<nions;n++)
	{	
	if (ion[n].drflag==0)   //We have no DR for this ion.
		{
		dr_coeffs[n]=0.0;
//		printf ("DDDDDD ion %i has no DR coefficients\n",n);
		}
	else
		{	
		x += one->vol * xplasma->ne * xplasma->density[n] * dr_coeffs[n] * meanke;
//		printf ("DDDDDD ion %i has DR lum of %e\n",n, one->vol * xplasma->ne * xplasma->density[n] * dr_coeffs[n] * meanke);
		}
	}
return (x);
}
