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
	02jul	ksl	Removed all references to the wind cell.
                                                                                                   
 ************************************************************************/

double
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
