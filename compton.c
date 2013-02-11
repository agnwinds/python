


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  These are the routines to handle compton scattering.

  Description:	There are two routines in here:
		kappa_compton (xplasma,freq) this calculates the opacity in the cell xplasma due to compton cooling

  Arguments: 		


  Returns:

  Notes:

  History:
2011	nsh	Coded as part of the effort to include compton scattering in August 2011 at Southampton.

 ************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "atomic.h"
#include "python.h"

/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  kappa_compton computes the opacity in the cell due to compton heating.

  Description:	

  Arguments:   xplasma - pointer to the current plasma cell
		freq - the frequency of the current photon packet being tracked

  Returns:   kappa - the compton opacity for the cell.

  Notes:   This implements the equation kappa=(sigmaT*ne*freq)/(me*c^2)
            Initally it is called by radiation

  History:
2011	nsh	Coded as part of the effort to include compton scattering in August 2011 at Southampton.

 ************************************************************************/


	double kappa_comp(xplasma,freq)
	PlasmaPtr xplasma;   // Pointer to current plasma cell
	double freq;  // Frequency of the current photon being tracked
{
	double x;  // The opacity of the cell by the time we return it.
	x=(THOMPSON * H)/(MELEC * C * C); //Calculate the constant
	x*= xplasma->ne * freq; //Multiply by cell electron density and frequency of the packet.
	return (x);
}

/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  total_comp computes the cooling in the cell due to compton cooling.

  Description:	

  Arguments:   one - pointer to the current wind cell
		t_e - electron temperature of the current cell

  Returns:   the compton luminosity for the cell.

  Notes:   This implements the equation C=16piVsigmatJ(kTe/me)

  History:
2011	nsh	Coded as part of the effort to include compton scattering in August 2011 at Southampton.

 ************************************************************************/

	double total_comp(one,t_e)
	WindPtr one; 	// Pointer to the current wind cell - we need the cell volume, this is not in the plasma structure
	double t_e;   	//Current electron temperature of the cell
{
	double x;   	//The returned variable
	int nplasma; 	//The cell number in the plasma array
	PlasmaPtr xplasma;   //pointer to the relevant cell in the plasma structure


	nplasma=one->nplasma;  //Get the correct plasma cell related to this wind cell
	xplasma=&plasmamain[nplasma];   //copy the plasma structure for that cell to local variable

	x=16. * PI * THOMPSON * BOLTZMANN / (MELEC * C * C); //Keep all the constants together
        x*=xplasma->ne*one->vol * xplasma->j * t_e; //multiply by the volume (from wind) and j (from plasma) and t_e


	return(x);
}

	













