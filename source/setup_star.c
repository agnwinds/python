#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


/***********************************************************
             University of Southampton

Synopsis:
  get_stellar_params sets rstar, mstar, tstar as well
  as secondary parameters based on user inputs

Arguments:

Returns:


Description:

Notes:

History:
	1502  JM 	Moved here from main()

**************************************************************/

double
get_stellar_params ()
{

  /* Describe the basic binary star system */

  geo.mstar /= MSOL;		// Convert to MSOL for ease of data entry
  rddoub ("Central_object.mass(msol)", &geo.mstar);
  geo.mstar *= MSOL;

  /* If a BH we want geo.rstar to be at least as large as the last stable orbit for
   * a non-rotating BH
   */

  if (geo.system_type == SYSTEM_TYPE_AGN)
    {
      geo.rstar = 6. * G * geo.mstar / (C * C);	//correction - ISCO is 6x Rg NSH 121025
    }

  rddoub ("Central_object.radius(cm)", &geo.rstar);


  geo.r_agn = geo.rstar;	/* At present just set geo.r_agn to geo.rstar */
  geo.rstar_sq = geo.rstar * geo.rstar;
  if (geo.system_type != SYSTEM_TYPE_AGN)
    {
      rdint ("Central_object.radiation(y=1)", &geo.star_radiation);
      get_spectype (geo.star_radiation,
		    //"Rad_type_for_star(0=bb,1=models)_to_make_wind",
		    "Central_object.rad_type_to_make_wind(0=bb,1=models)",
		    &geo.star_ion_spectype);

      if (geo.star_radiation)
	rddoub ("Central_object.temp", &geo.tstar_init);
    }
  else
    {
      geo.star_radiation = 0;
      geo.tstar_init = 0;
    }

  /* tstar_init and lum_star_init refer to values without the effects of backscattering */

  geo.tstar = geo.tstar_init;

  geo.lum_star = geo.lum_star_init =
    4 * PI * geo.rstar * geo.rstar * STEFAN_BOLTZMANN * pow (geo.tstar, 4.);


  /* Describe the secondary if that is required */

  if (geo.system_type == SYSTEM_TYPE_BINARY)	/* It's a binary system */
    {

      geo.m_sec /= MSOL;	// Convert units for ease of data entry
      rddoub ("msec(msol)", &geo.m_sec);
      geo.m_sec *= MSOL;

      geo.period /= 3600.;	// Convert units to hours for easy of data entry
      rddoub ("period(hr)", &geo.period);
      geo.period *= 3600.;	// Put back to cgs immediately
    }

  return (geo.lum_star_init);
}


