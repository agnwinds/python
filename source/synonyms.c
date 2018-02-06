


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

  This is a routine which is intended to help with program updates that
  change the name of a keyword in a parameter file.  It is called from
  rdpar.string_process_from_file

  Description:	

  The routine simply matches the string in new_question to a string in
  the array new_names, below.  If it finds a match, the old_name is 
  returned in old_question.



  Arguments:  	


  Returns:

  	0 if there was no match
	1 if there was a match of a name in the parameter file to
	  one of the new_names, and in old_question the name of
	  the old_value

  Notes:

  To add another variable to the list, just record the new_name and the
  old_name in the arrays below, and increase the number of names by 1

  Do not include the material that is in paren, that is for

  xxx(many_choices) -->  xxxx

  The routine is completely hardwired as wrtten though clearly this 
  could be changed.
  	


  History:
	16sept	ksl	Coded as part of an effort to make python more 
			robust to changes in parameter file names

 ************************************************************************/

#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "log.h"


#define	LINELEN 132

/* Note that there can be multiple old names that are the same due
 * to some old input formats that were not syntaciticaly perfect
 * See #319
 */


char *old_names[] =
  { "mstar", "rstar", "Disk.illumination.treatment", "disk.type",
  "Disk_radiation", "Rad_type_for_disk", "disk.mdot", "T_profile_file",
    "disk.radmax",
  "stellar_wind_mdot", "stellar.wind.radmin", "stellar.wind_vbase",
    "stellar.wind.v_infinity", "stellar.wind.acceleration_exponent",
    "spectrum_wavemin","spectrum_wavemax","no_observers","angle",
    "phase","live.or.die","spec.type","mstar","rstar","Star_radiation",
    "tstar","Rad_type_for_star","Rad_type_for_star",
    "Rad_type_for_disk","Rad_type_for_bl","Boundary_layer_radiation",
    "Rad_type_for_bl","t_bl","lum_bl","homologous_boundary_mdot","msec",
    "period"
};

char *new_names[] = { "Central.object.mass", "Central.object.radius",
  "Surface.reflection.or.absorption", "Disk.type", "Disk.radiation",
    "Disk.rad_type_to_make_wind",
  "Disk.mdot", "Disk.T_profile_file", "Disk.radmax",
  "Stellar_wind.mdot", "Stellar_wind.radmin", "Stellar_wind.vbase",
    "Stellar_wind.v_infinity", "Stellar_wind.acceleration_exponent",
    "Spectrum.wavemin","Spectrum.wavemax","Spectrum.no_observers",
    "Spectrum.angle","Spectrum.orbit_phase","Spectrum.live_or_die",
    "Spectrum.type","Central_object.mass","Central_object.radius",
    "Central_object.radiation","Central_object.temp","Central_object.rad_type_to_make_wind",
    "Central_object.rad_type_in_final_spectrum",
    "Disk.rad_type_in_final_spectrum","Boundary_layer.rad_type_in_final_spectrum",
    "Boundary_layer.radiation","Boundary_layer.rad_type_to_make_wind","Boundary_layer.temp",
    "Boundary_layer.luminosity","homologous.boundary_mdot","Binary.mass_sec",
    "Binary.period"
};

int nunber_of_names = 36;



int
check_synonyms (new_question, old_question)
     char new_question[], old_question[];
{
  int n;
  char *line;
  char firstword[LINELEN];
  int nwords, wordlength;
  char *ccc, *index ();


// Strip off any extra bits in the new question
  line = new_question;
  strcpy (firstword, "");
  nwords = sscanf (line, "%s ", firstword);
  wordlength = strlen (firstword);


  if (nwords == 0)
    {
      return (0);
    }


/* Strip off everthing prior to open paren, leaving the bare parameter name that was passed as the new
 * question
 */

  if ((ccc = index (firstword, '(')) != NULL)
    {
      wordlength = (int) (ccc - firstword);
      if (wordlength == 0)
	return (0);
    }


/* firstword is the bare parameter name for the current way the parameter is expressed elsewhere.
 * We must find this in the list of new_names
 */


  for (n = 0; n < nunber_of_names; n++)
    {
      if (strncmp (new_names[n], firstword, wordlength) == 0)
	{
	  Log
	    ("Matched keyword %s in .pf file to %s in current python version\n",
	     new_question, old_names[n]);
	  strcpy (old_question, old_names[n]);
	  return (1);
	}

    }

  return (0);

}
