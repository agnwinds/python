/* diag.c contains a number of routines relating to extra diagnositics
   and related file i/o. 
   JM 1410 -- moved a few routines here relating to extra diagnostics
*/

/***********************************************************
                                       Southampton

 Synopsis:
	 open_diagfile sets up diagnostic files for use in python when the extra.diagnostics flag is 

Arguments:

	PhotPtr p;	the photon
	double ds	the distance the photon has travelled in the cell

Returns:
	Always returns 0.  .
 
Description:	 
Notes:
	

History:

	12jun 	nsh	72 Added lines to set up a file for outputting photons in given cells. The cells are read in from a file called diag_cells.dat. The existance of this file defined wether the diagnostic takes place.
	13jul	jm	changed print statements to logs and made more descriptive
	14oct   jm  changed to reflect move to 'advanced mode', see #111 and #120

**************************************************************/




#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "atomic.h"
#include "python.h"



/***********************************************************
             University of Southampton

Synopsis:
  get_standard_care_factors provides more control over how the program is
  run

Arguments:

Returns:

Description:

Notes:
    ksl - It is not obvious that much recent thought has been
    given to the choices that are here.  The fractional distance
    that a photon travel is intended to make sure the velocity
    along the line of sight can be approximated linearly.  If 
    a photon travels too far in an azimuthal direction the sense
    of the velocity can change and this prevensts this

    The lowest ion density contributing to photoionization is used
    to determine what ions one has to calculate the photoionzation
    xsection for.  The lower this density; the more that have to be
    calculated, and as a result the slower the program.

    Keeping photoionizaion during final spectrum allows one to
    check the contribution of photoabsorption.

History:
  1802  ksl  Moved here  to consolidate some of the inputs

**************************************************************/
int
get_standard_care_factors ()
{
  int istandard;
  istandard = 1;
  SMAX_FRAC = 0.5;
  DENSITY_PHOT_MIN = 1.e-10;

  /* 141116 - ksl - Made care factors and advanced command as this is clearly somethng that is diagnostic */

  if (modes.iadvanced)
    {
      rdint ("@Diag.use_standard_care_factors(1=yes)", &istandard);

      if (!istandard)
	{
	  rddoub ("@Diag.fractional_distance_photon_may_travel", &SMAX_FRAC);
	  rddoub ("@Diag.lowest_ion_density_for_photoabs", &DENSITY_PHOT_MIN);
	  rdint ("@Diag.keep_photoabs_in_final_spectra(1=yes)",
		 &modes.keep_photoabs);
	}
    }
  return (0);
}

/***********************************************************
				University of Southampton

Synopsis:
	get_extra_diagnostics reads in extra diagnostics if 
	the user has asked for them. It uses rd_int() in rdpar.c
	to get the modes, which should always be one or zero.

Arguments:	
    none	

Returns:
    modifies the modes structure to turn on whichever modes
    are asked for.
 
Description:	
	
Notes:
    see #111 and #120

History:
    1410 -- JM -- Coded
**************************************************************/

int
get_extra_diagnostics ()
{
  if (modes.iadvanced == 0)
    Error ("Getting extra_diagnostics but advanced mode is off!\n");

  Log ("get_extra_diagnostics: Getting extra diagnostics as requested...\n");

  /* read the options. */
  rdint ("@Diag.save_cell_statistics", &modes.save_cell_stats);
  rdint ("@Diag.ispymode", &modes.ispy);
  rdint ("@Diag.keep_ioncycle_windsaves", &modes.keep_ioncycle_windsaves);
  rdint ("@Diag.make_ioncycle_tables", &modes.make_tables);
  rdint ("@Diag.save_extract_photons", &modes.save_extract_photons);
  rdint ("@Diag.print_dvds_info", &modes.print_dvds_info);
  rdint ("@Diag.track_resonant_scatters", &modes.track_resonant_scatters);

  if (modes.save_cell_stats || modes.ispy
      || modes.save_extract_photons | modes.track_resonant_scatters)
    {
      modes.extra_diagnostics = 1;
    }
  else
    {
      modes.extra_diagnostics = 0;
    }

  return 0;
}


/* the next few variables indicate whether certain diagnostic files have
 * already been opened
 */

int eplinit = 0;
int pstatinit = 0;		/*To say if we have checked to see if we need to log photons */
FILE *epltptr;			/* Extra diagnostics file */


int
init_extra_diagnostics ()
{
  FILE *cellfile;		/*File that may or may not exist, pointing to cells we want to write out photon stats for */
  int cell;			/*Temporary storage of cell to use */

  if (eplinit == 0 && modes.save_extract_photons)
    {
      epltptr = fopen ("python.ext", "w");
      eplinit = 1;
    }

  ncstat = 0;			/*Zero the counter for the number of cells to be tracked */
  if (pstatinit == 0 && modes.save_cell_stats)	/* Check we havent already done this */
    {
      cellfile = fopen ("diag_cells.dat", "r");	/*This is the file containing cells to track */
      if (cellfile != NULL)	/*If there actually *is* a file read it */
	{
	  while (fscanf (cellfile, "%d", &cell) == 1)	/*If the line contains only one integer number read it in, otherwise quit reading */
	    {
	      Log
		("open_diagfile: Cell diagnostics - we have a cell - %i, ncstat=%i, NCSTAT=%i\n",
		 cell, ncstat, NCSTAT);
	      if (-1 < cell && cell < geo.nplasma && ncstat < NCSTAT)	/*if the cells are real */
		{
		  Log
		    ("open_diagfile: Cell numbers have been accepted as real.\n");
		  ncell_stats[ncstat] = cell;
		  ncstat = ncstat + 1;
		}
	      else
		{
		  Error
		    ("open_diagfile: %i is an unacceptable cell number for photon tracking\n",
		     cell);
		}
	    }
	  fclose (cellfile);
	  pstatptr = fopen ("cell_phot_stats.dat", "w");
	}
      else
	{
	  Log
	    ("open_diagfile: We have no file of cells to track, so we wont be doing any cell tracking\n");
	}
      pstatinit = 1;		/* We have initialised this routine */
    }

  return (0);
}





/***********************************************************
                Southampton University

Synopsis: 
	save_photon_stats prints photon statistics to a file

Arguments:	
	One 		WindPtr for the cell
	p 			Photon pointer
	ds 			ds travelled

Returns:
 
Description:
   the loop below is if the user requires extra diagnostics and
   has provided a file diag_cells.dat to store photons stats for cells they have specified

Notes:
   Moved here to save duplicating code between bf_estimators_increment and radiation.

History:
   1410 JM 		Coding began
   1410 JM      Moved here from python.c	
 
**************************************************************/



int
save_photon_stats (one, p, ds, w_ave)
     WindPtr one;
     PhotPtr p;
     double ds, w_ave;
{
  int i;

  /* JM -- 1310 -- the loop below is if the user requires extra diagnostics and
     has provided a file diag_cells.dat to store photons stats for cells they have specified
   */

  for (i = 0; i < ncstat; i++)
    {
      /* check if the cell is in the specified list - ncell_stats is global variable */
      if (one->nplasma == ncell_stats[i])
	{
	  fprintf (pstatptr,
		   "PHOTON_DETAILS cycle %3d n_photon %d freq %8.3e  w %8.3e ave_w %8.3e ds %8.3e nscat %d plasma cell %3d wind cell %3d\n",
		   geo.wcycle, p->np, p->freq, p->w, w_ave, ds, p->nscat,
		   one->nplasma, one->nwind);
	}
    }
  return (0);
}

int
save_extract_photons (n,p, pp, v)
    int n;
     PhotPtr p, pp;
     double *v;
     {
       fprintf (epltptr,
		"EXTRACT %3d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %7.2f %7.2f \n",
		n, p->x[0], p->x[1], p->x[2], v[0], v[1], v[2],
		p->lmn[0], p->lmn[1], p->lmn[2], pp->lmn[0], pp->lmn[1],
		pp->lmn[2], 2.997925e18 / p->freq, 2.997925e18 / pp->freq);

       return(0);
     }
