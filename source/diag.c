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

int eplinit = 0;
int pstatinit = 0;              /*To say if we have checked to see if we need to log photons */

int
open_diagfile ()
{
  FILE *cellfile;               /*File that may or may not exist, pointing to cells we want to write out photon stats for */
  int cell;                     /*Temporary storage of cell to use */

  if (eplinit == 0 && modes.save_extract_photons)
  {
    epltptr = fopen ("python.ext", "w");
    eplinit = 1;
  }

  ncstat = 0;                   /*Zero the counter for the number of cells to be tracked */
  if (pstatinit == 0 && modes.save_cell_stats)  /* Check we havent already done this */
  {
    cellfile = fopen ("diag_cells.dat", "r");   /*This is the file containing cells to track */
    if (cellfile != NULL)       /*If there actually *is* a file read it */
    {
      while (fscanf (cellfile, "%d", &cell) == 1)       /*If the line contains only one integer number read it in, otherwise quit reading */
      {
        Log ("open_diagfile: Cell diagnostics - we have a cell - %i, ncstat=%i, NCSTAT=%i\n", cell, ncstat, NCSTAT);
        if (-1 < cell && cell < geo.nplasma && ncstat < NCSTAT) /*if the cells are real */
        {
          Log ("open_diagfile: Cell numbers have been accepted as real.\n");
          ncell_stats[ncstat] = cell;
          ncstat = ncstat + 1;
        }
        else
        {
          Error ("open_diagfile: %i is an unacceptable cell number for photon tracking\n", cell);
        }
      }
      fclose (cellfile);
      pstatptr = fopen ("cell_phot_stats.dat", "w");
    }
    else
    {
      Log ("open_diagfile: We have no file of cells to track, so we wont be doing any cell tracking\n");
    }
    pstatinit = 1;              /* We have initialised this routine */
  }

  // diag_on_off = 0;            // 0=off everything else is on
  return (0);
}






/***********************************************************
				University of Southampton

Synopsis:
	get_extra_diagnostics reads in extra diagnostics if 
	the user has asked for them. It uses rd_extra() in rdpar.c
	to get the actual lines, and then checks them against 
	hardwired options via a number of if loops.

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
  {
    Error ("Trying to get extra_diagnostics but advanced mode is off!\n");
    return 0;
  }
  
  Log ("get_extra_diagnostics: Getting extra diagnostics as requested...\n");

  rdint_extra ("save_cell_statistics", &modes.save_cell_stats);
  rdint_extra ("track_resonant_scatters", &modes.track_resonant_scatters);
  rdint_extra ("keep_ioncycle_windsaves", &modes.keep_ioncycle_windsaves);
  rdint_extra ("save_extract_photons", &modes.save_extract_photons);
  rdint_extra ("print_windrad_summary", &modes.print_windrad_summary);
  rdint_extra ("print_dvds_info", &modes.save_cell_stats);
  rdint_extra ("save_cell_statistics", &modes.save_cell_stats);
  rdint_extra ("ispymode", &modes.ispy);

  return 0;
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
               geo.wcycle, p->np, p->freq, p->w, w_ave, ds, p->nscat, one->nplasma, one->nwind);
    }
  }
  return (0);
}
