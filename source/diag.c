
/***********************************************************
                          Space Telescope Science Insittute

 Synopsis:
	 diag.c  consolidates an approach to providing extra diagnositcs 
     for debugging python, Some of the diagnostics allow one to change
     defaults in python while others allow traking of and contains initialization routines
     that control which diagnosics are a

Arguments:


Returns:
	Always returns 0.  .
 
Description:	 

    The first few functions are associated with gather information about what
    one wants to track

    The basic idea for the extra diagnositics routines are as follows:

    There is one file to which all of the diagnostics are written.
    Which diagnositica are written depends on how varaious of the modes
    variables are set.

    All of the various diagnostic routines should be line oriented and begin
    with a unique identifier to allow easy parsing of the outputs.


    To add a new diagnostic one should 
    
    (0) check that one of the existing diagnostics is insuffient for your needs 
    (1) add a routine here that writes to the diagnostic file.  
    (2) add a rdpar statement in get_extra_diagnostics
    (3) Assure that your new routine has the appropriate doxygne description, 
    (4) Add a yaml file for the new rdpare statemnt), and 
    (5) call the routine with the approriate condidtions from wherever you wish.

    As written, you must put your routine here.  This is done on purpose
    so that we will not proliferate extra diagnostics around the various
    subroutines.


Notes:

    The unified approach makes most sense when one wants similar information
    spread out over several routines.  There will be places where one
    wishes to print out detailed information within a certain routine.  
    A current example is print_dvds_info, which is totally contained within 
    gradv.c
	

History:
    180102 - Adopted the approach above based, as discussed in #338


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
  rdint ("@Diag.keep_ioncycle_windsaves", &modes.keep_ioncycle_windsaves);
  rdint ("@Diag.make_ioncycle_tables", &modes.make_tables);
  rdint ("@Diag.save_photons", &modes.save_photons);
  rdint ("@Diag.save_extract_photons", &modes.save_extract_photons);
  rdint ("@Diag.print_dvds_info", &modes.print_dvds_info);
  rdint ("@Diag.track_resonant_scatters", &modes.track_resonant_scatters);

  if (modes.save_cell_stats || modes.save_photons
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



/***********************************************************
				University of Southampton

Synopsis:
	init_extra_diagnositics opons a file for writing out 
    extra diagnostics and in some cases reads a file
    that specifies in which cells ones wants diagnostics

Arguments:	

    none	

Returns:
 
Description:	
	
Notes:
    see #111 and #120

    The diagnostic filename is hardwired

History:
    1410 -- JM -- Coded
**************************************************************/


int eplinit = 0;
int pstatinit = 0;		/*To say if we have checked to see if we need to log photons */
FILE *epltptr;			/* Extra diagnostics file */


int
init_extra_diagnostics ()
{
  FILE *cellfile;		/*File that may or may not exist, pointing to cells we want to write out photon stats for */
  int cell;			/*Temporary storage of cell to use */

  if (eplinit == 0 && modes.extra_diagnostics)
    {
      epltptr = fopen ("python.ext.txt", "w");
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


/***********************************************************
                Space Telescope Science Insittue

Synopsis: 
	save_extract_photons saves informations about phtoons in
    a particulare wavelength gange

Arguments:	
Returns:
 
Description:

Notes:
   Moved here to save duplicating code between bf_estimators_increment and radiation.

History:
   1802 ksl Functionality moved from extract.c to consolidiate how extra
            diagnostics were carried out.
 
**************************************************************/

int
save_extract_photons (n, p, pp, v)
     int n;
     PhotPtr p, pp;
     double *v;
{
  fprintf (epltptr,
	   "EXTRACT %3d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %7.2f %7.2f \n",
	   n, p->x[0], p->x[1], p->x[2], v[0], v[1], v[2],
	   p->lmn[0], p->lmn[1], p->lmn[2], pp->lmn[0], pp->lmn[1],
	   pp->lmn[2], 2.997925e18 / p->freq, 2.997925e18 / pp->freq);

  return (0);
}

/***********************************************************
             Space Telescope Science Institute

Synopsis: 
	save_photon

Arguments:	
	p 			Photon pointer
	comment 	An arbitrary comment

Returns:
 
Description:
   This is a diagnostic that allows one to print out information about
   a photon at any time.  It was written as part of the effort to
   debub imports for the fu_ori project.

Notes:

History:
   !802 ksl Coded
 
**************************************************************/
int save_photon_number=0;

int
save_photons (p, comment)
     PhotPtr p;
     char comment[];
{
    save_photon_number+=1;
    if (save_photon_number>100000) 
        return(0);

fprintf (epltptr,
"PHOTON %3d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %3d %3d %s \n",
p->np, p->x[0], p->x[1], p->x[2], p->lmn[0], p->lmn[1],
p->lmn[2], p->grid, p->istat,comment);
	   
return(0);
}



/***********************************************************
             Space Telescope Science Institute

Synopsis: 
	track_scatters

Arguments:	
	p 			Photon pointer
    nplasma     the cell in which the photon resides
	comment 	An arbitrary comment

Returns:
 
Description:
   Another way to track photons. This one is focussed on what happens
   during scatterin events.

Notes:

History:
   1802 ksl Coded, to move this from transphot to better 
            reflect our unified scheme for writing extra
            diagnostics.
 
**************************************************************/


int 
track_scatters(p,nplasma,comment)
    PhotPtr p;
    int nplasma;
    char *comment;
{

 fprintf (epltptr, "Scattter %i %.2e %.2e %.2e  %i %e %e %i %s\n", p->np, p->x[0], p->x[1], p->x[2], p->grid,p->freq,p->w, nplasma,comment);

 return(0);
}
