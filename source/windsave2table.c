
/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	windsave2table writes key variables in a wind save file to an astropy table    
		as calculated by python.  This is the main routine.

Arguments:		

	py_wind  windsave_root



Returns:
 
Description:	
	

	
Notes:



History:
	150428	ksl	Adpated from routines in py_wind.c

**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"
#include "python.h"



int
main (argc, argv)
     int argc;
     char *argv[];
{


  int i;
  int ochoice;
  char c;

  char root[LINELENGTH], input[LINELENGTH], wspecfile[LINELENGTH],
    specfile[LINELENGTH];
  char outputfile[LINELENGTH];
  char windradfile[LINELENGTH], windsavefile[LINELENGTH];
  char parameter_file[LINELENGTH];
  char photfile[LINELENGTH];
  double xdoub;
  int interactive;


  // py_wind uses rdpar, but only in an interactive mode. As a result 
  // there is no associated .pf file

  interactive = 1;		/* Default to the standard operating mofe for py_wind */
  strcpy (parameter_file, "NONE");

  /* Next command stops Debug statements printing out in py_wind */
  Log_set_verbosity (3);

  if (argc == 1)
    {
      printf ("Root for wind file :");
      fgets (input, LINELENGTH, stdin);
      get_root (root, input);
	}
  else
  {
      strcpy (input, argv[argc - 1]);
      get_root (root, input);
    }


  printf ("Reading data from file %s\n", root);

  /* Now create the names of all the files which will be written */

  strcpy (windsavefile, root);
  strcpy (outputfile, root);

  strcat (windsavefile, ".wind_save");
  strcat (outputfile,".txt");


/* Read in the wind file */

/* Note that wind_read allocates the space for the WindPtr array.  The
reason that this is done in wind_read is because wind_read also reads
the geo structure.  The geo struc contains the dimensions of the wind
array and so until it is read the space for the WindPtr structure cannot
be allocated.  Howver, one cannot assign w in a subroutine.  Therefore
the old call, namely wind_read(w,windsavefile, will not work.  The statement
w=wmain is a way around this since wmain is external and therefore can
be assigned.  Then w can be set to this ptr and all is restored. The
use of w is endemic in the program. and it is always called through main.
I did not change this now.  Though it could be done.  02apr ksl */

  if (wind_read (windsavefile) < 0)
    {
      Error ("py_wind: Could not open %s", windsavefile);
      exit (0);
    }

/* aaa is used to store variable for writing to files for the purpose of plotting*/
  aaa = calloc (sizeof (xdoub), NDIM2);

  printf ("Read wind_file %s\n", windsavefile);

  get_atomic_data (geo.atomic_filename);

  printf ("Read Atomic data from %s\n", geo.atomic_filename);


  ochoice=1;
  xcomplete_file_summary (wmain, root, ochoice);


  create_master_table(root);
  return (0);
}



int
create_master_table(rootname)
	char rootname[];
{
	double xdoub;

	int i,ii,jj;
  FILE *fopen (), *fptr;

  fptr=fopen("test.txt","w");

	get_t_e();
  	xt_r = calloc (sizeof (xdoub), NDIM2);
	get_t_r(&xt_r);

      if (geo.coord_type == SPHERICAL)
	{
    
    /* JM 1411 -- added header for reading with astropy ascii module */
    fprintf (fptr, "r var inwind i\n");



	  for (i = 0; i < NDIM; i++)
	    {
	      fprintf (fptr, "%8.2e %8.2e %3d %3d \n", wmain[i].r, aaa[i],
		       wmain[i].inwind, i);
	    }
	}
      else
	{

    /* JM 1411 -- added header for reading with astropy ascii module */
    fprintf (fptr, "x z var inwind i j\n");

	  for (i = 0; i < NDIM2; i++)
	    {
	      wind_n_to_ij (i, &ii, &jj);
	      fprintf (fptr, "%8.4e %8.4e %8.2e %3d %3d %3d\n",
		       wmain[i].xcen[0], wmain[i].xcen[2], aaa[i],
		       wmain[i].inwind, ii, jj);
	    }
	}
}


/**************************************************************************


  Synopsis:  

  complete_file summary produces a standardised set of output files 
  from the wind_save_file

  Description:	

  At present(1409), we we print temperatures, electron densities,
  convergence information and ion denisities and fractions for
  C III, IV, V, 
  N IV, V, VI
  O V, VI, VII
  Si III, IV

  Arguments:  

  Returns:

  Notes:

  History:

  1409	ksl	Eliminated several files that were of limited use
  		to users.  The deleted files are printed out if
		  DEBUG is set.
  1411  JM debug is now deprecated so replaced with FULL_ION_SUMMARTY

 ************************************************************************/

#define FULL_ION_SUMMARY 0  


int
xcomplete_file_summary (w, root, ochoice)
     WindPtr w;
     char root[];
     int ochoice;
{
  xtemp_summary (w, root, ochoice);

  temp_rad (w, root, ochoice);
  electron_summary (w, root, ochoice);
  convergence_summary (w, root, ochoice);
  xion_summary (w, 6, 3, 0, root, ochoice);
  xion_summary (w, 6, 4, 0, root, ochoice);
  xion_summary (w, 6, 5, 0, root, ochoice);
  xion_summary (w, 7, 4, 0, root, ochoice);
  xion_summary (w, 7, 5, 0, root, ochoice);
  xion_summary (w, 7, 6, 0, root, ochoice);
  xion_summary (w, 8, 4, 0, root, ochoice);
  xion_summary (w, 8, 5, 0, root, ochoice);
  xion_summary (w, 8, 6, 0, root, ochoice);
  xion_summary (w, 8, 7, 0, root, ochoice);
  xion_summary (w, 14, 3, 0, root, ochoice);
  xion_summary (w, 14, 4, 0, root, ochoice);
  xion_summary (w, 14, 5, 0, root, ochoice);

  xion_summary (w, 6, 3, 1, root, ochoice);
  xion_summary (w, 6, 4, 1, root, ochoice);
  xion_summary (w, 6, 5, 1, root, ochoice);
  xion_summary (w, 7, 4, 1, root, ochoice);
  xion_summary (w, 7, 5, 1, root, ochoice);
  xion_summary (w, 7, 6, 1, root, ochoice);
  xion_summary (w, 8, 4, 1, root, ochoice);
  xion_summary (w, 8, 5, 1, root, ochoice);
  xion_summary (w, 8, 6, 1, root, ochoice);
  xion_summary (w, 8, 7, 1, root, ochoice);
  xion_summary (w, 14, 3, 1, root, ochoice);
  xion_summary (w, 14, 4, 1, root, ochoice);
  xion_summary (w, 14, 5, 1, root, ochoice);

#if FULL_ION_SUMMARY
  xion_summary (w, 6, 3, 2, root, ochoice);
  xion_summary (w, 6, 4, 2, root, ochoice);
  xion_summary (w, 6, 5, 2, root, ochoice);
  xion_summary (w, 7, 4, 2, root, ochoice);
  xion_summary (w, 7, 5, 2, root, ochoice);
  xion_summary (w, 7, 6, 2, root, ochoice);
  xion_summary (w, 8, 4, 2, root, ochoice);
  xion_summary (w, 8, 5, 2, root, ochoice);
  xion_summary (w, 8, 6, 2, root, ochoice);
  xion_summary (w, 8, 7, 2, root, ochoice);
  xion_summary (w, 14, 3, 2, root, ochoice);
  xion_summary (w, 14, 4, 2, root, ochoice);
  xion_summary (w, 14, 5, 2, root, ochoice);


  xion_summary (w, 6, 3, 3, root, ochoice);
  xion_summary (w, 6, 4, 3, root, ochoice);
  xion_summary (w, 6, 5, 3, root, ochoice);
  xion_summary (w, 7, 4, 3, root, ochoice);
  xion_summary (w, 7, 5, 3, root, ochoice);
  xion_summary (w, 7, 6, 3, root, ochoice);
  xion_summary (w, 8, 4, 3, root, ochoice);
  xion_summary (w, 8, 5, 3, root, ochoice);
  xion_summary (w, 8, 6, 3, root, ochoice);
  xion_summary (w, 8, 7, 3, root, ochoice);
  xion_summary (w, 14, 3, 3, root, ochoice);
  xion_summary (w, 14, 4, 3, root, ochoice);
  xion_summary (w, 14, 5, 3, root, ochoice);
#endif

  return (0);
}


int
xion_summary (w, element, istate, iswitch, rootname, ochoice)
     WindPtr w;
     int element, istate;
     int iswitch;
     char rootname[];
     int ochoice;
{
  int nion, nelem;
  int n;
  char choice[LINELENGTH], iname[LINELENGTH];
  char name[LINELENGTH];
  char filename[LINELENGTH];
  double x;
  int nplasma;



/* Find the CIV ion */
  nion = 0;
  while (nion < nions
	 && !(ion[nion].z == element && ion[nion].istate == istate))
    nion++;
  if (nion == nions)
    {
      Log ("Error--element %d ion %d not found in define_wind\n", element,
	   istate);
      return (-1);
    }
  nelem = 0;
  while (nelem < nelements && ele[nelem].z != element)
    nelem++;

  strcpy (name, "");


  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      nplasma = w[n].nplasma;
      if (w[n].vol > 0.0 && plasmamain[nplasma].ne > 1.0)
	{
	  if (iswitch == 0)
	    {
  sprintf (name, "Element %d (%s) ion %d fractions\n", element,
	   ele[nelem].name, istate);
	      aaa[n] = plasmamain[nplasma].density[nion];
	      aaa[n] /=
		((plasmamain[nplasma].density[0] +
		  plasmamain[nplasma].density[1]) * ele[nelem].abun);
	    }
	  else if (iswitch == 1)
	    {
  sprintf (name, "Element %d (%s) ion %d density\n", element,
	   ele[nelem].name, istate);
	      aaa[n] = plasmamain[nplasma].density[nion];
	    }
	  else if (iswitch == 2)
	    {
  sprintf (name, "Element %d (%s) ion %d  #scatters\n", element,
	   ele[nelem].name, istate);
	      aaa[n] = plasmamain[nplasma].scatters[nion];
	    }
	  else if (iswitch == 3)
	    {
  sprintf (name, "Element %d (%s) ion %d scattered flux\n", element,
	   ele[nelem].name, istate);
	      aaa[n] = plasmamain[nplasma].xscatters[nion];
	    }
	  else
	    {
	      Error ("xion_summary : Unknown switch %d \n", iswitch);
	      exit (0);
	    }
	}
    }





  /* Store the appropriate values in a place where it does not matter */
  if (ochoice)
    {
      for (n = 0; n < NDIM2; n++)
	{
	  nplasma = w[n].nplasma;
	  if (w[n].vol > 0.0 && plasmamain[nplasma].ne > 1.0)
	    {
	      if (iswitch == 0)
		x /=
		  ((plasmamain[nplasma].density[0] +
		    plasmamain[nplasma].density[1]) * ele[nelem].abun);
	      else if (iswitch == 1)
		{
		  x = plasmamain[nplasma].density[nion];
		  x = log10 (x);
		}
	      else if (iswitch == 2)
		{
		  x = plasmamain[nplasma].scatters[nion];
		}
	      else if (iswitch == 3)
		{
		  x = plasmamain[nplasma].xscatters[nion];
		}
	      else
		{
		  Error ("xion_summary : Unknown switch %d \n", iswitch);
		  exit (0);
		}


	    }
	  else
	    x = 0.0;
	  w[n].x[1] = x;

	}

      strcpy (filename, rootname);
      if (iswitch == 0)
	strcpy (choice, ".ion");
      else if (iswitch == 1)
	strcpy (choice, ".ionc");
      else if (iswitch == 2)
	strcpy (choice, ".ions");
      else if (iswitch == 3)
	strcpy (choice, ".iona");
      else
	{
	  Error ("xion_summary : Unknown switch %d \n", iswitch);
	  exit (0);
	}

      strcat (choice, ele[nelem].name);
      sprintf (iname, "%d", istate);
      strcat (choice, iname);

      strcat (filename, choice);
      write_table (filename, ochoice);
    }

  return (0);
}



#define ODIM 256
int init_write_table = 1;
float aout[ODIM][ODIM];

int
write_table (filename, choice)
     char filename[];
     int choice;
{
  //Dynamical allocation is allowed, although I generally avoid it -- 05apr ksl
  float rin[NDIM], zin[MDIM];
  float r, z;
  float rmin, rmax, zmin, zmax;
  int ii, jj;
  FILE *fopen (), *fptr;
  char outfile[LINELENGTH];
  double length ();
  double xx[3];
  int i;
  int nn, nnn[4], nelem;
  double frac[4];

  if (choice == 0)		// Then we don't really want output files
    {
      return (0);
    }

// Open the appropriate file

  strcpy (outfile, filename);
  strcat (outfile, ".dat");	// Add standard extension to filenames, i.e. the one used by tecplot
  fptr = fopen (outfile, "w");


  /* Write out the header information for the file */
  fprintf (fptr, "# TITLE= \"%s\"\n", outfile);
  if (geo.coord_type == SPHERICAL)
    {
      fprintf (fptr, "# Coord_Sys SPHERICAL\n");
    }
  else if (geo.coord_type == CYLIND)
    {
      fprintf (fptr, "# Coord_Sys CYLIND\n");
    }
  else if (geo.coord_type == RTHETA)
    {
      fprintf (fptr, "# Coord_Sys RTHETA\n");
    }
  else if (geo.coord_type == CYLVAR)
    {
      fprintf (fptr, "# Coord_Sys CYLVAR\n");
    }
  else
    {
      Error ("write_table: Unknown coordinate system type: %d\n",
	     geo.coord_type);
    }


// Put the r and z coord. grid  into easier to understand arrays
  for (i = 0; i < NDIM; i++)
    {
      zin[i] = wmain[i].x[2];
      rin[i] = wmain[i * MDIM].x[0];
    }


/* In py_wind the filenames are set to start with z if the output coordinates
are linear, and x otherwise.  This is not particularly transparent ?? ksl */

  if (choice == 1)
    {				// print out the original array elements

      if (geo.coord_type == SPHERICAL)
	{
    
    /* JM 1411 -- added header for reading with astropy ascii module */
    fprintf (fptr, "r var inwind i\n");

	  for (i = 0; i < NDIM; i++)
	    {
	      fprintf (fptr, "%8.2e %8.2e %3d %3d \n", wmain[i].r, aaa[i],
		       wmain[i].inwind, i);
	    }
	}
      else
	{

    /* JM 1411 -- added header for reading with astropy ascii module */
    fprintf (fptr, "x z var inwind i j\n");

	  for (i = 0; i < NDIM2; i++)
	    {
	      wind_n_to_ij (i, &ii, &jj);
	      fprintf (fptr, "%8.4e %8.4e %8.2e %3d %3d %3d\n",
		       wmain[i].xcen[0], wmain[i].xcen[2], aaa[i],
		       wmain[i].inwind, ii, jj);
	    }
	}

    }
  else if (choice == 2)
    {				// Then regrid to a linear grid
/* Wind array is organized so that it increments up along z and then over */

/* Want linearly spaced output array */
      rmin = zmin = 0;
      rmax = zmax = geo.rmax;

      xx[1] = 0.0;
      for (ii = 0; ii < ODIM; ii++)
	{
	  xx[0] = r = rmin + ii * (rmax - rmin) / (ODIM - 1);
	  for (jj = 0; jj < ODIM; jj++)
	    {
	      xx[2] = z = zmin + jj * (zmax - zmin) / (ODIM - 1);
	      //OLD 70B if (where_in_wind (xx) == 0)
	      if (where_in_wind (xx) >= 0)
		{		// Then the position is in the wind region
		  coord_fraction (0, xx, nnn, frac, &nelem);
		  for (nn = 0; nn < nelem; nn++)
		    {
		      aout[ii][jj] += aaa[nnn[nn]] * frac[nn];
		    }

		}
	    }
	}

      /* Now print out the array */
      // fprintf (fptr, "VARIABLES= \"X\" \"Z\" \"Var1\" \n");
      // fprintf (fptr, "ZONE I=%d J=%d DATAPACKING=POINT\n", ODIM, ODIM);
      fprintf (fptr, "# Resampled outputs\n");
      for (jj = 0; jj < ODIM; jj++)
	{
    /* JM 1411 -- added header for reading with astropy ascii module */
    fprintf (fptr, "r z var\n");

	  z = zmin + (zmax - zmin) * jj / (ODIM - 1);
	  for (ii = 0; ii < ODIM; ii++)
	    {
	      r = rmin + (rmax - rmin) * ii / (ODIM - 1);
	      fprintf (fptr, "%8.2e %8.2e %8.2e\n", r, z, aout[ii][jj]);
	    }
	}
      /*Finished printing out the array */
    }
  else
    {
      Error ("write_table: Unknown choice %d\n", choice);
    }

  fclose (fptr);
  return (0);
}

/**************************************************************************


  Synopsis:  
	A summary of the electron temperatures

  Description:	

  Arguments:  

  Returns:

  Notes:

  History:

 ************************************************************************/
int
xtemp_summary (w, rootname, ochoice)
     WindPtr w;
     char rootname[];
     int ochoice;
{
  int n;
  char filename[LINELENGTH];
  int nplasma;

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (w[n].vol > 0.0)
	{
	  nplasma = w[n].nplasma;
	  aaa[n] = plasmamain[nplasma].t_e;
	}
    }

  if (ochoice)
    {
      strcpy (filename, rootname);
      strcat (filename, ".te");
      write_table (filename, ochoice);
    }
  return (0);

}

/**************************************************************************


  Synopsis:  
	A summary of the electron temperatures

  Description:	
  	Put t_e into the extenal array aaa

  Arguments:  

  Returns:

  Notes:
  	Getting any simple variable from the plama structure should
	follow this template.

  History:
  	150429 ksl Adapted from te_summary in py_wind

 ************************************************************************/
int
get_t_e ()
{
  int n;
  int nplasma;

  for (n = 0; n < NDIM2; n++)
    {
      aaa[n] = 0;
      if (wmain[n].vol > 0.0)
	{
	  nplasma = wmain[n].nplasma;
	  aaa[n] = plasmamain[nplasma].t_e;
	}
    }

  return (0);

}
/**************************************************************************


  Synopsis:  
	A summary of the radiation temperatures

  Description:	
  	Put t_r into the external array aaa

  Arguments:  

  Returns:

  Notes:
  	Getting any simple variable from the plama structure should
	follow this template.

  History:
  	150429 ksl Adapted from te_summary in py_wind

 ************************************************************************/
int
get_t_r (x)
	float x[]
{
  int n;
  int nplasma;

  for (n = 0; n < NDIM2; n++)
    {
      x[n] = 0;
      if (wmain[n].vol > 0.0)
	{
	  nplasma = wmain[n].nplasma;
	  x[n] = plasmamain[nplasma].t_r;
	}
    }

  return (0);

}
