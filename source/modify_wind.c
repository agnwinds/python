
/***********************************************************/
/** @file  modify_wind.c
 * @author ksl
 * @date   June, 2020
 *
 * @brief  Routines to modfify a wind structure to for example
 * change densities of certain ions
 *
 *###Notes###
 * At present this is just a prototype of a more general 
 * purpose routine to paint sections of a windsavefile
 * with new values.  At present the code just supports
 * changing densities, although it would be fairly 
 * straighfoward to extend to propeties like T_e
 *
 * At present the changes have to be hardcoded into 
 * main, and the program lacks an interface to for
 * example files of the tupe produced by windsave2table
 * 
 * The intent is to make the structure of this program
 * like the inverse of windsave2table.
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


char inroot[LINELENGTH], outroot[LINELENGTH];


/**********************************************************/
/**
 * @brief      parses the command line options
 *
 * @param [in]  int  argc   the number of command line arguments
 * @param [in]  char *  argv[]   The command line arguments
 *
 *
 * ###Notes###
 *
 * The general purpose of each of the command line options
 * should be fairly obvious from reading the code.
 *
 *
 * Although this routine uses the standard Log and Error commands
 * the diag files have not been created yet and so this information
 * is really simply written to the terminal.
 *
 **********************************************************/

int
xparse_command_line (argc, argv)
     int argc;
     char *argv[];
{
  int j = 0;
  int i;
  char dummy[LINELENGTH];
  int mkdir ();
  char *fgets_rc;


  sprintf (outroot, "%s", "new");


  if (argc == 1)
  {
    printf ("Parameter file name (e.g. my_model.pf, or just my_model):");
    fgets_rc = fgets (dummy, LINELENGTH, stdin);
    if (!fgets_rc)
    {
      printf ("Input rootname is NULL or invalid\n");
      exit (1);
    }
    get_root (inroot, dummy);
  }
  else
  {

    for (i = 1; i < argc; i++)
    {
      if (strcmp (argv[i], "-out_root") == 0)
      {
        if (sscanf (argv[i + 1], "%s", dummy) != 1)
        {
          printf ("python: Expected out_root after -out_root switch\n");
          exit (0);
        }

        get_root (outroot, dummy);
        i++;
        j = i;

      }
      else if (strcmp (argv[i], "--dry-run") == 0)
      {
        modes.quit_after_inputs = 1;
        j = i;
      }
      else if (strncmp (argv[i], "-", 1) == 0)
      {
        printf ("python: Unknown switch %s\n", argv[i]);
        exit (0);
      }
    }

    /* The last command line variable is always the windsave file */

    if (j + 1 == argc)
    {
      printf ("All of the command line has been consumed without specifying a parameter file name, so exiting\n");
      exit (1);
    }
    strcpy (dummy, argv[argc - 1]);
    get_root (inroot, dummy);

  }

  return (0);
}



/**********************************************************/
/**
 * @brief      the main routine which carries out the effort
 *
 * @param [in]  int  argc   the number of command line arguments
 * @param [in]  char *  argv[]   The command line arguments
 *
 *
 * ###Notes###
 *
 * This routine oversees the effort.  The basic steps are
 *
 * - parse the command line to get the names of files
 * - read the old windsave file
 * - read the densities from in this case H
 * - modify the densities
 * - write out the new windsave file
 *
 *
 **********************************************************/


int
main (argc, argv)
     int argc;
     char *argv[];
{

  double *den;
  char name[LINELENGTH];        /* file name extension */
  char infile[LINELENGTH], outfile[LINELENGTH];
  int i;
  int put_ion ();

  Log_set_verbosity (3);
  xparse_command_line (argc, argv);


  sprintf (infile, "%s.wind_save", inroot);
  sprintf (outfile, "%s.wind_save", outroot);

  printf ("Reading %s and writing to %s\n", infile, outfile);

  wind_read (infile);

  den = get_ion (0, 1, 1, 1, name);

  /* Having gotten the densities, modify them */

  printf ("%d\n", zdom[0].ndim2);


  for (i = 0; i < zdom[0].ndim2; i++)
  {
    printf ("%e\n", den[i]);
    den[i] *= 1e-6;
  }

  /* Now paint the densities into the appropriatle
     places in the plasma structure
   */


  put_ion (0, 1, 1, den);



  wind_save (outfile);


  printf ("gotcha %s\n", files.root);

  exit (0);

}


/**********************************************************/
/**
 * @brief      update a specific ion with new densities
 *
 * @param [in]  int ndom   the domain number
 * @param [in]  int element the element(z) for which 
 * @param [in]  int istate  the ion number
 * @param [in]  double den  an array with the new densities
 *
 *
 * ###Notes###
 *
 *
 **********************************************************/

int
put_ion (ndom, element, istate, den)
     int ndom, element, istate;
     double *den;
{
  int i, n;
  int nion;
  int nstart, ndim2;
  int nplasma;



  for (i = 0; i < zdom[0].ndim2; i++)
  {
    printf ("%e\n", den[i]);
  }
  nstart = zdom[ndom].nstart;
  ndim2 = zdom[ndom].ndim2;


  /* Find the ion */

  nion = 0;
  while (nion < nions && !(ion[nion].z == element && ion[nion].istate == istate))
    nion++;

  printf ("Found %d\n", nion);


  for (n = 0; n < ndim2; n++)
  {
    nplasma = wmain[nstart + n].nplasma;
    plasmamain[nplasma].density[nion] = den[n];
  }

  return (0);

}
