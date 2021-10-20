
/***********************************************************/
/** @file  inspect_wind.c
 * @author ksl
 * @date   October, 2021
 *
 * @brief  Routines to inspect variables in a wind structure 
 *
 *###Notes###
 * This is intended just as a diagnostic routine 
 * so that one can print out whatever variables in
 * a windstrucutre one wants in order to diagnose
 * a problem.  It was written so that we could inspect
 * some of the macro atom variables in paralell mode
 * in diagnosing issue #898 and #910, but anyon 
 * should change it so that other problems might 
 * be addressed.
 *
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "python.h"
//OLD #include "import.h"


char inroot[LINELENGTH], outroot[LINELENGTH], model_file[LINELENGTH];
int model_flag, ksl_flag, cmf2obs_flag, obs2cmf_flag;

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

  model_flag = ksl_flag = obs2cmf_flag = cmf2obs_flag = 0;

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
      if (strcmp (argv[i], "-model_file") == 0)
      {
        if (sscanf (argv[i + 1], "%s", dummy) != 1)
        {
          printf ("python: Expected a model file containing density, velocity and temperature after -model_file switch\n");
          exit (0);
        }
        get_root (model_file, dummy);
        i++;
        j = i;
        printf ("got a model file %s\n", model_file);
        model_flag = 1;
      }
      else if (strcmp (argv[i], "-ksl") == 0)
      {
        printf ("Carrying out a simple hard wired ion modification\n");
        ksl_flag = 1;
      }
      else if (strcmp (argv[i], "--dry-run") == 0)
      {
        modes.quit_after_inputs = 1;
        j = i;
      }
      else if (strcmp (argv[i], "-cmf") == 0)
      {
        obs2cmf_flag = 1;
      }
      else if (strcmp (argv[i], "-obs") == 0)
      {
        cmf2obs_flag = 1;
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

  char infile[LINELENGTH], outfile[LINELENGTH];
  int n, i;
  FILE *fptr, *fopen ();


  xparse_command_line (argc, argv);

  sprintf (infile, "%s.wind_save", inroot);
  sprintf (outfile, "%s.txt", inroot);

  printf ("Reading %s and writing to %s\n", infile, outfile);

  wind_read (infile);


  fptr = fopen (outfile, "w");

  fprintf (fptr, "Results for %s\n", infile);

  fprintf (fptr, "Size  %d  %d %d \n", size_Jbar_est, size_gamma_est, size_alpha_est);

  for (n = 0; n < NPLASMA; n++)
  {
    fprintf (fptr, "jbar ");
    fprintf (fptr, "%d ", n);
    for (i = 0; i < 10; i++)
    {
      fprintf (fptr, "%8.2e ", macromain[n].jbar[i]);
    }
    fprintf (fptr, "\n");
  }



  for (n = 0; n < NPLASMA; n++)
  {
    fprintf (fptr, "jbar_old ");
    fprintf (fptr, "%d ", n);
    for (i = 0; i < 10; i++)
    {
      fprintf (fptr, "%8.2e ", macromain[n].jbar_old[i]);
    }
    fprintf (fptr, "\n");
  }



  fclose (fptr);


  exit (0);

}
