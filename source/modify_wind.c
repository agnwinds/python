
/***********************************************************/
/** @file  modify_wind.c
 * @author ksl
 * @date   June, 2020     
 *
 * @brief  Routines to modfify a wind structure to for example
 * change densities of certain ions
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


char outroot[LINELENGTH];


/**********************************************************/
/** 
 * @brief      parses the command line options
 *
 * @param [in]  int  argc   the number of command line arguments
 * @param [in]  char *  argv[]   The command line arguments
 * @return      restart_stat   1 if restarting a previous model,
 * 0 in all other cases.
 *
 * Python has a fairly rich set of command line options, which
 * are parsed by this routine
 *
 * The routine also creates the diag folder, which is where most
 * of the log files are written 
 *
 * ###Notes###
 *
 * The general purpose of each of the command line options
 * should be fairly obvious from reading the code.
 *
 * If changes to the command line interface are made they should
 * be described in the routine help 
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
  int restart_stat, verbosity, max_errors, i;
  int j = 0;
  char dummy[LINELENGTH];
  int mkdir ();
  double time_max;
  char *fgets_rc;
  double x;

  restart_stat = 0;

  if (argc == 1)
  {
    printf ("Parameter file name (e.g. my_model.pf, or just my_model):");
    fgets_rc = fgets (dummy, LINELENGTH, stdin);
    if (!fgets_rc)
    {
      printf ("Input rootname is NULL or invalid\n");
      exit (1);
    }
    get_root (files.root, dummy);
    strcpy (files.diag, files.root);
    strcat (files.diag, ".diag");
  }
  else
  {

    for (i = 1; i < argc; i++)
    {
      if (strcmp (argv[i], "-out_root") == 0)
      {
        if (sscanf (argv[i + 1], "%s", outroot) != 1)
        {
          printf ("python: Expected out_root after -out_root switch\n");
          exit (0);
        }
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
    get_root (files.root, dummy);

  }

  return (0);
}



int
main (argc, argv)
     int argc;
     char *argv[];
{

    xparse_command_line (argc, argv);

    wind_read("star.wind_save");

    wind_save("foo.wind_save");


    printf("gotcha %s\n",files.root);



}
