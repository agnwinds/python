
/***********************************************************/
/** @file   setup_line_transfer.c
 * @author ksl
 * @date   October, 2018
 * @brief  Get the parameters needed to describe line transer
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>

#include "atomic.h"
#include "python.h"


/**********************************************************/
/** 
 * @brief       Get the line transfer mode for the wind and 
 * several other variables related to line transfer
 *
 * @param [in] None
 * @return  0 
 *
 * This rontinues simply gets the line tranfer mode for
 * all componensts of the wind.  After logging this
 * information, the routine also reads in the atomic 
 * data, and the variable that establishes in the 
 * non-macro mode whether the wind is allowed to radiate
 *
 *
 * ###Notes###
 * 
 * Rather than having inputs for several aspects of line transfer
 * the choices with regard to line_transfer_mode are used to
 * define several variables, geo.rt_mode, geo.line_mode, and geo.scatter_mode.
 * within the routine.  
***********************************************************/



int
get_line_transfer_mode ()
{
  int rc;                       // Return code from running Setup_Py_Dir
  struct stat file_stat;        // Used to check the atomic data exists
  char answer[LINELENGTH];

  int user_line_mode = 0;

  /* 181010-ksl-There a several line transfer modes which are diagnostic in nature which 
     cannot be reached by rdchoice easily, except as numbers.  
     * We need a better way to deal with these. */

  strcpy (answer, "thermal_trapping");
  user_line_mode =
    rdchoice
    ("Line_transfer(pure_abs,pure_scat,sing_scat,escape_prob,thermal_trapping,macro_atoms_escape_prob,macro_atoms_thermal_trapping)",
     "0,1,2,3,5,6,7", answer);

  /* JM 1406 -- geo.rt_mode and geo.macro_simple control different things. geo.rt_mode controls the radiative
     transfer and whether or not you are going to use the indivisible packet constraint, so you can have all simple 
     ions, all macro-atoms or a mix of the two. geo.macro_simple just means one can turn off the full macro atom 
     treatment and treat everything as 2-level simple ions inside the macro atom formalism */

  /* Set the default scattering and RT mode */
  geo.scatter_mode = SCATTER_MODE_ISOTROPIC;    // isotropic
  geo.rt_mode = RT_MODE_2LEVEL; // Not macro atom (SS)


  if (user_line_mode == 0)
  {
    Log ("Line_transfer mode:  Simple, pure absorption, no scattering\n");
    geo.line_mode = user_line_mode;
  }
  else if (user_line_mode == 1)
  {
    Log ("Line_transfer mode:  Simple, pure scattering, no absoprtion\n");
    geo.line_mode = user_line_mode;
  }
  else if (user_line_mode == 2)
  {
    Log ("Line_transfer mode:  Simple, single scattering, with absorption, but without escape probality treatment\n");
    geo.line_mode = user_line_mode;
  }
  else if (user_line_mode == 3)
  {
    Log ("Line_transfer mode:  Simple, isotropic scattering, escape probabilities\n");
    geo.line_mode = user_line_mode;
  }
  else if (user_line_mode == 5)
  {
    Log ("Line_transfer mode:  Simple, thermal trapping, Single scattering \n");
    geo.scatter_mode = SCATTER_MODE_THERMAL;    // Thermal trapping model
    geo.line_mode = 3;
    geo.rt_mode = RT_MODE_2LEVEL;       // Not macro atom (SS) 
  }
  else if (user_line_mode == 6)
  {
    Log ("Line_transfer mode:  macro atoms, isotropic scattering  \n");
    geo.scatter_mode = SCATTER_MODE_ISOTROPIC;  // isotropic
    geo.line_mode = 3;
    geo.rt_mode = RT_MODE_MACRO;        // Identify macro atom treatment (SS)
    geo.macro_simple = 0;       // We don't want the all simple case (SS)
  }
  else if (user_line_mode == 7)
  {
    Log ("Line_transfer mode:  macro atoms, anisotropic  scattering  \n");
    geo.scatter_mode = SCATTER_MODE_THERMAL;    // thermal trapping
    geo.line_mode = 3;
    geo.rt_mode = RT_MODE_MACRO;        // Identify macro atom treatment (SS)
    geo.macro_simple = 0;       // We don't want the all simple case (SS)
  }
  else if (user_line_mode == 8)
  {
    Log ("Line_transfer mode: simple macro atoms, isotropic  scattering  \n");
    geo.scatter_mode = SCATTER_MODE_ISOTROPIC;  // isotropic
    geo.line_mode = 3;
    geo.rt_mode = RT_MODE_MACRO;        // Identify macro atom treatment i.e. indivisible packets
    geo.macro_simple = 1;       // This is for test runs with all simple ions (SS)
  }
  else if (user_line_mode == 9)
  {
    Log ("Line_transfer mode: simple macro atoms, anisotropic  scattering  \n");
    geo.scatter_mode = SCATTER_MODE_THERMAL;    // thermal trapping
    geo.line_mode = 3;
    geo.rt_mode = RT_MODE_MACRO;        // Identify macro atom treatment i.e. indivisible packets
    geo.macro_simple = 1;       // This is for test runs with all simple ions (SS)
  }
  else
  {
    Error ("Unknown line_transfer mode %d\n");
    line_transfer_help_message ();
    Exit (0);
  }

  /* With the macro atom approach we won't want to generate photon 
     bundles in the wind so switch it off here. If not in macro atom mode
     ask the user whether we want the wind to radiate
   */
  if (geo.rt_mode == RT_MODE_MACRO)
  {
    geo.wind_radiation = FALSE;
  }
  else
  {
    strcpy (answer, "yes");
    geo.wind_radiation = rdchoice ("Wind.radiation(yes,no)", "1,0", answer);

  }

  /* Note the only other variable read in in this section is that for the atomic data,
     which can only be specified for a new model, becuase of various structures
     that have been allocated.  
   */

  if (geo.run_type == RUN_TYPE_NEW)
  {


    /* read in the atomic data */
    rdstr ("Atomic_data", geo.atomic_filename);

    /* read a variable which controls whether to save a summary of atomic data
       this is defined in atomic.h, rather than the modes structure */

    if (modes.iadvanced)
    {


      strcpy (answer, "no");
      write_atomicdata = rdchoice ("@Diag.write_atomicdata(yes,no)", "1,0", answer);
      if (write_atomicdata)
        Log ("You have opted to save a summary of the atomic data\n");
    }

    /*
     * Check that geo.atomic_filename exists - i.e. that the directory is readable
     * and in the directory Python is being executed from. If it isn't - then
     * try to run Setup_Py_Dir. If both fail, then warn the user and exit Python
     */

    if (stat (geo.atomic_filename, &file_stat))
    {
      Log ("Unable to open atomic masterfile %s\n", geo.atomic_filename);
      Log ("Running Setup_Py_Dir to try and fix the situation\n");
      rc = system ("Setup_Py_Dir");
      if (rc)
      {
        Error ("Unable to open %s and run Setup_Py_Dir\n", geo.atomic_filename);
        Exit (1);
      }
    }

    get_atomic_data (geo.atomic_filename);

    /* throw a fatal error if there are macro-atom levels but rt_mode is non macro */
    if (nlevels_macro > 0 && geo.rt_mode != RT_MODE_MACRO)
    {
      Error ("Fatal error: you specified macro-atom data but standard line transfer. Not supported.\n");
      Log ("Try changing either to a simple data set or to macro-atom line transfer\n");
      Exit (1);
    }
  }

  return (0);
}


/**********************************************************/
/** 
 * @brief      print out a help message about line transfer modes
 *
 * @return     Always returns 0
 **********************************************************/

int
line_transfer_help_message ()
{
  char *some_help;

  some_help = "\
\n\
Available line transfer modes and descriptions are: \n\
\n\
  pure_abs Pure Absorption\n\
  pure_scat Pure Scattering\n\
  sing_scat Single Scattering\n\
  escape_prob Escape Probabilities, isotropic scattering\n\
  thermal_trapping Escape Probabilities, anisotropic scattering\n\
  macro_atoms Indivisible energy packets / macro-atoms, isotropic scattering\n\
  macro_atoms_thermal_trapping Indivisible energy packets / macro-atoms, anisotropic scattering\n\
  8(deprecated) Indivisible energy packets, force all simple-atoms, anisotropic scattering\n\
  9(deprecated) Indivisible energy packets, force all simple-atoms, anisotropic scattering\n\
\n\
  Standard mode is thermal_trapping for runs involving weight reduction and no macro-atoms\n\
  Standard macro-atom mode is macro_atoms_thermal_trapping\n\
\n\
See this web address for more information: https://github.com/agnwinds/python/wiki/Line-Transfer-and-Scattering\n\
\n\
\n\
";                              // End of string to provide one with help

  Log ("%s\n", some_help);

  return (0);
}
