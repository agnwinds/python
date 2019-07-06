
/***********************************************************/
/** @file   setup_domains.c
 * @author ksl
 * @date   January, 2018
 * @brief  Get the parameters needed to describe winds
 *
 * File containing several routines that collectively
 * define the components to a wind in python. Each component
 * of the wind is said to be a domain, and the information
 * for each domain is stored in the elements of zdom
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "python.h"



/**********************************************************/
/** 
 * @brief       Get inputs that describe a particular component of the wind
 *
 * @param [in] ndom  The number (begining with 0) of this particular domain
 * @return  0 
 *
 * Sets up the one domain, which includes defining wind type, e.g whether
 * it is a shell, or a biconical flow, or an imported model, as well
 * as the type of coordinate system and its simensions.
 *
 * If the wind is to be imported from a file, it is imported in this
 * routine.
 *
 * ###Notes###
 * 1801 -   Refactored into this file in 1801.  Updates in the
 *          fall of 17 were made to allow for importing models.
 *          Note that cyl_var coordinates are not currently 
 *          working
***********************************************************/


int
get_domain_params (ndom)
     int ndom;
{
  int input_int;
  char answer[LINELENGTH];

  if (ndom >= geo.ndomain)
  {
    Error ("Trying to get grid params for a non-existent domain!\n");
    Exit (0);
  }


  strcpy (answer, "SV");
  zdom[ndom].wind_type = rdchoice ("Wind.type(SV,star,hydro,corona,kwd,homologous,yso,shell,imported)", "0,1,3,4,5,6,7,9,10", answer);


  strcat (zdom[ndom].name, "Wind");


  input_int = 1;


  /* Define the coordinate system for the grid and allocate memory for the wind structure */
  strcpy (answer, "cylindrical");
  zdom[ndom].coord_type = rdchoice ("Wind.coord_system(spherical,cylindrical,polar,cyl_var)", "0,1,2,3", answer);

  if (zdom[ndom].wind_type == IMPORT)
  {
    import_wind (ndom);
  }
  else
  {
    rdint ("Wind.dim.in.x_or_r.direction", &zdom[ndom].ndim);
    if (zdom[ndom].coord_type)
    {
      rdint ("Wind.dim.in.z_or_theta.direction", &zdom[ndom].mdim);
      if (zdom[ndom].mdim < 4)
      {
        Error ("python: domain mdim must be at least 4 to allow for boundaries\n");
        Exit (0);
      }
    }
    else
      zdom[ndom].mdim = 1;

  }

/* Check that NDIM_MAX is greater than NDIM and MDIM.  */

  if ((zdom[ndom].ndim > NDIM_MAX) || (zdom[ndom].mdim > NDIM_MAX))
  {
    Error ("NDIM_MAX %d is less than NDIM %d or MDIM %d. Fix in python.h and recompile\n", NDIM_MAX, zdom[ndom].ndim, zdom[ndom].mdim);
    Exit (0);
  }


  /* If we are in advanced then allow the user to modify scale lengths */
  if (modes.iadvanced)
  {
    strcpy (answer, "no");
    modes.adjust_grid = rdchoice ("@Diag.adjust_grid(yes,no)", "1,0", answer);

    if (modes.adjust_grid)
    {
      Log ("You have opted to adjust the grid scale lengths\n");
      rddoub ("@geo.xlog_scale", &zdom[ndom].xlog_scale);
      if (zdom[ndom].coord_type != SPHERICAL)
        rddoub ("@geo.zlog_scale", &zdom[ndom].zlog_scale);
    }
  }

  zdom[ndom].ndim2 = zdom[ndom].ndim * zdom[ndom].mdim;


  return (0);
}




/**********************************************************/
/**   
 * @brief       Get detailed particular component of the wind
 *
 * @param [in] ndom  The number (begining with 0) of this particular domain
 * @return  0 
 *
 * Continues the setup of a single domain begun in get_domain_params.    
 *
 * Much of this routine is a steering routine that calls other
 * subroutines depending on the type of wind, e.g sv, for this
 * particular wind domain.
 *
 *
 * ###Notes###
 * 1801 -   Refactored into this file in 1801.  

***********************************************************/

int
get_wind_params (ndom)
     int ndom;
{

  /* Now get parameters that are specific to a given wind model

     Note: When one adds a new model, the only things that should be read in and modified
     are parameters in geo.  This is in order to preserve the ability to continue a calculation
     with the same basic wind geometry, without reading in all of the input parameters.  
   */

  zdom[ndom].rmax = 0;

  if (zdom[ndom].wind_type == STAR)
  {
    get_stellar_wind_params (ndom);
  }
  else if (zdom[ndom].wind_type == SV)
  {
    get_sv_wind_params (ndom);
  }
  else if (zdom[ndom].wind_type == HYDRO)
  {
    get_hydro_wind_params (ndom);
  }
  else if (zdom[ndom].wind_type == CORONA)
  {
    get_corona_params (ndom);
  }
  else if (zdom[ndom].wind_type == KNIGGE)
  {
    get_knigge_wind_params (ndom);
  }
  else if (zdom[ndom].wind_type == HOMOLOGOUS)
  {
    get_homologous_params (ndom);
  }
  else if (zdom[ndom].wind_type == YSO)
  {
    get_yso_wind_params (ndom);
  }
  else if (zdom[ndom].wind_type == SHELL)       //NSH 18/2/11 This is a new wind type to produce a thin shell.
  {
    get_shell_wind_params (ndom);
  }
  else if (zdom[ndom].wind_type == IMPORT)      //Read in the wind model.
  {
    get_import_wind_params (ndom);
  }
  else
  {
    Error ("get_wind_parameters: Unknown wind type %d\n", zdom[ndom].wind_type);
    Exit (0);
  }

  /* For many models, zdom[ndom].rmax is defined within the get_parameter file, e.g
   * for stellar wind, but for others zdom[ndom].rmax can be set independently of
   * the details of the model for thesse, we clean up at the end. Ultimately
   * we may want to push all of this into ghe various get_..._params routines
   * but for now we just check if zdom[ndom].rmax has been set already and if not
   * we ask a generic quesiton here.  Because we ultimately may want to change this
   * we throw an error at this point
   *
   * These are models known not to have zdom[ndom].rmax defined as part of the model
   *
   * sv, knigge, 
   *
   * Others currently have it defined within get_whatever_params
   *
   * stellar, homologous, shell, corona
   */

  if (zdom[ndom].rmax == 0)
  {
    Error ("get_wind_params: zdom[ndom].rmax 0 for wind type %d\n", zdom[ndom].wind_type);


    zdom[ndom].rmax = 1e12;

    if (geo.system_type == SYSTEM_TYPE_AGN)
    {
      zdom[ndom].rmax = 50. * geo.r_agn;
    }


    rddoub ("Wind.radmax(cm)", &zdom[ndom].rmax);
  }

  if (zdom[ndom].rmax <= zdom[ndom].rmin)
  {
    Error ("get_wind_parameters: rmax (%10.4e) less than or equal to rmin %10.4e in domain %d\n", zdom[ndom].rmax, zdom[ndom].rmin, ndom);
    Exit (0);
  }

  zdom[ndom].twind = 40000;
  rddoub ("Wind.t.init", &zdom[ndom].twind);


  /* Next lines are to assure that we have the largest possible value of the 
   * sphere surrounding the system
   * JM 1710 -- if this is the first domain, then initialise geo.rmax see #305
   */
  if ((ndom == 0) || (zdom[ndom].rmax > geo.rmax))
  {
    geo.rmax = zdom[ndom].rmax;
  }
  geo.rmax_sq = geo.rmax * geo.rmax;

  /* Get the filling factor of the wind */

  zdom[ndom].fill = 1.;

  /* JM 1606 -- the filling factor is now specified on a domain by domain basis. See #212
     XXX allows any domain to be allowed a filling factor but this should be modified when
     we know what we are doing with inputs for multiple domains. Could create confusion */

  rddoub ("Wind.filling_factor(1=smooth,<1=clumped)", &zdom[ndom].fill);

  return (0);
}




/**********************************************************/
/** 
 * @brief      sets up the windcones for each domain
 *
 * @return     Always returns 0
 *
 * @details
 * 
 * This routine takes input variables, the minimim and maximum
 * radius of the wind at the disk, and the flow angles bounding
 * the wind, and calculates the variables used to define the cones
 * that bound each wind domain involving biconical flows.
 *
 * ### Notes ###
 *
 * The routine cycles through all of the existing domains, and
 * uses variables which have been read in or entered previously.
 *
 * The input variables that are used are typicall, wind_rho_min, wind_rho_max
 * and wind_thetamin and max.  They are defined in routines like,
 * get_sv_parameters.
 *
 * The angles thetamin and
 * thetamax are all defined from the z axis, so that an angle of 0
 * is a flow that is perpindicular to to the disk and one that is
 * close to 90 degrees will be parallel to the plane of the disk
 *
 * The routine files a structure of type cone, that contain, an element
 * z,  the place where the windcone intercepts the z axis, and an element
 * dzdr is the slope of the cone
 *
 * Wind cones are defined azimuthally around the z axis.
 *
 * For models that are not biconical flows, the windcones are set 
 * to include the entire domain.
 * 
 *
 **********************************************************/

int
setup_windcone ()
{
  int ndom;

  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {

    if (zdom[ndom].wind_thetamin > 0.0)
    {
      zdom[ndom].windcone[0].dzdr = 1. / tan (zdom[ndom].wind_thetamin);
      zdom[ndom].windcone[0].z = (-zdom[ndom].wind_rho_min / tan (zdom[ndom].wind_thetamin));
    }
    else
    {
      zdom[ndom].windcone[0].dzdr = VERY_BIG;
      zdom[ndom].windcone[0].z = -VERY_BIG;;
    }


    if (zdom[ndom].wind_thetamax > 0.0)
    {
      zdom[ndom].windcone[1].dzdr = 1. / tan (zdom[ndom].wind_thetamax);
      zdom[ndom].windcone[1].z = (-zdom[ndom].wind_rho_max / tan (zdom[ndom].wind_thetamax));
    }
    else
    {
      zdom[ndom].windcone[1].dzdr = VERY_BIG;
      zdom[ndom].windcone[1].z = -VERY_BIG;;
    }
  }
  return (0);
}
