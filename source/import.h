/* ************************************************************************** */
/**
 * @file  import.h
 * @author EJP
 * @date   Feb 2019
 *
 * @brief    Global structures for containing imported models
 *
 * @details
 * ************************************************************************** */

#define DEFAULT_IMPORT_TEMPERATURE 40000

struct
{
  int ncell;
  int *element;
  double *r;
  double *v_r;
  double *mass_rho;
  double *t_r;

} import_model_1d;

/*
 * i is the row number
 * j is the column number
 */

struct
{
  int ndim, mdim, ncell;
  int *i, *j, *inwind;
  double *x, *z, *r, *theta;
  double *v_x, *v_y, *v_z;
  double *mass_rho, *t_r;
  double *wind_x, *wind_z, *wind_midx, *wind_midz;
} import_model_2d;

