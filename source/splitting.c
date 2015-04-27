
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

int importance_map_init_spherical ( 
	WindPtr w,
	double r_imp_max,
	double r_thresh
	)
{
	int i, n;

	for(n=0; n<NDIM; n++)
	{
		if(w[n].xcen[0] < r_thresh)
		{
			w[n].importance = 1.0;
		}
		else
		{
			w[n].importance = r_imp_max; 
		}
	}

	return(0);
}

int importance_map_init_cylindrical (
	WindPtr w
	)
{
	int i, j, n;

	for(i=0; i<NDIM; i++)
	{
		for(j=0; j<MDIM; j++)
		{
			wind_ij_to_n(i, j, &n);
		}
	}

	return(0);
}

int importance_map_init (
	WindPtr w,
	double r_imp_max,
	double r_thresh_r,
	double r_thresh_theta
	)
{
	if(geo.vr_ionisation && geo.vr_spectrum)
		printf("Initialising importance maps for ionisation and spectrum cycles.\n");
	else if (geo.vr_ionisation)
		printf("Initialising importance maps for ionisation cycles.\n");
	else if (geo.vr_spectrum)
		printf("Initialising importance maps for spectrum cycles.\n");

	if(geo.coord_type == SPHERICAL)
	{
		importance_map_init_spherical (w, r_imp_max, r_thresh_r);
	}

	return(0);
}