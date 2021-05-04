/** ************************************************************************* */
/**
 * @file     py_optical_depth.h
 * @author   Edward Parkinson
 * @date     February 2021
 *
 * @brief    Globals used during optical depth diagnostics.
 *
 * @details
 *
 * Do not include this in any file other than py_optical_depth_sub.c, otherwise
 * there will be double definition errors. If you want to include it in
 * multiple, place MAXDIFF and N_PI_EDGES into py_optical_depth_sub.c to avoid
 * the errors. The reason they are in this header file, and why it exists, is
 * to keep things a bit cleaner in py_optical_depth_sub.c
 *
 * ************************************************************************** */

void do_optical_depth_diagnostics(void);

#define N_FREQ_BINS 10000

// Inclination angles structure

typedef struct SightLines_s
{
  char name[50];
  double lmn[3];
} SightLines_t;

SightLines_t *INCLINATION_ANGLES;
int N_INCLINATION_ANGLES;       // The number of inclination angles

// PI edges structure

typedef struct PIEdges_s
{
  char name[50];
  double freq;
} PIEdges_t;

#define N_PI_EDGES (int) (sizeof PHOTOION_EDGES_TO_MEASURE / sizeof *PHOTOION_EDGES_TO_MEASURE)     // The number of optical depths for the simple calculation

// Control the column density extracted

enum COLUMN_DENSITY
{
  COLUMN_MODE_RHO,
  COLUMN_MODE_ION,
};

int COLUMN_MODE;
int COLUMN_MODE_ION_NUMBER;

// Control which domain to send photons from

int N_DOMAIN;

// Control how the optical depth integration is done

enum {
  RUN_MODE_OUTWARD,
  RUN_MODE_PHOTOSPHERE,
} MODE;

double TAU_DEPTH;
