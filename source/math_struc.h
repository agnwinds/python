/***************************CDF STRUCTURE*****************************/
/* This is the structure for storing cumulative distribution functions. The CDFs are
generated from a function which is usually only proportional to the probability density
function or from an array.  It is sometimes useful, e.g. in calculating the reweighting function to
have access to the proper normalization.  


*/


#define NCDF 30000              //The default size for these arrays.  This needs to be greater than
                                //the size of any model that is read in, hence larger than NWAVE_EXTRACT in models.h
#define FUNC_CDF  2000          //The size for CDFs made from functional form CDFs
#define ARRAY_PDF 1000          //The size for PDFs to be turned into CDFs from arrays


typedef struct Cdf
{
  double x[NCDF];               /* Positions for which the CDF is calculated */
  double y[NCDF];               /* The value of the CDF at x */
  double d[NCDF];               /* 57i -- the rate of change of the CDF at x */
  double limit1, limit2;        /* Limits (running from 0 to 1) that define a portion
                                   of the CDF to sample */
  double x1, x2;                /* limits if they exist on what is returned */
  double norm;                  /* The scaling factor which would renormalize the CDF */
  int ncdf;                     /* Size of this CDF */
}
 *CdfPtr, cdf_dummy;


 struct basis
{
  double a[3][3];

};

