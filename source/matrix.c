/* ************************************************************************************************************************************** */
/**
 *  @file matrix.c
 *  @author Edward J. Parkinson (e.parkinson@soton.ac.uk)
 *  @data August 2023
 *
 *  @brief Functions for solving matrix problems, using both GSL and cuSOLVER.
 *
 *  ************************************************************************************************************************************* */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef CUDA_ON
#include <cuda.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#endif

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "atomic.h"
#include "python.h"

#ifdef CUDA_ON

/**
 *  @brief Check the return status of a CUDA function
 *
 *  @param [in] status  the status to check
 */
#define CUDA_CHECK(status)                                                                                             \
  do {                                                                                                                 \
    cudaError_t err = status;                                                                                          \
    if (err != cudaSuccess) {                                                                                          \
      Error("CUDA Error: %s\n", cudaGetErrorString(err));                                                              \
      Exit(EXIT_FAILURE);                                                                                              \
    }                                                                                                                  \
  } while (0)

/**
 *  @brief Check the return status of a cuSOLVER function
 *
 *  @param [in] status  the status to check
 */
#define CUSOLVER_CHECK(status)                                                                                         \
  do {                                                                                                                 \
    cusolverStatus_t err = status;                                                                                     \
    if (err != CUSOLVER_STATUS_SUCCESS) {                                                                              \
      Error("cuSolver Error: %d\n", err);                                                                              \
      Exit(EXIT_FAILURE);                                                                                              \
    }                                                                                                                  \
  } while (0)

#endif

/* ************************************************************************************************************************************** */
/**
 * @brief  Solve the linear system A x = b, for the vector x
 *
 * @param  [in]  a_matrix  a square matrix on the LHS
 * @param  [in]  b_vector  the B resultant vector
 * @param  [in]  size  the number of rows (and columns) in the square matrix matrix and vectors
 * @param  [out] x_vector  the x vector on the RHS
 *
 * @return an integer representing the error state
 *
 * @details
 * Performs LU decomposition to solve for x in the linear system A x = b. The calculation is perform, in serial, on
 * the CPU using GSL.
 *
 *  ************************************************************************************************************************************* */

#ifdef CUDA_ON

cusolverDnHandle_t cu_handle;

static int
gpu_solve_linear_system (double *a_matrix, double *b_vector, int size, double *x_vector)
{
  int error = EXIT_SUCCESS;

  return error;
}

#endif

/* ************************************************************************************************************************************** */
/**
 * @brief  Solve the linear system A x = b, for the vector x
 *
 * @param  [in]  a_matrix - a square matrix on the LHS
 * @param  [in]  b_vector - the B resultant vector
 * @param  [in]  size - the number of rows (and columns) in the square matrix matrix and vectors
 * @param  [out] x_vector - the x vector on the RHS
 * @param  [in]  nplasma - the index of the plasma cell we are working on
 *
 * @return an integer representing the error state
 *
 * @details
 * Performs LU decomposition to solve for x in the linear system A x = b. The calculation is perform, in serial, on
 * the CPU using GSL.
 *
 * ### Notes ###
 *
 * `nplasma` is only used to indicate a cell number, if the calculation fails.
 *
 *  ************************************************************************************************************************************* */

static int
cpu_solve_linear_system (double *a_matrix, double *b_vector, int size, double *x_vector, int nplasma)
{
  int mm, ierr, s;
  /* s is the 'sign' of the permutation - is has the value -1^n where n is the number of
     permutations. We dont use it anywhere, but in principle it can be used to refine the
     solution via gsl_linalg_LU_refine */
  double test_val;
  double lndet;
  int n_error = 0;

  gsl_permutation *p;
  gsl_matrix_view m;
  gsl_vector_view b;
  gsl_vector *test_vector, *populations;
  gsl_matrix *test_matrix;
  gsl_error_handler_t *handler;

  /* Turn off gsl error handling so that the code does not abort on error */

  handler = gsl_set_error_handler_off ();

  ierr = EXIT_SUCCESS;
  test_val = 0.0;

  /* create gsl matrix/vector views of the arrays of rates.
   * This is the structure that gsl uses do define an array.
   * It contains not only the data but the dimensions, etc.*/
  m = gsl_matrix_view_array (a_matrix, size, size);

  /* these are used for testing the solution below */
  test_matrix = gsl_matrix_alloc (size, size);
  test_vector = gsl_vector_alloc (size);

  gsl_matrix_memcpy (test_matrix, &m.matrix);   // create copy for testing


  /* gsl_vector_view_array creates the structure that gsl uses to define a vector
   * It contains the data and the dimension, and other information about where
   * the vector is stored in memory etc.
   */
  b = gsl_vector_view_array (b_vector, size);

  /* the populations vector will be a gsl vector which stores populations */
  populations = gsl_vector_alloc (size);


  /* permuations are special structures that contain integers 0 to size-1, which can
   * be manipulated */

  p = gsl_permutation_alloc (size);

  /* Decomposes m into its LU Components.  Store the L part in m and the
   * U part in s and p is modified.
   */

  ierr = gsl_linalg_LU_decomp (&m.matrix, p, &s);

  if (ierr)
  {
    Error ("Solve_matrix: gsl_linalg_LU_decomp failure %d for cell %i \n", ierr, nplasma);
    Exit (0);

  }

  ierr = gsl_linalg_LU_solve (&m.matrix, p, &b.vector, populations);

  if (ierr)
  {
    lndet = gsl_linalg_LU_lndet (&m.matrix);    // get the determinant to report to user
    Error ("Solve_matrix: gsl_linalg_LU_solve failure (%d %.3e) for cell %i \n", ierr, lndet, nplasma);

    return (4);

  }

  gsl_permutation_free (p);

  /* Check that the populations vector we have just created really is a solution to
     the matrix equation */

  /* The following line does the matrix multiplication test_vector = 1.0 * test_matrix * populations The CblasNoTrans
     statement just says we do not do anything to test_matrix, and the 0.0 means we do not add a second matrix to the result
     If the solution has worked, then test_vector should be equal to b_temp */

  ierr = gsl_blas_dgemv (CblasNoTrans, 1.0, test_matrix, populations, 0.0, test_vector);

  if (ierr != 0)
  {
    Error ("Solve_matrix: bad return (%d) when testing matrix solution to rate equations in plamsa cell.\n", ierr, nplasma);
  }

  /* now cycle through and check the solution to y = m * populations really is (1, 0, 0 ... 0) */

  for (mm = 0; mm < size; mm++)
  {

    /* get the element of the vector we want to check */
    test_val = gsl_vector_get (test_vector, mm);

    /* b_vector is (1,0,0,0..) when we do matom rates. test_val is normally something like
       1e-16 if it's supposed to be 0. We have a different error check if b_vector[mm] is 0 */


    if (b_vector[mm] > 0.0)
    {
      if (n_error == 0 && fabs ((test_val - b_vector[mm]) / test_val) > EPSILON)
      {
        Error ("Solve_matrix: test solution fails relative error for row %i %e != %e frac_error %3 in plasma cell %d\n", mm, test_val,
               b_vector[mm], fabs ((test_val - b_vector[mm]) / test_val), nplasma);
        ierr = 2;
        n_error += 1;
      }
    }
    else if ((n_error) == 0 && fabs (test_val - b_vector[mm]) > EPSILON)        // if b_vector is 0, check absolute error

    {
      Error ("Solve_matrix: test solution fails absolute error for row %i %e != %e in plasma cell %d\n", mm, test_val, b_vector[mm],
             nplasma);
      ierr = 3;
      n_error += 1;
    }
  }

  if (n_error > 1)
  {
    Error ("Solve_matrix: There were %d row errors in all for plasma cell %d\n", n_error, nplasma);
  }


  /* copy the populations to a normal array */
  for (mm = 0; mm < size; mm++)
    x_vector[mm] = gsl_vector_get (populations, mm);

  /* free memory */
  gsl_vector_free (test_vector);
  gsl_matrix_free (test_matrix);
  gsl_vector_free (populations);

  gsl_set_error_handler (handler);

  return (ierr);
}

/* ************************************************************************************************************************************** */
/**
 * @brief Solve the linear system A x = b, for the vector x
 *
 * @param  [in]  a_matrix - a square matrix on the LHS
 * @param  [in]  b_vector - the B resultant vector
 * @param  [in]  size - the number of rows (and columns) in the square matrix matrix and vectors
 * @param  [out] x_vector - the x vector on the RHS
 * @param  [in]  nplasma - the index of the plasma cell we are working on
 *
 * @return an integer representing the error state
 *
 * @details
 * Performs LU decomposition to solve for x in the linear system A x = b.
 *
 * This function is a wrapper around two others: `cpu_solve_linear_system` and `gpu_solve_linear_system`.
 *
 * ### Notes ###
 *
 * `nplasma` is only used to indicate a cell number, if the calculation fails.
 *
 *  ************************************************************************************************************************************* */
int
solve_matrix (double *a_matrix, double *b_matrix, int size, double *x_matrix, int nplasma)
{
  int error = EXIT_SUCCESS;

#ifdef CUDA_ON
  error = gpu_solve_linear_system (a_matrix, b_matrix, size, x_matrix, nplasma);
#else
  error = cpu_solve_linear_system (a_matrix, b_matrix, size, x_matrix, nplasma);
#endif

  return error;
}
