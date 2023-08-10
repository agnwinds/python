/* ****************************************************************************************************************** */
/**
 *  @file matrix_gpu.cu
 *  @author Edward J. Parkinson (e.parkinson@soton.ac.uk)
 *  @date August 2023
 *
 *  @brief
 *
 *  ***************************************************************************************************************** */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "log.h"

#if CUDA_ON

#include <cuda.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>

// #include "matrix.h"

cusolverDnHandle_t cusolver_handle = NULL;

/* ****************************************************************************************************************** */
/**
 *  @brief Check the return status of a CUDA function
 *
 *  @param [in] status  the status to check
 *
 *  ***************************************************************************************************************** */

#define CUDA_CHECK(status)                                                                                             \
  do {                                                                                                                 \
    cudaError_t err = status;                                                                                          \
    if (err != cudaSuccess) {                                                                                          \
      printf((char *)"CUDA Error: %s\n", cudaGetErrorString(err));                                                              \
      return(EXIT_FAILURE);                                                                                              \
    }                                                                                                                  \
  } while (0)

/* ****************************************************************************************************************** */
/**
 *  @brief Check the return status of a cuSOLVER function
 *
 *  @param [in] status  the status to check
 *
 *  ***************************************************************************************************************** */

#define CUSOLVER_CHECK(status)                                                                                         \
  do {                                                                                                                 \
    cusolverStatus_t err = status;                                                                                     \
    if (err != CUSOLVER_STATUS_SUCCESS) {                                                                              \
      printf((char *)"cuSolver Error: %d\n", err);                                                                              \
      return(EXIT_FAILURE);                                                                                              \
    }                                                                                                                  \
  } while (0)

/* ****************************************************************************************************************** */
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
 *  ***************************************************************************************************************** */

extern "C" int
gpu_solve_linear_system (double *a_matrix, double *b_vector, int size, double *x_vector)
{
  if (cusolver_handle == NULL)
  {
    CUSOLVER_CHECK (cusolverDnCreate (&cusolver_handle));
  }

  printf ((char *) "We are in the GPU function\n");

  // Device variables
  double *d_A, *d_b;
  int *devInfo;
  int lwork;
  double *d_work;

  // Allocate memory on the GPU
  cudaMalloc ((void **) &d_A, size * size * sizeof (double));
  cudaMalloc ((void **) &d_b, size * sizeof (double));
  cudaMalloc ((void **) &devInfo, sizeof (int));

  // Transfer data to the GPU
  cudaMemcpy (d_A, a_matrix, size * size * sizeof (double), cudaMemcpyHostToDevice);
  cudaMemcpy (d_b, b_vector, size * sizeof (double), cudaMemcpyHostToDevice);

  // Perform LU factorization
  cusolverDnDgetrf_bufferSize (cusolver_handle, size, size, d_A, size, &lwork);
  cudaMalloc ((void **) &d_work, lwork * sizeof (double));

  int *d_pivot;                 // device array of pivoting sequence
  cudaMalloc ((void **) &d_pivot, size * sizeof (int));

  cusolverDnDgetrf (cusolver_handle, size, size, d_A, size, d_work, d_pivot, devInfo);

  // Solve the linear system
  // cusolverDnDgetrs (cusolver_handle, CUBLAS_OP_N, size, 1, d_A, size, d_pivot, d_b, size, devInfo);
  cusolverDnDgetrs (cusolver_handle, CUBLAS_OP_T, size, 1, d_A, size, d_pivot, d_b, size, devInfo);

  // Transfer the solution back to the host
  cudaMemcpy (x_vector, d_b, size * sizeof (double), cudaMemcpyDeviceToHost);

  // Clean up
  cudaFree (d_A);
  cudaFree (d_b);
  cudaFree (d_work);
  cudaFree (d_pivot);
  cusolverDnDestroy (cusolver_handle);

  return EXIT_SUCCESS;
}

/* ****************************************************************************************************************** */
/**
 * @brief
 *
 * @param  [in]  a_matrix
 * @param  [out] a_inverse
 * @param  [in]  num_rows
 *
 * @return an integer representing the error state
 *
 * @details
 *
 *  ***************************************************************************************************************** */

extern "C" int
gpu_invert_matrix (double *matrix, double *inverse, int num_rows)
{
  if (cusolver_handle == NULL)
  {
    CUSOLVER_CHECK (cusolverDnCreate (&cusolver_handle));
  }

  return EXIT_SUCCESS;
}

#endif
