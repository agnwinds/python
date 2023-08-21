/* ****************************************************************************************************************** */
/**
 *  @file matrix_gpu.cu
 *  @author Edward J. Parkinson (e.parkinson@soton.ac.uk)
 *  @date August 2023
 *
 *  @brief
 *
 *  ***************************************************************************************************************** */

#include <stdlib.h>
#include <stdio.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>

/* NVCC is a C++ compiler at heart, so anything we re-use from regular C source
   has to be defined here with `extern "C"` to tell the compiler that the
   function has been compiled by a C compiler (and does some computer science
   stuff to make linking possible) */

extern "C" int Exit (int error_code);
extern "C" int Error (const char *format, ...);
extern "C" int Log (const char *format, ...);
extern "C" int Debug (const char *format, ...);

/* `cusolver_handle` is a variable used to interact with the cuSolver/CUDA
    runtime and is used to initialise and clean up the resources required for
    both runtimes */

static cusolverDnHandle_t CUSOLVER_HANDLE = NULL;

/* ****************************************************************************************************************** */
/**
 *  @brief  Convert a cuSolver status into an error string
 *
 *  @param [in] status The cuSolver status to get the error string for
 *
 *  @return  The error string
 *
 *  @details
 *
 *  This function aims to emulate the beahaviour (and naming style) of cudaGetErrorString. The cuSolver library does not
 *  have a function like this, so we had to create it ourselves.
 *
 *  ***************************************************************************************************************** */

extern "C" const char *
cusolver_get_error_string (int error)
{
  cusolverStatus_t status = (cusolverStatus_t) error;

  switch (status)
  {
  case CUSOLVER_STATUS_SUCCESS:
    return "The operation completed successfully";
  case CUSOLVER_STATUS_NOT_INITIALIZED:
    return "The cuSolver library was not initialized";
  case CUSOLVER_STATUS_ALLOC_FAILED:
    return "Resource allocation failed inside the cuSolver library";
  case CUSOLVER_STATUS_INVALID_VALUE:
    return "An unsupported value or parameter was passed to the function (e.g. negative vector size)";
  case CUSOLVER_STATUS_ARCH_MISMATCH:
    return "The function requires a feature abscent from the device architecture";
  case CUSOLVER_STATUS_EXECUTION_FAILED:
    return "The GPU program failed to executed";
  case CUSOLVER_STATUS_INTERNAL_ERROR:
    return "An internal cuSolver operation failed";
  case CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
    return "The matrix type is not supported by this function";
  default:
    return "There was an unknown failure";
  }
}

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
      Error("CUDA Error (%d): %s (%s:%d)\n", err, cudaGetErrorString(err), __FILE__, __LINE__);                        \
      return err;                                                                                                      \
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
      Error("cuSolver Error (%d): %s (%s:%d)\n", err, cusolver_get_error_string(err), __FILE__, __LINE__);     \
      return err;                                                                                                      \
    }                                                                                                                  \
  } while (0)

/* ****************************************************************************************************************** */
/**
 * @brief Intialise CUDA and libraries.
 *
 * @return error status
 *
 * @details
 *
 * This function initialises CUDA and cuSolver (using the global cusolver_handle). If MPI is enabled, a diagnostic
 * message is printed stating if the MPI library is CUDA-aware or not.
 *
 *  ***************************************************************************************************************** */

extern "C" int
cuda_init (void)
{
  if (CUSOLVER_HANDLE)
    return EXIT_SUCCESS;

#ifdef MPI_ON
#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
  Log ("This MPI library has CUDA-aware support.\n");
#elif defined(MPIX_CUDA_AWARE_SUPPORT) && !MPIX_CUDA_AWARE_SUPPORT
  Log ("This MPI library does not have CUDA-aware support.\n");
#else
  Log ("This MPI library cannot determine if there is CUDA-aware support.\n");
#endif
#endif

  CUSOLVER_CHECK (cusolverDnCreate (&CUSOLVER_HANDLE));
  Debug ("cuSolver initialised\n");

  return EXIT_SUCCESS;
}

/* ****************************************************************************************************************** */
/**
 * @brief Clean up CUDA
 *
 * @return error status
 *
 * @details
 *
 * Destroyed the cuSolver handle and does any other cleaning up required.
 *
 *  ***************************************************************************************************************** */

extern "C" int
cuda_finish (void)
{
  if (!CUSOLVER_HANDLE)
    return EXIT_SUCCESS;

  CUSOLVER_CHECK (cusolverDnDestroy (CUSOLVER_HANDLE));
  CUSOLVER_HANDLE = NULL;
  Debug ("cuSolver destroyed\n");

  return EXIT_SUCCESS;
}

/* ****************************************************************************************************************** */
/**
 * @brief Transpose a square row major matrix to column major
 *
 * @param [in] row_major the input matrix
 * @param [out] column_major the output matrix
 * @param [in] matrix_size the size of the matrix
 *
 * @details
 *
 *  ***************************************************************************************************************** */

__global__ void
tranpose_row_to_column_major (double *row_major, double *column_major, int matrix_size)
{
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int idy = blockIdx.y * blockDim.y + threadIdx.y;

  if (idx < matrix_size && idy < matrix_size)
  {
    column_major[idx * matrix_size + idy] = row_major[idy * matrix_size + idx];
  }
}

/* ****************************************************************************************************************** */
/**
 * @brief Transpose a square column major matrix to row major
 *
 * @param [in] column_major the input matrix
 * @param [out] row_major the output matrix
 * @param [in] matrix_size the size of the matrix
 *
 * @details
 *
 *  ***************************************************************************************************************** */

__global__ void
tranpose_column_to_row_major (double *column_major, double *row_major, int matrix_size)
{
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int idy = blockIdx.y * blockDim.y + threadIdx.y;

  if (idx < matrix_size && idy < matrix_size)
  {
    row_major[idy * matrix_size + idx] = column_major[idx * matrix_size + idy];
  }
}

/* ****************************************************************************************************************** */
/**
 * @brief Create an identity matrix of size `size`
 *
 * @param [out] d_identity A pointer to device memory of `size` * `size` elements
 * @param [in] size  The size of the identity matrix (rows and columns)
 *
 * @details
 *
 *  ***************************************************************************************************************** */

__global__ void
create_identity_matrix (double *d_identity, int size)
{
  const int col = threadIdx.x + blockIdx.x * blockDim.x;
  const int row = threadIdx.y + blockIdx.y * blockDim.y;

  if (row < size && col < size)
  {
    d_identity[row * size + col] = (row == col) ? 1.0 : 0.0;
  }
}

/* ****************************************************************************************************************** */
/**
 * @brief  Solve the linear system A x = b, for the vector x
 *
 * @param  [in]  a_matrix  a square matrix on the LHS
 * @param  [in]  b_vector  the B vector on the RHS
 * @param  [in]  size  the number of rows (and columns) in the square matrix matrix and vectors
 * @param  [out] x_vector  the x vector (solution)
 *
 * @return an integer representing the error state
 *
 * @details
 * Performs LU decomposition to solve for x in the linear system A x = b. The calculation is perform, in serial, on
 * the CPU using GSL.
 *
 * The matrix `a_matrix` is transposed before being processed by cuSolver. This is because cuSolver (I think, it doesn't
 * seem to be documented but cuBLAS, based on the Fortran BLAS library, is column major) stores matrices in
 * column-major order, whereas arrays are typically row-major.
 *
 *  ***************************************************************************************************************** */

extern "C" int
gpu_solve_matrix (double *a_matrix, double *b_vector, int matrix_size, double *x_vector)
{
  int *devInfo;
  int lwork;
  int *d_pivot;                 /* device array of pivoting sequence */
  double *d_A, *d_Arow, *d_b;   /* d_A_row is row-major (from C) */
  double *d_work;               /* cuSolver needs a "workspace" to do stuff, which we have to allocate manually */

  /* Allocate memory on the GPU (device) to store the matrices/vectors */
  CUDA_CHECK (cudaMalloc ((void **) &d_A, matrix_size * matrix_size * sizeof (double)));
  CUDA_CHECK (cudaMalloc ((void **) &d_Arow, matrix_size * matrix_size * sizeof (double)));
  CUDA_CHECK (cudaMalloc ((void **) &d_b, matrix_size * sizeof (double)));
  CUDA_CHECK (cudaMalloc ((void **) &devInfo, sizeof (int)));
  CUDA_CHECK (cudaMalloc ((void **) &d_pivot, matrix_size * sizeof (int)));

  /* Copy the matrix and vector to the device memory */
  CUDA_CHECK (cudaMemcpy (d_Arow, a_matrix, matrix_size * matrix_size * sizeof (double), cudaMemcpyHostToDevice));
  CUDA_CHECK (cudaMemcpy (d_b, b_vector, matrix_size * sizeof (double), cudaMemcpyHostToDevice));

  dim3 blockDim (16, 16);
  dim3 gridDim ((matrix_size + blockDim.x - 1) / blockDim.x, (matrix_size + blockDim.y - 1) / blockDim.y);
  tranpose_row_to_column_major <<< gridDim, blockDim >>> (d_Arow, d_A, matrix_size);
  CUDA_CHECK (cudaFree (d_Arow));

  /* XXXX_bufferSize is used to compute the size of the workspace we need, and depends on the size of the linear
     system being solved */
  CUSOLVER_CHECK (cusolverDnDgetrf_bufferSize (CUSOLVER_HANDLE, matrix_size, matrix_size, d_A, matrix_size, &lwork));
  CUDA_CHECK (cudaMalloc ((void **) &d_work, lwork * sizeof (double)));

  /* Perform LU factorization and solve the linear system. The vector d_b is not used in `getrs` (the solver), but
     it's the same size of the solution vector so we'll re-use that. d_b is then copied back to host memory (CPU RAM) */
  CUSOLVER_CHECK (cusolverDnDgetrf (CUSOLVER_HANDLE, matrix_size, matrix_size, d_A, matrix_size, d_work, d_pivot, devInfo));
  CUSOLVER_CHECK (cusolverDnDgetrs (CUSOLVER_HANDLE, CUBLAS_OP_N, matrix_size, 1, d_A, matrix_size, d_pivot, d_b, matrix_size, devInfo));

  /* don't need to tranpose vectors */
  CUDA_CHECK (cudaMemcpy (x_vector, d_b, matrix_size * sizeof (double), cudaMemcpyDeviceToHost));

  CUDA_CHECK (cudaFree (d_A));
  CUDA_CHECK (cudaFree (d_b));
  CUDA_CHECK (cudaFree (d_work));
  CUDA_CHECK (cudaFree (d_pivot));

  return EXIT_SUCCESS;
}

/* ****************************************************************************************************************** */
/**
 * @brief  Compute the inverse of a square matrix
 *
 * @param  [in]  matrix  the matrix to inverse
 * @param  [out] inverse_matrix  the inverse of `matrix`
 * @param  [in]  matrix_size  the size of `matrix` and `inverse_matrix`
 *
 * @return an integer representing the error state
 *
 * @details
 *
 * This function calculates the inverse of `matrix` through LU decomposition. We are essentially solving the system
 * A A^-1 = I, where I is the indentity matrix. In this case, A = S = LU. We compute it this way because we can solve
 * the matrix inverse (A A^-1 = I) directly with forward/backward substitution using suSolver.
 *
 * The matrix `matrix` is transposed before being processed by cuSolver. This is because cuSolver (I think, it doesn't
 * seem to be documented but cuBLAS, based on the Fortran BLAS library, is column major) stores matrices in
 * column-major order, whereas arrays are typically row-major.
 *
 *  ***************************************************************************************************************** */

extern "C" int
gpu_invert_matrix (double *matrix, double *inverse_matrix, int matrix_size)
{
  int *d_pivot;
  int work_size;
  int *dev_info;
  double *d_matrix;
  double *d_matrix_row;
  double *d_identity;
  double *d_identity_row;
  double *d_workspace;

  /* Allocate memory on the GPU (device) to store the matrices/vectors */
  CUDA_CHECK (cudaMalloc ((void **) &d_matrix_row, matrix_size * matrix_size * sizeof (double)));
  CUDA_CHECK (cudaMalloc ((void **) &d_matrix, matrix_size * matrix_size * sizeof (double)));
  CUDA_CHECK (cudaMalloc ((void **) &d_pivot, matrix_size * sizeof (int)));
  CUDA_CHECK (cudaMalloc ((void **) &d_identity, matrix_size * matrix_size * sizeof (double)));
  CUDA_CHECK (cudaMalloc ((void **) &d_identity_row, matrix_size * matrix_size * sizeof (double)));
  CUDA_CHECK (cudaMalloc ((void **) &dev_info, sizeof (int)));

  /* XXXX_bufferSize is used to compute the size of the workspace we need, and depends on the size of the matrix */
  CUSOLVER_CHECK (cusolverDnDgetrf_bufferSize (CUSOLVER_HANDLE, matrix_size, matrix_size, d_matrix, matrix_size, &work_size));
  CUDA_CHECK (cudaMalloc ((void **) &d_workspace, work_size * sizeof (double)));

  /* Copy matrix from host (CPU) to device memory (d_matrix) */
  dim3 blockDim (16, 16);
  dim3 gridDim ((matrix_size + blockDim.x - 1) / blockDim.x, (matrix_size + blockDim.y - 1) / blockDim.y);
  CUDA_CHECK (cudaMemcpy (d_matrix_row, matrix, matrix_size * matrix_size * sizeof (double), cudaMemcpyHostToDevice));
  tranpose_row_to_column_major <<< gridDim, blockDim >>> (d_matrix_row, d_matrix, matrix_size);
  CUDA_CHECK (cudaFree (d_matrix_row));

  /* We'll use a CUDA kernel to create an indentity matrix, which we'll use to solve for the inverse. The
     syntax is a bit strange. We can imagine this bit as being similar to OpenMP, as CUDA is shared memory. We're
     sending work to a bunch of GPU cores, which will construct their own element in the identity matrix */
  create_identity_matrix <<< gridDim, blockDim >>> (d_identity, matrix_size);

  /* We first need to facotrise the matrix to get the pivot indcies, for getrs. The function getrs is solves a linear
     system to solve for the inverse matrix. The inverse matrix is placed back into `d_identity` */
  CUSOLVER_CHECK (cusolverDnDgetrf (CUSOLVER_HANDLE, matrix_size, matrix_size, d_matrix, matrix_size, d_workspace, d_pivot, dev_info));
  CUSOLVER_CHECK (cusolverDnDgetrs
                  (CUSOLVER_HANDLE, CUBLAS_OP_N, matrix_size, matrix_size, d_matrix, matrix_size, d_pivot, d_identity, matrix_size,
                   dev_info));

  /* Copy the inverse matrix from host to device */
  tranpose_column_to_row_major <<< gridDim, blockDim >>> (d_identity, d_identity_row, matrix_size);
  CUDA_CHECK (cudaMemcpy (inverse_matrix, d_identity_row, matrix_size * matrix_size * sizeof (double), cudaMemcpyDeviceToHost));

  CUDA_CHECK (cudaFree (d_matrix));
  CUDA_CHECK (cudaFree (d_identity));
  CUDA_CHECK (cudaFree (d_workspace));
  CUDA_CHECK (cudaFree (dev_info));

  return EXIT_SUCCESS;
}
