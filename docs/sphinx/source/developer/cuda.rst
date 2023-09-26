Matrix Acceleration using CUDA
###############################

Python can use CUDA/GPUs to accelerate solving linear systems using the cuSOLVER library, which is part of the NVIDIA
CUDA Toolkit.

*This pilot study into using GPUs in Python was conducted as an HPC RSE project at the University of Southampton.*

When should you use CUDA?
=========================

Given the pilot study nature of this, the current matrix acceleration implementation (September 2023) is very simple.
Therefore, in most small to mid-sized models, using GPU acceleration will not improve model performance. As of writing,
CUDA is only being used to accelerate matrix calculations, so there will be no performance improvements other than in
models which are bogged down by any matrix caulations; e.g. such as when calculating the ionization state for a large
number of ions, or models with lots of macro atom levels. Even then, you may only see modest improvements in model
performance if the majority of time is spent transporting and scattering photons.

You should only prefer to use CUDA when matrices are large enough to warrant GPU acceleration. The size of where this is
important is tricky to know, as it is hardware dependent. If you are using an old CPU, then you are likely to see
improvements from matrices as small as 200 x 200. Otherwise, you may need to reach up to matrix sizes of 500 x 500 (or
larger!) before there is any tangible benefit.

This situation could be improved further with smarter memory management. Most time is spent transferring data back and
forth between the GPU and CPU. As an example, consider the matrix ionization state calculation. Currently only solving
the linear system to calculate the ionization state has been ported to the GPU, which means each iteration toward a
converged ionization state requires memory to be copied to and from the GPU. If you could instead port the entire
process onto the GPU, you no longer need to copy memory between iterations and would save a significant amount of time.

Requirements
============

To be able to use the CUDA matrix acceleration in Python, you need to meet to the following requirements,

- A CUDA-capable NVIDIA GPU
- NVIDIA CUDA Toolkit
- NVIDIA GPU drivers
- A supported operating system (Windows or Linux) with a gcc compiler and toolchain

NVIDIA provides a list of CUDA-enabled GPUs `here <https://developer.nvidia.com/cuda-gpus>`_. Whilst the GeForce series
of NVIDIA GPUs are more affordable and generally *good enough*, from a purely raw computation standpoint, NVIDIA's
workstation and data center GPUs are more suited due differences (and additional) in hardware which is not included in
the gaming focused GeForce GPUs.

Installing the CUDA toolkit
---------------------------

You can either install the NVIDIA CUDA Toolkit either through an installer downloaded from NVIDIA or through your
package manager on Linux. It should be noted that the CUDA Toolkit *does not* come with NVIDIA drivers, so you need to
make sure you install these as well. The NVIDIA CUDA Toolkit is available at `https://developer.nvidia.com/cuda-downloads
<https://developer.nvidia.com/cuda-downloads>`_ and NVIDIA drivers at `https://www.nvidia.co.uk/Download/index.aspx
<https://www.nvidia.co.uk/Download/index.aspx?lang=en-uk>`_.

On Ubuntu 22.04, you can install both the toolkit and NVIDIA's proprietary drivers through :code:`apt`,

.. code:: bash

    sudo apt install nvidia-cuda-toolkit nvidia-driver-535

How to Enable and Run CUDA
==========================

Compilation
-----------

To enable CUDA, we need to compile Python with the CUDA code and libraries. There are several ways to enable CUDA, with
the most simple being to run the configure script in the root directory with the argument :code:`--with-cuda`,

.. code:: bash

    [$PYTHON] $ ./configure --with-cuda

    Configuring Makefile for Python radiative transfer code
    Checking for mpicc...yes
    Checking for gcc...yes
    Checking for nvcc...yes
    Preparing Makefile...Done.

If the NVIDIA CUDA Toolkit is found, you will see output, like above, saying the CUDA compiler :code:`nvcc` was found.
What essentially happens when you run :code:`code` is that a value for the variable :code:`NVCC` is set in the Makefile
in Python's source directory. If you re-run :code:`configure` without :code:`--with-cuda`, then :code:`NVCC` will be
unset and CUDA will not be used. You can easily enable or disable CUDA *"on the fly"* by modifying this variable
without running tone configure script, e.g.,

.. code:: bash

    [$PYTHON/source] $ make clean
    [$PYTHON/source] $ make python NVCC=nvcc

:code:`make clean` has to be run whenever you switch CUDA on or off, due to code conditionally compiling depending on if
CUDA is enabled or not.

Running
-------

Running Python is the exact same, even when running in parallel with MPI. If you are running on a HPC system, then you
will need to request GPU resources. For example, on Iridis at the University of Southampton, a Slurm submission script
may look like,

.. code:: bash

    #!/bin/bash

    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=40
    #SBATCH --time=06:00:00
    #SBATCH --partition=gpu

    module load openmpi/4.1.5/gcc

    mpirun -n $SLURM_NTASKS py model.pf

If you run with CUDA enabled and no GPU resources, Python will exit at the beginning of the program when it cannot
initialise the CUDA environment.

Implementation
==============

This part of the documentation covers the important implementation details of the matrix acceleration. For the most
part, cuSolver can be treated as a library just like GSL where we write wrapper functions around the functionality of
GSL to solve a problem.

- writing wrapper functions
- conditional compilation
- extern "C" linkage

Basics
------

.. code:: cpp
    :caption: A CUDA kernel to transpose a matrix from row to column major

    __global__ void
    tranpose_row_to_column_major(double *row_major, double *column_major, int matrix_size)
    {
        const int idx = blockIdx.x * blockDim.x + threadIdx.x;
        const int idy = blockIdx.y * blockDim.y + threadIdx.y;

        if (idx < matrix_size && idy < matrix_size) {
            column_major[idx * matrix_size + idy] = row_major[idy * matrix_size + idx];
        }
    }

.. code:: c
    :caption: A useful macro for error checking cuSolver returns

    #define CUSOLVER_CHECK(status)                                                                                     \
        do {                                                                                                           \
            cusolverStatus_t err = status;                                                                             \
            if (err != CUSOLVER_STATUS_SUCCESS) {                                                                      \
                Error("cuSolver Error (%d): %s (%s:%d)\n", err, cusolver_get_error_string(err), __FILE__, __LINE__);   \
                return err;                                                                                            \
            }                                                                                                          \
        } while (0)

Structure
---------

To separate CPU and GPU code, it is convention to put CUDA code into :code:`.cu` files. Therefore most of the code
associated with the GPU accelerated matrix code are in the file :code:`$PYTHON/source/matrix_gpu.cu` with the header
file :code:`$PYTHON/source/matrix_gpu.h` which defines function prototypes. The GSL matrix operations are now kept in
the file :code:`$PYTHON/source/matrix_cpu.c` and function prototypes are kept in :code:`$PYTHON/source/templates.h`.

.. code:: c
    :caption: The wrapper function which calls the appropriate matrix solver

    int solve_matrix(double *a_matrix, double *b_vector, int matrix_size, double *x_vector)
    {
        int error;

    #ifdef CUDA_ON
        error = gpu_solve_matrix(...);  /* CUDA implementation */
    #else
        error = cpu_solve_matrix(...);  /* GSL implementation */
    #endif

        return error;
    }

Compiling and Linking
---------------------

- modifying the Makefile
- the CUDA compiler (nvcc)
- additional libraries, lcudart, lcusolver
- linking object code

.. code:: bash
    :caption: A brief overview on how to compile and link C and CUDA code

    # Define compilers
    CC = mpicc
    NVCC = nvcc

    # Define C and CUDA libraries
    C_LIBS = -lgsl -lgslcblas -lm
    CUDA_LIBS = -lcudart -lcusolver

    # Define flags for C and CUDA compilers
    C_FLAGS = -O3 -DCUDA_ON -DMPI_ON -I../includes -L../libs
    CUDA_FLAGS = -O3 -DCUDA_ON

    # Compile CUDA source to object code
    $(NVCC) $(CUDA_FLAGS) $(CUDA_SOURCE) -c -o $(CUDA_OBJECTS)

    # Compile the C code
    $(CC) $(C_FLAGS) $(C_SOURCE) -c -o $(C_OBJECTS)

    # Link the CUDA and C object code together
    $(CC) $(CUDA_OBJECTS) $(C_OBJECTS) -o python $(CUDA_LIBS) $(C_LIBS)


.. code:: bash
    :caption: The recipe in the Python Makefile to build CUDA objects

    # Recipe to create CUDA object code. If NVCC is blank, then nothing happens
    # in thie recipe
    $(CUDA_OBJECTS): $(CUDA_SOURCE)
    ifneq ($(CUDA_FLAG),)
        $(NVCC) $(NVCC_FLAGS) -DCUDA_ON -I$(INCLUDE) -c $< -o $@
    endif

Code Example
------------

Using cuSolver is like using any other CPU library. Below is an example of using the cuSolverDn (dense matrix) library
to perform LU decomposition and to solve a linear system.

.. code:: c
    :caption: A pedagogical (and in-complete) examples of a cuSolver implementation for solving a linear system

    #include <stdlib.h>
    #include <cuSolverDn.h>

    extern "C" int  /* extern "C" has to be used to make it available to the C run time */
    gpu_solve_matrix(double *a_matrix, double *b_vector, int matrix_size, double *x_vector)
    {
        /* First of all, allocate memory on the GPU and copy data from the CPU to the
        GPU. This is part of the code which takes the most time. */
        allocate_memory_for_gpu();
        copy_data_to_gpu();

        /* cuSolver and cuBLAS are both ports of Fortran libraries, which expect arrays to
        be in column-major format and we therefore need to transpose our row-major arrays */
        transpose_row_to_column_major<<<grid_dim, block_dim>>>(d_matrix_row, d_matrix_col, matrix_size);

        /* Perform LU decomposition. Variables prefixed with d_ are kept in GPU memory where we
        allocated space for them in `allocate_memory_for_gpu` */
        CUSOLVER_CHECK(cusolverDnDgetrf(
            CUSOLVER_HANDLE, matrix_size, matrix_size, d_matrix_col, matrix_size, d_workspace, d_pivot, d_info
        ));

        /* Solve the linear system A x = b. The final solution is returned in the variable d_v_vector */
        CUSOLVER_CHECK(cusolverDnDgetrs(
            CUSOLVER_HANDLE, CUSOLVER_OP_N, matrix_size, matrix_size, d_matrix_col, matrix_size, d_pivot,
            d_b_vector, matrix_size, d_info
        ));

        /* We now have to copy d_b_vector back to the CPU, so we can use that value in
        the rest of Python */
        copy_data_to_cpu();

        return EXIT_SUCCESS;
    }


.. code:: c
    :caption: Fundamentally, the API hasn't changed at all

    #include <stdlib.h>
    #include "python.h"

    double *populations = malloc(nions * sizeof(*populations));
    double *ion_density = malloc(nions * sizeof(*populations));
    double *rate_matrix = malloc(nions * nions * sizeof(*populations));

    /* The wrapper function is named the same as the original GSL implementation
       and accepts the same arguments */
    int error = solve_matrix(
        rate_matrix, ion_density, nions, populations, xplasma->nplasma
    );

    /* One user difference is that error handling is more robust now, and there
       is a function to convert error codes into error messages */
    if (error != EXIT_SUCCESS) {
        Error(
            "Error whilst solving for ion populations: %d (%d)\n",
            get_matrix_error_string(error), error
        );
    }

