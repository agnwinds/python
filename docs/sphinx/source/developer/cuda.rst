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

Writing CUDA
------------

- .cu files (C++)
- writing wrapper functions
- conditional compilation
- extern "C" linkage

Compiling and Linking CUDA
--------------------------

- modifying the Makefile
- the CUDA compiler (nvcc)
- additional libraries, lcudart, lcusolver
- linking object code
