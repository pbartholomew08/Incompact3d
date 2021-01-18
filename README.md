# xcompact3d_developers
Developers repository for xcompact3d
=======
<a name="logo"/>
<div align="center">
<a href="https://www.incompact3d.com" target="_blank">
<img src="https://www.incompact3d.com/uploads/5/8/7/2/58724623/incompact3d.png" alt="Logo" width="300" height="100"></img>
</a>
</div>

[![Build Status](https://travis-ci.org/xcompact3d/Incompact3d.svg?branch=master)](https://travis-ci.org/xcompact3d/Incompact3d)
![Build Status](https://github.com/xcompact3d/Incompact3d/workflows/Build/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/xcompact3d/badge/?version=latest)](https://xcompact3d.readthedocs.io/en/latest/?badge=latest)

## The Incompact3d code

Incompact3d is a Fortran-MPI based, finite difference high-performance code for solving the incompressible Navier-Stokes equation and as many as you need scalar transport equations.

This repository contains a major upgrade in the Incompact3d code. The new version is faster, more flexible and user friendly.

The main homepage for Incompact3d can be found at [incompact3d.com](https://www.incompact3d.com/). You can find a list of the work related to the code.

This is the GitHub repository of Incompact3d source code, including instructions for running and compiling Incompact3d, below.

## Resources

- **Homepage:** <https://www.incompact3d.com/>
- **Binaries:** <https://www.incompact3d.com/download.html>
- **Documentation:** <https://www.incompact3d.com/docs>
- **Discussion section (FORUM):** <https://github.com/xcompact3d/Incompact3d/issues>
- **Git clone URL-SSH:** <git@github.com:xcompact3d/Incompact3d.git>
- **Git clone URL-HTTPS:** <https://github.com/xcompact3d/Incompact3d.git>

New users and developers are welcome to join.

### External Resources

- [**Twitter**](https://twitter.com/incompact3d)

### Main improvements/changes ###

- [**Wiki**](https://github.com/xcompact3d/Incompact3d/wiki/New-Features)

### Benchmark test cases for code comparison ###

We estabilished a solid and easy way to run a range of benchmark test cases to verify the code. Incompact3d works now on a flow configuration specific file. You must choose a case and set it on the 'Makefile' and recompile. The following cases are set to match the parameters for cases of reference articles obtained with different codes.

|Code| Flow configuration             | BC File         | Reference | Dataset |
|:----------------:|:----------------:|:----------------:|:----------------:|:----------------:|
|1| Taylor-Green Vortices        | TGV              |[Beck et al. (2014)](https://link.springer.com/article/10.1007/s00162-011-0253-7)
|2| Periodic Channel            | Channel-flow     |[Moser, Kim & Mansour (1999)](https://www.researchgate.net/publication/243777258_Direct_Numerical_Simulation_of_Turbulent_Channel_Flow_up_to_Re590)|[Dataset](http://turbulence.ices.utexas.edu/data/MKM/chan180/profiles/)
|3| Flow over a Cylinder         | Cylinder         |[Mittal and Balachandar (1995)](https://www.researchgate.net/publication/252073966_Effect_of_three-dimensionality_on_the_lift_and_drag_of_nominal_two-dimensional_cylinders)
|4| Periodic Hill                | Periodic-Hill    |[Breuer et al. (2009)](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&ved=2ahUKEwiR9fjTx97cAhWCC5AKHbItAwEQFjAAegQIAhAC&url=https%3A%2F%2Fcfd.spbstu.ru%2Fagarbaruk%2Fdoc%2F2009_Breuer-et-al._Flow-over-periodic-hills.-Numerical-and-experimental-study-in-a-wide-range-of-Reynolds-numbers.pdf&usg=AOvVaw0pL7mF2TaPEbY0a8dOK-fN)
|5| Gravity Current              | Lock-exchange    |[Necker et al. (2002)](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.483.5774&rep=rep1&type=pdf)
|6| Boundary Layer               | TBL              |


## New compiling FLAGS
If the flags are not specified in the Makefile, the compile ignore the sections related to each flag. For example, if you do not need IBM in your simulation, do not compile the code with -DIBM,

   -DDOUBLE_PREC - use double-precision
   -DSAVE_SINGLE - save 3D data in single-precision
   -DDEBG        - debuggin incompact3d.f90
   -DIBM         - enable IBM calls
   -DPOST        - enable statistics processing
   -DVISU        - enable visu.f90
   -DVISUEXTRA   - enable extra options visu.f90
   -DELES        - enable explicit LES modelling
   -DSTRETCHING  - enable mesh stretching in y direction

**Note:** In order to compile the code with the apropiate flags you must enter the -D$FLAG, i.e., -DDOUBLE_PREC

You can do an export of the variable in order to the Makefile knows what BC you want to compiler
       export FLOW_TYPE=Channel-flow
or uncomment line 19 from the Makefile

## Source Download and Compilation

First, make sure you have all the [required dependencies](#required-build-tools-and-external-libraries) installed.
Then, acquire the source code by cloning the git repository:

    git clone git@github.com:xcompact3d/Incompact3d.git

(If you are behind a firewall, you may need to use the `https` protocol instead of the `git` protocol:

    git config --global url."https://".insteadOf git@

Be sure to also configure your system to use the appropriate proxy settings, e.g. by setting the `https_proxy` and `http_proxy` variables.)

By default you will be building the latest unstable version of Incompact3d. However, most users should use the most recent stable version of Incompact3d, which is currently the `3.0` series of releases. You can get this version by changing to the Incompact3d directory and running

    git checkout v3.0

Now run `make` to build the `Incompact3d` executable. To perform a parallel build, use `make -j N` and supply the maximum number of concurrent processes. (See [Platform Specific Build Notes] for details.)
This takes a while, but only has to be done once. If the defaults in the build do not work for you, and you need to set specific make parameters, you can save them in `Make.user`. The build will automatically check for the existence of `Makefile` and use it if it exists.
Building Incompact3d requires very little of disk space and virtual memory.

**Note:** The compiling process

Once it is built, you can run the `Incompact3d` executable using its full path in the directory created above (the `Incompact3d` directory).

Now you should be able to run Incompact3d like this:

    mpirun -n 4 ./incompact3d

If everything works correctly, you will see a Incompact3d banner and an interactive prompt into which you can enter expressions for evaluation. (Errors related to libraries might be caused by old, incompatible libraries sitting around in your PATH. In this case, try moving the `Incompact3d` directory earlier in the PATH).

### Updating an existing source tree

If you have previously downloaded `incompact3d` using `git clone`, you can update the
existing source tree using `git pull` rather than starting anew:

    cd incompact3d
    git pull && make

Assuming that you had made no changes to the source tree that will conflict
with upstream updates, these commands will trigger a build to update to the
latest version.

#### General troubleshooting

1. Over time, the base library may accumulate enough changes such that the
   bootstrapping process in building the system image will fail. If this
   happens, the build may fail with an error like

   ```sh
    *** This error is usually fixed by running 'make clean'. If the error persists, try 'make cleanall' ***
   ```
   As described, running `make clean && make` is usually sufficient. Occasionally, the stronger cleanup done by `make cleanall` is needed.


## Platform-Specific Notes

### Linux

#### General

* GCC version 4.7 or later is required to compile the code.

We recommended that you remove the limits of the environment (e.g. in `.bash_profile`)

    ulimit -c unlimited
    ulimit -s unlimited

You must exit and re-login from your terminal for the change to take effect

    ulimit -a

#### Architecture Customization

Incompact3d can be compiled for a non-generic architecture by configuring the `ARCH` Makefile variable. See the appropriate section of `Makefile` for additional customization options, such as `MARCH` and `CPU_TARGET`. You can also set `march=native` for a maximum-performance build customized for the current machine CPU.

#### Ubuntu

In order to compile and execute Incompact3d in the latest Ubuntu version please install the following packages:

    sudo apt install gfortran libopenmpi-dev

#### Fedora/RHEL/CentOS

On RHEL/CentOS 6 systems:

    sudo dnf install gcc-gfortran

### OS X

You need to have the current Xcode command line utilities installed: run `xcode-select --install` in the terminal.
You will need to rerun this terminal command after each OS X update, otherwise you may run into errors involving missing libraries or headers. You will also need a 64-bit gfortran to compile the code. The gfortran-4.7 (and newer) compilers in Brew, Fink, and MacPorts work for building Incompact3d. On current systems, we recommend that you install the command line tools as described above. Older systems do not have a separate command line tools package from Apple, and will require a full Xcode install. On these, you will need at least Xcode 4.3.3. In Xcode prior to v5.0, you can alternatively go to Preferences -> Downloads and select the Command Line Utilities.


### Intel compilers and Math Kernel Library (MKL)

To build Incompact3d using the Intel compilers, and link the [MKL] FFTW libraries, first make sure you have a recent version
of the Intel Parallel Studio XE Cluster Edition (version 15 or later). Qualified students, classroom educators and open-source contributors can obtain a license and download the latest version on https://software.intel.com/en-us/parallel-studio-xe/choose-download/student-linux-fortran

After installation, for a 64-bit architecture, the environment should be set up as follows:

    # bash
    source /path/to/intel/bin/compilervars.sh intel64

You also need to activate the shared memory copy [LMT mechanism](https://software.intel.com/en-us/mpi-developer-reference-linux-shared-memory-control)

    export I_MPI_SHM_LMT=shm

The code is compiled with the mpiifort command and the appropriate flags are set in the `Makefile` file.
If you eish to compile the code with Intel Compiler use:

      LCL = local# local,lad,sdu,archer
      IVER = 17# 15,16,17,18
      CMP = intel
      FFT = mkl

 If you have a non-standard installation the supercomputer you're using, you have to edit and create a custom flag.

## GitHub Configuration

To add you SSH key to your GitHub account please follow the steps https://help.github.com/articles/adding-a-new-ssh-key-to-your-github-account/ or just copy the content of the id_rsa.pub file to your clipboard, go to Personal settings and add a new SSH key:

    cat ~/.ssh/id_rsa.pub
