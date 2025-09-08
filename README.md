# b-spline-two-e
A program that computes matrix elements for one- and two-electron atoms in the non-relativistic dipole approximation (including the singlet approximation)
The radial part of the wavefunctions are represented using B-splines.

## Requirements
The currently supported compilers are
- gfortran 10+
  
The program needs the following external libraries
- The GNU Scientific Library (GSL) https://www.gnu.org/software/gsl/
- Intel MKL (can be installed as part of their Fortran essentials package)
- The Fortran standard (with support for ilp64) library https://github.com/fortran-lang/stdlib
- ARPACK-NG (Note, ARPACK-NG should be built with ilp64 support, see their README on instructions how) https://github.com/opencollab/arpack-ng
- FEAST (optional, currently not used for any of the main programs) https://www.feast-solver.org/index.htm

## How to build the executables
Once you have cloned the repository, do the following:
- Enter the top-level directory of the project.
- Create the build directory and enter it
  ```
  mkdir build
  cd build  
  ```
- Make sure that the MKL environment variables are set
  ```
  source /path/to/MKL/setvars.sh
  ```
- Run CMake (with default build type)
  ```
  cmake ..
  ```
- Compile the executables by running make
  ```
  make
  ```

A shell script for building the program and its prerequisites on AWS-Linux 2023 can be found in ```utils/install_bs2e_aws_linux_2023.sh```

## Turn on debug build
In the build directory run cmake with the debug-flag set
```
cmake -DCMAKE_BUILD_TYPE=debug ..
```
This adds runtime checks, debug symbols, and sets the optimization level to -O0.
If you additionaly want to enable the address sanitizer add
```
-DCMAKE_USE_ASAN=1
```
when you run CMake.

## Currently available programs
  - ```basis_setup``` For setting up the basis and computing matrix elements.
  - ```diag``` For finding targeted eigenstates of the Hamiltonian.
  - ```quasi``` For finding complex quasienergies (Floquet, rotating-frame, or static-field).
  - ```feval``` Evaluate one-particle wavefunction or the electron density from a given state vector.

### How to run 
There are example input files located in the ```examples``` directory
To run e.g. ```basis_setup```, run 
```
/path/to/build/bin/basis_setup /path/to/basis_setup.in
```
Note that the output directory needs to be specfied in the input file, and it should exist before running the program (otherwise you will get an error).
The output of the program will be stored in a numbered folder created inside the output directory specified in the input file 
(except in the case of ```feval``` where a ```function_eval``` directory is created in the directory where the vectors are stored).
