# b-spline-two-e
A program that computes matrix elements for one- and two-electron atoms in the non-relativistic and dipole approximations (including the singlet approximation).
The radial part of the wavefunctions are represented using B-splines. The program utilizes OpenMP for parallel processing.

## Requirements
The currently supported compilers are
- ```gfortran``` 10+, ```ifx```
  
The program needs the following external libraries
- The GNU Scientific Library (GSL) https://www.gnu.org/software/gsl/
- Intel MKL (can be installed as part of their Fortran essentials package)
- The Fortran standard library (with support for ilp64) https://github.com/fortran-lang/stdlib
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
- To build with Ninja use
  ```
  cmake .. -G Ninja
  ninja
  ```

A shell script for building the program and its prerequisites on AWS-Linux 2023 can be found in ```utils/install_bs2e_aws_linux_2023.sh```

## Configure the build
These are options that can be set when running CMake in order to configure your build.

### Setting the Fortran compiler
A specific compiler and version can be specified using ```-DCMAKE_Fortran_COMPILER=compiler```
Currently, only ```gfortran``` and ```ifx``` are supported.

### Setting the build type
There are currently four build types available, and they can be set via the ```-DCMAKE_BUILD_TYPE=Type``` flag.
The current options are:
- ```Debug``` Optimization level: -O0. Runtime checks, debug symbols, and backtrace enabled.
- ```Sanitize``` Optimization level: -O0. Address sanitizer, debug symbols, and backtrace enabled. 
- ```Release``` Optimization level: -O3. For ```gfortran``` also enables ```-march=native -mtune=native```, and for ```ifx``` enables ```-xHost```.
- ```RelWithDebInfo``` Same as ```Release```, but with debug symbols and backtrace enabled.

```Release``` is the default build type, and is used if ```-DCMAKE_BUILD_TYPE``` is not specified.

### Debug PARDISO
To enable more checks for PARDISO and printing of information add:
```
-DUSE_DEBUG_PARDISO=ON
```

### Build with FEAST
If you want to enable FEAST set
```
-DUSE_FEAST=ON
```

## Currently available programs
  - ```basis_setup``` For setting up the basis and computing matrix elements.
  - ```diag``` For finding targeted eigenstates of the Hamiltonian.
  - ```quasi``` For finding complex quasienergies (Floquet, rotating-frame, or static-field).
  - ```feval``` Evaluate one-particle wavefunction or the electron density from a given state vector.
  - ```time_prop``` Perform time propagation to solve the TDSE.
  - ```RKH``` Find rotating frame Siegert states and adiabatic potentials
  
### How to run 
There are example input files located in the ```examples``` directory.
To run e.g. ```basis_setup```, run 
```
/path/to/build/bin/basis_setup /path/to/basis_setup.in
```
Note that the output directory needs to be specfied in the input file, and it should exist before running the program (otherwise you will get an error).
The output of the program will be stored in a numbered folder created inside the output directory specified in the input file 
(except in the case of ```feval``` where a ```function_eval``` directory is created in the directory where the vectors are stored).

The number of threads used in OpenMP sections of the program can be controlled by setting the ```OMP_NUM_THREADS``` environment variable.
