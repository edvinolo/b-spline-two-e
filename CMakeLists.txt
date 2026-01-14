# Figure out which minimum version should be used
cmake_minimum_required(VERSION 3.10)

# If you want to print variables
include(CMakePrintHelpers)

# Project information
project(bs2e)
set(PROJECT_VERSION "0.1.0")
set(LANGUAGES Fortran)
enable_language(Fortran)

# Binary output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# Add include directorires
include_directories("${CMAKE_BINARY_DIR}/src/tools")
include_directories("${CMAKE_BINARY_DIR}/src/mat_els")
include_directories("${CMAKE_BINARY_DIR}/src/quasienergies")
include_directories("${CMAKE_BINARY_DIR}/src/diagonalization")
include_directories("${CMAKE_BINARY_DIR}/src/function_eval")
include_directories("${CMAKE_BINARY_DIR}/src/time_dep")
include_directories("${CMAKE_BINARY_DIR}/src/RKH_frame")

# Make sure color ouput is used by the build system for more readable error messages and warnings
set(CMAKE_COLOR_DIAGNOSTICS ON)

# Set allowed config types
set(ALLOWED_CONFIGS Debug Sanitize Release RelWithDebInfo)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build" FORCE
    )
endif ()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${ALLOWED_CONFIGS})
message(STATUS "Selected CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

# Set compiler flags for compilers and build types
add_library(project_options INTERFACE)

# Set compiler flags that should always be on
target_compile_options(project_options INTERFACE
    # GNU Fortran
    $<$<Fortran_COMPILER_ID:GNU>:
        -Wall
        -cpp
        -fdefault-integer-8>

    # Intel ifx
    $<$<Fortran_COMPILER_ID:IntelLLVM>:
        -warn all
        -fpp
        -i8>
)

# Set flags for compilation choices
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(project_options INTERFACE
        $<$<Fortran_COMPILER_ID:GNU>: -O0 -fcheck=all -g -fbacktrace>
        $<$<Fortran_COMPILER_ID:IntelLLVM>: -O0 -g -check all -traceback>
    )
elseif(CMAKE_BUILD_TYPE STREQUAL "Sanitize")
    target_compile_options(project_options INTERFACE
        $<$<Fortran_COMPILER_ID:GNU>: -O0 -g -fsanitize=address -fno-omit-frame-pointer -fbacktrace>
        $<$<Fortran_COMPILER_ID:IntelLLVM>: -O0 -g -fsanitize=address -fno-omit-frame-pointer -traceback>
    )
    target_link_options(project_options INTERFACE
        -fsanitize=address
    )
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(project_options INTERFACE
        $<$<Fortran_COMPILER_ID:GNU>: -O3 -march=native -mtune=native>
        $<$<Fortran_COMPILER_ID:IntelLLVM>: -O3 -xHost>
    )
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    target_compile_options(project_options INTERFACE
        $<$<Fortran_COMPILER_ID:GNU>: -O3 -march=native -mtune=native -g -fbacktrace>
        $<$<Fortran_COMPILER_ID:IntelLLVM>: -O3 -xHost -g -traceback>
    )
else()
    message(FATAL_ERROR
        "Unknown CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}\n"
        "CMAKE_BUILD_TYPE should be one of: ${ALLOWED_CONFIGS}"
    )
endif()

# Make sure OpenMP flags are set
find_package(OpenMP REQUIRED)
target_link_libraries(project_options INTERFACE OpenMP::OpenMP_Fortran)

option(USE_DEBUG_PARDISO "Enable debugging mode of PARDISO" OFF)
if (USE_DEBUG_PARDISO)
    add_compile_definitions(DEBUG_PARDISO)
endif()

# Tell MKL to use 32-bit integer interface
# set(MKL_INTERFACE lp64)

# Tell MKL to use GNU OpenMP
# set(MKL_THREADING gnu_thread)

# Find required packages
find_package(GSL REQUIRED)
find_package(MKL CONFIG REQUIRED)
find_package(fortran_stdlib REQUIRED)
find_package(arpackng REQUIRED)

# Find BLAS and LAPACK, use MKL versions
set(BLA_VENDOR Intel10_64ilp)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Check if FEAST is available
option(USE_FEAST "Enable FEAST" OFF)
if (USE_FEAST)
    find_library(feast NAMES feast PATHS "$ENV{FEASTROOT}/lib/x64/")
    if (feast)
        add_compile_definitions(WITH_FEAST)
    endif()
endif()

# Add subdirectories
add_subdirectory(src/tools)
add_subdirectory(src/mat_els)
add_subdirectory(src/quasienergies)
add_subdirectory(src/diagonalization)
add_subdirectory(src/function_eval)
add_subdirectory(src/time_dep)
add_subdirectory(src/RKH_frame)

# Add executables
add_executable(basis_setup src/apps/main_basis_setup.f90)
add_executable(quasi src/apps/main_quasienergies.f90)
add_executable(diag src/apps/main_diagonalization.f90)
add_executable(feval src/apps/main_eval_rad.f90)
add_executable(time_prop src/apps/main_time_prop.f90)
add_executable(RKH src/apps/main_RKH.f90)

# Link libraries to executable
target_link_libraries(basis_setup PRIVATE mat_els tools project_options)
target_link_libraries(quasi PRIVATE quasienergies tools project_options)
target_link_libraries(diag PRIVATE diagonalization tools project_options)
target_link_libraries(feval PRIVATE function_eval tools project_options)
target_link_libraries(time_prop PRIVATE time_dep quasienergies tools project_options)
target_link_libraries(RKH PRIVATE RKH_frame tools project_options)

# Install executable in binary direcetory
set(executables basis_setup quasi diag feval time_prop RKH)
install(TARGETS ${executables} DESTINATION bin)