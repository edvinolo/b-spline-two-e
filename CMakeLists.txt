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

# Set compiler flags that should always be on
add_compile_options(-Wall -fopenmp)
add_link_options(-fopenmp)

# Set flags for compilation choices
if (CMAKE_BUILD_TYPE STREQUAL "debug")
    add_compile_options(-O0 -fcheck=all -g)
else()
    add_compile_options(-O3 -march=native)
endif()

# Force MKL to use 32-bit interface
set(MKL_INTERFACE lp64)

# Find required packages
find_package(GSL REQUIRED)
find_package(MKL CONFIG REQUIRED)
find_package(fortran_stdlib REQUIRED)
find_library(feast NAMES feast PATHS "$ENV{FEASTROOT}/lib/x64/" REQUIRED)

# Add subdirectories
add_subdirectory(src/mat_els)
add_subdirectory(src/tools)

# Add executables
add_executable(basis_setup src/apps/main_basis_setup.f90)

# Link libraries to executable
target_link_libraries(basis_setup mat_els tools)

# Install executable in binary direcetory
set(executables basis_setup)
install(TARGETS ${executables} DESTINATION bin)
