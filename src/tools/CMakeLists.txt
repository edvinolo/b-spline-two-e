set(tool_files
    "block_tools.f90"
    "bspline_tools.f90"
    "CAP_tools.f90"
    "constants_tools.f90"
    "dir_tools.f90"
    "eig_tools.f90"
    "GMRES_tools.f90"
    "grid_tools.f90"
    "ILU0_tools.f90"
    "input_tools.f90"
    "Jacobi_tools.f90"
    "kind_tools.f90"
    "orbital_tools.f90"
    "PARDISO_tools.f90"
    "precond_tools.f90"
    "quad_tools.f90"
    "sparse_array_tools.f90"
    "wigner_tools.f90"
)

add_library(tools ${tool_files})

target_link_libraries(tools PRIVATE fortran_stdlib::fortran_stdlib)
target_link_libraries(tools PRIVATE MKL::MKL)
target_link_libraries(tools PRIVATE GSL::gsl)
target_link_libraries(tools PRIVATE ARPACK::ARPACK)
target_link_libraries(tools PRIVATE project_options)

if (feast)
    target_link_libraries(tools PRIVATE ${feast})
endif()

install(TARGETS tools DESTINATION lib)