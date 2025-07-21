set(quasienergies_files
    "essential_states.f90"
    "quasi_calcs.f90"
    "quasi_circ.f90"
    "quasi_floquet.f90"
    "static.f90"
)

add_library(quasienergies ${quasienergies_files})

target_link_libraries(quasienergies PRIVATE fortran_stdlib::fortran_stdlib)
target_link_libraries(quasienergies PRIVATE tools)
target_link_libraries(quasienergies PRIVATE diagonalization)

install(TARGETS quasienergies DESTINATION lib)