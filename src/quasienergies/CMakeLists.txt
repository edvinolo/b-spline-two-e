set(quasienergies_files
    "quasi_calcs.f90"
    "quasi_circ.f90"
    "quasi_floquet.f90"
)

add_library(quasienergies ${quasienergies_files})

target_link_libraries(quasienergies PRIVATE fortran_stdlib::fortran_stdlib)
target_link_libraries(quasienergies PRIVATE tools)

install(TARGETS quasienergies DESTINATION lib)