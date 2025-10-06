set(time_dep_files
    "fields.f90"
)

add_library(time_dep ${time_dep_files})

# target_link_libraries(time_dep PRIVATE fortran_stdlib::fortran_stdlib)
target_link_libraries(time_dep PRIVATE tools)
target_link_libraries(time_dep PRIVATE quasienergies)

install(TARGETS time_dep DESTINATION lib)