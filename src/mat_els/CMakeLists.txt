set(mat_els_files
    "hamiltonian.f90"
    "mat_els.f90"
    "orbitals.f90"
    "potentials.f90"
)

add_library(mat_els ${mat_els_files})

target_link_libraries(mat_els tools)

install(TARGETS mat_els DESTINATION lib)