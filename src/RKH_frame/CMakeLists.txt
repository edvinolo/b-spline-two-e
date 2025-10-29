set(RKH_frame_files
    "adiabatic_potentials.f90"
    "RKH_hamiltonian.f90"
)

add_library(RKH_frame ${RKH_frame_files})

target_link_libraries(RKH_frame PRIVATE fortran_stdlib::fortran_stdlib)
target_link_libraries(RKH_frame PRIVATE tools)
target_link_libraries(RKH_frame PRIVATE function_eval)

install(TARGETS RKH_frame DESTINATION lib)