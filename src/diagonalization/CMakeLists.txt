set(diagonalization_files
    "diagonalization.f90"
)

add_library(diagonalization ${diagonalization_files})

target_link_libraries(diagonalization PRIVATE tools)
target_link_libraries(diagonalization PRIVATE project_options)

install(TARGETS diagonalization DESTINATION lib)