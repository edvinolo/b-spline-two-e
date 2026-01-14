set(function_eval_files
    "function_eval.f90"
)

add_library(function_eval ${function_eval_files})

target_link_libraries(function_eval PRIVATE fortran_stdlib::fortran_stdlib)
target_link_libraries(function_eval PRIVATE tools)
target_link_libraries(function_eval PRIVATE mat_els)
target_link_libraries(function_eval PRIVATE project_options)

install(TARGETS function_eval DESTINATION lib)