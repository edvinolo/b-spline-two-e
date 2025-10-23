program eval_rad
    use kind_tools
    use input_tools
    use dir_tools
    use grid_tools
    use bspline_tools
    use block_tools
    use orbital_tools
    use function_eval
    use omp_lib, only: omp_get_wtime
    implicit none

    character(len=:), allocatable :: input_file
    character(len=:), allocatable :: eval_res_dir

    type(basis) :: bas
    type(b_spline) :: splines

    real(dp) :: t_prog_start,t_prog_end

    t_prog_start = omp_get_wtime()

    call set_eval_defaults()
    input_file = get_input_file()
    call get_eval_input(input_file)
    call get_basis_input(eval_root_dir//"/basis_input.dat")

    call bas%load(eval_root_dir)
    call splines%load(eval_root_dir)

    eval_res_dir = eval_root_dir//'function_eval/'
    if (.not.is_dir(eval_res_dir)) call mkdir(eval_res_dir)
    call write_eval_input(eval_res_dir)

    if (two_el) then
        call eval_from_states_2p(eval_res_dir,eval_root_dir,bas,splines)
    else
        call eval_from_states_1p(eval_res_dir,eval_root_dir,bas,splines)
    end if

    t_prog_end = omp_get_wtime()
    write(6,*) "Total time for function eval (s): ", t_prog_end - t_prog_start
    write(6,*) "Total time for function eval (min): ", (t_prog_end - t_prog_start)/60.d0
    write(6,*) "Total time for function eval (h): ", (t_prog_end - t_prog_start)/(60.d0*60.d0)
end program eval_rad