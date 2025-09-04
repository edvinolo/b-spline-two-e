program eval_rad
    use kind_tools
    use input_tools
    use dir_tools
    use grid_tools
    use bspline_tools
    use block_tools
    use orbital_tools
    use function_eval
    implicit none

    character(len=:), allocatable :: input_file
    character(len=:), allocatable :: eval_res_dir

    type(basis) :: bas
    type(b_spline) :: splines

    call set_eval_defaults()
    input_file = get_input_file()
    call get_eval_input(input_file)
    call get_basis_input(eval_root_dir//"/basis_input.dat")

    call bas%load(eval_root_dir)
    call splines%load(eval_root_dir)

    eval_res_dir = eval_root_dir//'function_eval/'
    if (.not.is_dir(eval_res_dir)) call mkdir(eval_res_dir)

    if (two_el) then
        call eval_from_states_2p(eval_res_dir,eval_root_dir,bas,splines)
    else
        call eval_from_states_1p(eval_res_dir,eval_root_dir,bas,splines)
    end if

end program eval_rad