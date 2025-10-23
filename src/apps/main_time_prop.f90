program main_time_prop
    use kind_tools
    use constants_tools
    use input_tools
    use dir_tools
    use bspline_tools, only: b_spline
    use orbital_tools, only: basis
    use sparse_array_tools, only: CSR_matrix
    use block_tools, only: block_CS, block_diag_CS
    use essential_states, only: setup_essential_states, find_essential_states
    use quasi_calcs, only: red_basis, reorder_blocks
    use fields, only: pulse
    use time_propagation, only: init_pulses, init_time_grid, init_CN_propagator, init_psi, propagate, write_results, ess_states
    use omp_lib, only: omp_get_wtime
    implicit none

    character(len=:), allocatable :: input_file
    character(len=:), allocatable :: time_prop_res_dir

    type(basis) :: bas
    type(b_spline) :: splines
    type(block_diag_CS) :: H_0_block,S_block
    type(block_CS) :: dip(-1:1)

    real(dp) :: t_prog_start,t_prog_end

    t_prog_start = omp_get_wtime()

    ! Read the input, and the input file used for the basis_setup program
    call set_time_prop_defaults()
    input_file = get_input_file()
    call get_time_prop_input(input_file)
    call get_basis_input(basis_dir//"/basis_input.dat")

    ! Load basis and matrix elements
    call bas%load(basis_dir)
    call splines%load(basis_dir)
    call H_0_block%load(basis_dir//"H_diag.dat")
    call S_block%load(basis_dir//"S_diag.dat")
    call dip(-1)%load(basis_dir//"D_-1.dat")
    call dip(0)%load(basis_dir//"D_0.dat")
    call dip(1)%load(basis_dir//"D_1.dat")

    ! Create directory for results and copy inputs
    call make_res_dir(time_prop_root_dir,time_prop_res_dir)
    call write_basis_input(time_prop_res_dir)
    call write_time_prop_input(time_prop_res_dir)

    if (use_essential_states) then
        call setup_essential_states(ess_states,n_ess,target_blocks,bas)
        call find_essential_states(ess_states,n_ess,targets,H_0_block,S_block)
    end if

    if (reduce_basis) then
        if (use_essential_states) then
            call red_basis(bas,H_0_block,S_block,dip,ess_states)
        else
            call red_basis(bas,H_0_block,S_block,dip)
        end if
    else
        precond_blocks = precond_blocks_in
    end if

    if ((block_precond).and.(block_precond_type=='PQ')) then
        if (use_essential_states) then
            call reorder_blocks(bas,H_0_block,S_block,dip,ess_states)
        else
            call reorder_blocks(bas,H_0_block,S_block,dip)
        end if
    end if

    ! Store information about the actual basis that was used
    call bas%store(time_prop_res_dir)
    call splines%store(time_prop_res_dir)

    ! Prepare the time propagation
    call init_pulses()
    call init_time_grid()
    call init_CN_propagator(H_0_block,S_block,dip)
    call init_psi(bas)

    ! Propagate TDSE
    call propagate()

    ! Store results
    call write_results(time_prop_res_dir)

    t_prog_end = omp_get_wtime()
    write(6,*) "Total time for time_prop (s): ", t_prog_end - t_prog_start
    write(6,*) "Total time for time_prop (min): ", (t_prog_end - t_prog_start)/60.d0
    write(6,*) "Total time for time_prop (h): ", (t_prog_end - t_prog_start)/(60.d0*60.d0)

    write(6,*) ""
    write(6,*) "Output was written to ", time_prop_res_dir
end program main_time_prop