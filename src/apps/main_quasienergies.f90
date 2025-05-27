program main_quasienergies
    use kind_tools
    use constants_tools
    use input_tools
    use dir_tools
    use orbital_tools, only: basis
    use sparse_array_tools, only: CSR_matrix
    use block_tools, only: block_CS, block_diag_CS
    use quasi_floquet
    use quasi_circ
    use quasi_calcs
    implicit none

    character(len=:), allocatable :: input_file
    character(len=:), allocatable :: quasi_res_dir

    type(basis) :: bas
    type(block_diag_CS) :: H_0_block,S_block
    type(block_CS) :: dip(-1:1)
    type(block_CS) :: D

    ! Read the input, and the input file used for the basis_setup program
    call set_quasi_defaults()
    input_file = get_input_file()
    call get_quasi_input(input_file)
    call get_basis_input(basis_dir//"/basis_input.dat")

    ! Load basis and matrix elements
    call bas%load(basis_dir)
    call H_0_block%load(basis_dir//"H_diag.dat")
    call S_block%load(basis_dir//"S_diag.dat")
    call dip(-1)%load(basis_dir//"D_-1.dat")
    call dip(0)%load(basis_dir//"D_0.dat")
    call dip(1)%load(basis_dir//"D_1.dat")

    ! Create directory for results and copy inputs
    call make_res_dir(quasi_root_dir,quasi_res_dir)
    call write_basis_input(quasi_res_dir)
    call write_quasi_input(quasi_res_dir)

    ! Setup the structure of the interaction matrix D (to be multiplied by field strength to get true interaction V)
    if (z_pol) then
        !call setup_dip_floquet(dip(0),D,gauge,n_blocks,bas)
    else
        call setup_dip_circ(dip,D,gauge)
    end if

    if (reduce_basis) then
        call red_basis(bas,H_0_block,S_block,D)
    else
        precond_blocks = precond_blocks_in
    end if

    if (block_precond) then
        call reorder_blocks(bas,H_0_block,S_block,D)
    end if

    ! Store information about the actual basis that was used
    call bas%store(quasi_res_dir)

    if (calc_type == 'omega') then
        call omega_scan(H_0_block,S_block,D,bas,gs_energy,quasi_res_dir)
    else if (calc_type == 'intensity') then
        call intensity_scan(H_0_block,S_block,D,bas,gs_energy,quasi_res_dir)
    else if (calc_type == 'omega_follow') then
        call omega_follow(H_0_block,S_block,D,bas,gs_energy,quasi_res_dir)
    else if (calc_type == 'intensity_follow') then
        call intensity_follow(H_0_block,S_block,D,bas,gs_energy,quasi_res_dir)
    end if

    write(6,*) ""
    write(6,*) "Output was written to ", quasi_res_dir
end program main_quasienergies