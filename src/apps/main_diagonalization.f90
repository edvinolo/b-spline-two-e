program main_diagonalization
    use kind_tools
    use constants_tools
    use input_tools
    use dir_tools
    use bspline_tools, only: b_spline
    use orbital_tools, only: basis
    use sparse_array_tools, only: CSR_matrix
    use block_tools, only: block_CS, block_diag_CS
    use diagonalization
    implicit none

    character(len=:), allocatable :: input_file
    character(len=:), allocatable :: diag_res_dir

    type(basis) :: bas
    type(b_spline) :: splines
    type(block_diag_CS) :: H_block,S_block

    ! Read the input, and the input file used for the basis_setup program
    call set_diag_defaults()
    input_file = get_input_file()
    call get_diag_input(input_file)
    call get_basis_input(basis_dir//"/basis_input.dat")

    ! Load basis and matrix elements
    call bas%load(basis_dir)
    call splines%load(basis_dir)
    call H_block%load(basis_dir//"H_diag.dat")
    call S_block%load(basis_dir//"S_diag.dat")

    ! Create directory for results and copy inputs
    call make_res_dir(diag_root_dir,diag_res_dir)
    call write_basis_input(diag_res_dir)
    call write_diag_input(diag_res_dir)

    ! Store information about the actual basis that was used
    call bas%store(diag_res_dir)
    call splines%store(diag_res_dir)

    ! Find the targeted energies for the requested symmetries
    call do_diag(H_block,S_block,bas,diag_res_dir)

    write(6,*) ""
    write(6,*) "Output was written to ", diag_res_dir
end program main_diagonalization