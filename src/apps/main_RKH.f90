program main_RKH
    use kind_tools
    use constants_tools
    use input_tools
    use dir_tools
    use bspline_tools, only: b_spline
    use orbital_tools, only: basis
    use CAP_tools, only: CAP
    use block_tools, only: block_CS, block_diag_CS
    ! use essential_states, only: setup_essential_states, setup_Floquet_essential_states, find_essential_states
    use adiabatic_potentials, only: setup_U_r,compute_adiabatic_potentials,write_adiabatic_potentials
    use RKH_hamiltonian, only: setup_shift_Coul_Ham, compute_RKH_energies
    implicit none

    character(len=:), allocatable :: input_file
    character(len=:), allocatable :: RKH_res_dir

    type(basis) :: bas
    type(b_spline) :: splines
    type(CAP) :: CAP_c
    ! type(block_diag_CS) :: H_0_block
    type(block_diag_CS) :: S_block
    type(block_CS) :: H_block
    ! type(block_CS) :: dip(-1:1)
    ! type(block_CS) :: D

    ! Read the input, and the input file used for the basis_setup program
    call set_RKH_defaults()
    input_file = get_input_file()
    call get_RKH_input(input_file)
    call get_basis_input(basis_dir//"/basis_input.dat")

    ! Load basis and matrix elements
    call bas%load(basis_dir)
    call splines%load(basis_dir)
    ! call H_0_block%load(basis_dir//"H_diag.dat")
    call S_block%load(basis_dir//"S_diag.dat")
    ! call dip(-1)%load(basis_dir//"D_-1.dat")
    ! call dip(0)%load(basis_dir//"D_0.dat")
    ! call dip(1)%load(basis_dir//"D_1.dat")

    ! Create directory for results and copy inputs
    call make_res_dir(RKH_root_dir,RKH_res_dir)
    call write_basis_input(RKH_res_dir)
    call write_RKH_input(RKH_res_dir)

    ! if (use_essential_states) then
    !     if (z_pol.and.(calc_type/='static')) then
    !         call setup_Floquet_essential_states(ess_states,n_ess,target_blocks,n_blocks,target_Floquet_blocks,bas)
    !     else
    !         call setup_essential_states(ess_states,n_ess,target_blocks,bas)
    !     end if
    !     call find_essential_states(ess_states,n_ess,targets,H_0_block,S_block)
    ! end if

    ! if (reduce_basis) then
    !     if (use_essential_states) then
    !         call red_basis(bas,H_0_block,S_block,dip,ess_states)
    !     else
    !         call red_basis(bas,H_0_block,S_block,dip)
    !     end if
    ! else
    !     precond_blocks = precond_blocks_in
    ! end if

    ! if ((block_precond).and.(block_precond_type=='PQ')) then
    !     if (use_essential_states) then
    !         call reorder_blocks(bas,H_0_block,S_block,dip,ess_states)
    !     else
    !         call reorder_blocks(bas,H_0_block,S_block,dip)
    !     end if
    ! end if

    ! Store information about the actual basis that was used
    call bas%store(RKH_res_dir)
    call splines%store(RKH_res_dir)

    call setup_U_r(bas,RKH_intensity,RKH_omega)
    call compute_adiabatic_potentials()
    call write_adiabatic_potentials(RKH_res_dir)

    call CAP_c%init(CAP_order,CAP_r_0,CAP_eta)
    call setup_shift_Coul_Ham(H_block,S_block,RKH_omega,gs_energy,bas,CAP_c,k_GL,splines,B_subset)
    call compute_RKH_energies(H_block,S_block,gs_energy,bas,n_quasi,splines,RKH_res_dir)


    write(6,*) ""
    write(6,*) "Output was written to ", RKH_res_dir
end program main_RKH